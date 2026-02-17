 
 ############################################################## 
 
 ###   FiCUS: FItting the stellar Continuum of Uv Spectra   ### 
 
 ############################################################## 
 
 ### ficus.py ### 

""" "ficus.py" -> main script. 
     
     This code fits the observed input SPECTRUM (in F_\lambda units) with a linear combination 
     of single-stellar population MODELS, and returns the best-fit light-fractions of each model 
     (X_i) and attenuation parameter (E_BV, assuming a uniform screen of DUST). 
     Additionally, a bunch of secondary SED parameters are calculated. 
     
     - To be RUN in console as (Python version: 3.7.4, or later):
         > python3.7 ficus-path/FiCUS/ficus.py SPEC-NAME REDSHIFT 
     
     - The INPUT ".fits" file must contain:
         WAVE     > observed-frame wavelength array, in \AA, 
         FLUX     > spectral flux-density, in F_\lambda units (erg/s/cm2/AA), 
         FLUX_ERR > 1\sigma error on the spectral flux-density, 
         MASK     > mask array (0 = masked, 1 = un-masked), 
       
       (see "FiCUS/examples/CDFS017345.fits"). 
     
     - The user parameters are defined in the CONFIGURATION "ficus.ini" file:
         ... ...
         ssp_models > pick your preferred stellar LIBRARY 
                      [sb99 (Starburst99, Leitherer et al. 2011; ApJS, 189, 2) | 
                       bpass (BPASSv2.2.1, Eldridge et al. 2017; PASA, 34, E058) | 
                       sb99stripped (Gotberg et al. 2023; ApJ, 959, 2)], 
                       
         neb_mode   > activate or desactivate NEBULAR CONTIUUM contribution
                      [on | off], 
         
         zarray     > specify a set of METALLICITIES 
                      [0001,0002,0008,0014 (standing for 1/10, 1/5, 1/2 and 1 Z_\sun)], 
                      
         att_law    > choose the DUST attenuation law 
                      [r16 (Reddy et al. 2016; ApJ, 828, 2), 
                       smc (Prevot et al. 1994; A&A, 132, 389-392)], 
         
         wave_range > rest-frame WAVELENGTH range to be considered in the fit 
                      [e.g., 1200.,1920. (\lambda_min, \lambda_max; in \AA)], 
         
         wave_norm  > rest-frame wavelength interval for spectral NORMALIZATION 
                      [e.g., 1350.,1370. (\lambda_min, \lambda_max; in \AA)], 
         
         r_obs      > instrumental RESOLUTION of the input spectra 
                      [e.g., 600.; as R = (\Delta \lambda) / \lambda], 
         
         nsim       > number of Monte-Carlo (MC) ITERATIONS 
                      [e.g., 100.], 
         
         plot_mode  > activate or desactivate PLOTTING mode 
                      [yes | no].
         ... ...
     
     (see "FiCUS/examples/CDFS017345.ini"). 
     
     - List of OUTPUT ".txt" files:
         *_ficus_fit.txt     > best-fit reduced chi-squared [chi^2], average light-weighted stellar 
                               metallicity [Z(Zo)] and age [Age(Myr)], light-fractions [X_i] and
                               dust attenuation parameter [E_BV(mag.)], with errors, 
         
         *_ficus_par.txt     > secondary SED parameters, with errors, 
         
         *_ficus_SED.txt     > original spectrum and best-fit stellar continuum, with errors, 
         
         *_ficus_SEDfull.txt > same as "*_ficus_SED.txt", but NOT restricted to "wave_range", 
         
         *.pdf               > ".pdf" file with figure (only when "plot_mode = yes"), 
       
       (* = SPEC-NAME of the INPUT file). 
       
""" 

# ----------------------------------- #
#  Python3.7 packages and libraries:  #
# ----------------------------------- #

""" basics, Python3.7... 
""" 
import datetime
import numpy as np
import os
import sys

""" basics, astropy...
""" 
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, Column
from astropy import convolution

""" basics, matplotlib... 
""" 
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#matplotlib.use('Qt5Agg') #for MacOS graphical outputs
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#matplotlib.rcParams['text.usetex']=True; #for LaTeX format

""" lmfit() - see https://lmfit.github.io/lmfit-py/...
""" 
import lmfit
from lmfit import Parameters, minimize
from scipy.optimize import least_squares

""" configuration file (.ini)...
""" 
import configparser

""" external modules and functions...
""" 
from ficus_scripts import *

# -------------------------------------- #
#  General variables and working paths:  #
# -------------------------------------- #

""" set the date, path and create output directories:
""" 
now = datetime.datetime.now();
path = os.getcwd();
try:
    os.mkdir(path+'/outputs/%s_ficus_OutputFiles' %(now.strftime('%Y%m%d')));
    print(path+'/outputs/%s_ficus_OutputFiles created...' %(now.strftime('%Y%m%d')))
except OSError:
    print(path+'/outputs/%s_ficus_OutputFiles already exists...' %(now.strftime('%Y%m%d')))

""" initialize the seed:
""" 
np.random.seed(1234567);

""" define speed-of-light (km/s):
""" 
c = 2.99e5;

# ------------------------------------------ #
#  Configuration file and input parameters:  #
# ------------------------------------------ #

""" Command-line input arguments:
""" 
spec_name, z_spec = str(sys.argv[1]), float(sys.argv[2]);

""" Read the CONFIGURATION file (.ini) and initialize the input parameters:
""" 
config = configparser.ConfigParser();
config.read('ficus.ini');

ssp_models = config['ficus']['ssp_models'];
neb_mode = config['ficus']['neb_mode'];
Zarray = config['ficus']['Zarray'].split(',');
att_law = config['ficus']['att_law'];
wave_range = np.float_(config['ficus']['wave_range'].split(','));
wave_norm = np.float_(config['ficus']['wave_norm'].split(','));
R_obs = float(config['ficus']['r_obs']);
nsim = int(np.float_(config['ficus']['nsim']));
plot_mode = config['ficus']['plot_mode'];

# ------------------------------- #
#  Fit and results (FiCUS core):  #
# ------------------------------- #

def ficus(path, spec_name, ssp_models, neb_mode, Zarray, att_law, wave_range, z_spec, wave_norm, R_obs, nsim, plot_mode):
    
    """ Load normalized SPECTRUM and rest-frame WAVELENGTH arrays (in \AA); 
        apply MASK to spectrum: ficus_scripts > load_spec()
    """ 
    wave, flux_norm, err_norm, norm_factor, mask_array, normID = load_spec(path, spec_name, wave_range, z_spec, wave_norm);
    
    """ Load MODELS: ficus_scripts > load_ssp_bases()
    """ 
    wl_model, models_array = load_models(ssp_models, neb_mode, Zarray, wave_norm);
    wl_full, models_full = load_models(ssp_models, neb_mode, Zarray, wave_norm, fullSED=True);
    
    """ Set MODELS resolution:
    """ 
    R_mod = np.nanmean(wave_norm) / dwave(wl_model, np.nanmean(wave_norm));

    """ Convolve MODELS to intrumental resolution; or convolve 
        DATA to theoretical resolution: ficus_scripts > model_conv()
    """ 
    custom_lib, custom_data, custom_error = model_conv(R_mod, R_obs, wave_norm, wave, flux_norm, err_norm, wl_model, models_array);
    
    """ Apply MASK to models:
    """ 
    custom_lib_orig = custom_lib.copy();
    for lib in custom_lib:
        lib[mask_array==0.] = 0.;
    
    """ Initialize MC-chains:
    """ 
    params_array = np.zeros((10*len(Zarray)+1+3,nsim));
    obs_spec, sed_model, sed_hRmodel, sed_fullmod = [], [], [], [];
    SEDparams = np.zeros((26,nsim));
    
    # ------------------------------------- #
    #  Run MC iterations: lmfit.minimize()  #
    # ------------------------------------- #
    
    for ns in np.arange(0,nsim,1):
        
        """ Apply MASK to FLUX and ERROR arrays and sample randomly...
        """ 
        if ns == 0:
            flux_normSIM = custom_data.copy();
            err_normSIM = custom_error.copy();
        else:
            flux_normSIM = np.random.normal(custom_data, custom_error);
            err_normSIM = custom_error.copy();
        
        obs_spec.append([wave,flux_normSIM*norm_factor]);
        
        flux_normSIM[mask_array==0.] = 0.;
        err_normSIM[mask_array==0.] = 999e6;
        
        """ Run fitting algorithm...
        """ 
        params = Parameters();
        for n in range(10*len(Zarray)):
            params.add('X%s' %(n), value=0.1, min=0., max=10.);
        params.add('ebv', value=0.1, min=0., max=0.5);
        
        fit_ = minimize(residuals, params=params, args=(wave, custom_lib, att_law, flux_normSIM, err_normSIM), method='leastsq', nan_policy='omit', scale_covar=True, reduce_fcn='neglogcauchy');
        
        # reduced chi-2
        chi2 = fit_.redchi;
        
        """ Average light-weighted AGES and METALLICITIES...
        """ 
        param_array = np.array(list(fit_.params.valuesdict().values()))[0:-1];
        
        # luminosity-weighted metallicity (Zo)
        Z_dict = {'0001': 1/10, '0002': 1/5, '0008': 1/2, '0014': 1., '0040': 2.};
        
        Z_set = [];
        for Z in range(len(Zarray)):
            Z_set.append(np.ones(10)*Z_dict[str(Zarray[Z])]);
        Z_w = np.sum(np.array(Z_set).flatten()*param_array[0:10*len(Zarray)])/np.sum(param_array[0:10*len(Zarray)]);
        
        # luminosity-weighted age (Myr)
        age_array = np.array([1., 2., 3., 4., 5., 8., 10., 15., 20., 40.]*len(Zarray)).flatten();
        age_w = np.sum(age_array*param_array[0:10*len(Zarray)])/np.sum(param_array[0:10*len(Zarray)]);
        
        """ SED parameters...
        """ 
        params_array[:,ns] = np.hstack([chi2, Z_w, age_w, np.array(list(fit_.params.valuesdict().values()))]);
        
        if att_law == 'r16':
            kl, klfull = R16(wave), R16full(wl_full);
        elif atta_law == 'smc':
            kl, klfull = SMC(wave), SMC(wl_full);
        else:
            print('/!\ Warning /!\: %s not a valid extinction law' %att_law)
        
        hRspec = ssp2obs(fit_.params, wave, custom_lib_orig, att_law);
        sed_model.append([wave, hRspec*10**(0.4*kl*params_array[-1,ns])*norm_factor, hRspec*norm_factor]);
        
        lRspec = ssp2obs(fit_.params, wl_full, models_full, att_law, fullSED=True);
        sed_fullmod.append([wl_full, lRspec*10**(0.4*klfull*params_array[-1,ns])*norm_factor, lRspec*norm_factor]);
        
        SEDparams[:,ns] = sed_params(wl_full, np.array(sed_fullmod)[-1,2], np.array(sed_fullmod)[-1,1], z_spec);
    
    """ Store OUTPUT parameters, with errors:
    """ 
    params_output = np.c_[params_array[:,0], np.nanstd(params_array, axis=1)];
    SEDparams_output = np.c_[SEDparams[:,0], np.nanstd(SEDparams, axis=1)];
    
    """ Save MC results into (.npy) files:
    """ 
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_lightfracs.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(params_array));
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_params.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(SEDparams));
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_obs.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(obs_spec));
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_SPEC.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(sed_model));
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_SED.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(sed_fullmod));
    
    # ------------------- #
    #  Plotting options:  #
    # ------------------- #
    
    if plot_mode == 'yes':
        
        chi2_def, Zw_def, agew_def, ebv_def = params_output[0,:], params_output[1,:], params_output[2,:], params_output[-1,:];
        agesws_def = params_output[3:10*len(Zarray)+3,0]/np.sum(params_output[3:10*len(Zarray)+3,0]);
        
        """ Save PLOT in (.pdf) format: ficus_scripts > ficus_plot()
        """ 
        pdf_file = PdfPages(path+'/outputs/%s_ficus_OutputFiles/%s.pdf' %(now.strftime('%Y%m%d'), spec_name));
        ficus_plot(pdf_file, spec_name, z_spec, wave, custom_data, custom_error, normID, mask_array, att_law, ebv_def, np.array(sed_model)[0,2]/norm_factor, agew_def, Zw_def, agesws_def, age_array, Zarray, Z_set, chi2_def);
        pdf_file.close();
    
    # ---------------------------------------- #
    #  Save data and results to (.txt) files:  #
    # ---------------------------------------- #

    """ (1) Best-fit parameters (light-fractions), with errors:
    """ 
    col_names = np.array(['param.name', 'param.value', 'param.error']);
    age_labels = ['1Myr', '2Myr', '3Myr', '4Myr', '5Myr', '8Myr', '10Myr', '15Myr', '20Myr', '40Myr'];
    
    Xi_labels = [];
    for Z in Zarray:
        X_array = [Z+'Z-'+age for age in age_labels];
        Xi_labels.append(X_array);
    Xi_labels = np.array(Xi_labels).flatten();
    
    tab = np.c_[Xi_labels, params_output[3:-1,0], params_output[3:-1,1]];
    tab_comments=''' 
 # %s
 # ---------------------------------------------------------------------------
 #    %s - stellar continuum (light-fractions) - ficus.py                   
 # ---------------------------------------------------------------------------
 # 
 ### inputs ### 
 # spec_name   --> %s
 # z_spec      --> %s
 # ssp_models  --> %s
 # neb_mode    --> %s
 # Zarray      --> %s
 # att_igm     --> %s
 # wave_range  --> %s
 # wave_norm   --> %s
 # r_obs       --> %s
 # nsim        --> %s
 # plot_mode   --> %s  
 #
             ''' %(now.strftime('%Y/%m/%d %H:%M:%S'), spec_name, spec_name, z_spec, ssp_models, neb_mode, Zarray, att_law, wave_range, wave_norm, R_obs, nsim, plot_mode);
    
    ascii.write(tab, path+'/outputs/%s_ficus_OutputFiles/%s_ficus_lightfracs.txt' %(now.strftime('%Y%m%d'), spec_name), format='commented_header', names=col_names, comment=tab_comments, overwrite=True);
    show = ascii.read(path+'/outputs/%s_ficus_OutputFiles/%s_ficus_lightfracs.txt' %(now.strftime('%Y%m%d'), spec_name));
    
    formats = defaultdict(str);
    formats['param.name'] = '%s';
    for a in col_names[1::]:
        formats[a] = '%5.4f';
    ascii.write(show, path+'/outputs/%s_ficus_OutputFiles/%s_ficus_lightfracs.txt' %(now.strftime('%Y%m%d'), spec_name), 
                formats=formats, overwrite=True);

    """ (2) Secondary SED-derived parameters, with errors:
    """ 
    param = np.hstack([params_output[0,0], params_output[1,0], params_output[2,0], params_output[-1,0], SEDparams_output[:,0]]);
    param_err = np.hstack([params_output[0,1], params_output[1,1], params_output[2,1], params_output[-1,1], SEDparams_output[:,1]]);
    
    row_names = np.array(['chi2_nu', 'Z', 'Age', 'EBVuv', 'f1100', 'f1100int', 'f1500', 'f1500int', 'M1500', 'M1500int',
                          'beta1200', 'beta1200int', 'beta1550', 'beta1550int', 'beta2000', 'beta2000int', 
                          'f500f1500int', 'f700f1500int', 'f900f1500', 'f900f1500int', 
                          'f900f1100', 'f900f1100int', 'f1100f1500',  'f1100f1500int',
                          'QH', 'IHbeta', 'xiion', 'IHbImod', 'fLyCmod', 'fLyCmodINT']);
    col_names = np.array(['param.name', 'param.value', 'param.error']);
    
    tab = np.c_[row_names, param, param_err];
    tab_comments=''' 
 # %s
 # ---------------------------------------------------------------------------
 #    %s - stellar continuum (SED-derived parameters) - ficus.py                   
 # ---------------------------------------------------------------------------
 # 
 ### inputs ### 
 # spec_name   --> %s
 # z_spec      --> %s
 # ssp_models  --> %s
 # neb_mode    --> %s
 # Zarray      --> %s
 # att_igm     --> %s
 # wave_range  --> %s
 # wave_norm   --> %s
 # r_obs       --> %s
 # nsim        --> %s
 # plot_mode   --> %s
 # 
 ### outputs ###   
 # ---------------------------------------------------------------------------
 # Column              Units                      Description                      
 # ---------------------------------------------------------------------------
 # chi2_nu                                        reduced \chi^2 value, 
 # Z                  (Zo, solar)                 light-weighted stellar metallicity, 
 # Age                (Myr)                       light-weighted stellar age, 
 # EBVuv              (mag.)                      dust-attenuation parameter (B-V color excess), 
 # f1100              (1e-18 erg/s/cm2/AA)        observed flux density modeled at 1100\AA, 
 # f1100int           (1e-18 erg/s/cm2/AA)        intrinsic flux density modeled at 1100\AA, 
 # f1500              (1e-18 erg/s/cm2/AA)        observed flux density modeled at 1500\AA, 
 # f1500int           (1e-18 erg/s/cm2/AA)        intrinsic flux density modeled at 1500\AA, 
 # M1500              (AB)                        derived observed absolute AB magnitude at 1500\AA, 
 # M1500int           (AB)                        derived intrinsic absolute AB magnitude at 1500\AA, 
 # beta1200                                       observed UV beta-slope modeled around 1200\AA, 
 # beta1200int                                    intrinsic UV beta-slope modeled around 1200\AA, 
 # beta1550                                       observed UV beta-slope modeled around 1550\AA, 
 # beta1550int                                    intrinsic UV beta-slope modeled around 1550\AA, 
 # beta2000                                       observed UV beta-slope modeled around 2000\AA, 
 # beta2000int                                    intrinsic UV beta-slope modeled around 2000\AA, 
 # f500f1500int                                   intrinsic 500-to-1500\AA flux ratio (in F\lambda), 
 # f700f1500int                                   intrinsic 700-to-1500\AA flux ratio (in F\lambda), 
 # f900f1500                                      observed 900-to-1500\AA flux ratio (in F\lambda), 
 # f900f1500int                                   intrinsic 900-to-1500\AA flux ratio (in F\lambda), 
 # f900f1100                                      observed 900-to-1100\AA flux ratio (in F\lambda), 
 # f900f1100int                                   intrinsic 900-to-1100\AA flux ratio (in F\lambda), 
 # f1100f1500                                     observed 1100-to-1500\AA flux ratio (in F\lambda), 
 # f1100f1500int                                  intrinsic 1100-to-1500\AA flux ratio (in F\lambda), 
 # QH                 (1e+54 1/s)                 intrinsic ionizing photon flux Q(H), 
 # IHbeta             (1e-15 erg/s/cm2)           modeled H\beta flux, 
 # xiion              (log10 Hz/erg)              intrinsic ionizing photon production efficiency, 
 # IHbImod            (AA)                        modeled H\beta versus LyC flux ratio, 
 # fLyCmod            (1e-18 erg/s/cm2/AA)        modeled LyC flux at LyC window, 
 # fLyCmodINT         (1e-18 erg/s/cm2/AA)        modeled intrinsic LyC flux at LyC window, 
 #
             ''' %(now.strftime('%Y/%m/%d %H:%M:%S'), spec_name, spec_name, z_spec, ssp_models, neb_mode, Zarray, att_law, wave_range, wave_norm, R_obs, nsim, plot_mode);
    
    ascii.write(tab, path+'/outputs/%s_ficus_OutputFiles/%s_ficus_params.txt' %(now.strftime('%Y%m%d'), spec_name), format='commented_header', names=col_names, comment=tab_comments, overwrite=True);
    show = ascii.read(path+'/outputs/%s_ficus_OutputFiles/%s_ficus_params.txt' %(now.strftime('%Y%m%d'), spec_name));
    
    formats = defaultdict(str);
    formats['param.name'] = '%s';
    for a in col_names[1::]:
        formats[a] = '%5.4f';
    ascii.write(show, path+'/outputs/%s_ficus_OutputFiles/%s_ficus_params.txt' %(now.strftime('%Y%m%d'), spec_name), 
                formats=formats, overwrite=True);
    
    """ (3) Observed and best-fit stellar spectra, with errors:
    """ 
    np.savetxt(path+'/outputs/%s_ficus_OutputFiles/%s_ficus_SPEC.txt' %(now.strftime('%Y%m%d'), spec_name), np.c_[wave, flux_norm, err_norm, np.array(sed_model)[0,2], np.nanstd(np.array(sed_model)[:,2], axis=0)], 
fmt='  '.join(['%+15.5e'] + ['%+15.5e'] + ['%+15.5e'] + ['%+15.5e'] + ['%+15.5e']),
header = '  '.join(['%-15s'] + ['%-15s'] + ['%-15s'] + ['%-5s'] + ['%-15s']) %('# wave(A)', 'Flux.OBS', 'Flux.ERR', 'Flux.SPEC', 'Flux.SPECerr'),
delimiter='\t', 
comments= '# %s' %now.strftime('%Y/%m/%d %H:%M:%S') + '\n' + '#\n'
+    '# ---------------------------------------------------------------------------\n'
+    '#    %s - stellar continuum spectrum - ficus.py\n' %spec_name
+    '# ---------------------------------------------------------------------------\n'
+ '# \n'
+ '### inputs ### \n'
+ '# spec_name   --> %s\n' %spec_name
+ '# z_spec      --> %s\n' %z_spec
+ '# ssp_models  --> %s\n' %ssp_models
+ '# neb_mode    --> %s\n' %neb_mode
+ '# Zarray      --> %s\n' %Zarray
+ '# att_law     --> %s\n' %att_law
+ '# wave_range  --> %s\n' %wave_range
+ '# wave_norm   --> %s\n' %wave_norm
+ '# r_obs       --> %s\n' %R_obs
+ '# nsim        --> %s\n' %nsim
+ '# plot_mode   --> %s\n' %plot_mode
+ '# \n');

    """ (4) Full best-fit SED (dust-attenuated and dust-free), with errors:
    """ 
    np.savetxt(path+'/outputs/%s_ficus_OutputFiles/%s_ficus_SED.txt' %(now.strftime('%Y%m%d'), spec_name), np.c_[wl_full, np.array(sed_fullmod)[0,2], np.nanstd(np.array(sed_fullmod)[:,2], axis=0), np.array(sed_fullmod)[0,1], np.nanstd(np.array(sed_fullmod)[:,1], axis=0)], 
fmt='  '.join(['%+15.5e'] + ['%+15.5e'] + ['%+15.5e'] + ['%+15.5e'] + ['%+15.5e']),
header = '  '.join(['%-15s'] + ['%-15s'] + ['%-15s'] + ['%-15s'] + ['%-15s']) %('# wave(A)', 'Flux.SED', 'Flux.SEDerr', 'Flux.SEDINT', 'Flux.SEDINTerr'),
delimiter='\t', 
comments= '# %s' %now.strftime('%Y/%m/%d %H:%M:%S') + '\n' + '#\n'
+    '# ---------------------------------------------------------------------------\n'
+    '#    %s - stellar continuum SED - ficus.py\n' %spec_name
+    '# ---------------------------------------------------------------------------\n'
+ '# \n'
+ '### inputs ### \n'
+ '# spec_name   --> %s\n' %spec_name
+ '# z_spec      --> %s\n' %z_spec
+ '# ssp_models  --> %s\n' %ssp_models
+ '# neb_mode    --> %s\n' %neb_mode
+ '# Zarray      --> %s\n' %Zarray
+ '# att_law     --> %s\n' %att_law
+ '# wave_range  --> %s\n' %wave_range
+ '# wave_norm   --> %s\n' %wave_norm
+ '# r_obs       --> %s\n' %R_obs
+ '# nsim        --> %s\n' %nsim
+ '# plot_mode   --> %s\n' %plot_mode
+ '# \n');

    print('   ')
    return print(' # done!')

# ----------------------------- #
#  Deafault "__main__" header:  #
# ----------------------------- #

if __name__ == "__main__":
    
    print(r"""
  ███████╗██╗ ██████╗██╗   ██╗███████╗                                
  ██╔════╝██║██╔════╝██║   ██║██╔════╝                                
  █████╗  ██║██║     ██║   ██║███████╗                                
  ██╔══╝  ██║██║     ██║   ██║╚════██║                                
  ██║     ██║╚██████╗╚██████╔╝███████║                                
  ╚═╝     ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝       
     FItting the stellar Continuum of Uv Spectra                        
    """)
    
    print('   ')
    print(' ### Running FiCUS (ficus.py) for %s.fits ...' %spec_name)
    
    to_print = ''' 
 ### inputs ### 
 # spec_name   --> %s
 # z_spec      --> %s
 # ssp_models  --> %s
 # neb_mode    --> %s
 # Zarray      --> %s
 # att_law     --> %s
 # wave_range  --> %s
 # wave_norm   --> %s
 # r_obs       --> %s
 # nsim        --> %s
 # plot_mode   --> %s
               ''' %(spec_name, z_spec, ssp_models, neb_mode, Zarray, att_law, wave_range, wave_norm, R_obs, nsim, plot_mode);
    
    print(to_print)
    ficus(path, spec_name, ssp_models, neb_mode, Zarray, att_law, wave_range, z_spec, wave_norm, R_obs, nsim, plot_mode);
    print('   ')

# EOF
