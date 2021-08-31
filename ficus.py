 
 ######################################################################
 
 ###   FiCUS: FItting the stellar Continuum of Uv Spectra (FiCUS)   ### 
 
 ######################################################################

 
 ### ficus.py ###

""" "ficus.py" -> main script. 
     Fit the observed input spectra with the best combination of SSPs models, and return 
     the mean light-fractions (X_i) and attenuation parameter (E_BV, uniform screen of dust). 
     Additionally, a bunch of different SED-derived parameters are calculated. 
"""

# ----------------------------------- #
#  Python3.7 packages and libraries:  #
# ----------------------------------- #

""" basics, Python3.7... 
"""
import numpy as np
import os
import sys
import datetime

""" basics, astropy...
"""
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
from astropy import convolution

""" basics, matplotlib... 
"""
import matplotlib
#matplotlib.use('Qt5Agg') #for MacOS graphical outputs
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages

""" lmfit()...
"""
from scipy.optimize import least_squares
from lmfit import Parameters, minimize

""" to handle with configuration files (.ini)...
"""
import configparser

""" extrenal "ficus.py" modules and functions...
"""
from ficus_functions import *


# -------------------------------------- #
#  General variables and working paths:  #
# -------------------------------------- #

""" set the date, time, path and create output directories:
"""
now = datetime.datetime.now();
path = os.getcwd();
try:
    os.mkdir(path+'/outputs/%s_ficus_OutputFiles' %(now.strftime('%Y%m%d')));
    print(path+'/outputs/%s_ficus_OutputFiles created...' %(now.strftime('%Y%m%d')))
except OSError:
    print(path+'/outputs/%s_ficus_OutputFiles already exists...' %(now.strftime('%Y%m%d')))

""" set the seed:
"""
np.random.seed(1234567);

""" speed-of-light (km/s):
"""
c = 2.99e5;


# ------------------------------------------ #
#  Configuration file and input parameters:  #
# ------------------------------------------ #

""" Read the configuration file (.ini) and initialize the input parameters:
"""
config = configparser.ConfigParser();
config.read('ficus.ini');

spec_name = config['ficus']['spec_name'];
plot_mode = config['ficus']['plot_mode'];
ssp_models = config['ficus']['ssp_models'];
Zarray = config['ficus']['Zarray'].split(',');
att_law = config['ficus']['att_law'];
att_igm = config['ficus']['att_igm'];
wave_range = np.float_(config['ficus']['wave_range'].split(','));
z_spec = np.float_(config['ficus']['z_spec']);
wave_R = np.float_(config['ficus']['wave_R']);
wave_norm = np.float_(config['ficus']['wave_norm'].split(','));
nsim = np.int(np.float_(config['ficus']['nsim']));


# ------------------------------- #
#  Fit and results (FiCUS core):  #
# ------------------------------- #

# ### Load the data and SSPs models, run the fitting algorithm
def ficus(path, spec_name, plot_mode, ssp_models, Zarray, att_law, att_igm, wave_range, z_spec, wave_R, wave_norm, nsim):
    """ Load the normalized spectrum and rest-frame wavelength arrays, 
        apply mask-array to spectrum: load_spec()
    """
    wave, flux_norm, err_norm, norm_factor, mask_array, normID = load_spec(path, spec_name, wave_range, z_spec, wave_R, wave_norm);
    
    """ Load SSPs bases: load_ssp_bases()
    """
    wl_model, models_array = load_ssp_bases(ssp_models, Zarray, wave_norm);
    wl_full, models_full = load_ssp_bases(ssp_models, Zarray, wave_norm, fullSED=True);
    
    
    """ Give the resolution of the data and the models:
    """ 
    if ssp_models == 'sb99':
        if wave_norm[0] <= 1200.:
            R_mod = 2500.;
        else:
            R_mod = 4000.;
    if ssp_models == 'bpass':
        R_mod = 1000.; ### CHECK!
    
    #R_obs = specR(wave, wave_R);
    R_obs = 1100.;
    
    """ Convolve models to intrumental resolution; otherwise convolve 
        data to theoretical resolution: spec2R_conv()
    """
    custom_lib, custom_data, custom_error = spec2R_conv(R_mod, R_obs, wave_R, z_spec, wave, flux_norm, err_norm, wl_model, models_array);
    
    """ Apply mask-array to models:
    """ 
    custom_lib_orig = custom_lib.copy();
    for lib in custom_lib:
        lib[mask_array==0.] = 0.;
    
    
    """ Initialize MC-chains:
    """ 
    params_array = np.zeros((10*len(Zarray)+1+3,nsim));
    obs_spec, sed_model, sed_fullmod = [], [], [];
    SEDparams = np.zeros((25,nsim));
    
    
    # ------------------------------------------------------------ #
    #  Run the (nsim) MC iterations of the fit: lmfit.minimize():  #
    # ------------------------------------------------------------ #
    
    for ns in np.arange(0,nsim,1):
        """ Apply mask to flux and error arrays...
        """ 
        if ns == 0:
            flux_normSIM = custom_data.copy();
            err_normSIM = custom_error.copy();
        else:
            flux_normSIM = np.random.normal(custom_data,custom_error);
            err_normSIM = custom_error.copy();
        
        flux_normSIM[mask_array==0.] = 0.;
        err_normSIM[mask_array==0.] = 999.e6;
        
        """ Run fitting algorithm...
        """ 
        params = Parameters();
        for n in range(10*len(Zarray)):
            params.add('X%s' %(n), value=0.1, min=0., max=10.);
        params.add('ebv', value=0.1, min=1e-6, max=1.);
        
        fit_ = minimize(residuals, params=params, args=(wave, custom_lib, att_law, flux_normSIM, err_normSIM), method='least_squares', nan_policy='omit', scale_covar=True, reduce_fcn='neglogcauchy');
        
        # reduced chi-2
        chi2 = fit_.chisqr/(len(mask_array[mask_array==1.])-len(params));
        
        """ Light-weighted parameters...
        """
        # param_array
        param_array = np.array(list(fit_.params.valuesdict().values()));
        # luminosity-weighted metallicity (Zo)
        Z_dict = {'001':1/20., '004':1/5., '008':2/5., '02':1., '04':2.};
        Z_set = [];
        for Z in range(len(Zarray)):
            Z_set.append(np.ones(10)*Z_dict[str(Zarray[Z])]);
        Z_w = np.sum(np.array(Z_set).flatten()*param_array[0:10*len(Zarray)])/np.sum(param_array[0:10*len(Zarray)]);
        # luminosity weighted age (Myr)
        age_array = np.array([1., 2., 3., 4., 5., 8., 10., 15., 20., 40.]*len(Zarray)).flatten();
        age_w = np.sum(age_array*param_array[0:10*len(Zarray)])/np.sum(param_array[0:10*len(Zarray)]);
        
        """ Finally, prepare the resulting parameters and SEDs...
        """ 
        params_array[:,ns] = np.hstack([chi2, Z_w, age_w, list(fit_.params.valuesdict().values())]);
        obs_spec.append([wave,flux_norm]);
        
        if att_law == 'reddy16':
            kl = 2.191 + 0.974/(wave*1e-4);
        else:
            kl = 0.17684325 + 1.44022523/(wave*1e-4) + 0.08560055/(wave*1e-4)**2;
        hRspec = ssp2obs(fit_.params, wave, custom_lib_orig, att_law);
        sed_model.append([wave,hRspec*10**(0.4*kl*params_array[-1,ns]),hRspec]);
        
        if att_law == 'reddy16':
            kl = R16(wl_full);
        else:
            kl = 0.17684325 + 1.44022523/(wl_full*1e-4) + 0.08560055/(wl_full*1e-4)**2;
        lRspec = ssp2obs(fit_.params, wl_full, models_full, att_law, fullSED=True);
        sed_fullmod.append([wl_full,lRspec*10**(0.4*kl*params_array[-1,ns]),lRspec]);
        
        """ Other secondary SED-derived parameters...
        """
        SEDparams[:,ns] = sed_params(wl_full,np.array(sed_fullmod)[-1,2],np.array(sed_fullmod)[-1,1],z_spec,norm_factor);
        
        
    """ Obtaining the median and standard deviation of the output parameters:
    """ 
    params_output = np.c_[params_array[:,0], np.nanstd(params_array, axis=1)];
    SEDparams_output = np.c_[SEDparams[:,0], np.nanstd(SEDparams, axis=1)];
    
    """ Save MC results into (.npy) files:
    """
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_fit.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(params_array));
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_par.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(SEDparams));
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_obs.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(obs_spec));
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_SED.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(sed_model));
    np.save(path+'/outputs/%s_ficus_OutputFiles/%s_ficusMC_SEDfull.npy' %(now.strftime('%Y%m%d'), spec_name), np.array(sed_fullmod));
    
    
    # ------------------- #
    #  Plotting options:  #
    # ------------------- #
    
    if plot_mode == 'yes':
        cont = np.array(sed_model)[0,2];
        chi2_def = params_output[0,:];
        Zw_def = params_output[1,:];
        agew_def = params_output[2,:];
        agesws_def = params_output[3:10*len(Zarray)+3,0]/np.sum(params_output[3:10*len(Zarray)+3,0]);
        ebv_def = params_output[-1,:];
        
        """ Save plot in (.pdf) formtat: ficus_plot()
        """
        pdf_file = PdfPages(path+'/outputs/%s_ficus_OutputFiles/%s.pdf' %(now.strftime('%Y%m%d'), spec_name));
        ficus_plot(pdf_file, spec_name, z_spec, wave, custom_data, custom_error, normID, mask_array, ebv_def, cont, agew_def, Zw_def, agesws_def, age_array, Zarray, Z_set, chi2_def);
        pdf_file.close();
    
    
    # ------------------------------------------ #
    #  Save data and results into (.txt) files:  #
    # ------------------------------------------ #
    
    """ (1) Derived basic parameters from the fit (light-fractions and E_BV):
    """
    np.savetxt(path+'/outputs/%s_ficus_OutputFiles/%s_ficus_fit.txt' %(now.strftime('%Y%m%d'), spec_name), params_output, 
fmt='  '.join(['%+12.5e'] + ['%+12.5e']),
header = '  '.join(['%-12s']+ ['%-12s']) %('# param', 'param.err'), 
delimiter='\t', 
comments= '### Alberto Saldana-Lopez (Obs. Geneva - UniGE)\n' + '# ' + '%s' %(now.strftime('%Y/%m/%d %H:%M:%S')) + '\n' + '#\n'
+    '# --------------------------------------------------------------------------------\n'
+    '#    %s - stellar continuum SED (light-fractions) - ficus.py\n' %(spec_name)
+    '# --------------------------------------------------------------------------------\n'
+ '# \n'
+ '### inputs ### \n'
+ '# spec_name   --> %s\n' %(spec_name)
+ '# plot_mode   --> %s\n' %(plot_mode)
+ '# ssp_models  --> %s\n' %(ssp_models)
+ '# Zarray      --> %s\n' %(Zarray)
+ '# att_law     --> %s\n' %(att_law)
+ '# att_igm     --> %s\n' %(att_igm)
+ '# wave_range  --> %s\n' %(wave_range)
+ '# z_spec      --> %s\n' %(z_spec)
+ '# wave_R      --> %s\n' %(wave_R)
+ '# wave_norm   --> %s\n' %(wave_norm)
+ '# nsim        --> %s\n' %(nsim)
+ '# \n'
+ '### outputs ### \n'
+ '# norm_factor --> %s\n' %(norm_factor)
+ '# params      --> chi^2, Z(Zo), Age(Myr), X_i + E_BV(mag.) \n'
+ '# \n'
+ '# \n');
    
    
    """ (2) Secondary SED-derived parameters:
    """
    param = np.hstack([params_output[0,0], params_output[1,0], params_output[2,0], params_output[-1,0], SEDparams_output[:,0]]);
    param_err = np.hstack([params_output[0,1], params_output[1,1], params_output[2,1], params_output[-1,1], SEDparams_output[:,1]]);
    
    row_names = np.array(['chi2_nu', 'Z', 'Age', 'EBVuv', 'f1100', 'f1100int', 'f1500', 'f1500int', 'M1500',  'M1500int',
                          'beta1200', 'beta1200int', 'beta1550', 'beta1550int', 'beta2000', 'beta2000int', 
                          'f500f1500int', 'f700f1500int', 'f900f1500int', 'f1100f1500',  'f1100f1500int', 
                          'f900f1100', 'f900f1100int', 'Q_H', 'I_Hbeta', 'xi-ion', 'IHbImod', 'fLyCmod', 'fLyCmodINT']);
    col_names = np.array(['param.name', 'param.value', 'param.error']);
    
    tab = np.c_[row_names, param, param_err];
    tab_comments=''' 
 # Alberto Saldana-Lopez (Obs. Geneva - UniGE)
 # %s
 # --------------------------------------------------------------------------------
 #    %s - stellar continuum SED (derived parameters) - ficus.py                   
 # --------------------------------------------------------------------------------
 # 
 ### inputs ### 
 # spec_name   --> %s
 # plot_mode   --> %s
 # ssp_models  --> %s
 # Zarray      --> %s
 # att_law     --> %s
 # att_igm     --> %s
 # wave_range  --> %s
 # z_spec      --> %s
 # wave_R      --> %s
 # wave_norm   --> %s
 # nsim        --> %s
 # 
 ### outputs ### 
 # norm_factor --> %s
 # params      -->  
 # --------------------------------------------------------------------------------
 # Column              Units                      Description                      
 # --------------------------------------------------------------------------------
 #
 # chi2_nu                                        reduced -chi- squared value, 
 # Z                  (Zo, solar)                 light-weighted stellar metallicity, 
 # Age                (Myr)                       light-weighted stellar age, 
 # EBVuv              (mag.)                      dust-attenuation parameter (B-V color excess assuming a Uniform Foreground Screen), 
 # f1100              (1e-17 erg/s/cm2/AA)        observed flux density modelled at 1100\AA, 
 # f1100int           (1e-17 erg/s/cm2/AA)        intrinsic flux density modelled at 1100\AA, 
 # f1500              (1e-17 erg/s/cm2/AA)        observed flux density modelled at 1500\AA, 
 # f1500int           (1e-17 erg/s/cm2/AA)        intrinsic flux density modelled at 1500\AA, 
 # M1500              (mag.)                      derived observed absolute AB magnitude at 1500\AA, 
 # M1500int           (mag.)                      derived intrinsic absolute AB magnitude at 1500\AA, 
 # beta1200                                       observed UV beta slope modelled around 1200\AA, 
 # beta1200int                                    intrinsic UV beta slope modelled around 1200\AA, 
 # beta1550                                       observed UV beta slope modelled around 1550\AA, 
 # beta1550int                                    intrinsic UV beta slope modelled around 1550\AA, 
 # beta2000                                       observed UV beta slope modelled around 2000\AA, 
 # beta2000int                                    intrinsic UV beta slope modelled around 2000\AA, 
 # f500f1500int                                   intrinsic 500-to-1500\AA flux ratio (in F\lambda), 
 # f700f1500int                                   intrinsic 700-to-1500\AA flux ratio (in F\lambda), 
 # f900f1500int                                   intrinsic 900-to-1500\AA flux ratio (in F\lambda), 
 # f1100f1500                                     observed 1100-to-1500\AA flux ratio (in F\lambda), 
 # f1100f1500int                                  intrinsic 1100-to-1500\AA flux ratio (in F\lambda), 
 # f900f1100                                      observed 900-to-1100\AA flux ratio (in F\lambda), 
 # f900f1100int                                   intrinsic 900-to-1100\AA flux ratio (in F\lambda), 
 # Q_H                (1/s)                       ionizing photon flux Q(H), 
 # xi-ion             (log10 Hz/erg)              ionizing photon production efficiency, 
 # I_Hbeta            (1e-15 erg/s/cm2)           modelled H\beta flux, 
 # IHbImod            (AA)                        modelled H\beta versus LyC flux ratio, 
 # fLyCmod            (1e-17 erg/s/cm2/AA)        modelled LyC flux at LyC window, 
 # fLyCmodINT         (1e-17 erg/s/cm2/AA)        modelled intrinsic LyC flux at LyC window, 
 #
             ''' %(now.strftime('%Y/%m/%d %H:%M:%S'), spec_name, spec_name, plot_mode, ssp_models, Zarray, att_law, att_igm, wave_range, z_spec, wave_R, wave_norm, nsim, norm_factor);
    
    ascii.write(tab, path+'/outputs/%s_ficus_OutputFiles/%s_ficus_par.txt' %(now.strftime('%Y%m%d'), spec_name), format='commented_header', names=col_names, comment=tab_comments, overwrite=True);
    
    show = ascii.read(path+'/outputs/%s_ficus_OutputFiles/%s_ficus_par.txt' %(now.strftime('%Y%m%d'), spec_name));
    
    formats = defaultdict(str)
    formats['param.name'] = '%s';
    for a in col_names[1::]:
        formats[a] = '%5.4f'; 
    ascii.write(show, path+'/outputs/%s_ficus_OutputFiles/%s_ficus_par.txt' %(now.strftime('%Y%m%d'), spec_name), 
                formats=formats, overwrite=True);
    
    
    
    """ (3) Oringal spectrum and fitted stellar continuum (SED), with errors:
    """
    np.savetxt(path+'/outputs/%s_ficus_OutputFiles/%s_ficus_SED.txt' %(now.strftime('%Y%m%d'), spec_name), np.c_[wave, flux_norm, err_norm, np.array(sed_model)[0,2], np.nanstd(np.array(sed_model)[:,2], axis=0)], 
fmt='  '.join(['%+12.5e'] + ['%+12.5e'] + ['%+12.5e'] + ['%+12.5e'] + ['%+12.5e']),
header = '  '.join(['%-12s'] + ['%-12s'] + ['%-12s'] + ['%-12s'] + ['%-12s']) %('# wave(A)', 'Flux.OBS', 'Flux.ERR', 'Flux.FIT', 'Flux.FITerr'),
delimiter='\t', 
comments= '### Alberto Saldana-Lopez (Obs. Geneva - UniGE)\n' + '# ' + '%s' %(now.strftime('%Y/%m/%d %H:%M:%S')) + '\n' + '#\n'
+    '# --------------------------------------------------------------------------------\n'
+    '#    %s - stellar continuum SED - ficus.py\n' %(spec_name)
+    '# --------------------------------------------------------------------------------\n'
+ '# \n'
+ '### inputs ### \n'
+ '# spec_name   --> %s\n' %(spec_name)
+ '# plot_mode   --> %s\n' %(plot_mode)
+ '# ssp_models  --> %s\n' %(ssp_models)
+ '# Zarray      --> %s\n' %(Zarray)
+ '# att_law     --> %s\n' %(att_law)
+ '# att_igm     --> %s\n' %(att_igm)
+ '# wave_range  --> %s\n' %(wave_range)
+ '# z_spec      --> %s\n' %(z_spec)
+ '# wave_R      --> %s\n' %(wave_R)
+ '# wave_norm   --> %s\n' %(wave_norm)
+ '# nsim        --> %s\n' %(nsim)
+ '# \n'
+ '### outputs ### \n'
+ '# norm_factor --> %s\n' %(norm_factor)
+ '# arrays      --> wave(A), Flux.OBS, Flux.ERR, Flux.FIT, Flux.FITerr \n'
+ '# \n'
+ '# \n');
    
    
    """ (4) Full SED (dust-attenuated and dust-free), with errors:
    """ 
    np.savetxt(path+'/outputs/%s_ficus_OutputFiles/%s_ficus_SEDfull.txt' %(now.strftime('%Y%m%d'), spec_name), np.c_[wl_full, np.array(sed_fullmod)[0,2], np.nanstd(np.array(sed_fullmod)[:,2], axis=0), np.array(sed_fullmod)[0,1], np.nanstd(np.array(sed_fullmod)[:,1], axis=0)], 
fmt='  '.join(['%+12.5e'] + ['%+12.5e'] + ['%+12.5e'] + ['%+12.5e'] + ['%+12.5e']),
header = '  '.join(['%-12s'] + ['%-12s'] + ['%-12s'] + ['%-12s'] + ['%-12s']) %('# wave(A)', 'Flux.MOD', 'Flux.MODerr', 'Flux.INT', 'Flux.INTerr'),
delimiter='\t', 
comments= '### Alberto Saldana-Lopez (Obs. Geneva - UniGE)\n' + '# ' + '%s' %(now.strftime('%Y/%m/%d %H:%M:%S')) + '\n' + '#\n'
+    '# --------------------------------------------------------------------------------\n'
+    '#    %s - stellar continuum SED (full) - ficus.py\n' %(spec_name)
+    '# --------------------------------------------------------------------------------\n'
+ '# \n'
+ '### inputs ### \n'
+ '# spec_name   --> %s\n' %(spec_name)
+ '# plot_mode   --> %s\n' %(plot_mode)
+ '# ssp_models  --> %s\n' %(ssp_models)
+ '# Zarray      --> %s\n' %(Zarray)
+ '# att_law     --> %s\n' %(att_law)
+ '# att_igm     --> %s\n' %(att_igm)
+ '# wave_range  --> %s\n' %(wave_range)
+ '# z_spec      --> %s\n' %(z_spec)
+ '# wave_R      --> %s\n' %(wave_R)
+ '# wave_norm   --> %s\n' %(wave_norm)
+ '# nsim        --> %s\n' %(nsim)
+ '# \n'
+ '### outputs ### \n'
+ '# norm_factor --> %s\n' %(norm_factor)
+ '# arrays      --> wave(A), Flux.MOD, Flux.MODerr, Flux.INT, Flux.INTerr \n'
+ '# \n'
+ '# \n');
    
    print('   ')
    return print(' # done!')



# ----------------------------- #
#  Deafault "__main__" header:  #
# ----------------------------- #

if __name__ == "__main__":
    print('   ')
    print(' ###################################################################### ')
    print('   ')
    print(' ###   FiCUS: FItting the stellar Continuum of Uv Spectra (FiCUS)   ### ')
    print('   ')
    print(' ###################################################################### ')
    print('   ')
    print('   ')
    print(' ### Running FiCUS (ficus.py) for %s.fits ...' %(spec_name))
    
    to_print = ''' 
 ### inputs ### 
 # spec_name   --> %s
 # plot_mode   --> %s
 # ssp_models  --> %s
 # Zarray      --> %s
 # att_law     --> %s
 # att_igm     --> %s
 # wave_range  --> %s
 # z_spec      --> %s
 # wave_R      --> %s
 # wave_norm   --> %s
 # nsim        --> %s
               ''' %(spec_name, plot_mode, ssp_models, Zarray, att_law, att_igm, wave_range, z_spec, wave_R, wave_norm, nsim);
    
    print(to_print)
    ficus(path, spec_name, plot_mode, ssp_models, Zarray, att_law, att_igm, wave_range, z_spec, wave_R, wave_norm, nsim);
    print('   ')

# EOF
