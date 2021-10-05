 
 ### ficus_functions.py ###

""" "ficus_functions.py" -> secondary script. 
     After being called, all the analysis and plotting functions are imported into "ficus.py".
"""

import numpy as np
import scipy as sci
from scipy.integrate import trapz
import os
import sys
import datetime

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
from astropy import convolution
from astropy.cosmology import Planck15 as cosmo

import matplotlib
#matplotlib.use('Qt5Agg') #for MacOS graphical outputs
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages

c = 2.99e5; #km/s



# --------------------------------- #
#  basics for UV spectral analysis  #
# --------------------------------- #

# ### Specific spectral Resolution (R)
def specR(wl_spec, wave_R):
    """ Compute the resolution from the spectra at a
        given input wavelength, wave_R: 'R'
    """
    wave_i = np.argmin(np.abs(wl_spec-wave_R));
    dwave = np.abs(wl_spec[wave_i] - wl_spec[wave_i+1]);
    R = wave_R / dwave;
    
    return R


# ### Velocity dispersion for a given Resolution (R)
def R2vdisp(R):
    """ Calculate the velocity dispersion from R
        (in km/s): 'vdisp'
    """
    vdisp = c / R;
    
    return vdisp


# ### Apply Gaussian convolution to a given spectrum
def spec_convolution(spec_array, sig, xsize):
    """ Convolve a spectral-like array with a Gaussian kernel of a given 
        stdandard deviation 'sig', and size 'xsize': 'spec_conv'
    """
    kernel = convolution.Gaussian1DKernel(stddev=sig, x_size=xsize);
    spec_conv = convolution.convolve(spec_array, kernel, boundary='extend', normalize_kernel=True);
    
    return spec_conv


# ### Reddy et al. 2015 attenuation law
def Reddy2015(lam, rv=False):
    """ https://ui.adsabs.harvard.edu/abs/2015ApJ...806..259R/abstract
    """ 
    RV = 2.505;
    mulam = lam/1e4;
    klam = np.where((mulam<0.60) & (mulam>=0.15), 
                 -5.726 + 4.004/mulam - 0.525/mulam**2 + 0.029/mulam**3 + RV, 0);
    klam = np.where((mulam<2.85) & (mulam>=0.6), 
                 -2.672 - 0.010/mulam + 1.532/mulam**2 - 0.412/mulam**3 + RV, klam);
    if rv:
        return klam, RV;
    else:
        return klam;
    
    
# ### Reddy et al. 2016 attenuation law
def R16(lam,rv=False):
    """ https://ui.adsabs.harvard.edu/abs/2016ApJ...828..107R/abstract
    """ 
    mulam = lam/1.e4;
    klam, RV = Reddy2015(lam, rv=True);
    klam = np.where((mulam<=0.15), 2.191+0.974/mulam, klam);
    #klam = np.where(mulam<0.095, 0., klam);
    klam = np.where(mulam>2.2, 0., klam);
    
    if rv:
        return klam, RV;
    else:
        return klam;



# --------------------------------------------- #
#  functions for handling with data and models  #
# --------------------------------------------- #

# ### Read and prepare the data
def load_spec(path, spec_name, wave_range, z_spec, wave_R, wave_norm):
    """ Retrieve 'WAVE' and 'FLUX' from file
    """
    spec_file = fits.open(path+'/inputs/%s.fits' %(str(spec_name)), mmap=True);
    specTable = Table(spec_file[1].data);
    spec_wave = np.array(specTable['WAVE']);    #\AA
    spec_flux = np.array(specTable['FLUX']);    #flux-density units (normalized or not)
    spec_err = np.array(specTable['FLUX_ERR']); #flux-density units (normalized or not)
    spec_mask = np.array(specTable['MASK']);    #mask-array (0,1 array-like)
    
    """ Define the fitting window and transform to 
        rest-frame wavelength
    """
    I = wave_range[0]*(1.+z_spec);
    F = wave_range[1]*(1.+z_spec);
    if spec_wave[0] > I:
        I = spec_wave[0];
    if spec_wave[-1] < F:
        F = spec_wave[-1];
    fit_region = np.nonzero((spec_wave>I)&(spec_wave<F));
    
    wave = spec_wave[fit_region]/(1.+z_spec);
    flux = spec_flux[fit_region];
    err = spec_err[fit_region];
    mask_array = spec_mask[fit_region];
    
    """ Normalize to 'wave_norm' interval
    """
    normID = np.nonzero((wave>=wave_norm[0])&(wave<=wave_norm[1]))[0];
    norm_factor = np.median(flux[normID]);
    flux_norm = flux/norm_factor;
    err_norm = err/norm_factor;
    
    return wave, flux_norm, err_norm, norm_factor, mask_array, normID


# ### Load SSPs bases
def load_ssp_bases(ssp_models, Zarray, wave_norm, fullSED=False):
    """ Load all the 'Zarray'-SSPs files and normalize them to 
        the 'wave_norm' range. Return s set of normalized 
        SSPs models: 'models_array'
    """
    models_array = [];
    
    models_files = [];
    if ssp_models == 'sb99' and fullSED == True:
        [models_files.append(filename) for filename in os.listdir('inputs/ssp_bases/lowres_sb99/')];
    else:
        [models_files.append(filename) for filename in os.listdir('inputs/ssp_bases/')];
    models_files = np.array(models_files);
    
    for i in np.arange(0,len(Zarray),1):
        file_id = [];
        
        [file_id.append(str(filename.endswith(ssp_models+Zarray[i]+'.fits'))) for filename in models_files];
        file_id = np.array(file_id);
        
        ifile = np.nonzero(file_id=='True')[0];
        if  (ssp_models == 'sb99') and (fullSED == True):
            ssp_file = fits.open('inputs/ssp_bases/lowres_sb99/%s' %(models_files[ifile][0]), mmap=True)[1].data;
            
            ssp = Table(ssp_file);
            norm_factorID = np.nonzero((ssp['WAVE'][0]>=wave_norm[0])&(ssp['WAVE'][0]<=wave_norm[1]))[0];
            for j in range(len(ssp['FLUX'][0])):
                models_array.append(ssp['FLUX'][0][j]/np.median(ssp['FLUX'][0][j][norm_factorID]));
    
        else:
            ssp_file = fits.open('inputs/ssp_bases/%s' %(models_files[ifile][0]), mmap=True)[1].data;
        
            ssp = Table(ssp_file);
            norm_factorID = np.nonzero((ssp['WAVE'][0]>=wave_norm[0])&(ssp['WAVE'][0]<=wave_norm[1]))[0];
            for j in range(len(ssp['FLUX'][0])):
                models_array.append(ssp['FLUX'][0][j]/np.median(ssp['FLUX'][0][j][norm_factorID]));
    
    wl_model = np.array(ssp['WAVE'][0]);
    models_array = np.array(models_array);
        
    return wl_model, models_array



# ----------------------------------- #
#  functions for the fitting routine  #
# ----------------------------------- #

# ### Interpolate and convolve the models to the input instrumental Resolution (R), or viceversa
def spec2R_conv(R_mod, R_obs, wave_R, z_obs, wl_obs, data_array, error_array, wl_model, models_array):
    """ If the observed velocity dispersion 'vdisp_obs' is larger than the modelled one 'vdisp_mod',
        then the models have to be convolved to the instrumental resolution: 'custom_lib'
        Otherwise, the data need to be degraded to the theoretical resolution of 
        the SSP models: 'custom_data', 'custom_error'
    """
    vdisp_mod = R2vdisp(R_mod);
    vdisp_obs = R2vdisp(R_obs);
    
    if vdisp_obs >= vdisp_mod:
        vdisp_add = np.sqrt(vdisp_obs**2 - vdisp_mod**2);
        
        sig = (vdisp_add/vdisp_mod);
        n = np.ceil(3.034854259*sig);
        xsize = int(2*n) + 1;
        
        custom_data, custom_error = data_array, error_array;

        custom_lib = np.zeros((len(models_array),len(wl_obs)));
        for i in np.arange(0,len(models_array),1):
            model_conv = spec_convolution(models_array[i,:], sig, xsize);
            custom_lib[i,:] = np.interp(wl_obs, wl_model, model_conv);
    
    else:
        vdisp_add = np.sqrt(vdisp_mod**2 - vdisp_obs**2);
        
        sig = (vdisp_add/vdisp_obs);
        n = np.ceil(3.034854259*sig);
        xsize = int(2*n) + 1;
        
        custom_data = spec_convolution(data_array, sig, xsize);
        custom_error = np.sqrt(spec_convolution(error_array**2, sig, xsize));
        
        custom_lib = np.zeros((len(models_array),len(wl_obs)));
        for i in np.arange(0,len(models_array),1):
            model_interp = np.interp(wl_obs, wl_model, models_array[i,:]);
            custom_lib[i,:] = model_interp;
    
    """ Now, all data and SSP-libraries have been interpolated 
        into the same wavelength-grid
    """
    return custom_lib, custom_data, custom_error


# ### SSPs weighted combination and attenuation by dust
def ssp2obs(params, wl_model, models, att_law, fullSED=False):
    """ Build the linear dust-attenuated combination of the convolved models
        and return the normalized modelled spectrum: 'spec_model'
    """
    light_fracs = [];
    for n in range(len(params)-1):
        light_fracs.append(params['X%s' %(n)].value);
    light_fracs = np.array(light_fracs);
    ebv = params['ebv'].value;
    
    light_fracs_new = np.hstack([light_fracs, np.zeros(len(models)-len(light_fracs))]);
    spec_model = np.zeros(len(wl_model));
    for i in np.arange(0,len(light_fracs_new),1):
        
        if att_law == 'reddy16':
            if fullSED == True:
                klam = R16(wl_model);
            else:
                klam = 2.191 + 0.974/(wl_model*1e-4);
        else:
            klam = 0.17684325 + 1.44022523/(wl_model*1e-4) + 0.08560055/(wl_model*1e-4)**2;
        
        spec_model = (spec_model + (models[i]*light_fracs_new[i])*10**(-0.4*klam*ebv));
    
    return spec_model


# ### Residuals for fitting algorithm
def residuals(params, wl_model, models, att_law, data, data_err, fullSED=False):
    """ Make use of the ssp2obs() function to return the minimization
        parameter that feeds into lmfit.minimize(): 'chi2'
    """
    resids = data - ssp2obs(params, wl_model, models, att_law);
    chi2 = resids/data_err;
    return chi2


    
# ------------------------ #
#  SED-derived parameters  #
# ------------------------ #

# ### Mean flux through a uniform band-pass (squared filter)
def flam_x(wave,flam,lc,wth):
    """ Flux in F_\lambda units: 'flam_x'
    """
    l1_lim = lc - wth;
    l2_lim = lc + wth;
    w_int = wave[sci.where((wave>=l1_lim)&(wave<=l2_lim))];
    f_int = flam[sci.where((wave>=l1_lim)&(wave<=l2_lim))];
    flam_x = trapz(f_int,w_int) / (2*wth);
    
    return flam_x


# ### Mean flux through a uniform band-pass (squared filter)
def fnu_x(wave,flam,lc,wth):
    """ Flux in F_\nu units: 'fnu_x'
    """
    flam = flam_x(wave,flam,lc,wth);
    fnu_x = lc**2/(c*1e3*1e10)*flam;
    
    return fnu_x


# ### UV beta-slope
def beta_x(wave,flam,lc1,lc2):
    """ Compute the slope between two spectral-regions
        (in log-space): 'beta_x'
    """
    flux1 = flam_x(wave,flam,lc1,20.);
    flux2 = flam_x(wave,flam,lc2,20.);
    if ((flux1*flux2) != 0.):
        beta_x = (sci.log10(flux1)-sci.log10(flux2)) / (sci.log10(lc1)-sci.log10(lc2));
    else:
        beta_x = -99.;
    
    return beta_x


# ### Flux-density to absolute magnitude (AB) conversion
def fnu2MAB(z,fnu):
    """ Flux in F_\nu units, absolute magnitude
        in AB system: 'MAB'
    """
    dL = cosmo.luminosity_distance(z).value;
    if fnu != 0.:
        mAB = -2.5*np.log10(fnu*1e23) + 8.9;
        MAB = mAB - 5.*np.log10(dL*1e6/10.) + 2.5*np.log10(1.+z);
    else:
        MAB = -99.;

    return MAB


# ### Ionizing photon-flux (Q(H))
def QH_IHb(wave,flam):
    """ Ionizing photons-flux (1/s): 'q_h', 
        and H\beta integrated flux (case B, 1e-15 erg/s/cm2): 'IHb'
    """ 
    n_lam = flam * wave/(6.626e-9*2.998);
    eh = 109678.758;
    wedge_h   = 1.e8 / eh;
    
    n_int = n_lam[sci.where(wave<=wedge_h)];
    w_int = wave[sci.where(wave<=wedge_h)];
    
    q_h  = trapz(n_int,w_int);
    IHb = 4.76e-13 * q_h;
    
    return q_h, IHb


# ### Ionizing photon production efficiency (xi_ion)
def xiion(wave,flam):
    """ Ionizing photon production efficiency in 
        units of (log10 Hz/erg): 'xi_ion'
    """ 
    q_h = QH_IHb(wave,flam)[0];
    f_1500 = fnu_x(wave,flam,1500.,20.);
    xi_ion = q_h / f_1500;
    
    return np.log10(xi_ion)


# ### Compilation of ALL the SED derived parameters
def sed_params(wl_obs,flam,flamINT,z,nfactor):
    wave, fMOD, fINT = wl_obs, flam*nfactor, flamINT*nfactor;
    
    """ Fluxes and AB magnitudes: 
    """
    f1100int, f1100obs, f1500int, f1500obs, M1500int, M1500obs = flam_x(wave,fINT,1100.,20.), flam_x(wave,fMOD,1100.,20.), flam_x(wave,fINT,1500.,20.), flam_x(wave,fMOD,1500.,20.), fnu2MAB(z,fnu_x(wave,fINT,1500.,20.)), fnu2MAB(z,fnu_x(wave,fMOD,1500.,20.));
    
    """ UV beta-slopes: 
    """
    beta1200int, beta1200obs, beta1550int, beta1550obs, beta2000int, beta2000obs = beta_x(wave,fINT,1050.,1350.), beta_x(wave,fMOD,1050.,1350.), beta_x(wave,fINT,1300.,1800.), beta_x(wave,fMOD,1300.,1800.), beta_x(wave,fINT,1800.,2200.), beta_x(wave,fMOD,1800.,2200.);
    
    """ UV ionizing vs. non-ionizing flux ratios: 
    """
    f500f1500int, f700f1500int, f900f1500int, f1100f1500int, f1100f1500obs, f900f1100int, f900f1100obs = flam_x(wave,fINT,500.,20.)/f1500int, flam_x(wave,fINT,700.,20.)/f1500int, flam_x(wave,fINT,900.,20.)/f1500int, flam_x(wave,fINT,1100.,20.)/f1500int, flam_x(wave,fMOD,1100.,20.)/f1500obs, flam_x(wave,fINT,900.,20.)/f1100int, flam_x(wave,fMOD,900.,20.)/f1100obs;
    
    """ LyC flux (dust-attenuated and dust-free):
    """
    fLyCmod = flam_x(wave,fMOD,890.,10.);
    fLyCmodINT = flam_x(wave,fINT,890.,10.);
    
    """ Q(H), xi_ion and Hb flux: 
    """
    qh, IHb, xi_ion = QH_IHb(wave,fMOD)[0], QH_IHb(wave,fMOD)[1], xiion(wave,fINT);
    IHbImod = IHb / fLyCmodINT; #AA
    
    return np.hstack([f1100int/1e-17, f1100obs/1e-17, f1500int/1e-17, f1500obs/1e-17, M1500int, M1500obs, beta1200int, beta1200obs, beta1550int, beta1550obs, beta2000int, beta2000obs, f500f1500int, f700f1500int, f900f1500int, f1100f1500int, f1100f1500obs, f900f1100int, f900f1100obs, qh, IHb/1e-15, IHbImod, xi_ion, fLyCmod/1e-17, fLyCmodINT/1e-17])
    

# --------------------------------- #
#  plotting function: ficus_plot()  #
# --------------------------------- #

def ficus_plot(pdf_file, spec_name, z_spec, wave, flux_norm, err_norm, normID, mask_array, ebv, cont, agew_def, Zw_def, agesws_def, age_array, Zarray, Z_set, chi2_def):
    fig = plt.figure(figsize=(15,9));
    gs = fig.add_gridspec(2, 2);
    
    """ Plot-1: observed spectra and fitted continuum (SEDs)
    """
    matplotlib.rcParams['ytick.right'] = True;
    ax1 = fig.add_subplot(gs[0, :]);
    ax1.step(wave, flux_norm, 'k-', where='mid', lw=2, alpha=0.8, solid_capstyle='round', label=r'$\mathrm{data}$');
    ax1.step(wave, err_norm, '-', color='orange', where='mid', lw=2, alpha=0.8, solid_capstyle='round', label=r'$\mathrm{error}$');
    ax1.plot(wave, cont, 'r-', lw=3., alpha=1., solid_capstyle='round', label=r'$\mathrm{S99~best~fit}$');
    
    for a in np.split(np.array(wave), np.nonzero(mask_array==1)[0]):
        if len(a) > 1:
            ax1.fill_between([a[1], a[-1]], -0.2, np.max(flux_norm[normID])*2. ,alpha=0.25, color='gray');
    
    select_lines = {'Names': ['OVI', 'CIII', 'NV', 'OV', 'SiIV', 'CIV', 'HeII'], 'WAVES': [1037., 1175., 1240., 1371., 1400., 1550., 1640.]};
    for nid, wlid in zip(select_lines['Names'],select_lines['WAVES']):
        if wave[0] <= wlid < wave[-1]:
            ax1.plot([wlid, wlid], [np.max(flux_norm[normID])*1.4, np.max(flux_norm[normID])*1.5], ls='-', color='indigo', lw=2.5);
            ax1.text(x=wlid, y=np.max(flux_norm[normID])*1.275, s=r'$\mathrm{%s}$' %(nid), fontsize=14, color='indigo', ha='center', va='center');
    
    select_lines = {'Names': ['Ly6', 'Ly5', 'Ly \delta', 'Ly \gamma', 'SiII', 'Ly {\\beta}', 'SiII', 'Ly {\\alpha}', 'SiII', 'OI+SiII', 'CII', 'FeII', 'AlII'], 'WAVES': [930., 938., 950., 973., 1020., 1026., 1191.5, 1216., 1260., 1302., 1334., 1608., 1670.]};
    for nid, wlid in zip(select_lines['Names'],select_lines['WAVES']):
        if wave[0] <= wlid < wave[-1]:
                ax1.plot([wlid, wlid], [np.max(flux_norm[normID])*1.4, np.max(flux_norm[normID])*1.5], ls='-', color='k', lw=2.);
                ax1.text(x=wlid, y=np.max(flux_norm[normID])*1.175, s=r'$\mathrm{%s}$' %(nid), fontsize=12, color='k', ha='center', va='center');
            
    #ax1.plot([1216./(1+z_spec), 1216./(1+z_spec)], [0.075, 0.2], ls='-', color='k', lw=2.);
    #ax1.text(x=1216./(1+z_spec)-6., y=0.25, s=r'$\mathrm{geoLy\alpha}$', fontsize=14, color='k');
    #ax1.plot([1302./(1+z_spec), 1302./(1+z_spec)], [0.075, 0.2], ls='-', color='k', lw=2.);
    #ax1.text(x=1302./(1+z_spec)-6., y=0.25, s=r'$\mathrm{geoOI}$', fontsize=14, color='k');
    
    ax1.set_xlim(wave[0], wave[-1]);
    ax1.set_ylim(0.,np.max(flux_norm[normID])*1.5);

    ax1.set_xlabel(r'$\mathrm{Rest-frame~wavelength~(\AA)}$', fontsize = 22, labelpad=12);
    ax1.set_ylabel(r'$\mathrm{Normalized~flux}$', fontsize = 22, labelpad=12);
    ax1.set_title(r'$\mathrm{%s~(z = %s) - stellar ~continuum}$' %(spec_name.replace('_', '~'), "{0:.5f}".format(z_spec)), fontsize = 22, pad=18);
    ax1.legend(loc='lower right', fontsize = 15, edgecolor = 'k', framealpha = 1., ncol=3);
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(100.));
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(25.));
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(1.));
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.25));
    matplotlib.rcParams['ytick.right'] = False;
    
    
    """ Plot-2: light-fractions (histos.)
    """
    ax2 = fig.add_subplot(gs[1, 0]);
    Z_dict = {'001':1/20., '004':1/5., '008':2/5., '02':1., '04':2.};
    
    cmap = plt.get_cmap('jet', 201);
    norm = matplotlib.colors.LogNorm(vmin=0.04, vmax=2.1);
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm);
    pxl = np.logspace(np.log10(0.04), np.log10(2.), 201);
    cmap_new = sm.get_cmap();
    new_cmaplist = [cmap_new(np.argmin(np.abs(pxl-i))+1) for i in [1/20., 1/5., 2/5., 1., 2.]];
    
    bottom = np.zeros(10);
    for Z in range(len(Zarray)):
        idx = np.nonzero(np.array(list(Z_dict.values()))==Z_dict[str(Zarray[Z])])[0];
        ax2.bar(age_array[10*Z:10*Z+10], agesws_def[10*Z:10*Z+10]*100., width=0.75, yerr=0., bottom=bottom, align='center', color=new_cmaplist[idx[0]], edgecolor='k', linewidth=1.5, ecolor='k');
        bottom = bottom + agesws_def[10*Z:10*Z+10]*100.;
    
    ax2.axvline(agew_def[0], ls='--', dashes=[6, 2], color='k', linewidth=1.25, solid_capstyle='round');  
    ax2.set_xlim(age_array[0]-1., age_array[-1]+1.);
    ax2.set_ylim(0., bottom.max()+5.);
    ax2.set_xlabel(r'$\mathrm{Stellar~age~(Myr)}$', fontsize = 18, labelpad=12);
    ax2.set_ylabel(r'$\mathrm{Light~fraction~(\%)}$', fontsize = 18, labelpad=12);
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(5.));
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(1.));
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(10.));
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(5.));
    
    inset_ax = inset_axes(ax2,height='6%', width='50%', loc='upper right');
    cab = plt.colorbar(sm, cax=inset_ax, ticks = [1/20., 1/5., 2/5., 1., 2.], orientation='horizontal');
    cab.set_ticklabels([r'$1/20$', r'$1/5$', r'$2/5$', r'$1$', r'$2$']);
    inset_ax.set_xlabel(r'$\mathrm{Z_*~(Z_{\odot})}$', fontsize = 16, labelpad=18, rotation=0);
    for axis in ['top','bottom','left','right']:
        inset_ax.spines[axis].set_linewidth(.75);
    inset_ax.tick_params(labelsize = 16, axis = 'x', pad=10, direction='inout', length = 10, width = .75);
    inset_ax.tick_params(labelsize = 0, axis = 'y', pad=10, direction='inout', length = 0, width = 0);
    inset_ax.tick_params(labelsize = 0, axis = 'x', which = 'minor', pad=10, direction='inout', length = 0, width = 0);
    inset_ax.tick_params(labelsize = 0, axis = 'y', which = 'minor', pad=10, direction='inout', length = 5, width = 0);
    
    
    """ Plot-3: summary chart with the results of the fit
    """
    ax3 = fig.add_subplot(gs[1, 1]);
    ax3.axis('off');
    ax3.text(x=0.01, y=0.625, s=r'$\mathrm{-~\chi_{\nu}^{2}=%s}$' %("{0:.2f}".format(chi2_def[0])), fontsize=18, color='k', transform=ax3.transAxes);
    ax3.text(x=0.01, y=0.425, s=r'$-~E\mathrm{_{B-V}~(mag.)=%s ~\pm~ %s}$' %("{0:.3f}".format(ebv[0]), "{0:.3f}".format(ebv[1])), fontsize=18, color='k', transform=ax3.transAxes);
    ax3.text(x=0.01, y=0.225, s=r'$\mathrm{-~Age~(Myr)=%s ~\pm~ %s}$' %("{0:.2f}".format(agew_def[0]), "{0:.2f}".format(agew_def[1])), fontsize=18, color='k', transform=ax3.transAxes);
    ax3.text(x=0.01, y=0.025, s=r'$\mathrm{-~Z~(Z_{\odot})=%s ~\pm~ %s}$' %("{0:.2f}".format(Zw_def[0]), "{0:.2f}".format(Zw_def[1])), fontsize=18, color='k', transform=ax3.transAxes);
    props = dict(boxstyle='round', facecolor='wheat', alpha=1., edgecolor = 'k', linewidth=1.);
    ax3.text(x=0.5, y=0.85, s=r'$\mathrm{Fitting~parameters}$', 
            fontsize=20, ha='center', va='center', transform=ax3.transAxes, bbox=props, zorder=3);
    
    
    for ax in [ax1, ax2, ax3]:
        for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(.75);
        ax.tick_params(labelsize = 22, axis = 'x', pad=10, direction='inout', length = 10, width = .75);
        ax.tick_params(labelsize = 22, axis = 'y', pad=10, direction='inout', length = 10, width = .75);
        ax.tick_params(labelsize = 22, axis = 'x', which = 'minor', pad=8, direction='inout', length = 5, width = .75);
        ax.tick_params(labelsize = 22, axis = 'y', which = 'minor', pad=8, direction='inout', length = 5, width = .75);
    
    #fig.tight_layout();
    fig.subplots_adjust(hspace=0.45);
    fig.subplots_adjust(wspace=0.20);
    fig.savefig(pdf_file, format='pdf', bbox_inches='tight');
    
    print('   ')
    return print(' ### plotting...   ')

#EOF
