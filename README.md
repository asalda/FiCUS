# FiCUS (FItting the stellar Continuum of Uv Spectra)
- [Description](https://github.com/asalda/FiCUS/README.md#description)
- [Installation](https://github.com/asalda/FiCUS/README.md#install)
- [Input and Configuration files](https://github.com/asalda/FiCUS/README.md#the-input-and-configuration-files)
- [Running FiCUS](https://github.com/asalda/FiCUS/README.md#running-ficus)
- [Outputs and Plots](https://github.com/asalda/FiCUS/README.md#outputs0-and-plots)


## Description
`FiCUS` is a customized `Python` script to fit the stellar continuum of extragalactic ultraviolet (UV) spectra. In short, it takes observed-frame wavelength, flux density (with errors) and user-defined mask arrays, and returns an estimation of the galaxy light-weighted stellar age, metallicity and dust extinction, as well as other secondary Spectral Energy Distribution (SED) parameters. 

The code was presented in [Saldana-Lopez et al. 2022b](https://ui.adsabs.harvard.edu/abs/2022arXiv221101351S/abstract), but the methodology was first described and tested in [Chisholm et al. 2019](https://ui.adsabs.harvard.edu/abs/2022arXiv221101351S/abstract). A previous version of the code has been extensively used in other papers such as [Gazagnes et al. 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..29G/abstract), [Gazagnes et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...639A..85G/abstract) and [Saldana-Lopez et al. 2022a](https://ui.adsabs.harvard.edu/abs/2022A%26A...663A..59S/abstract). 

The UV stellar continuum modeling $F_{\lambda}^{\star}$ is achieved by fitting observed spectra with a linear combination of single-burst stellar population (SSP) theoretical models: `Starburst99`[^1] ([Leitherer et al. 2010](https://ui.adsabs.harvard.edu/abs/2010ApJS..189..309L/abstract)) or `BPASSv2.2.1`[^2] ([Eldridge et al. 2017](https://ui.adsabs.harvard.edu/abs/2017PASA...34...58E/abstract)). These models assume a initial mass function (IMF) with a high-(low-)mass exponent of 2.3 (1.3), and a high-mass cutoff at $100 M_{\odot}$. The models include five different metallicities (0.05, 0.2, 0.4, 1 and $2 Z_{\odot}$) and ten ages for each metallicity (1, 2, 3, 4, 5, 8, 10, 15, 20 and 40 Myr). A nebular continuum was added to every model by self-consistently processing the original SSPs through the `Cloudy v17.0` code[^3] ([Ferland et al. 2017](https://ui.adsabs.harvard.edu/abs/2017RMxAA..53..385F/abstract)), assuming similar gas-phase and stellar metallicities, an ionization parameter of $\log(U)=-2.5$, and a volume hydrogen density of $n_H = 100 cm^{-3}$. Adopting a simple geometry where _all_ the light is attenuated by a uniform foreground slab of gas surrounding the galaxy, this results in: 

$$ F_{\lambda}^{\star} = 10^{-0.4 k_{\lambda} E_{B-V}} \sum_{i,j} X_{ij} F_{\lambda}^{ij} $$

$$ i \equiv 1, 2, 3, 4, 5, 8, 10, 15, 20, 40 Myr $$

$$ i \equiv 0.05, 0.2, 0.4, 1, 2 Z_{\odot} $$

... where $F_{\lambda}^{ij}$ represents the corresponding model at the i-_th_ age and j-_th_ metallicity, and the $X_{ij}$ linear coefficients determine the weight of every model within the fit. $k_{\lambda}$ is given by the dust-attenuation law (either [Reddy et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...828..107R/abstract) or SMC, [Prevot et al. 1984](https://ui.adsabs.harvard.edu/abs/1984A%26A...132..389P/abstract)), and $E_{B-V}$ is the so-called dust-attenuation parameter (in magnitudes). 

Finally, the best fit is chosen via a non-linear $\chi^2$ minimization algorithm with respect to the observed data (`lmfit` package[^4], see [Newville et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ascl.soft06014N/abstract)), and the errors are derived in a Monte-Carlo (MC) way, varying the observed pixel fluxes by a Gaussian distribution whose mean is zero and standard deviation is the $1 \sigma$ error of the flux at the same pixel, and re-fitting the continuum over a certain number of iterations.

$$ --- $$

- The code is composed of two `.py` files:
  - ```ficus.py``` is the main script. It reads the INPUT file provided by the user, and performs the fit (see [Running FiCUS](https://github.com/asalda/FiCUS/README.md#running-ficus)) according to the options enclosed in the CONFIGURATION file (see [Input and Configuration files](https://github.com/asalda/FiCUS/README.md#the-input-and-configuration-files)). Apart from the best-fit parameters, it creates the OUTPUT files and figures (see [Outputs and Plots](https://github.com/asalda/FiCUS/README.md#outputs0-and-plots)). 
  - ```ficus_functions.py``` is a secondary script. After being called, all the functions are imported into the main script. This file includes pre-defined scripts for spectral analysis, loading INPUT files and handling wityh data and modeld, as well as functions for the fitting routine, SED parameters calculations and plotting. 


## Installation
`FiCUS` is written in `Python` language. Before installing the code, the current working environment must be equipped with the following basic packages, that otherwise can be easily installed/updated using [pip](https://pypi.org/project/pip/), [conda](https://docs.conda.io/en/latest/) or any other package manager:
```
> python 3.7.4 
> matplotlib 3.1.1 
> astropy 4.0 
> lmfit 1.0.0 
```

... or later versions. Once the previous dependencies are fulfilled, `FiCUS` can be cloned from this repository using [git](https://git-scm.com/), by just plugging into the terminal the following command:
```
> git clone https://github.com/asalda/FiCUS/ficus.git
```


## Input and Configuration files
- The INPUT `.fits` file must contain, at least, four columns in the first extension after an empty Primary() HDU, with the following names: 
  - 'WAVE'.......... observed-frame wavelength array, in \AA, 
  - 'FLUX'........... spectral flux-density, in $F_{\lambda}$ units (erg/s/cm2/AA), 
  - 'FLUX_ERR'... $1 \sigma$ error on the spectral flux-density, 
  - 'MASK'.......... mask array (0 = masked, 1 = un-masked).
  
  The INPUT file can, for example, inherit the name (`SPEC-NAME`) of the spectrum to be fitted, and must always be placed into the [/FiCUS/inputs/](https://github.com/asalda/FiCUS/inputs/) folder beforehand. The 'MASK' extension of the INPUT file is an binary-array of the same length as 'WAVE', and indicates whether the 'FLUX' and 'FLUX_ERR' values at a certain wavelength $\lambda_i$ will be considered (1 = un-masked) or excluded (0 = masked) from the fit. 

- The CONFIGURATION `.ini` file contains all the input parameters that feed the ```ficus.py``` code. Its format is as follows, with the options (and descrption) for each input parameter spcified in brackets: 
  ```
  [FiCUS]
   plot_mode  > activate or desactivate PLOTTING mode 
                [yes // no], 
   ssp_models > pick your preferred stellar LIBRARY 
                [sb99 (Starburst99) // bpass (BPASSv2.2.1)], 
   zarray     > specify a set of METALLICITIES 
                [001,004,008,02,04 (standing for 1/20, 1/5, 2/5, 1 and 2 Z_\sun)], 
   att_law    > choose the DUST attenuation law 
                [r16 (Reddy et al. 2016) // smc (Prevot et al. 1994)], 
   wave_range > rest-frame WAVELENGTH range to be considered in the fit 
                [e.g., 1200.,1920. (\lambda_min, \lambda_max; in \AA)], 
   wave_norm  > rest-frame wavelength interval for NORMALIZATION of the spectrum
                [e.g., 1350.,1370. (\lambda_min, \lambda_max; in \AA)], 
   r_obs      > instrumental RESOLUTION of the input spectra 
                [e.g., 600.; as R = (\Delta \lambda) / \lambda], 
   nsim       > number of Monte-Carlo (MC) ITERATIONS 
                [e.g., 100.]. 
  ```

Examples of the INPUT and CONFIGURATION files can be found at the [/FiCUS/examples/](https://github.com/asalda/FiCUS/examples/) dedicated folder: [example.fits](https://github.com/asalda/FiCUS/examples/example.fits) and [example.ini](https://github.com/asalda/FiCUS/examples/example.ini), respectively. 


## Running FiCUS
Given the name of the INPUT file (`SPEC-NAME`) and the redshift of the source (`REDSHIFT`), the code can be run in console as a normal `.py` script:
```
> python3.7 ficus-path/FiCUS/ficus.py SPEC-NAME REDSHIFT
```

... where `ficus-path` corresponds to the local path in which the code was initially placed (see [Installation](https://github.com/asalda/FiCUS/README.md#install)). The code can also work within a jupyter-notebook environment (`.ipynb`) using the magic command `%run`:
```
> import os
> os.chdir(path-to-ficus/FiCUS/);
> %run -i ficus.py
```

When the code runs successfully, the terminal will print the name of the INPUT file, as well as the input parameters selected in the CONFIGURATION file. Succesfull fits will prompt the quote `# done!` when the process is finished. A typical console interface can be like this:
```
 ############################################################## 
 
 ###   FiCUS: Fitting the stellar Continuum of Uv Spectra   ### 
 
 ############################################################## 
   
 ### Running FiCUS (ficus.py) for example.fits ...
 
 ### inputs ### 
 # spec_name   --> example
 # z_spec      --> 3.6052
 # plot_mode   --> yes
 # ssp_models  --> sb99
 # Zarray      --> ['001', '004', '008', '02']
 # att_law     --> r16
 # wave_range  --> [1200. 1920.]
 # wave_norm   --> [1350. 1370.]
 # r_obs       --> 600.0
 # nsim        --> 100
 
 ### plotting...   
   
 # done!
```
... in which the [example.fits](https://github.com/asalda/FiCUS/examples/example.fits) spectrum at z = 3.6052 (VANDELS-ID: `CDFS017345`) was fitted using the `Starburst99` stellar library and a set of x4 metallicities. The dust attenuates the stellar continuum following [Reddy et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...828..107R/abstract) prescription. The wavelength range considered in the fit is $1200-1920$ angstroms, and the output SEDs are normalized to $1350-1370$ angstroms. The VANDELS resolution is $R = 600$ (example taken from [Saldana-Lopez et al. 2022b](https://ui.adsabs.harvard.edu/abs/2022arXiv221101351S/abstract)).

## Outputs and Plots
If the fit goes well, `FiCUS` generates different OUTPUT files

If `plot_mode == yes`, a `# plotting!` message will appear as soon as the new figure is created. 


[^1]: https://www.stsci.edu/science/starburst99/docs/default.htm
[^2]: https://bpass.auckland.ac.nz/
[^3]: https://trac.nublado.org/
[^4]: https://lmfit.github.io/lmfit-py/
