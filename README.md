# FiCUS (FItting the stellar Continuum of Uv Spectra)
- [Description](https://github.com/asalda/FiCUS/main/README.md#description)
- [Installation](https://github.com/asalda/FiCUS/main/README.md#install)
- [Input and Configuration files](https://github.com/asalda/FiCUS/main/README.md#the-input-and-configuration-files)
- [Running FiCUS](https://github.com/asalda/FiCUS/main/README.md#running-ficus)
- [Outputs and Plots](https://github.com/asalda/FiCUS/main/README.md#outputs0-and-plots)
- [Tutorials](https://github.com/asalda/FiCUS/main/README.md#tutorials)
- [Links and Notes](https://github.com/asalda/FiCUS/main/README.md#links-and-notes)


## Description
`FiCUS` is a customized `Python` script to fit the stellar continuum of extragalactic ultraviolet (UV) spectra. In short, it takes observed-frame wavelength, flux density (with errors) and user-defined mask arrays, and returns an estimation of the galaxy light-weighted stellar age, metallicity and dust extinction, as well as other secondary Spectral Energy Distribution (SED) parameters. 

The code was presented in [Saldana-Lopez et al. 2022b](https://ui.adsabs.harvard.edu/abs/2022arXiv221101351S/abstract), but the methodology was first described and tested in [Chisholm et al. 2019](https://ui.adsabs.harvard.edu/abs/2022arXiv221101351S/abstract). A previous version of the code has been extensively used in other papers such as [Gazagnes et al. 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..29G/abstract), [Gazagnes et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...639A..85G/abstract) and [Saldana-Lopez et al. 2022a](https://ui.adsabs.harvard.edu/abs/2022A%26A...663A..59S/abstract). 

The UV stellar continuum modeling $F_{\lambda}^{\star}$ is achieved by fitting observed spectra with a linear combination of single-burst stellar population (SSP) theoretical models: STARBURST99 ([Leitherer et al. 2010](https://ui.adsabs.harvard.edu/abs/2010ApJS..189..309L/abstract)) or BPASS ([Eldridge et al. 2017](https://ui.adsabs.harvard.edu/abs/2017PASA...34...58E/abstract)). These models assume a initial mass function (IMF) with a high-(low-)mass exponent of 2.3 (1.3), and a high-mass cutoff at 100 $M_{\odot}$. The models include five different metallicities (0.05, 0.2, 0.4, 1 and 2 $Z_{\odot}$) and ten ages for each metallicity (1, 2, 3, 4, 5, 8, 10, 15, 20 and 40 Myr). A nebular continuum was added to every model by self-consistently processing the original SSPs through the `Cloudy v17.0` code[^1] [(Ferland et al. 2017)](https://ui.adsabs.harvard.edu/abs/2017RMxAA..53..385F/abstract), assuming similar gas-phase and stellar metallicities, an ionization parameter of $\log(U)=-2.5$, and a volume hydrogen density of $n_H = 100 cm^{-3}$. Adopting a simple geometry where _all_ the light is attenuated by a uniform foreground slab of gas surrounding the galaxy, this results in: 

$$ F_{\lambda}^{\star} = 10^{-0.4 k_{\lambda} E_{B-V}} \sum_{i,j} X_{ij} F_{\lambda}^{ij} $$

$$ i \equiv 1, 2, 3, 4, 5, 8, 10, 15, 20, 40 Myr $$

$$ i \equiv 0.05, 0.2, 0.4, 1, 2 Z_{\odot} $$

where $F_{\lambda}^{ij}$ represents the corresponding model at the i-_th_ age and j-_th_ metallicity, and the $X_{ij}$ linear coefficients determine the weight of every model within the fit. $k_{\lambda}$ is given by the dust-attenuation law (either [Reddy et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...828..107R/abstract) or SMC, [Prevot et al. 1984](https://ui.adsabs.harvard.edu/abs/1984A%26A...132..389P/abstract)), and $E_{B-V}$ is the so-called dust-attenuation parameter (in magnitudes). 

Finally, the best fit is chosen via a non-linear $\chi^2$ minimization algorithm with respect to the observed data (`lmfit` package[^2], see [Newville et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ascl.soft06014N/abstract)), and the errors are derived in a Monte-Carlo (MC) way, varying the observed pixel fluxes by a Gaussian distribution whose mean is zero and standard deviation is the 1 $\sigma$ error of the flux at the same pixel, and re-fitting the continuum over a certain number of iterations.

$$ --- $$

- The code is structured in two `.py` files:
  - ```ficus.py``` is the main script. It reads the INPUT file provided by the user, and performs the fit (see [Running FiCUS](https://github.com/asalda/FiCUS/edit/main/README.md#running-ficus)) according to the options enclosed in the CONFIGURATION file (see [Input and Configuration files](https://github.com/asalda/FiCUS/edit/main/README.md#the-input-and-configuration-files)). Apart from the best-fit parameters, it creates the OUTPUT files and figures (see [Outputs and Plots](https://github.com/asalda/FiCUS/edit/main/README.md#outputs0-and-plots)). 
  - ```ficus_functions.py``` is a secondary script. After being called, all the functions are imported into `ficus.py`. This file includes pre-defined scripts for spectral analysis, loading INPUT files and handling wityh data and modeld, as well as functions for the fitting routine, SED parameters calculations and plotting. 


## Installation
`FiCUS` is written in `Python` language. Before installing the code, the current working environment must be equipped with the following versions and basic packages, that otherwise can be easily installed/updated using [pip](https://pypi.org/project/pip/), [conda](https://docs.conda.io/en/latest/) or any other package manager:
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
- The INPUT file is a `.fits` file that must contain, at least, the following extensions and names: 
  - 'WAVE'.......... observed-frame wavelength array, in \AA, 
  - 'FLUX'........... spectral flux-density array, in F_\lambda units (e.g., erg/s/cm2/AA), 
  - 'FLUX_ERR'... 1 $\sigma$ error on the spectral flux-density, in F_\lambda units, 
  - 'MASK'.......... mask array (0 = masked, 1 = un-masked).
  
  This INPUT file can, for example, inherit the `NAME` of the spectrum to be fitted, and must always be placed into the `ficus-path/inputs/` folder. We note that 'WAVE', 'FLUX' and 'FLUX_ERR'. The binary-array 'MASK' match the length of the WAVE array, and 

- The CONFIGURATION file


Dedicated examples of the INPUT (`example.fits`) and CONFIGURATION files (`example.ini`) can be found at [./examples/](https://github.com/asalda/FiCUS/main/examples/) folder.


## Running FiCUS
Given the name of the INPUT file (`NAME`) and the redshift of the source (`REDSHIFT`), the code can be run in console as a normal `.py` script:
```
> python3.7 ficus-path/ficus.py NAME REDSHIFT
```

The code can also work within a jupyter-notebook environment (`.ipynb`) using the magical command `%run`:
```
> import os
> os.chdir(ficus-path);
> %run -i ficus.py
```


## Outputs and Plots



## Tutorials

`/!\  ... work in progress ...  /!\`


## Links and Notes
[1]: https://trac.nublado.org/



[2]: https://lmfit.github.io/lmfit-py/

#### ```ficus.py``` -> main script. 
 Fit the observed input spectra with the best combination of SSPs models, and return 
 the mean light-fractions (X_i) and attenuation parameter (E_BV, uniform screen of dust). 
 Additionally, a bunch of different SED-derived parameters are calculated. 
 
 The input spectra have to be placed in the ```/inputs/``` folder as a .fits file with 
 a similar format than in the ```/examples/example.fits``` provided example (i.e. four .fits columns 
 in the first extension after an empty Primary() HDU: 'WAVE, FLUX, FLUX_ERR, MASK'). 
 
 The masks are included in the same file as an additional binary-array (0 = masked, 1 = un-masked).
 Outputs are saved in .txt/.npy/.pdf format in the ```/outputs/``` folder.

 To be run in console as:
 ```
 > python3.7 your-ficus-path/ficus.py
 ```
 
 To be run into a jupyter-notebook (.ipynb) as:
 ```
 > import os
 > os.chdir(your-ficus-path);
 > %run -i your-ficus-path/ficus.py
 ```


#### ```ficus_functions.py``` -> secondary script. 
 After being called, all the analysis and plotting functions are imported into ```ficus.py```.


#### ```ficus.ini``` -> configuration file. 
 It contains all the input parameters that feed the ```ficus.py``` code.
 It is structured as in the ```/examples/example.ini``` default example:
 
 ```> less /examples/example.ini
   [ficus]
   spec_name = example
   plot_mode = yes
   ssp_models = sb99
   zarray = 001,004,008,02
   att_law = reddy16
   att_igm = no
   wave_range = 925.,1375.
   z_spec = 0.33314
   wave_r = 1000.
   wave_norm = 1070.,1110.
   nsim = 100.
   
   #EOF
 ```
 
 [^1]: https://trac.nublado.org/
 [^2]: https://lmfit.github.io/lmfit-py/
