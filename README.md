# FiCUS

## FItting the stellar Continuum of Uv Spectra


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
