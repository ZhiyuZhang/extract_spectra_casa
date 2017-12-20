# extract_spectra_casa
extract spectra and plot/explore redshifts within CASA

extract.py is for extracting spectral files from datacubes, using CASA. 

The output are ascii files.


You also need to convert the ascii file into class format. 

readcasa.class is for the conversion. 


zd.class and zd.hlp are for CLASS/Gildas: 

https://www.iram.fr/IRAMFR/GILDAS/

To install, you can use macport which is the easiest way for Mac. 
Note that ! is the comment symbol in class. 

To run it, simply type: 

> class ( in your terminal) 

Then 

> zd redshift velocity_resolution 
example: 

>  zd 2.34  100  ! (km/s) 


updated 20 Dec. 2017 - zy

