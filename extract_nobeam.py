## extract spectra from a datacube using a small aperture 
# use: change the file name in the filelist,  the location of the aperture, and the diameter of the aperture 
# k is the number of the desired datacube in the filelist. 

#-----------------------------------------------------------
# Setup: 
#
# 1). install astroph in casa
#  https://docs.astropy.org/en/stable/install.html 
#  

# 2). install the package supports 
#  
#  from setuptools.command import easy_install
#  subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'astropy'])
#  subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'photutils']) 
#  /Applications/CASA.app/Contents/Frameworks/Python.framework/Versions/3.6/bin/pip3 install astropy 
#  sys.path.append("/Users/zhiyu/Library/Python/3.6/lib/python/site-packages")

#-----------------------------------------------------------

# usage: 

# change the filenames in filelist and, corresponding locations, and select your circular Aperture diameter, and k (which file )
# then
# execfile('extract.py') 
#
# --written by Zhiyu Zhang  pmozhang@gmail.com 
# last update: 26 Apr. 2017  

import numpy as np
import copy
import astropy.units as u
import matplotlib.pyplot as plt
from   astropy.io          import fits
from   astropy             import wcs
from   photutils           import aperture_photometry, CircularAperture
from   matplotlib.patches  import Ellipse

os.system("rm -rf *fit_beam*")

# -----------   list of the files     ---------- 
filelist=[ 's18bt_lsb.cube.image',   # 0 
           's18bt_usb.cube.image',   # 1 
          ]

location=['12h10m32.557s, 39d24m20.989s', # 0
          '12h10m32.557s, 39d24m20.989s', # 1
           ] 

AperDiameter = 3 # arcsec 
k            = 1 # the number in the filelist and associated location 


#----------------------- do not change below -------------------


# ------read header file --------
myhead    =  imhead(filelist[k],mode='list')
channels  =  myhead['perplanebeams']['nChannels'] #  #channels  =  SpecExtrCube.shape[2] The same 

# ---------- corresponding psf files  ----------
psffilelist  = copy.copy(filelist)
beamfilelist = copy.copy(filelist)


for i in range(0, len(filelist)):
     psffilelist[k]  = filelist[k].replace("image", "psf")
     beamfilelist[k] = filelist[k].replace("image", "fit_beam")
# ----------------------------------------------


# --------- the regions for spectral subtraction ------- 
#AperDiameter      =  3   # diameter of the aperture for the spectral extraction (in arcsec)
region            =  'circle[['+location[k]+'], '+str(AperDiameter/2.)+'arcsec]'
region_beam_fit   =  'circle[['+str(myhead['shape'][0]/2.)+'pix,' +str(myhead['shape'][1]/2.)+'pix],' +str(AperDiameter*2.)+'arcsec]'

# region is the aperture region for spectral extraction 
# region_beam_fit is the region for fitting the beam shape for the .psf files 


# Make circular apertures with diameter of xxx arcsec, centralised in (R.A. Dec. J2000)

# --------- extract spectra from the defined region ----
xval         =  imval(imagename=filelist[k], region=region)
SpecExtrCube =  xval['data']
SpecExtrMask =  xval['mask']


# --------------------------------
# get the extracted spectral cube from the selected region  
# This is a 3-D cube, with a dimension something like (15x15x127) 
# So, the spectral values are in box, rather than in circles or in Ellipses
# Xval goes with a mask array to turn the box into circles or Ellipses 
# Xval['data'] is the cube
# Xval['mask'] is the mask cube 
# -------------------------------

#--------------------------------------------
# apply mask to the extracted spectra 
SpecExtrCube[np.where(SpecExtrMask==False)] = np.NaN
# from now on, SpecExtrCube is the cube of the mask-applied extracted spectra 
#--------------------------------------------

# -------------------------------------------
# Obtain the average flux density within the masked region, using the mask-applied extracted spectra 
Flux_jy_p_mean = np.nanmean(SpecExtrCube,axis=(0,1))
# -------------------------------------------


# ------read header file --------
myhead    =  imhead(filelist[k],mode='list')
channels  =  myhead['perplanebeams']['nChannels'] #  #channels  =  SpecExtrCube.shape[2] The same 



# -initialise bmaj and bmin array for each channel 
bmaj      =  np.arange(channels)*1.0 
bmin      =  np.arange(channels)*1.0
bpa       =  np.arange(channels)*1.0


for i in range(0, channels):
#   print(i)
    bmaj[i] = myhead['perplanebeams']['*'+str(i)]['major']['value']
    bmin[i] = myhead['perplanebeams']['*'+str(i)]['minor']['value']
    bpa[i]  = myhead['perplanebeams']['*'+str(i)]['positionangle']['value']




f0       = float(myhead['crval4'])          # reference freq in Hz
df       = float(myhead['cdelt4'])          # channel width in Hz
i0       =       myhead['crpix4']           # reference pixel
freqspec = ((np.arange(channels))*df + f0)  # channel array of the observing frequencies in Hz 
xdat     =  freqspec /1E9                   # Convert from Hz to GHz



# --- Read and assign beam values to the arrays 
#     Obtain beamsize/shape of each channel varies.

#get pixel size in arcsec  
pixelsize    = np.abs((myhead['cdelt1']* u.rad).to(u.arcsec).value)
# beam area calculated from a 2-D Gaussian  in arcsec^2  
MyBeamArea     = np.pi*  bmaj * bmin / ( 4 * np.log(2) ) 
# How many pixels in one beam 
MyBeamPixels   = MyBeamArea / (pixelsize**2 )
#------------------

# Area of the aperture for line extraction in arcsec^2 (Area = pi * r^2 ) 
MyAperArea     =  np.pi*(AperDiameter/2.)**2
# How many pixels in one aperture 
MyAperPixels   = MyAperArea/ pixelsize**2 
#  ---------------------------------

# ----- Flux density in the masked aperture -----
Flux           = Flux_jy_p_mean * MyAperPixels / MyBeamPixels 
# -----------------------------------------------

plt.clf()
ax1        =  plt.subplot(111) 
ax1.plot(xdat,Flux,drawstyle='steps-mid')

plt.show()
plt.savefig('test.pdf')

# --------------------- fitting 
#
#
imfit(imagename=psffilelist[k], region=region_beam_fit,model=beamfilelist[k])

exportfits(imagename=beamfilelist[k],fitsimage=beamfilelist[k]+'.fits')


BeamHeader    = imstat(beamfilelist[k]) 
image         =  fits.open(beamfilelist[k]+'.fits')[0]
header        =  image.header
image         =  image.data



target_location_pix  = (BeamHeader['maxpos'][0], BeamHeader['maxpos'][1])


f, (ax1, ax2)  = plt.subplots(2, sharex=False, sharey=False)

Aper_Corr = image[0,:,0,0]


for channel in range(0, BeamHeader['trc'][3]+1):
    beam       = image[0,channel,:]
    radii      = np.arange(1,12,0.5)
    growth     = []
    fluxes     = []

    for radius in radii:
        radius_pix           = radius / abs(header['CDELT1'] * 3600)               ## in pixel
        apertures            = CircularAperture(target_location_pix, r=radius_pix)
        phot_table           = aperture_photometry(beam, apertures)
        flux_aper            = phot_table['aperture_sum'].data
#       print(flux_aper)
        fluxes.append(flux_aper)

    mean   = np.mean(fluxes[-5:-1])
    growth = fluxes / mean
    bad    = np.isnan(growth)
    growth[bad] = 1
    Aper_Corr[channel] = np.interp(AperDiameter/2., radii, growth[:,0])

    ax1.plot(radii,growth)


Flux_corrected= Flux/Aper_Corr



ax1.set_xlabel('Aperture radii [\"]')
ax1.set_ylabel('Power fraction')
ax1.set_title(r'Curve of growth for the clean beams (2-D Gaussian fitted)')



ax2.set_aspect('equal')
ax2.set_xlim(30,100)
ax2.set_ylim(3,18)

j=channels*0.1 
ells = Ellipse(xy= [40,12], width=bmaj[j], height=bmin[j], angle=bpa[j])
ax2.add_artist(ells)

j=channels* 0.25
ells = Ellipse(xy= [50,12], width=bmaj[j], height=bmin[j], angle=bpa[j])
ax2.add_artist(ells)

j=channels * 0.45
ells = Ellipse(xy= [60,12], width=bmaj[j], height=bmin[j], angle=bpa[j])
ax2.add_artist(ells)

j=channels* 0.6
ells = Ellipse(xy= [70,12], width=bmaj[j], height=bmin[j], angle=bpa[j])
ax2.add_artist(ells)

j=channels* 0.75 
ells = Ellipse(xy= [80,12], width=bmaj[j], height=bmin[j], angle=bpa[j])
ax2.add_artist(ells)

j=channels* 0.9 
ells = Ellipse(xy= [90,12], width=bmaj[j], height=bmin[j], angle=bpa[j])
ax2.add_artist(ells)


ax2.set_xlabel(r'$\Delta$ R.A  ["]')
ax2.set_ylabel(r'$\Delta$ Dec. ["]')

plt.savefig(filelist[k]+'aper_correction.pdf')

plt.clf()
ax1        =  plt.subplot(111) 
ax1.plot(xdat,np.array(Flux_corrected),drawstyle='steps-mid')

plt.show()
plt.savefig(filelist[k]+'corrected.pdf')

import sys
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
header="this is a test"
data = Table([xdat,Flux_corrected], names=['#freq(GHz)', 'flux(Jy)'] )
ascii.write(data, filelist[k]+'.dat', format='csv', comment=header)



