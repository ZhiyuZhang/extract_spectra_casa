"""
This script creates the plot `line.pdf`, a stepped spectral line plot with GHz
along the lower x-axis and km/s along the upper x-axis


USAGE: 
    zd(input_file, redshift, velo_resolution) 


Example: 
$ipython
In [1]: from zd import *
In [2]: zd('test_spectrum.txt',3,100)

# modified py NHH
In [2]: zd('path_to_data',3,100)


Written by Zhiyu Zhang
pmozhang@gmail.com
1st update 29 Dec. 2017 

"""

from astropy.io import ascii
from resample import *
from astropy.cosmology import WMAP9 as cosmo

import sys,os
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt





def zd(input_spectra_file, redshift, velo_resolution): 
#def zd(input_spectra_file, redshift, velo_resolution): 

    linedata       = ascii.read('line.dat')
    kall            = 17 # the number in the filelist and associated location 

    
    filelist=['spw0.image.dat',   # 0 
          'spw1.image.dat',   # 1
          'spw2.image.dat',   # 2
          'spw3.image.dat',   # 3
          '276_278GHz.image.dat',   # 4 
          '277_279GHz.image.dat',   # 5 
          '289_291GHz.image.dat',
          '290_292GHz.image.dat',
          '312_314GHz_high.image.dat',
          '313_315GHz_high.image.dat',
          '316_328_high.image.dat',
          '326_328GHz.image.dat',
          '327_329GHz.image.dat',
          '338_340GHz.image.dat',
          '339_341GHz.image.dat',
          'spw03.image.dat',
          'spw12.image.dat',
          'spw_3_15.625MHz_uvlin.image.dat'#17 not use for now
          ]


    for k in range(kall): 

        path_to_data =  input_spectra_file+filelist[k]
        '''
        data           = np.loadtxt(path_to_data)
        data           = np.loadtxt(input_spectra_file)    
        frequency      = floatdata[:,0] * u.GHz
        flux           = data[:,1] * u.Jy
        '''
        readfile = open(path_to_data, 'r')
        frequency = []
        flux = []
        i = 0
        for line in readfile:
            if i > 0:
                line=line.rstrip('\n')
                element=line.split(',')
                
                frequency += [float(element[0])] #* u.GHz
                flux += [float(element[1])]
            i +=1
        frequency = np.array(frequency)
        flux = np.array(flux)
        readfile.close()

        '''
        restfreq       = frequency[0]
        freq_to_vel    = u.doppler_relativistic(restfreq)
        velocity       = frequency.to(u.km / u.s, equivalencies = freq_to_vel)
    
        new_resolution = velo_resolution 
        new_velocity   = np.arange(np.min(velocity.value), np.max(velocity.value), new_resolution)
        new_flux       = resample(velocity.value, flux.value, new_velocity)
        '''
    
        plt.clf()
        f, ax1 = plt.subplots(1, sharex=True, sharey=False)
        #l = plt.step(velocity, flux, 'b-',         lw=0.2,  where='mid', label='original') 
        #l = plt.step(new_velocity, new_flux, 'r-', lw=0.5,  where='mid', label='resampled',alpha=0.5,) 
        l = plt.step(frequency, flux, 'r-', lw=0.5,  where='mid', label='resampled',alpha=0.5,) 
        ax1.legend(loc=2, borderaxespad=0.)
 
        '''
        vmin = np.min(velocity.value)
        vmax = np.max(velocity.value)    
        '''
        vmin = np.min(frequency)
        vmax = np.max(frequency)
        ax1.set_xlim(vmin,vmax)

        
        for num, name in enumerate(linedata['line']):
            line_freq = linedata[num][0] / (1 + redshift) /1e3 #* from u.MHz to u.GHz
            '''
            line_velo = line_freq.to(u.km / u.s, equivalencies = freq_to_vel).value
            if line_velo < vmax and line_velo > vmin:  
                ax1.axvline(line_velo, ymax= 1, color='b', linewidth=0.2)
                ax1.annotate(name, (line_velo, 0.005), rotation='vertical', horizontalalignment='center', verticalalignment='bottom', size=3,)
                ax1.set_xlabel('Velocity ['+str(velocity.unit)+']')
                ax1.set_ylabel('Flux ['+str(flux.unit)+']')
            '''
            if line_freq < vmax and line_freq > vmin:  
                ax1.axvline(line_freq, ymax= 1, color='b', linewidth=0.2)
                ax1.annotate(name, (line_freq, 0.005), rotation='vertical', horizontalalignment='center', verticalalignment='bottom', size=3,)
                ax1.set_xlabel('frequency [GHz]')
                ax1.set_ylabel('Flux [Jy]')
  
        
        
        f.savefig(path_to_data+'z%.2f.pdf' % redshift, bbox_inches='tight', pad_inches=0.1)
        #os.system('open test.pdf')
        
    return 







 
