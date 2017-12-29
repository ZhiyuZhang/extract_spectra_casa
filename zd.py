"""
This script creates the plot `line.pdf`, a stepped spectral line plot with GHz
along the lower x-axis and km/s along the upper x-axis


USAGE: 
    zd(input_file, redshift, velo_resolution) 


Example: 

In [1]: from zd import *
In [2]: zd('test_spectrum.txt',3,100)


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

    data           = np.loadtxt(input_spectra_file)
    linedata       = ascii.read('line.dat')
    
    frequency      = data[:,0] * u.GHz
    flux           = data[:,1] * u.Jy
    
    restfreq       = frequency[0]
    freq_to_vel    = u.doppler_relativistic(restfreq)
    velocity       = frequency.to(u.km / u.s, equivalencies = freq_to_vel)
    
    new_resolution = velo_resolution 
    new_velocity   = np.arange(np.min(velocity.value), np.max(velocity.value), new_resolution)
    new_flux       = resample(velocity.value, flux.value, new_velocity)
    
    
    plt.clf()
    f, ax1 = plt.subplots(1, sharex=True, sharey=False)
    l = plt.step(velocity, flux, 'b-',         lw=0.2,  where='mid', label='original') 
    l = plt.step(new_velocity, new_flux, 'r-', lw=0.5,  where='mid', label='resampled',alpha=0.5,) 
    ax1.legend(loc=2, borderaxespad=0.)
    
    vmin = np.min(velocity.value)
    vmax = np.max(velocity.value)
    
    ax1.set_xlim(np.min(velocity.value),np.max(velocity.value))
    
    for num, name in enumerate(linedata['line']):
        line_freq = linedata[num][0] / (1 + redshift) * u.MHz
        line_velo = line_freq.to(u.km / u.s, equivalencies = freq_to_vel).value
        if line_velo < vmax and line_velo > vmin:  
            ax1.axvline(line_velo, ymax= 1, color='b', linewidth=0.2)
            ax1.annotate(name, (line_velo, 0.005), rotation='vertical', horizontalalignment='center', verticalalignment='bottom', size=3,)
            ax1.set_xlabel('Velocity ['+str(velocity.unit)+']')
            ax1.set_ylabel('Flux ['+str(flux.unit)+']')
    
    f.savefig('test.pdf', bbox_inches='tight', pad_inches=0.1)
    os.system('open test.pdf')
    return 







 
