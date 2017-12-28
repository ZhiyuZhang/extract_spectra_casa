"""
This script creates the plot `line.pdf`, a stepped spectral line plot with GHz
along the lower x-axis and km/s along the upper x-axis
"""


from astropy.io import fits
from astropy.io import ascii
from astropy.cosmology import WMAP9 as cosmo

import copy
import re
import sys

from resample import *


import matplotlib
import matplotlib.gridspec
import matplotlib.ticker
import numpy as np

import astropy.units as u
import matplotlib.pyplot as plt



from astropy.table import Table, Column, MaskedColumn


import numpy as n
import scipy.interpolate
import scipy.ndimage



data = np.loadtxt('test_spectrum.txt')

frequency = data[:,0] * u.GHz 
flux      = data[:,1] * u.Jy 

restfreq      = frequency[0]         
freq_to_vel   = u.doppler_relativistic(restfreq)
velocity      = frequency.to(u.km / u.s, equivalencies=freq_to_vel)

new_resolution = 250
new_velocity   = np.arange(np.min(velocity.value), np.max(velocity.value), new_resolution)

new_flux     = resample(velocity.value, flux.value, new_velocity) 


plt.clf()
f, ax1 = plt.subplots(1, sharex=True, sharey=False)
l = plt.step(velocity, flux, 'b-',         lw=0.2,  where='mid', label='ori') 
l = plt.step(new_velocity, new_flux, 'r-', lw=0.2,  where='mid', label='new') 

ax1.legend(loc=2, borderaxespad=0.)


f.savefig('test.pdf', bbox_inches='tight', pad_inches=0.1)




 
