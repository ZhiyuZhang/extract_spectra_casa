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

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    print('ndims:  ', ndims)
    if len( newdims ) != ndims:
        print("[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions.")
        return None
    newdims = n.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print("Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported.")



data = np.loadtxt('test_spectrum.txt')

frequency = data[:,0] * u.GHz 
flux      = data[:,1] * u.Jy 

plt.clf()
f, ax1 = plt.subplots(1, sharex=True, sharey=False)
l = plt.step(frequency, flux, 'b-', linewidth=1, where='mid', label='mid') 
ax1.legend(loc=7, borderaxespad=0.)


f.savefig('test.pdf', bbox_inches='tight', pad_inches=0.1)




#     table = fits.open('stack_all_base_cut-scale-to-z-H2O.fits')[1].data
# 
#     # change unit from Jy to W*u.m ** -2 vs. Hz
#     table.flux  = (table.flux  * u.Jy).to(u.W * u.cm ** -2 / u.Hz)
# 
#     #--------------------------------------------------------
#     #    flux in  W / m** 2 /Hz
#     #    get luminosity density (spectral luminosity) in W/Hz
#     #--------------------------------------------------------
#     z                       = 2.32591
#     luminosity_distance     = cosmo.luminosity_distance(z)
#     luminosity_distance     = luminosity_distance.to(u.cm)
#     table.flux              = table.flux  * 4 * np.pi * luminosity_distance.value ** 2
#     table.flux              = (table.flux * u.W / u.Hz ).to(u.Lsun / u.Hz)
#     spectrum                = specu.Spectrum1D(disp=table.wave * u.GHz, flux=table.flux * u.Lsun / u.Hz)
# 
# 
# 
