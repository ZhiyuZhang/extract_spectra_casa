"""
Resample spectra  
Old Array
New Array 
Diff Array 

old array: 2-D array
frequency, flux 

new array: 2-D array
frequency, flux

Diff array: 2-D array
frequency, flux 

"""

import matplotlib.pyplot as plt
import numpy as np 
import sys

spectrum = np.loadtxt("data.txt")

old_x  = spectrum[:,0]
old_y  = spectrum[:,1]

new_x  = np.arange(6510., 9085, 5.) + 2.5
resampled  = np.empty(len(new_x))

# For each bin in the diff array, calculate the centers, edges and fluxes   

## initialisation for input and output arrays 
#old_x      = np.arange(0, 125)
#old_y      = np.random.rand(len(old_x))
#new_x      = np.arange(-3, 127, 3)
#resampled  = np.empty(len(new_x))

# one_width and new_width are the bin width arrays, with the same size as the two arrays. 
# currently they are regularly gridding, but one day they can be non-uniform size for each small channel (bin). 
old_width  = np.ones(len(old_x))* (old_x[1] - old_x[0])
new_width  = np.ones(len(new_x))* (new_x[1] - new_x[0])

# old_x_edge and new_x_edge are the edge values for each channel (bin).  

old_x_edge = old_x - old_width / 2
old_x_edge = np.append(old_x_edge,    old_x[-1] + old_width[-1] / 2 )
new_x_edge = new_x - new_width / 2
new_x_edge = np.append(new_x_edge,     new_x[-1] + new_width[-1] / 2 )


#def diff_array(old_x, old_width, old_y, old_err, new_x, new_start, new_width, new_err):

for i in range(len(new_x)):
    index = np.where( (old_x - (old_width[0] / 2) < new_width[0] / 2 + new_x[i]   ) & (new_x[i] - new_width[0] / 2 < old_x +  (old_width[0]) / 2))

    if index[0].size > 1:

        sub_array_width     = old_width[index]
        sub_array_y         = old_y[index]
        sub_weight          = np.ones(len(index[0]))
        start               = np.abs(old_x_edge[index][0]  - new_x_edge[i]  )
        end                 = np.abs(old_x_edge[index][-1] - new_x_edge[i+1]) 


        if start > 0:
            sub_weight[0]       = start / old_width[0]
            sub_array_width[0]  = start
        if end > 0:
            sub_weight[-1]      = end / old_width[0]
            sub_array_width[-1] = end 
        resampled[i]        = np.sum(sub_array_y * sub_weight) / np.sum(sub_array_width)

        print( " ~~~~~ ") 
        print('sub_array_y:    ', sub_array_y)
        print('sub_weight:     ', sub_weight)
        print('sub_array_width ', sub_array_width)
        print( " ~~~~~ ") 

    elif index[0].size == 1:  
        resampled[i]        = old_y[index]
    else:
        resampled[i]        = np.nan


plt.clf()
f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
l                = ax1.step(old_x, old_y,where='mid', color="blue",     lw=0.3)
l                = ax1.step(new_x, resampled, where='mid', color="red", lw=0.3)

ax1.legend(loc=7, borderaxespad=0.)
l                = ax2.step(new_x, resampled, where='mid', color="blue", lw=0.5)
ax2.legend(loc=7, borderaxespad=0.)

ax1.set_xlim(6800, 9000)
plt.show()
#f.savefig('test.pdf', bbox_inches='tight', pad_inches=0.1)




#     bin_widths = np.zeros(wavelengths.shape[0])
# 
#     # This option makes the final entry in the left hand sides array the right hand side of the final bin
#     if make_rhs == "True":
#         bin_lhs = np.zeros(wavelengths.shape[0]+1)
#         #The first lhs position is assumed to be as far from the first central wavelength as the rhs of the first bin.
#         bin_lhs[0] = wavelengths[0] - (wavelengths[1]-wavelengths[0])/2
#         bin_widths[-1] = (wavelengths[-1] - wavelengths[-2])
#         bin_lhs[-1] = wavelengths[-1] + (wavelengths[-1]-wavelengths[-2])/2
#         bin_lhs[1:-1] = (wavelengths[1:] + wavelengths[:-1])/2
#         bin_widths[:-1] = bin_lhs[1:-1]-bin_lhs[:-2]
# 
#     # Otherwise just return the lhs positions of each bin
#     else:
#         bin_lhs = np.zeros(wavelengths.shape[0])
#         bin_lhs[0] = wavelengths[0] - (wavelengths[1]-wavelengths[0])/2
#         bin_widths[-1] = (wavelengths[-1] - wavelengths[-2])
#         bin_lhs[1:] = (wavelengths[1:] + wavelengths[:-1])/2
#         bin_widths[:-1] = bin_lhs[1:]-bin_lhs[:-1]
# 
#     return bin_lhs, bin_widths
# 
# 
# 
# # Function for performing spectral resampling on a spectrum or array of spectra.
# def spectres(spec_wavs, spec_fluxes, resampling, spec_errs=None):
# 
#     # Generate arrays of left hand side positions and widths for the old and new bins
#     filter_lhs, filter_widths = make_bins(resampling, make_rhs="True")
#     spec_lhs, spec_widths = make_bins(spec_wavs)
# 
# 
#     # Check that the range of wavelengths to be resampled onto falls within the initial sampling region
#     if filter_lhs[0] < spec_lhs[0] or filter_lhs[-1] > spec_lhs[-1]:
#         print "Spec_lhs, filter_lhs, filter_rhs, spec_rhs ", spec_lhs[0], filter_lhs[0], filter_lhs[-1], spec_lhs[-1]
#         sys.exit("spectres was passed a spectrum which did not cover the full wavelength range of the specified filter curve.")
#     
# 
#     #Generate output arrays to be populated
#     if spec_fluxes.ndim == 1:
#         resampled = np.zeros((resampling.shape[0]))
# 
#     elif spec_fluxes.ndim == 2:
#         resampled = np.zeros((len(resampling), spec_fluxes.shape[1]))
# 
#     if spec_errs is not None:
#         if spec_errs.shape != spec_fluxes.shape:
#             sys.exit("If specified, spec_errs must be the same shape as spec_fluxes.")
#         else:
#             resampled_errs = np.copy(resampled)
# 
#     start = 0
#     stop = 0
# 
#     # Calculate the new spectral flux and uncertainty values, loop over the new bins
#     for j in range(len(filter_lhs)-1):
# 
#         # Find the first old bin which is partially covered by the new bin
#         while spec_lhs[start+1] <= filter_lhs[j]:
#             start += 1
# 
#         # Find the last old bin which is partially covered by the new bin
#         while spec_lhs[stop+1] < filter_lhs[j+1]:
#             stop += 1
# 
#         if spec_fluxes.ndim == 1:
# 
#             # If the new bin falls entirely within one old bin the are the same the new flux and new error are the same as for that bin
#             if stop == start:
# 
#                 resampled[j] = spec_fluxes[start]
#                 if spec_errs is not None:
#                     resampled_errs[j] = spec_errs[start]
# 
#             # Otherwise multiply the first and last old bin widths by P_ij, all the ones in between have P_ij = 1 
#             else:
# 
#                 start_factor = (spec_lhs[start+1] - filter_lhs[j])/(spec_lhs[start+1] - spec_lhs[start])
#                 end_factor = (filter_lhs[j+1] - spec_lhs[stop])/(spec_lhs[stop+1] - spec_lhs[stop])
# 
#                 spec_widths[start] *= start_factor
#                 spec_widths[stop] *= end_factor
# 
#                 # Populate the resampled spectrum and uncertainty arrays
#                 resampled[j] = np.sum(spec_widths[start:stop+1]*spec_fluxes[start:stop+1])/np.sum(spec_widths[start:stop+1])
# 
#                 if spec_errs is not None:
#                     resampled_errs[j] = np.sqrt(np.sum((spec_widths[start:stop+1]*spec_errs[start:stop+1])**2))/np.sum(spec_widths[start:stop+1])
#                 
#                 # Put back the old bin widths to their initial values for later use
#                 spec_widths[start] /= start_factor
#                 spec_widths[stop] /= end_factor
# 
# 
#         # The same as above, except operates on each row of the array, resampling all of the input models
#         elif spec_fluxes.ndim == 2:
# 
#             if stop == start:
# 
#                 resampled[j, :] = spec_fluxes[start, :]
#                 if spec_errs is not None:
#                     resampled_errs[j, :] = spec_errs[start, :]
# 
#             else:
# 
#                 start_factor = (spec_lhs[start+1] - filter_lhs[j])/(spec_lhs[start+1] - spec_lhs[start])
#                 end_factor = (filter_lhs[j+1] - spec_lhs[stop])/(spec_lhs[stop+1] - spec_lhs[stop])
# 
#                 spec_widths[start] *= start_factor
#                 spec_widths[stop] *= end_factor
# 
#                 resampled[j, :] = np.sum(np.expand_dims(spec_widths[start:stop+1], axis=1)*spec_fluxes[start:stop+1, :], axis=0)/np.sum(spec_widths[start:stop+1])
#                 
#                 if spec_errs is not None:
#                     resampled_errs[j, :] = np.sqrt(np.sum((np.expand_dims(spec_widths[start:stop+1], axis=1)*spec_errs[start:stop+1])**2, axis=0))/np.sum(spec_widths[start:stop+1])
#                 
#                 spec_widths[start] /= start_factor
#                 spec_widths[stop] /= end_factor
# 
# 
#     # If errors were supplied return the resampled spectrum and error arrays
#     if spec_errs is not None:
#         return resampled, resampled_errs
# 
#     # Otherwise just return the resampled spectrum array
#     else: 
#         return resampled
# 
# 
