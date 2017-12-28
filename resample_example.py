"""
Resample spectra  
From Old Array to  New Array 

old array: 2-D array
frequency, flux 

new array: 2-D array
frequency, flux


"""

import matplotlib.pyplot as plt
import numpy as np 
import sys

#spectrum  = np.loadtxt("data.txt")
spectrum  = np.loadtxt("test_spectrum.txt")
old_x     = spectrum[:,0]
old_y     = spectrum[:,1]
#old_x     = old_x[::-1]
#old_y     = old_y[::-1]

new_x     = np.arange(80., 90, 0.1) 

## initialisation for input and output arrays 
#old_x      = np.arange(0, 125)
#old_y      = np.random.rand(len(old_x))
#new_x      = np.arange(-3, 127, 3)



def resample(old_x, old_y, new_x):
    if old_x[1] - old_x[0] < 0:
        old_x   = old_x[::-1]
        old_y   = old_y[::-1]
        reverse = True 
    else: 
        reverse = False

    # define output flux array, using the same size of the x axis
    new_y = np.empty(len(new_x))

    # one_width and new_width are the bin width arrays, with the same size as the two arrays. 
    # currently they are regularly gridding, but one day they can be non-uniform size for each small channel (bin). 
    old_width  = np.ones(len(old_x))* (old_x[1] - old_x[0])
    new_width  = np.ones(len(new_x))* (new_x[1] - new_x[0])
    
    # old_x_edge and new_x_edge are the edge values for each channel (bin).  
    
    old_x_edge = old_x - old_width / 2
    old_x_edge = np.append(old_x_edge,    old_x[-1] + old_width[-1] / 2 )
    new_x_edge = new_x - new_width / 2
    new_x_edge = np.append(new_x_edge,     new_x[-1] + new_width[-1] / 2 )
    
    
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
            new_y[i]        = np.sum(sub_array_y * sub_weight * old_width[0]) / np.sum(sub_array_width)
    
            print( " ~~~~~ ") 
            print('sub_array_y:    ', sub_array_y)
            print('sub_weight:     ', sub_weight)
            print('sub_array_width ', sub_array_width)
            print( " ~~~~~ ") 
    
        elif index[0].size == 1:  
            new_y[i]        = old_y[index]
        else:
            new_y[i]        = np.nan
    return new_y 
    
plt.clf()
f, ax1        = plt.subplots(1, sharex=True, sharey=False)
l             = ax1.step(old_x, old_y,                        where='mid', color="blue", lw=0.3)
l             = ax1.step(new_x, resample(old_x, old_y, new_x), where='mid', color="red", lw=0.3)
#ax1.set_xlim(6800, 9000)
ax1.legend(loc=7, borderaxespad=0.)
#plt.show()
f.savefig('test.pdf', bbox_inches='tight', pad_inches=0.1)
stop    
