from   resample import *
import numpy as np
import matplotlib.pylab as plt

old_x = np.arange(1,100,1)
old_y = old_x**2

new_x  = np.arange(1,100,10)
new_y = resample(old_x,old_y,new_x)

#plt.step(xx,yy)


plt.clf()
f, ax1 = plt.subplots(1, sharex=True, sharey=True)
#l      = ax1.step(x,y)
l      = ax1.step(old_x, old_y,  where='mid', color="blue", lw=0.3)
l      = ax1.step(new_x, new_y, where='mid', color="red", lw=0.3)
plt.show()

