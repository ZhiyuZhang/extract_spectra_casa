from resample import *
x= np.arange(100)
y= x**2
import numpy as np
import matplotlib.pylab as plt
plt.step(x,y)
#plt.show()


xx = np.arange(10)*10

yy =resample(x,y,xx)

plt.step(xx,yy)
plt.show()

