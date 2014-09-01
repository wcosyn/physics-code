import numpy as np
import pylab as pl
import sys

x,y,r,i = np.loadtxt(sys.argv[1],unpack=True)

fig = pl.figure()
ax  = fig.add_subplot(111)
CS  = ax.scatter(x,y,c=r)
pl.colorbar(CS)

pl.show()
