import sys
import pylab as pl
import numpy as np


b,z,res,err = np.loadtxt(sys.argv[1],unpack=True)

fig = pl.figure()
ax = fig.add_subplot(111)
CS = ax.scatter(b,z,c=res,s=20,linewidths=0.)
pl.colorbar(CS)
ax.set_xlabel("b")
ax.set_ylabel("z")

pl.show()
