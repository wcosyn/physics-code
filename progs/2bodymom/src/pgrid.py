import numpy as np
import pylab as pl
import sys
import matplotlib
f = file(sys.argv[1],"r")

i = 0
brange = None
zrange = None
data_r = None
error_r = None
data_i = None
error_i = None
for line in f:
    if not line.strip().startswith('#'):
        verts = line.split()
        if (brange==None):
            brange = np.linspace(*map(float,verts))
            print "brange is ", brange
            print "length is ", len(brange)
        elif (zrange==None):
            zrange = np.linspace(*map(float,verts))
            print "zrange is ", zrange
            print "length is ", len(zrange)
        else:
            if (data_r==None):
                data_r = np.zeros((len(brange),len(zrange)))
                data_i = np.zeros((len(brange),len(zrange)))
                error_r = np.zeros((len(brange),len(zrange)))
                error_i = np.zeros((len(brange),len(zrange)))
            print ("processing item nr {:d}, maps to indices {:d},{:d}".format(i,i/len(zrange),i%len(zrange)))
            data_r[i/len(zrange),i%len(zrange)] = float(verts[0]) 
            data_i[i/len(zrange),i%len(zrange)] = float(verts[1])
            error_r[i/len(zrange),i%len(zrange)] = float(verts[2]) 
            error_i[i/len(zrange),i%len(zrange)] = float(verts[3])
            i+=1


print "brange", brange
print "zrange", zrange
print data_r

fig = pl.figure()
fig.suptitle("real part")
ax  = fig.add_subplot(111)
CS  = ax.contourf(brange,zrange,np.transpose(data_r),100)
pl.colorbar(CS)

fig = pl.figure()
fig.suptitle("imag part")
ax  = fig.add_subplot(111)
CS  = ax.contourf(brange,zrange,np.transpose(data_i),100)
pl.colorbar(CS)

fig = pl.figure()
fig.suptitle("real error")
ax  = fig.add_subplot(111)
CS  = ax.contourf(brange,zrange,np.transpose(error_r),100)
pl.colorbar(CS)

fig = pl.figure()
fig.suptitle("imag error")
ax  = fig.add_subplot(111)
CS  = ax.contourf(brange,zrange,np.transpose(error_i),100)
pl.colorbar(CS)

pl.show()
