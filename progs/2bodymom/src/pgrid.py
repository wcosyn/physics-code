import numpy as np
import pylab as pl
import sys
import matplotlib
#################################
# set a bunch of plot parameters#
#################################

params = {'legend.fontsize' : 35,
          'legend.linewidth': 2,
          'axes.linewidth'  : 3.5,
          'axes.labelsize'  : 30,
          'xtick.major.pad' : 10,
          'xtick.major.width' : 3,
          'xtick.major.size'  : 6,
          'ytick.major.width' : 3,
          'ytick.major.size'  : 6,
          'xtick.labelsize'    : 28,
          'ytick.labelsize'    : 28,
          'text.usetex'        : True,
	  'font.size'          : 30,
	  'font'               : 'serif',
          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts}"]
pl.rcParams.update(params)
pl.rc("font",family="serif")


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
blength = np.max(brange)-np.min(brange)
zlength = np.max(zrange)-np.min(zrange)
A = blength*zlength
print "min/max of data_r  : ", (np.min(data_r),np.max(data_r))
print "min/max of data_i  : ", (np.min(data_i),np.max(data_i))
print "integral of data_r : ", np.mean(data_r)*A
print "integral of data_i : ", np.mean(data_i)*A

#print data_r

fig = pl.figure()
fig.suptitle("CX scattering probability")
fig.subplots_adjust(bottom=0.17,left=0.15)
ax  = fig.add_subplot(111)
ax.contour(brange,zrange,np.transpose(data_r),100,lw=0.1)
CS  = ax.contourf(brange,zrange,np.transpose(data_r),100)#,antialiased=False)
pl.colorbar(CS)
ax.set_xlabel(r"$\mathbf{b}$ [fm]")
ax.set_ylabel(r"$\mathbf{z}$ [fm]")

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
