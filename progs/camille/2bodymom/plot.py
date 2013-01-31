#!/usr/bin/python

import numpy as np
import pylab as pl
import sys

f = None
if (len(sys.argv) > 1):
	print "opening file ", sys.argv[1]
	f = file(sys.argv[1])
else:
	#f = file("rho_test.dat")
	f = file("rho_symm.dat")

def getvars(verts): # assumes each vert is of type "a=5"
	vard = dict()
	for v in verts:
		a,b = v.split("=")
		vard[a]=int(b)
	return vard

def getLchar(li):
	l=['s','p','d','f']
	return l[li]

def makeTitle(var):
	return r""+str(var["n1"])+getLchar(var["l1"])+str(var["j1"])+"/2"+" - "+str(var["n2"])+getLchar(var["l2"])+str(var["j2"])+"/2"

def makePlot(X,var):
	X = np.array(X) # for 2D slicing
	xind = np.nonzero(X[:,0])[0]
	yind = np.nonzero(X[:,1])[0]
	zind = np.nonzero(X[:,2])[0]
	
	if len(xind > 0):
		fig = pl.figure(figsize=(8,6))
		ax = fig.add_subplot(111)
		fig.subplots_adjust(bottom=0.2,left=0.2)
		#ax.plot(X[0:len(X)/3,0],X[0:len(X)/3,3])
		ax.plot(X[xind,0],X[xind,3])
		ax.ticklabel_format(style='sci',axis='y')
		#ax.set_xticks([-800,-400,0,400,800])
		ax.set_xlabel(r"$P_x [MeV/c]$")
		ax.set_xticks([-500,0,500])
		ax.set_ylabel(r"$\rho [fm^3]$")
		ax.set_title(makeTitle(var))
	if len(yind > 0):
		fig = pl.figure(figsize=(8,6))
		ax = fig.add_subplot(111)
		fig.subplots_adjust(bottom=0.2,left=0.2)
		#ax.plot(X[len(X)/3:2*len(X)/3,1],X[len(X)/3:2*len(X)/3,3])
		ax.plot(X[yind,1],X[yind,3])
		ax.ticklabel_format(style='sci',axis='y')
		#ax.set_xticks([-800,-400,0,400,800])	
		ax.set_xticks([-500,0,500])
		ax.set_xlabel(r"$P_y [MeV/c]$")
		ax.set_ylabel(r"$\rho [fm^3]$")
		ax.set_title(makeTitle(var))
	if len(zind > 0):
		fig = pl.figure(figsize=(8,6))
		fig.subplots_adjust(bottom=0.2,left=0.2)
		ax = fig.add_subplot(111)
		#ax.plot(X[2*len(X)/3:len(X),2],X[2*len(X)/3:len(X),3])
		ax.plot(X[zind,2],X[zind,3])
		ax.ticklabel_format(style='sci',axis='y')
		#ax.set_xticks([-800,-400,0,400,800])
		ax.set_xticks([-500,0,500])
		ax.set_xlabel(r"$P_z [MeV/c]$")
		ax.set_ylabel(r"$\rho [fm^3]$")
		ax.set_title(makeTitle(var))
X = None
firstLine = True
for line in f:
	if (line[0] == '*'):
		#print var
		if not(firstLine): # if not first line: X,var contains data, we can make a plot
			makePlot(X,var)
		line = line.lstrip('*').rstrip('\n')
		#print line
		verts = line.split('|')
		var = getvars(verts)
		X = [] # clear data
		firstLine=False
	else:
		verts = line.split() # split in spaces, first three colums are X,Y,Z of P third is two body dist
		X.append([float(verts[0]),float(verts[1]),float(verts[2]),float(verts[3])])

# now last X has not been plotted yet because file does not end with an variable line so call makePlot(X,var) one last time
makePlot(X,var)

pl.show()
