# -*- coding: utf-8 -*-

#inputs
#1. qindex = (0=2.5, 1=5, 2=10)
#2. symmetry (1=no mass terms in prop, 0=always mass terms in prop, -1=mass terms in prop like in SIDIS)
#3. fix (1=fixed resonance masses in prop terms, 0=mas terms depend on spectator integration loop)
#4. 1=fixed sigma values, 0=parametrized sigma values

import numpy as np
import pylab as p
import matplotlib.ticker as ticker
import matplotlib.font_manager
import sys

fig=p.figure(1,figsize=(14,10))
fig.subplots_adjust(left=0.08,right=0.98,bottom=0.08, top=0.98, wspace=0.18, hspace=0.22)


data= np.loadtxt(sys.argv[1])
plots=[None]*3
title=[r'$\sigma$=160 MeV',r'$\sigma$=180 MeV',r'$\sigma$=200 MeV']

for i in np.arange(0,3):
  plots[i], =p.plot(data[:,0],data[:,i+2]/data[:,1])

  
  
ax=p.gca()

#p.axis([0,1,0.9,2.])
p.text(0.9,-0.08,r'$p [MeV]$',horizontalalignment='center',
      verticalalignment='bottom',transform=ax.transAxes, fontsize=30)
p.text(-0.08,0.65,r'ratio',horizontalalignment='left',
	verticalalignment='center',rotation='vertical',transform=ax.transAxes, fontsize=30)

ticks=ax.xaxis.get_major_ticks()
for i in np.arange(1,len(ticks),2):
  ticks[i].label1On = False
ticks=ax.yaxis.get_major_ticks()
for i in np.arange(1,len(ticks),2):
  ticks[i].label1On = False
for tick in ax.xaxis.get_major_ticks():
  tick.label1.set_fontsize(20)
for tick in ax.yaxis.get_major_ticks():
  tick.label1.set_fontsize(20)

  
#p.subplot(len(W),len(ps),25)  
fontprop = matplotlib.font_manager.FontProperties(size=25)
p.legend(plots,title, prop=fontprop, bbox_to_anchor=(1., 1.), loc=1, borderaxespad=0.)
  
p.show()