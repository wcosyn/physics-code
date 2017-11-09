#script used to test kinematics for rho decay


import math
import numpy as np
import pylab as p
import matplotlib.ticker as ticker
import matplotlib.font_manager
import sys

p=float(sys.argv[1])

massrho=775.49
masspi=139.57
Gamma=145.
E=np.sqrt(massrho*massrho+p*p)
gamma=1/np.sqrt(1.-p*p/(E*E))
dil_decay=Gamma/gamma
print 1/dil_decay*197.
pcm=np.sqrt(massrho*massrho/4.-masspi*masspi)
Ecm = massrho/2.
for i in np.arange(0,20):  
  costh=1.-0.1*i
  pz=pcm*costh
  px=pcm*np.sqrt(1.-costh*costh)
  pzlab=gamma*(pz+p/E*Ecm)
  pzlab2=gamma*(-pz+p/E*Ecm)
  plab=np.sqrt(px*px+pzlab*pzlab)
  plab2=np.sqrt(px*px+pzlab2*pzlab2)
  theta1=np.arccos(pzlab/plab)*180/math.pi
  theta2=np.arccos(pzlab2/plab2)*180/math.pi
  print costh, theta1+theta2, theta1, theta2