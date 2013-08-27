# -*- coding: utf-8 -*-
import numpy as np
import math

Ein=3245
Eout=2255
the=28.6/180*math.pi
Q2=4*Ein*Eout*math.sin(the/2)*math.sin(the/2)
print Q2


Q2=1820000
w=990

#f = open('monaghan_dutta.sh', 'w')
#for i in np.arange(0, 33):
  #p=10*i
  #f.write('lua Monaghan.lua '+str(Q2)+' '+str(w)+' '+str(p)+'\n')
