# -*- coding: utf-8 -*-
import numpy as np
import math

Q2=1820000
w=990

f = open('monaghan_dutta.sh', 'w')
for i in np.arange(0, 33):
  p=10*i
  f.write('lua Monaghan.lua '+str(Q2)+' '+str(w)+' '+str(p)+'\n')
