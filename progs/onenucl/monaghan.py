# -*- coding: utf-8 -*-
import numpy as np
import math

Q2= np.float_([1.696,1.719,1.742,1.768,1.796,1.828,1.869,1.917,1.977,2.047,2.114,1.741,1.770,1.798,1.833,1.877,2.011,2.081,2.146,2.203])
Q2*=1.E06
w=np.float_([851.5,846.0,841.2,835.9,829.7,823.7,819.0,815.4,813.8,813.1,812.7,821.8,816.3,809.9,803.8,799.0,813.5,812.9,812.0,808.9])
pm=np.float_([190,210,230,250,270,290,310,330,350,370,390,350,370,390,410,430,360,380,400,420])
q=np.zeros(len(w))
pf=np.zeros(len(w))
thetaf=np.zeros(len(w))
thetam=np.zeros(len(w))
for i in np.arange(0,len(w)):
	q[i]=math.sqrt(w[i]*w[i]+Q2[i])
MAmin1=10255.12
MA=11177.93
MP=938.27231
for i in np.arange(0,len(w)):
	pf[i] = math.sqrt((w[i]+MA-math.sqrt(MAmin1*MAmin1+pm[i]*pm[i]))*(w[i]+MA-math.sqrt(MAmin1*MAmin1+pm[i]*pm[i]))-MP*MP)
print pf
for i in np.arange(0,len(w)):
	theta[i]=math.acos((pm[i]*pm[i]-pf[i]*pf[i]-q[i]*q[i])/(-2.*pf[i]*q[i]))*57.2957795
	print q[i], math.cos(theta[i])*pm[i]+

#f = open('monaghan.sh', 'w')
#for i in np.arange(0, len(w)):
  #f.write('./bin/kak '+str(Q2[i])+' '+str(w[i])+' '+str(pm[i])+'\n')
#f.close()

