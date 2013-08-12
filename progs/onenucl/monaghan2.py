# -*- coding: utf-8 -*-
import numpy as np
import math

Q2= np.float_([1.696,1.719,1.742,1.768,1.796,1.828,1.869,1.917,1.977,2.047,2.114,1.741,1.770,1.798,1.833,1.877,2.011,2.081,2.146,2.203])
Q2*=1.E06
w=np.float_([851.5,846.0,841.2,835.9,829.7,823.7,819.0,815.4,813.8,813.1,812.7,821.8,816.3,809.9,803.8,799.0,813.5,812.9,812.0,808.9])
pm=np.float_([190,210,230,250,270,290,310,330,350,370,390,350,370,390,410,430,360,380,400,420])
q=np.zeros(len(w))
pf=np.zeros(len(w))
t=np.zeros(len(w))
u=np.zeros(len(w))
beta=np.zeros(len(w))
gamma=np.zeros(len(w))
pzcm=np.zeros(len(w))
pmx2=np.zeros(len(w))
cthetaf=np.zeros(len(w))
cthetam=np.zeros(len(w))
cthetacm=np.zeros(len(w))
for i in np.arange(0,len(w)):
	q[i]=math.sqrt(w[i]*w[i]+Q2[i])
MAmin1=10255.12
MA=11177.93
MP=938.27231
for i in np.arange(0,len(w)):
	pf[i] = math.sqrt((w[i]+MA-math.sqrt(MAmin1*MAmin1+pm[i]*pm[i]))*(w[i]+MA-math.sqrt(MAmin1*MAmin1+pm[i]*pm[i]))-MP*MP)
#print pf
for i in np.arange(0,len(w)):
	cthetaf[i]=(pm[i]*pm[i]-pf[i]*pf[i]-q[i]*q[i])/(-2.*pf[i]*q[i])
	cthetam[i]=(pf[i]*pf[i]-pm[i]*pm[i]-q[i]*q[i])/(-2.*pm[i]*q[i])	
	t[i]=-Q2[i]+MAmin1*MAmin1-2.*w[i]*math.sqrt(MAmin1*MAmin1+pm[i]*pm[i])+2.*q[i]*pm[i]*cthetam[i]
	u[i]=-Q2[i]+MP*MP-2.*w[i]*math.sqrt(MP*MP+pf[i]*pf[i])+2.*q[i]*pf[i]*cthetaf[i]
	beta[i]=-q[i]/(MA+w[i])
	gamma[i]=1/math.sqrt(1-beta[i]*beta[i])
	pzcm[i]=gamma[i]*(pm[i]*cthetam[i]+beta[i]*math.sqrt(MAmin1*MAmin1+pm[i]*pm[i]))
	pmx2[i]=pm[i]*pm[i]*(1-cthetam[i]*cthetam[i])
	print q[i], cthetaf[i], cthetam[i], t[i], u[i], math.sqrt(pzcm[i]*pzcm[i]+pmx2[i])

f = open('monaghan_trillian.sh', 'w')
for i in np.arange(0, len(w)):
  f.write('lua Monaghan.lua '+str(Q2[i])+' '+str(w[i])+' '+str(pm[i])+'\n')
f.close()

