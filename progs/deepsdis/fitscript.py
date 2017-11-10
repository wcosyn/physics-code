import numpy as np
import sys
import math, copy, cmath
import shlex, subprocess
import pylab as p
import matplotlib.ticker as ticker
import matplotlib.font_manager
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from matplotlib import cm
#from pylab import *
#from scipy import *
from scipy import optimize


def apply(poly, x):
    """ Apply a value to a polynomial (simple). """ 
    sum = 0 
    for term in poly:
        sum += term[0] * (x ** term[1]) 
    return sum


def floatorcomplex(x):
  if not isinstance(x, complex):
      return float(x)
  else:
      return x
      
def simplify(poly):
  """ Simplify a polynomial (collect terms) """

  d = {}
  newlist = []

  # collect terms
  for term in poly:
      if d.has_key(term[1]):
	  d[term[1]] += term[0]
      else:
	  d[term[1]] = term[0]

  # remove terms with zero coefficient
  for k in d.keys():
      if d[k] == 0:
	  d.pop(k)

  keys = d.keys()
  keys.sort()
  keys.reverse()

  for exponent in keys:
      newlist.append([d[exponent], exponent])

  return newlist
 
def add(poly1, poly2):
    """ Add two polynomials """
 
    newlist = copy.deepcopy(poly1)
 
    newlist.extend(poly2)
 
    return simplify(newlist)


def multiply(poly1, poly2):
    """ Multiply two polynomials """
 
    newlist = []
 
    for term1 in poly1:
        for term2 in poly2:
            newlist.append([term1[0] * term2[0], term1[1] + term2[1]])
 
    return simplify(newlist)
 

def interpolate(inputs, outputs):
  """ Uses the Newton polynomial to compute an interpolation polynomial using the inputs and outputs """
  levels = []
  count = len(inputs)
  poly = []

  if count == 0:
      return []

  if count == 1:
      return [[outputs[0], 0]]

  levels.append(outputs)

  for i in xrange(count):
      if i == 0:
	  continue

      newlevel = []

      for j in xrange(count):
	  if j >= count - i:
	      continue

	  newlevel.append(floatorcomplex(levels[i - 1][j + 1] - levels[i - 1][j]) / (inputs[i + j] - inputs[j]))

      levels.append(newlevel)

  for i in xrange(count):
      in_poly = [[levels[i][0], 0]]

      for j in xrange(i):
	  in_poly = multiply(in_poly, [[1, 1], [-inputs[j], 0]])

      poly = add(poly, in_poly)

  return simplify(poly)


def diff(poly):
    """ Differentiate a polynomial. """
 
    newlist = copy.deepcopy(poly)
 
    for term in newlist:
        term[0] *= term[1]
        term[1] -= 1
 
    return simplify(newlist)
 

#fitset=int(sys.argv[2])  #how many free parameters in the fit, 4 has epsilon restricted to negative values
#fitfile=sys.argv[1] #file with scattering parameter fits to the data
numpoints=int(sys.argv[2]) #number of data points used to do the fit
Q2index=int(sys.argv[1]) #Q2 value index
costhetaindex=int(sys.argv[3]) #spectator costh index
signifier=sys.argv[4] #idintifier for output filenames
plot=int(sys.argv[5]) #yes or no?

costhetar=-0.9+costhetaindex*0.2
if(costhetar>1):
  print 'illegal costhetaindex [0-9]'
  sys.exit()
  

#filename='/home/wim/DeuteronDIS/fitcode/output/'+fitfile
#calc=np.loadtxt(filename)

#sigma=[['5.84321e+01','3.59042e+01','4.76730e+01','5.10952e+01','6.18136e+01'],
#['6.13945e+01','3.21320e+01','4.10246e+01','4.27569e+01','5.33368e+01']]


#Q=['','18','28']
Qcom=[0.93,1.665,3.38]
#W=['125','15','173','2','24']
Wcom=[1.175,1.475,1.725,2.025,2.44]
#ps=['300','340','390','460']
pscom=[0.078,0.093,0.110,0.135]


Q2=Qcom[Q2index]
print 'Computing for Q2 = '+str(Q2)+'GeV^2'
print 'Using fit from deeps analysis for sigma, fixed beta and epsilon'
print r'cos(\theta_r) = '+str(costhetar)
sigma=np.ones((2,5))
beta=np.ones((2,5))
eps=np.ones((2,5))
#sigma*=40.
#fix values of the parameters
beta*=8.
eps*=-0.5


Wpoints=np.zeros(numpoints)
Windex=np.zeros(numpoints)
pspoints=np.zeros(numpoints)
psindex=np.zeros(numpoints)
print 'Possible values for W (GeV):'+ str(Wcom)
print 'Possible values for ps (GeV) :'+ str(pscom)

#get points to which we will perform the fit
#successive points should be decreasing in W (go up in Bjorken x), and strict decreasing in ps (go down in offshellness of p_i)
#Bjorken x along the trajectory should be <0.5 because there we have good limits on the structure function
for i in np.arange(0,numpoints):
  print 'Enter W index for point '+str(i+1)+' of '+str(numpoints)+' [0-4]:'
  val=int(raw_input(''))
  Windex[i]=val
  Wpoints[i]=Wcom[val]
  if (val<0 or val>4): 
    print 'illegal value'
    sys.exit()  
  if (i!=0 and Wpoints[i]>Wpoints[i-1]):
    print 'indexes of W should be descending row'
    sys.exit()
  print 'Enter ps index for point '+str(i+1)+' of '+str(numpoints)+' [0-3]:'
  val=int(raw_input(''))
  psindex[i]=val
  pspoints[i]=pscom[val]
  if (val<0 or val>4): 
    print 'illegal value'
    sys.exit()
  if (i!=0 and psindex[i]>=psindex[i-1]):
    print 'indexes of ps should be strict descending row'
    sys.exit()




MASSD=1.8756
massi=0.93956536
massr=0.93827231

#these are all arrays!!  for all the considered points
Er=np.sqrt(pspoints*pspoints+massr*massr)
Einoff=MASSD-Er
massoff=np.sqrt(Einoff*Einoff-pspoints*pspoints)
tprime = massi*massi-massoff*massoff
xprime=Q2/(Wpoints*Wpoints-massoff*massoff+Q2)

print "ps:", pspoints
print "x':",xprime
print "t':",tprime
print "W:", Wpoints


print 'final x value' +str(xprime[numpoints-1])+'):'
#print 'Give x value (x must be larger than '+str(xprime[numpoints-1])+'):'
#xfinal=float(raw_input(''))
#if(xfinal<xprime[numpoints-1]):
  #print 'not a good x value'
  #sys.exit()
xfinal=Q2/(Wpoints*Wpoints-massi*massi+Q2)
print xfinal

#sys.exit()


for i in np.arange(0,numpoints):
  #find index for data file
  #calc the quantity that will converge to F2(x_ref) for t'->0
  #and put output in array
  #output is [xprime,W,pr,costhetar,alphar,t',F2extrpw ,F2extrFSI,F2ref,bonusdata,bonuserror]
  #alpha_r isn't used anywhere in the code, just for reference I guess (check difference over all points)
  commandline='/home/wim/Code/trunk/bin/bonusextrapolate 0 ' +str(Q2index)+' '+str(Windex[i])+' '+str(psindex[i])+' '+str(costhetaindex)+' 3 1 '+str(xfinal)
  #print commandline
  args=shlex.split(commandline)
  prun=subprocess.Popen(args, stdout=subprocess.PIPE)
  outarray=shlex.split(prun.communicate()[0])
  if i==0: 
    datacomparray=[outarray]
  else:
    datacomparray=np.append(datacomparray,[outarray],0)
  if(Q2index>0):
    commandline='/home/wim/Code/trunk/bin/bonusextrapolate 1 ' +str(Q2index)+' '+str(Windex[i])+' '+str(psindex[i])+' '+str(costhetaindex)+' 3 1 '+str(xprime[0])
    print commandline
    args=shlex.split(commandline)
    prun=subprocess.Popen(args, stdout=subprocess.PIPE)
    outarray=shlex.split(prun.communicate()[0])
    datacomparray=np.append(datacomparray,[outarray],0)
    

f = open('/home/wim/Calculations/Bonus/results/'+signifier+'.q'+str(Q2index)+'.'+str(costhetaindex)+'.data.dat', 'w')
for i in np.arange(0, numpoints):
  f.write(' '.join(n for n in datacomparray[i][:])+'\n')
if(Q2index>0):
  for i in np.arange(0, numpoints):
    f.write(' '.join(n for n in datacomparray[i+numpoints][:])+'\n')

f.close()



nufinal=Q2/(2.*massi*xfinal)
xprimefinal=Q2/(2.*(MASSD-massr)*nufinal)
Wpointsfinal=np.sqrt((MASSD-massr)*(MASSD-massr)-Q2+Q2/xprimefinal)



#for i in np.arange(1,10):
  #pspoints=np.append(pspoints,pspoints[numpoints-1]+i*0.01*(0.-pspoints[numpoints-1]))
  #Wpoints=np.append(Wpoints,Wpoints[numpoints-1]+i*0.01*(Wpointsfinal-Wpoints[numpoints-1]))

#pspoints=np.append(pspoints,0.)
#Wpoints=np.append(Wpoints,Wpointsfinal)
print "final W", Wpointsfinal

#here we try to find a trajectory 
print "interpolating..."

#xprime as a function of p_s
poly1=interpolate(pspoints,xprime)

#derivative value at last data point used
dpoly1=diff(poly1)
AA=apply(dpoly1,pspoints[numpoints-1])

#trajectory that connects last used data point kinematics with extrapolation point
poly2=[[(xprime[numpoints-1]-xfinal)/math.pow(pspoints[numpoints-1],2),2],[xprimefinal,0]]

#poly2=interpolate([pspoints[numpoints-1],0.],[Wpoints[numpoints-1],Wpointsfinal])

print "calculating trajectory...."
f = open('/home/wim/Calculations/Bonus/results/'+signifier+'.q'+str(Q2index)+'.'+str(costhetaindex)+'.traj.dat', 'w')

#trajectory part in between data points used
for j in np.arange(0,40):

  ps=pspoints[0]+j*0.025*(pspoints[numpoints-1]-pspoints[0])
  xxprime=apply(poly1,ps)
  #Wlower,Wupper, sigmalower, sigmaupper, betalower, betaupper, epslower, epsupper
  xEr=np.sqrt(ps*ps+massr*massr)
  xEinoff=MASSD-xEr
  xmassoff=np.sqrt(xEinoff*xEinoff-ps*ps)
  xtprime = massi*massi-xmassoff*xmassoff
  W=Wpoints[0]

  commandline='/home/wim/Code/trunk/bin/bonusextrapolate2 1 ' +str(Q2*1.E06)+' '+str(W*1.E03)+' '+str(ps*1.E03)+' '+str(costhetaindex)+' 3 1 '+str(xfinal)
  #print commandline
  args=shlex.split(commandline)
  prun=subprocess.Popen(args, stdout=subprocess.PIPE)
  outarray=shlex.split(prun.communicate()[0])
  f.write(' '.join(n for n in outarray))
  f.write('\n')
  
#trajectory part from last used data point to p_s=0
for j in np.arange(0,60):
  ps=pspoints[numpoints-1]+1./60*j*(-pspoints[numpoints-1])


  xxprime=apply(poly2,ps)
  #Wlower,Wupper, sigmalower, sigmaupper, betalower, betaupper, epslower, epsupper
  xEr=np.sqrt(ps*ps+massr*massr)
  xEinoff=MASSD-xEr
  xmassoff=np.sqrt(xEinoff*xEinoff-ps*ps)
  xtprime = massi*massi-xmassoff*xmassoff
  W=Wpoints[0]
  #print j, ps, xtprime, W, xxprime

  commandline='/home/wim/Code/trunk/bin/bonusextrapolate2 1 ' +str(Q2*1.E06)+' '+str(W*1.E03)+' '+str(ps*1.E03)+' '+str(costhetaindex)+' 3 1 '+str(xprime[0])
  args=shlex.split(commandline)
  prun=subprocess.Popen(args, stdout=subprocess.PIPE)
  outarray=shlex.split(prun.communicate()[0])
  f.write(' '.join(n for n in outarray))
  f.write('\n')
f.close()

print "Finished calculating trajectory"


print signifier, Q2index, costhetaindex
calcresult=np.loadtxt('/home/wim/Calculations/Bonus/results/'+signifier+'.q'+str(Q2index)+'.'+str(costhetaindex)+'.traj.dat') #trajectory calculations results
data=np.loadtxt('/home/wim/Calculations/Bonus/results/'+signifier+'.q'+str(Q2index)+'.'+str(costhetaindex)+'.data.dat')

print "Calculating quadratic fit to data using info from theoretical curve..."
##minimum of the curve does not depend on model, we use that here to fix tmin
#F2min=1.
#tmin=0
#for i in np.arange(40,len(calcresult[:,0])):
  #if(calcresult[i,7]<F2min):
    #F2min=calcresult[i,7]
    #tmin=calcresult[i,5]

#print "Minimum in calc fsi curve is for t' = ", tmin

commandlinepre='/home/wim/Code/trunk/bin/parabfit '
if(Q2index>0):
  commandlinepre+=str(numpoints*2)
else:
  commandlinepre+=str(numpoints)
for i in np.arange(0,numpoints):
  commandlinepre+=' '+str(data[i,5])+' '+str(data[i,9])+' '+str(data[i,10])
if(Q2index>0):
  for i in np.arange(0,numpoints):
    commandlinepre+=' '+str(data[numpoints+i,5])+' '+str(data[numpoints+i,9])+' '+str(data[numpoints+i,10])

print commandlinepre
#central data value
print "Fit 1 starting, calling Minuit..."  
prun=subprocess.Popen(shlex.split(commandlinepre), stdout=subprocess.PIPE)
fitoutput=shlex.split(prun.communicate()[0])
aa=float(fitoutput[-3])
bb=float(fitoutput[-2])
cc=float(fitoutput[-1])

parafit1=[[aa,2],[-2.*aa*bb,1],[aa*bb*bb+cc,0]]
print aa, -2.*aa*bb, aa*bb*bb+cc

print "making a parabole fit to data points in the low x region"
fitfunc=lambda p, x: p[0]*x*x+p[1]*x+p[2]
errfunc = lambda p, x, y: fitfunc(p, x) - y
p0 = [1., -1., 1.] # Initial guess for the parameters

if(Q2index>0):
  tprime2=np.zeros(numpoints*2)
  for i in np.arange(0,len(tprime)):
    tprime2[2*i]=tprime[i]
    tprime2[2*i+1]=tprime[i]
else:
  tprime2=tprime
  
#print tprime2, data[:,9]
#p1, success = optimize.leastsq(errfunc, p0[:], args=(tprime2, data[:,9]))
#parab1=[[p1[0],2],[p1[1],1],[p1[2],0]] #fit without constraints to data (doesn't work with the tmin restriction).


#lower=0
##only fit to the 0.3<xprime <0.5 region
#for i in np.arange(0,len(calcresult[:,0])):
  #if calcresult[i,0]>0.3:
    #lower=i
    #break
#else:
  #lower=len(calcresult[:,0])-1
#upper=lower
#for i in np.arange(lower,len(calcresult[:,0])):
  #if calcresult[i,0]>0.5:
    #upper=i
    #break
  
#print lower, upper
#if(upper==lower):
  #upper=len(calcresult[:,0])
#print "making a parabole fit to the plane wave calculation in the low x region"
#p2, success = optimize.leastsq(errfunc, p0[:], args=(calcresult[lower:upper,5],calcresult[lower:upper,6]))  
#parab2=[[p2[0],2],[p2[1],1],[p2[2],0]]
#print "making a parabole fit to the fsi calculation in the low x region"
#p3, success = optimize.leastsq(errfunc, p0[:], args=(calcresult[lower:upper,5],calcresult[lower:upper,7]))  
#parab3=[[p3[0],2],[p3[1],1],[p3[2],0]]

print '\n\n\n***********************************\nResult for F2 at Q2 = '+str(Q2)+' and x = '+str(xfinal)+':'
print 'middle: '+str(apply(parafit1,0.))
print 'CB value: '+str(calcresult[0,8])

if(plot):
  print "Preparing plot"
  mpl.rcParams['font.size'] = 12
  mpl.rcParams['axes.titlesize'] = 'larger'
  mpl.rcParams['xtick.labelsize'] = 'small'

  #fig=p.figure(1,figsize=(10,8))
  #ax = Axes3D(fig) 
  #ax.plot(calcresult[:,0], calcresult[:,5], calcresult[:,6], label='pw')
  #ax.plot(calcresult[:,0], calcresult[:,5], calcresult[:,7], label='fsi')
  #ax.scatter(data[:,0],data[:,5],data[:,9])
  ##ax.plot_surface(X,Y,ratio1,rstride=1, cstride=1, cmap=cm.Greys)
  #ax.set_ylabel(r'$t^\prime$ [GeV$^2$]')
  #ax.set_xlabel(r'$x^\prime$')
  #ax.set_zlabel('F2')

  #ax.legend()

  #p.show()
  


  

  #fig2=p.figure(1,figsize=(10,8))
  #title=['pw','fsi', 'fitpw','fitfsi', 'fitdata']
  #plots=[None]*5
  #plots[0], =p.plot(calcresult[:,5],calcresult[:,6])
  #plots[1], =p.plot(calcresult[:,5],calcresult[:,7])
  #p.errorbar(data[:,5],data[:,9],data[:,10], fmt='o',color='k')
  
  #xpar = np.linspace(0., calcresult[0,5], 100)
  #plots[2], =p.plot(xpar,apply(parab2,xpar),'--')
  #plots[3], =p.plot(xpar,apply(parab3,xpar),'--')
  #plots[4], =p.plot(xpar,apply(parab1,xpar),'--')
  ##p.plot([calcresult[lower,5],calcresult[upper,5]],[0.,0.],'o')
  #ax=p.gca()
  #fontprop = matplotlib.font_manager.FontProperties(size=12)
  #p.legend(plots,title, prop=fontprop, bbox_to_anchor=(0., 0.), loc=3, borderaxespad=0.)
  ##p.axis([0,6,0,0.65])
  #p.text(0.5,-0.1,'$t^\prime [\mathrm{GeV}^2]$',horizontalalignment='center',
	    #verticalalignment='bottom',transform=ax.transAxes, fontsize=20)
  #p.text(-0.12,0.14,r'$F_{2N}$',horizontalalignment='left',
	    #verticalalignment='center',rotation='vertical',transform=ax.transAxes, fontsize=20)
  #p.text(0.17,0.9,r'$x=$'+str(xfinal)+r'$, Q^2=$'+str(Q2)+r'GeV$^2$, $\cos\theta=$'+str(costhetar),horizontalalignment='left',
	    #verticalalignment='bottom',transform=ax.transAxes, fontsize=20)
  #ticks=ax.xaxis.get_major_ticks()
  #for i in np.arange(1,len(ticks),2):
    #ticks[i].label1On = False
  #for tick in ax.xaxis.get_major_ticks():
    #tick.label1.set_fontsize(14)
  #for tick in ax.yaxis.get_major_ticks():
    #tick.label1.set_fontsize(14)
  ##ax.annotate(r'$x^\prime=0.3$', xy=(calcresult[lower,5],0.0), xycoords='data',xytext=(calcresult[lower,5],-0.05), arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=12)
  ##ax.annotate(r'$x^\prime=0.5$', xy=(calcresult[upper,5],0.0), xycoords='data',xytext=(calcresult[upper,5],-0.05), arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=12)
  #p.show() 
  
  fig3=p.figure(1,figsize=(10,6))
  fig3.patch.set_alpha(0.0)
  ai=p.subplot(1,1,1)
  ai.patch.set_alpha(0.0)
  title=['IA','FSI', 'Pole extrapolation']
  plots=[None]*3
  plots[0], =p.plot(calcresult[:,5],calcresult[:,6],linewidth=2)
  plots[1], =p.plot(calcresult[:,5],calcresult[:,7],linewidth=2)
  p.errorbar(data[:,5],data[:,9],data[:,10], fmt='o',color='k')
  
  xpar = np.linspace(0., calcresult[0,5], 100)
  plots[2], =p.plot(xpar,apply(parafit1,xpar),'--',linewidth=2)
  #plots[3], =p.plot(xpar,apply(parafit2,xpar),'--')
  #plots[4], =p.plot(xpar,apply(parafit3,xpar),'--')
  ax=p.gca()
  fontprop = matplotlib.font_manager.FontProperties(size=18)
  leg=p.legend(plots,title, prop=fontprop, bbox_to_anchor=(0., 0.), loc=3, borderaxespad=0.)
  leg.get_frame().set_alpha(0)
  #p.axis([0,6,0,0.65]
  p.text(0.93,-0.13,'$t^\prime [\mathrm{GeV}^2]$',horizontalalignment='center',
	    verticalalignment='bottom',transform=ax.transAxes, fontsize=20)
  p.text(-0.14,0.8,r'$F_{2N}^{extr}$',horizontalalignment='left',
	    verticalalignment='center',rotation='vertical',transform=ax.transAxes, fontsize=20)
  p.text(0.17,0.9,r'$x=$'+str("%0.2f" % xprime[0])+r'$, Q^2=$'+str(Q2)+r'GeV$^2$,$\cos\theta=$'+str(costhetar),horizontalalignment='left',
	    verticalalignment='bottom',transform=ax.transAxes, fontsize=20)
  p.text(0.17,0.8,r'$F_{2N}=$'+str("%0.4f" % apply(parafit1,0.)),horizontalalignment='left',
	    verticalalignment='bottom',transform=ax.transAxes, fontsize=20)
  ax.tick_params(direction='in', pad=6)
  ticks=ax.xaxis.get_major_ticks()
  for i in np.arange(1,len(ticks),2):
    ticks[i].label1On = False
  ticks=ax.yaxis.get_major_ticks()
  for i in np.arange(1,len(ticks),2):
    ticks[i].label1On = False
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(18)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(18)
  #ax.annotate(r'$x^\prime=0.3$', xy=(calcresult[lower,5],0.0), xycoords='data',xytext=(calcresult[lower,5],-0.05), arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=12)
  #ax.annotate(r'$x^\prime=0.5$', xy=(calcresult[upper,5],0.0), xycoords='data',xytext=(calcresult[upper,5],-0.05), arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=12)
  p.show()   
