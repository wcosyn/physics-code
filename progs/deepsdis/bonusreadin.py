# script generates input to run all miniboone calculations on the hpc

# -*- coding: utf-8 -*-
import numpy as np
import pylab as p
import matplotlib.ticker as ticker
import sys
import math

data=np.loadtxt('bonusdata4.txt')

print "//Bonus data beam 4GeV"
print "//first index is Q^2={0.93,1.66,3.38}"
print "//second index is W={1.17,1.48,1.73,2.03,2.44}"
print "//third index is p_s[MeV]={78,93,110,135}"
print "//fourth index is costh={-0.9,-0.7,..,0.9} (10 values)"
print "//fifth index is [costh value, bonus ratio value, stat error, syst error]"

print "bonusdata4[3][5][4][10][4]={ "
for Q2 in np.arange(0,3):
  print "{ "
  for W in np.arange(0,5):
    print "{ "
    for ps in np.arange(0,4):
      print "{ "
      for costh in np.arange(0,10):
	print "{",data[ps*10+W*40+Q2*200+costh,0],",",data[ps*10+W*40+Q2*200+costh,1],",",data[ps*10+W*40+Q2*200+costh,2],",",data[ps*10+W*40+Q2*200+costh,3],"}"
	print "" if costh==9 else ","
      print "}" if ps==3 else "},"
    print "}" if W==4 else "},"
  print "}" if Q2==2 else "},"
print "}\n\n";
	  
print "//Bonus data beam 4GeV"
print "//first index is Q^2={1.66,3.38}"
print "//second index is W={1.17,1.48,1.73,2.03,2.44}"
print "//third index is p_s[MeV]={78,93,110,135}"
print "//fourth index is costh={-0.9,-0.7,..,0.9} (10 values)"
print "//fifth index is [costh value, bonus ratio value, stat error, syst error]"

data=data('bonusdata5.txt')

print "bonusdata4[2][5][4][10][4]={ "
for Q2 in np.arange(0,2):
  print "{ "
  for W in np.arange(0,5):
    print "{ "
    for ps in np.arange(0,4):
      print "{ "
      for costh in np.arange(0,10):
	print "{",data[ps*10+W*40+Q2*200+costh,0],",",data[ps*10+W*40+Q2*200+costh,1],",",data[ps*10+W*40+Q2*200+costh,2],",",data[ps*10+W*40+Q2*200+costh,3],"}"
	print "" if costh==9 else ","
      print "}" if ps==3 else "},"
    print "}" if W==4 else "},"
  print "}" if Q2==2 else "},"
print "}";
	  
      