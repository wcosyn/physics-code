#include "He3Cross.hpp"
#include <fstream>
#include <FourVector.h>
#include <TVector3.h>

using namespace std;

He3Cross::He3Cross(string wfname, string inputdir, string strucname, bool proton_): 
  proton(proton_),
  wavefunction(wfname,inputdir+"/he3"),
  strucfunc(proton,strucname){
  
}

He3Cross::~He3Cross(){
}

double He3Cross::getDensity(LightConeKin2to3 &kin){
  double result=0.;
  for(int mA=-1;mA<=1;mA+=2){
    for(int ms1=-1;ms1<=1;ms1+=2){
      for(int ms2=-1;ms2<=1;ms2+=2){
	for(int mt1=(proton?-1:1);mt1<=1;mt1+=2){
	  int mt2=1-mt1-(proton?1:-1);
	  TVector3 sp1(kin.getSp1_mu()[1]-kin.getA_mu()[1]/3.,
		       kin.getSp1_mu()[2]-kin.getA_mu()[2]/3.,
		       kin.getSp1_mu()[3]-kin.getA_mu()[3]/3.);
	  TVector3 sp2(kin.getSp2_mu()[1]-kin.getA_mu()[1]/3.,
		       kin.getSp2_mu()[2]-kin.getA_mu()[2]/3.,
		       kin.getSp2_mu()[3]-kin.getA_mu()[3]/3.);
	  
	  result+=norm(wavefunction.getWF(sp1,sp2,ms1,ms2,mt1,mt2,mA));		
	    
	}
      }
    }
  }
  return result*0.5; //He3 spin average
}

double He3Cross::getCross(LightConeKin2to3 &kin){
  return 2.*pow(kin.getYA()*ALPHA/kin.getQ2(),2.)*getDensity(kin)*strucfunc.getStructureLC(kin)
      *kin.getSp1_mu()[0]*kin.getSp1_mu()[0]*HBARC*HBARC*1.E25;  //conversion from MeV-8 to nb/GeV-6
}