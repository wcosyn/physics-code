#include "DeuteronCross.hpp"

DeuteronCross::DeuteronCross(TElectronKinematics elec, string name, bool proton, string strucname,
			     double sigmain, double betain, double epsilonin, double betaoffin, double lambdain, int offshellset):
massi(proton? MASSP:MASSN),
electron(elec),
momdistr(name,massi,offshellset,sigmain,betain,epsilonin,betaoffin,lambdain),
structure(elec,proton,strucname)
{
  
  
}


DeuteronCross::~DeuteronCross(){
  
}

double DeuteronCross::getavgCross(TKinematics2to2 &kin,int pw, double Einoff){
  
  double y=kin.GetWlab()/electron.GetBeamEnergy(kin); 
  double x=kin.GetQsquared()/(2.*massi*kin.GetWlab());
  double front=(4.*PI*ALPHA*ALPHA)/(x*kin.GetQsquared()*kin.GetQsquared())
		*(1-y-(x*x*y*y*massi*massi)/kin.GetQsquared());
  double dens=pw?momdistr.getMomDistrpw(kin,0.):momdistr.getMomDistrfsi(kin,0.);
  double Dstrucs=structure.getavgStructure(kin,Einoff);
  return front*dens*Dstrucs*sqrt(kin.GetPklab()*kin.GetPklab()+kin.GetMesonMass()*kin.GetMesonMass())*HBARC*HBARC*1.E19;
  
}
