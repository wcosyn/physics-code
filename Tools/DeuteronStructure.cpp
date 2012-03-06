#include "DeuteronStructure.hpp"
#include "NuclStructure.hpp"

#include <cmath>
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <cstdlib>



DeuteronStructure::DeuteronStructure(TElectronKinematics &el, int pr, string nm)
:electron(el),proton(pr),name(nm),massi(proton?MASSP:MASSN){
  
}


void DeuteronStructure::getStructureFunctions(TKinematics2to2 &kin, double &FL, double &FT, double &FTT, double &FTL,double Einoff) const{
  double alphaq=(kin.GetWlab()-kin.GetKlab())/(MASSD/2.);
  double alphai=(Einoff+kin.GetPklab()*kin.GetCosthklab())/(MASSD/2.);  //vec{pi}=-vec{pr} in PW!
  double pt=kin.GetPklab()*sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  double cosdelta=kin.GetWlab()/kin.GetKlab();
  double sindelta2=kin.GetQsquared()/(kin.GetKlab()*kin.GetKlab());
  double piq=(Einoff*kin.GetWlab()+kin.GetKlab()*kin.GetPklab()*kin.GetCosthklab());  //vec{pi}=-vec{pr} in PW!
  double nutilde=piq/massi;
  double wstar2 = Einoff*Einoff-kin.GetPklab()*kin.GetPklab(); //effective mass off-shell nucleon
  double xtilde=kin.GetQsquared()/(2*piq);
  double nuoffshell=(wstar2-massi*massi+2.*piq)/(2.*massi); //(m_i+qoffshell)^2=(p_i+q)^2
  double xoffshell=kin.GetQsquared()/(2.*massi*nuoffshell);
  double fm_sq=wstar2+2.*piq-kin.GetQsquared();//sqrt((p_i+q)^2)
  NuclStructure *strfunction;
  double F1,F2;
  if(!name.compare("SLAC")) strfunction=new NuclStructure(proton,kin.GetQsquared(),xoffshell,fm_sq,name);
  else if(!name.compare("CB")){
    if(fm_sq<25.E06) strfunction=new NuclStructure(proton,kin.GetQsquared(),xoffshell,fm_sq,name);
    else strfunction=new NuclStructure(proton,kin.GetQsquared(),xoffshell,fm_sq,"SLAC");
  }
  else if(!name.compare("Alekhin")) strfunction=new NuclStructure(proton,kin.GetQsquared(),xtilde,fm_sq,name);
  if(!(strfunction->getName().compare("SLAC"))){
    F2=strfunction->getF2();
    double R=0.18;
    F1=F2*2.*xtilde/(1+R)*(pow(alphai/alphaq+1/(2.*xtilde),2.)-pt*pt/(2.*kin.GetQsquared())*R);
  }
  else strfunction->getF(F1,F2);
  //cout << F1 << " " << F2 << " " << nutilde << " " << pt << " " << kin.GetCosthklab() << endl;
  FL=pow((alphai+alphaq*piq/kin.GetQsquared())*(1+cosdelta),2.)*kin.GetWlab()/nutilde*F2-kin.GetWlab()/massi*sindelta2*F1;
  FT=2.*F1+pt*pt/(massi*nutilde)*F2;
  FTT=kin.GetWlab()*pt*pt*sindelta2/(nutilde*massi*massi*2.)*F2;
  FTL=2.*(1.+cosdelta)*pt*kin.GetWlab()/(massi*nutilde)*(alphai+alphaq*piq/kin.GetQsquared())*F2;
  /*cout << F1 << " " << F2 << endl;
  cout << FL << " " << FT << " " << FTT << " " << FTL << endl;
  cout << (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT) << " " << F2*(pow(alphai/alphaq+1/(2.*xtilde),2.)+pt*pt/(2.*kin.GetQsquared()))*2.*xtilde*kin.GetWlab()/massi*(sindelta2+2.*electron.GetTan2HalfAngle(kin)/(1+0.18)) << endl;
  */delete strfunction;
  return;
}

double DeuteronStructure::getStructure(TKinematics2to2 &kin,double phir, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions(kin,FL, FT, FTT, FTL, Einoff);
  return (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT)
	+sqrt(kin.GetQsquared()/(kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*cos(phir)*FTL+cos(2.*phir)*FTT;
  
}


double DeuteronStructure::getavgStructure(TKinematics2to2 &kin, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions(kin,FL, FT, FTT, FTL, Einoff);
  return (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT);
  
}


double DeuteronStructure::getInclStructure(TKinematics2to2 &kin, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions(kin,FL, FT, FTT, FTL, Einoff);
  //cout << FL << " " << FT << " " << FTT << " " << FTL << " " << kin.GetQsquared() << " " << kin.GetKlab() << " " << kin.GetWlab() << endl;
  return FL+kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())*kin.GetWlab()/massi*FT;
  
}
