#include "DeuteronStructure.hpp"
#include "NuclStructure.hpp"

#include <cmath>
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <cstdlib>

using namespace std;


DeuteronStructure::DeuteronStructure(const bool pr, string nm)
:proton(pr),name(nm),massi(pr?MASSP:MASSN){
  
}


DeuteronStructure::DeuteronStructure(){
}

DeuteronStructure::DeuteronStructure(const DeuteronStructure& rhs){
  proton=rhs.proton;
  name=rhs.name;
  massi=rhs.massi;
}

DeuteronStructure& DeuteronStructure::operator=(const DeuteronStructure& rhs){
  if(this!=&rhs) { // avoid self-assignment
    proton=rhs.proton;
    name=rhs.name;
    massi=rhs.massi;
  }
  return *this;

  
}



void DeuteronStructure::getStructureFunctions(TKinematics2to2 &kin, double &FL, double &FT, double &FTT, double &FTL,double Einoff) const{
  double alphaq=(kin.GetWlab()-kin.GetKlab())/(MASSD/2.);
  double alphai=(Einoff+kin.GetPklab()*kin.GetCosthklab())/(MASSD/2.);  //vec{pi}=-vec{pr} in PW!
  double pt=kin.GetPklab()*sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  if(std::isnan(pt)) pt=0.;
  double cosdelta=kin.GetWlab()/kin.GetKlab();
  double sindelta2=kin.GetQsquared()/(kin.GetKlab()*kin.GetKlab());
  double piq=(Einoff*kin.GetWlab()+kin.GetKlab()*kin.GetPklab()*kin.GetCosthklab());  //vec{pi}=-vec{pr} in PW!
  double nutilde=piq/massi;
  double mi_off = Einoff*Einoff-kin.GetPklab()*kin.GetPklab(); //effective mass off-shell nucleon SQUARED
  double xtilde=kin.GetQsquared()/(2*piq);
  double nuoffshell=(mi_off-massi*massi+2.*piq)/(2.*massi); //(m_i+qoffshell)^2=(p_i+q)^2
  double xoffshell=kin.GetQsquared()/(2.*massi*nuoffshell); //xoffshell consistent with Q^2,m_i,W
  double W_sq=mi_off+2.*piq-kin.GetQsquared();//W^2=(p_i+q)^2
  NuclStructure strfunction(proton,kin.GetQsquared(),xoffshell,W_sq,name);
  double F1,F2;
  if(!(strfunction.getName().compare("SLAC"))){
    F2=strfunction.getF2();
    double R=0.18;
    F1=F2*2.*xtilde/(1+R)*(pow(alphai/alphaq+1/(2.*xtilde),2.)-pt*pt/(2.*kin.GetQsquared())*R);
  }
  else strfunction.getF(F1,F2);
  if(std::isnan(F2)){cout << "bla " << endl; FL=FT=FTT=FTL=F2; return;}
  if(std::isnan(F1)){FL=FT=FTT=FTL=F1; return;}
  //cout << F1 << " " << F2 << " " << nutilde << " " << pt << " " << kin.GetCosthklab() << endl;
  FL=pow((alphai+alphaq*piq/kin.GetQsquared())*(1+cosdelta),2.)*kin.GetWlab()/nutilde*F2-kin.GetWlab()/massi*sindelta2*F1;
  FT=2.*F1+pt*pt/(massi*nutilde)*F2;
  if(std::isnan(FT)) cout << F1 << " " << F2 << " " << pt << " " << piq << " " << kin.GetCosthklab() << endl;
  FTT=kin.GetWlab()*pt*pt*sindelta2/(nutilde*massi*massi*2.)*F2;
  FTL=2.*(1.+cosdelta)*pt*kin.GetWlab()/(massi*nutilde)*(alphai+alphaq*piq/kin.GetQsquared())*F2;
  /*cout << F1 << " " << F2 << endl;
  cout << FL << " " << FT << " " << FTT << " " << FTL << endl;
  cout << (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT) << " " << F2*(pow(alphai/alphaq+1/(2.*xtilde),2.)+pt*pt/(2.*kin.GetQsquared()))*2.*xtilde*kin.GetWlab()/massi*(sindelta2+2.*electron.GetTan2HalfAngle(kin)/(1+0.18)) << endl;
  */
  return;
}

void DeuteronStructure::getStructureFunctions(TKinematics2to2 &kin, double &FL, double &FT, double &FTT, 
					      double &FTL, double &F2, double Einoff) const{
  double alphaq=(kin.GetWlab()-kin.GetKlab())/(MASSD/2.);
  double alphai=(Einoff+kin.GetPklab()*kin.GetCosthklab())/(MASSD/2.);  //vec{pi}=-vec{pr} in PW!
  double pt=kin.GetPklab()*sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  if(std::isnan(pt)) pt=0.;
  double cosdelta=kin.GetWlab()/kin.GetKlab();
  double sindelta2=kin.GetQsquared()/(kin.GetKlab()*kin.GetKlab());
  double piq=(Einoff*kin.GetWlab()+kin.GetKlab()*kin.GetPklab()*kin.GetCosthklab());  //vec{pi}=-vec{pr} in PW!
  double nutilde=piq/massi;
  double mi_off = Einoff*Einoff-kin.GetPklab()*kin.GetPklab(); //effective mass off-shell nucleon
  double xtilde=kin.GetQsquared()/(2*piq);
  double nuoffshell=(mi_off-massi*massi+2.*piq)/(2.*massi); //(m_i+qoffshell)^2=(p_i+q)^2
  double xoffshell=kin.GetQsquared()/(2.*massi*nuoffshell); //xoffshell consistent with Q^2,m_i,W
  double W_sq=mi_off+2.*piq-kin.GetQsquared();//W^2=(p_i+q)^2
  NuclStructure strfunction(proton,kin.GetQsquared(),xoffshell,W_sq,name);
  double F1;
  if(!(strfunction.getName().compare("SLAC"))){
    F2=strfunction.getF2();
    double R=0.18;
    F1=F2*2.*xtilde/(1+R)*(pow(alphai/alphaq+1/(2.*xtilde),2.)-pt*pt/(2.*kin.GetQsquared())*R);
  }
  else strfunction.getF(F1,F2);
  if(std::isnan(F2)){cout << "bla " << endl; FL=FT=FTT=FTL=F2; return;}
  if(std::isnan(F1)){FL=FT=FTT=FTL=F1; return;}
  //cout << F1 << " " << F2 << " " << nutilde << " " << pt << " " << kin.GetCosthklab() << endl;
  FL=pow((alphai+alphaq*piq/kin.GetQsquared())*(1+cosdelta),2.)*kin.GetWlab()/nutilde*F2-kin.GetWlab()/massi*sindelta2*F1;
  FT=2.*F1+pt*pt/(massi*nutilde)*F2;
  if(std::isnan(FT)) cout << F1 << " " << F2 << " " << pt << " " << piq << " " << kin.GetCosthklab() << endl;
  FTT=kin.GetWlab()*pt*pt*sindelta2/(nutilde*massi*massi*2.)*F2;
  FTL=2.*(1.+cosdelta)*pt*kin.GetWlab()/(massi*nutilde)*(alphai+alphaq*piq/kin.GetQsquared())*F2;
//   cout << "VNA " << F1 << " " << F2 << endl;
//   cout << FL << " " << FT << " " << FTT << " " << FTL << endl;
//   cout << (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT) << " " << F2*(pow(alphai/alphaq+1/(2.*xtilde),2.)+pt*pt/(2.*kin.GetQsquared()))*2.*xtilde*kin.GetWlab()/massi*(sindelta2+2.*electron.GetTan2HalfAngle(kin)/(1+0.18)) << endl;
  
  return;
}


void DeuteronStructure::getStructureFunctions_off (TKinematics2to2 &kin, double &FL, double &FT, double &FTT, double &FTL,
						   double Wsq, double Einoff) const{
  double alphaq=(kin.GetWlab()-kin.GetKlab())/(MASSD/2.);
  double alphai=(Einoff+kin.GetPklab()*kin.GetCosthklab())/(MASSD/2.);  //vec{pi}=-vec{pr} in PW!
  double pt=kin.GetPklab()*sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  double cosdelta=kin.GetWlab()/kin.GetKlab();
  double sindelta2=kin.GetQsquared()/(kin.GetKlab()*kin.GetKlab());
  double piq=(Einoff*kin.GetWlab()+kin.GetKlab()*kin.GetPklab()*kin.GetCosthklab());  //vec{pi}=-vec{pr} in PW!
  double nutilde=piq/massi;
  double mi_off = Einoff*Einoff-kin.GetPklab()*kin.GetPklab(); //effective mass off-shell nucleon
  double xtilde=kin.GetQsquared()/(2*piq);
  NuclStructure strfunction(proton,kin.GetQsquared(),Wsq,1,name);
  double F1,F2;
  if(!(strfunction.getName().compare("SLAC"))){
    F2=strfunction.getF2();
    double R=0.18;
    F1=F2*2.*xtilde/(1+R)*(pow(alphai/alphaq+1/(2.*xtilde),2.)-pt*pt/(2.*kin.GetQsquared())*R);
  }
  else strfunction.getF(F1,F2);
  if(std::isnan(F2)){FL=FT=FTT=FTL=F2; return;}
  if(std::isnan(F1)){FL=FT=FTT=FTL=F1; return;}
  //cout << F1 << " " << F2 << " " << nutilde << " " << pt << " " << kin.GetCosthklab() << endl;
  FL=pow((alphai+alphaq*piq/kin.GetQsquared())*(1+cosdelta),2.)*kin.GetWlab()/nutilde*F2-kin.GetWlab()/massi*sindelta2*F1;
  FT=2.*F1+pt*pt/(massi*nutilde)*F2;
  FTT=kin.GetWlab()*pt*pt*sindelta2/(nutilde*massi*massi*2.)*F2;
  FTL=2.*(1.+cosdelta)*pt*kin.GetWlab()/(massi*nutilde)*(alphai+alphaq*piq/kin.GetQsquared())*F2;
  /*cout << F1 << " " << F2 << endl;
  cout << FL << " " << FT << " " << FTT << " " << FTL << endl;
  cout << (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT) << " " << F2*(pow(alphai/alphaq+1/(2.*xtilde),2.)+pt*pt/(2.*kin.GetQsquared()))*2.*xtilde*kin.GetWlab()/massi*(sindelta2+2.*electron.GetTan2HalfAngle(kin)/(1+0.18)) << endl;
  */
  return;
}

double DeuteronStructure::getStructure(TKinematics2to2 &kin, TElectronKinematics &electron, double phir, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions(kin,FL, FT, FTT, FTL, Einoff);
  if(std::isnan(FL)) return FL;
  return (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT)
	+sqrt(kin.GetQsquared()/(kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*cos(phir)*FTL+cos(2.*phir)*FTT;
  
}

void DeuteronStructure::getResponses(TKinematics2to2 &kin, TElectronKinematics &electron, 
				     double &ResLandT, double &ResTT, double &ResTL, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions(kin,FL, FT, FTT, FTL, Einoff);
  if(std::isnan(FL)) ResLandT=ResTT=ResTL=FL;
  else{ 
    ResLandT=(FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT);
    ResTL=sqrt(kin.GetQsquared()/(kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*FTL;
    ResTT=FTT;
  }
}

void DeuteronStructure::getResponses_off(TKinematics2to2 &kin, TElectronKinematics &electron, double Wsq,
				     double &ResLandT, double &ResTT, double &ResTL, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions_off(kin,FL, FT, FTT, FTL, Wsq, Einoff);
  if(std::isnan(FL)) ResLandT=ResTT=ResTL=FL;
  else{ 
    ResLandT=(FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT);
    ResTL=sqrt(kin.GetQsquared()/(kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*FTL;
    ResTT=FTT;
  }
}


double DeuteronStructure::getavgStructure(TKinematics2to2 &kin, TElectronKinematics &electron, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions(kin,FL, FT, FTT, FTL, Einoff);
  if(std::isnan(FL)) return FL;
  return (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT);
  
}

double DeuteronStructure::getavgStructureLC(LightConeKin2to2 &kin) const{
  double Einoff=kin.getA_mu()[0]-kin.getPs_mu()[0];
  double mi_off2 = Einoff*Einoff-kin.getPs_norm2(); //effective mass off-shell nucleon SQUARED
  double double_piq=kin.getQ2()/kin.getXN();
  double nuoffshell=(mi_off2-massi*massi+double_piq)/(2.*massi); //(m_i+qoffshell)^2=(p_i+q)^2
  double xoffshell=kin.getQ2()/(2.*massi*nuoffshell); //xoffshell consistent with Q^2,m_i,W
  double W_sq=mi_off2+double_piq-kin.getQ2();//W^2=(p_i+q)^2
  double alphaq=(kin.getQ_mu()[0]+kin.getQ_mu()[3])/sqrt(2.)/kin.getPA_plus()*2.;
  double x_i=kin.getQ2()/double_piq; // x'=Q^2/2pi*q
  NuclStructure strfunction(proton,kin.getQ2(),x_i,W_sq,name);
  double F1,F2;
  if(!(strfunction.getName().compare("SLAC"))){
    F2=strfunction.getF2();
    double R=0.18;
    F1=F2*2.*kin.getXN()/(1+R)*(pow(kin.getAlpha_i()/alphaq+1/(2.*kin.getXN()),2.)-kin.getPs_perp().Mag2()/(2.*kin.getQ2())*R);
  }
  else strfunction.getF(F1,F2);
  if(std::isnan(F2)){cout << "bla " << endl;  return F2;}
//   cout << x_i << " " << xoffshell << endl;
  return 2.*kin.getXN()*F2/kin.getQ2()*(-(mi_off2+kin.getQ2()/4./pow(kin.getXN(),2.))
	+(2*kin.getEpsilon()+1.)*pow(MASSD*MASSD+(1.-2.*kin.getZs()*kin.getXN())*kin.getQ2()/(2.*kin.getXA()*kin.getXN()),2.)
	  /(MASSD*MASSD+kin.getQ2()/kin.getXA()/kin.getXA()))+2.*(1-kin.getEpsilon())*F1;
  

}

double DeuteronStructure::getStructureLC(LightConeKin2to2 &kin) const{
  double Einoff=kin.getA_mu()[0]-kin.getPs_mu()[0];
  double mi_off2 = Einoff*Einoff-kin.getPs_norm2(); //effective mass off-shell nucleon SQUARED
  double double_piq=kin.getQ2()/kin.getXN();
  double nuoffshell=(mi_off2-massi*massi+double_piq)/(2.*massi); //(m_i+qoffshell)^2=(p_i+q)^2
  double xoffshell=kin.getQ2()/(2.*massi*nuoffshell); //xoffshell consistent with Q^2,m_i,W
  double W_sq=mi_off2+double_piq-kin.getQ2();//W^2=(p_i+q)^2
  double alphaq=(kin.getQ_mu()[0]+kin.getQ_mu()[3])/sqrt(2.)/kin.getPA_plus()*2.;
  double x_i=kin.getQ2()/double_piq; // x'=Q^2/2pi*q
  NuclStructure strfunction(proton,kin.getQ2(),x_i,W_sq,name);
  double F1,F2;
  if(!(strfunction.getName().compare("SLAC"))){
    F2=strfunction.getF2();
    double R=0.18;
    F1=F2*2.*kin.getXN()/(1+R)*(pow(kin.getAlpha_i()/alphaq+1/(2.*kin.getXN()),2.)-kin.getPs_perp().Mag2()/(2.*kin.getQ2())*R);
  }
  else strfunction.getF(F1,F2);
  if(std::isnan(F2)){cout << "bla " << endl;  return F2;}
  return F1+(1.-kin.getYN()-mi_off2*pow(kin.getYN()*kin.getXN(),2.)/kin.getQ2())/(pow(kin.getYN(),2.)*kin.getXN())*F2;
  

}

double DeuteronStructure::getavgPrefactor(TKinematics2to2 &kin, TElectronKinematics &electron, double Einoff) const{
  double FL, FT, FTT, FTL,F2;
  getStructureFunctions(kin,FL, FT, FTT, FTL, F2,Einoff);
  if(std::isnan(FL)) return FL;
  return (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT)/F2;
  
}

double DeuteronStructure::getavgStructure_off(TKinematics2to2 &kin, TElectronKinematics &electron, double Wsq, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions_off(kin,FL, FT, FTT, FTL, Wsq, Einoff);
  if(std::isnan(FL)) return FL;
  //cout << FL << " " << FT << " " << FTT << " " << FTL << " " << kin.GetQsquared() << " " << kin.GetKlab() << " " << kin.GetWlab() << endl;
  return (FL+(kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())+electron.GetTan2HalfAngle(kin))*kin.GetWlab()/massi*FT);
  
}


double DeuteronStructure::getInclStructure(TKinematics2to2 &kin, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions(kin,FL, FT, FTT, FTL, Einoff);
  if(std::isnan(FL)) {cout << "blaaa  " << endl; return FL;}
  //cout << FL << " " << FT << " " << FTT << " " << FTL << " " << kin.GetQsquared() << " " << kin.GetKlab() << " " << kin.GetWlab() << endl;
  return FL+kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())*kin.GetWlab()/massi*FT;
  
}

double DeuteronStructure::getInclStructure_off(TKinematics2to2 &kin, double Wsq, double Einoff) const{
  double FL, FT, FTT, FTL;
  getStructureFunctions_off(kin,FL, FT, FTT, FTL, Wsq, Einoff);
  if(std::isnan(FL)) return FL;
  //cout << FL << " " << FT << " " << FTT << " " << FTL << " " << kin.GetQsquared() << " " << kin.GetKlab() << " " << kin.GetWlab() << endl;
  return FL+kin.GetQsquared()/(2.*kin.GetKlab()*kin.GetKlab())*kin.GetWlab()/massi*FT;
  
}
