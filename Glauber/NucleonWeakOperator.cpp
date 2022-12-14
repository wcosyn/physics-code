#include "NucleonWeakOperator.hpp"

using namespace std;

const FourVector<GammaStructure> NucleonWeakOperator::gamma_5mu=FourVector<GammaStructure>(GammaStructure(0.,0.,0.,0.,0.,0.,1.),
											GammaStructure(0.,0.,0.,0.,0.,0.,0.,1.),
				      GammaStructure(0.,0.,0.,0.,0.,0.,0.,0.,1.),GammaStructure(0.,0.,0.,0.,0.,0.,0.,0.,0.,1.));
const GammaStructure NucleonWeakOperator::gamma_5=GammaStructure(0.,1.);


NucleonWeakOperator::NucleonWeakOperator(){
}

NucleonWeakOperator::NucleonWeakOperator(const double q2, const bool prot, const int para, const bool charge,
					 const double M_A_in, const double r_s2_in, const double mu_s_in, const double gA_s_in):
NucleonEMOperator(q2,prot,para),charged(charge),M_A(M_A_in),r_s2(r_s2_in),mu_s(mu_s_in),gA_s(gA_s_in),isopartner(q2,!prot,para){
  if(charged) tau=Q2/4./MASSP/MASSN;
  setGA_weak();
  setF1_weak();
  setF2_weak();
  setGE_weak();
  setGM_weak();
}

NucleonWeakOperator::~NucleonWeakOperator(){
}


NucleonWeakOperator::NucleonWeakOperator(const NucleonWeakOperator &rhs): NucleonEMOperator(rhs){
  charged=rhs.charged;
  M_A=rhs.M_A;
  r_s2=rhs.r_s2;
  mu_s=rhs.mu_s;
  gA_s=rhs.gA_s;
  isopartner=rhs.isopartner;
  GA_weak=rhs.GA_weak;
  F1_weak=rhs.F1_weak;
  F2_weak=rhs.F2_weak;
  GE_weak=rhs.GE_weak;
  GM_weak=rhs.GM_weak;
  
}

NucleonWeakOperator & NucleonWeakOperator::operator=(const NucleonWeakOperator &rhs){
  if(this!=&rhs) { // avoid self-assignment
    NucleonEMOperator::operator=(rhs);
    charged=rhs.charged;
    M_A=rhs.M_A;
    r_s2=rhs.r_s2;
    mu_s=rhs.mu_s;
    gA_s=rhs.gA_s;
    isopartner=rhs.isopartner;
    GA_weak=rhs.GA_weak;
    F1_weak=rhs.F1_weak;
    F2_weak=rhs.F2_weak;
    GE_weak=rhs.GE_weak;
    GM_weak=rhs.GM_weak;
  }
  return *this;

}




FourVector<GammaStructure> NucleonWeakOperator::getCC1_weak(const FourVector<double> &pi, const FourVector<double> &pf, 
						     const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  if(medium) return getGM_weak(r,medium,nucleus)*gamma_mu
		      -getF2_weak(r,medium,nucleus)/(2.*(proton?MASSP:MASSN))*(pi+pf)*Id;
  else return getGM_weak()*gamma_mu-getF2_weak()/(2.*(proton?MASSP:MASSN))*(pi+pf)*Id;
}
FourVector<GammaStructure> NucleonWeakOperator::getCC2_weak(const FourVector<double> &q, const double r, const int medium,
						     const MeanFieldNucleusThick &nucleus) const{
  if(medium) return gamma_mu*getF1_weak(r,medium,nucleus)
    +getF2_weak(r,medium,nucleus)/(4.*(proton?MASSP:MASSN))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
  else return gamma_mu*getF1_weak()+getF2_weak()/(4.*(proton?MASSP:MASSN))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
  
}
FourVector<GammaStructure> NucleonWeakOperator::getCC3_weak(const FourVector<double> &q, const FourVector<double> &pi,
						     const FourVector<double> &pf, const double r, const int medium,
						     const MeanFieldNucleusThick &nucleus) const{
  if(medium) return getF1_weak(r,medium,nucleus)/(2.*(proton?MASSP:MASSN))*(pi+pf)*Id 
    + getGM_weak(r,medium,nucleus)/(4.*(proton?MASSP:MASSN))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q)); 
  else return getF1_weak()/(2.*(proton?MASSP:MASSN))*(pi+pf)*Id 
    + getGM_weak()/(4.*(proton?MASSP:MASSN))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q)); 
}

FourVector<GammaStructure> NucleonWeakOperator::getCC_weak(const int current, const FourVector<double>& q, 
						    const FourVector<double>& pi, const FourVector<double>& pf, 
						    const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
switch(current){
    case(1):
      return getCC1_weak(pi,pf,r,medium,nucleus);
      break;
    case(2):
      return getCC2_weak(q,r,medium,nucleus);
      break;
    case(3):
      return getCC3_weak(q,pi,pf,r,medium,nucleus);
      break;
    default:
      cerr << "Current operator not supported " << current << endl;
      exit(1);
  }
}

FourVector<GammaStructure> NucleonWeakOperator::getCC_weak(const int current, const FourVector<double>& q, 
						    const FourVector<double>& pi, const FourVector<double>& pf) const{
switch(current){
    case(1):
      return getGM_weak()*gamma_mu-getF2_weak()/(2.*(proton?MASSn:MASSn))*(pi+pf)*Id;
      break;
    case(2):
      return gamma_mu*getF1_weak()+getF2_weak()/(4.*(proton?MASSn:MASSn))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
      break;
    case(3):
      return getF1_weak()/(2.*(proton?MASSn:MASSn))*(pi+pf)*Id 
    + getGM_weak()/(4.*(proton?MASSn:MASSn))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
      break;
    default:
      cerr << "Current operator not supported " << current << endl;
      exit(1);
  }
/*  
switch(current){
    case(1):
      return getGM_weak()*gamma_mu-getF2_weak()/(2.*(proton?MASSP:MASSN))*(pi+pf)*Id;
      break;
    case(2):
      return gamma_mu*getF1_weak()+getF2_weak()/(4.*(proton?MASSP:MASSN))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
      break;
    case(3):
      return getF1_weak()/(2.*(proton?MASSP:MASSN))*(pi+pf)*Id 
    + getGM_weak()/(4.*(proton?MASSP:MASSN))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
      break;
    default:
      cerr << "Current operator not supported " << current << endl;
      exit(1);
  }*/
  
}


FourVector<GammaStructure> NucleonWeakOperator::getAxial(const FourVector<double> &q) const{
 return +getGA_weak()*gamma_5mu-getGP_weak()/(2.*sqrt(MASSP*MASSN))*gamma_5*q; //-gamma^mu*gamma5G_A - q^mu gamma5 GP
}



double NucleonWeakOperator::Get_dipole_mass(const double Q2, const double M) const
{
  return pow(1+Q2/(M*M),-2);
}

void NucleonWeakOperator::setGA_weak(){
  double gA_null=1.2695;
  if(charged){
    GA_weak=gA_null*Get_dipole_mass(Q2,M_A);
  }
  else
    GA_weak=(-gA_s+(proton? 1.:-1.)*gA_null)*0.5*Get_dipole_mass(Q2,M_A);
  //G_P should not contribute for neutral currents (prop to lepton mass)
  GP_weak=4.*MASSP*MASSN*GA_weak/(MASSPI*MASSPI+Q2)*(1.-MASSPI*MASSPI/(M_A*M_A)); 

}

void NucleonWeakOperator::setF1_weak(){
  if(charged)
    F1_weak=(proton?1.:-1.)*(getF1()-isopartner.getF1());
  else
    F1_weak=0.5*(getF1()-isopartner.getF1())-2.*SIN2W*getF1()+0.5*r_s2*INVHBARC*INVHBARC/6.*Q2*Get_dipole_mass(Q2,1300.);
}
void NucleonWeakOperator::setF2_weak(){
  if(charged)
    F2_weak=(proton?1.:-1.)*(getF2()-isopartner.getF2());
  else
    F2_weak=0.5*(getF2()-isopartner.getF2())-2.*SIN2W*getF2()-0.5*mu_s*Get_dipole_mass(Q2,1260.); 
}

void NucleonWeakOperator::setGE_weak(){
  GE_weak=F1_weak-tau*F2_weak;
}
void NucleonWeakOperator::setGM_weak(){
  GM_weak=F1_weak+F2_weak;
}

double NucleonWeakOperator::getF1_weak(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  double F1mod=(tau*getGM()*getMmod(r,medium,nucleus)+getGE()*getEmod(r,medium,nucleus))/(1+tau);
  double isoF1mod=(tau*isopartner.getGM()*isopartner.getMmod(r,medium,nucleus)
		    +isopartner.getGE()*isopartner.getEmod(r,medium,nucleus))/(1+tau);
  if(charged)
    return (proton?1.:-1.)*(F1mod-isoF1mod);
  else
    return 0.5*(F1mod-isoF1mod)-2.*SIN2W*F1mod+0.5*r_s2*INVHBARC*INVHBARC/6.*Q2*Get_dipole_mass(Q2,1300.);
  
}

double NucleonWeakOperator::getF2_weak(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  double F2mod=(getGM()*getMmod(r,medium,nucleus)-getGE()*getEmod(r,medium,nucleus))/(1+tau);
  double isoF2mod=(isopartner.getGM()*isopartner.getMmod(r,medium,nucleus)
		    -isopartner.getGE()*isopartner.getEmod(r,medium,nucleus))/(1+tau);
  if(charged)
    return (proton?1.:-1.)*(F2mod-isoF2mod);
  else
    return 0.5*(F2mod-isoF2mod)-2.*SIN2W*F2mod-0.5*mu_s*Get_dipole_mass(Q2,1260.); 
}


double NucleonWeakOperator::getGE_weak(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  return getF1_weak(r,medium,nucleus)-tau*getF2_weak(r,medium,nucleus);
}

double NucleonWeakOperator::getGM_weak(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  return getF1_weak(r,medium,nucleus)+getF2_weak(r,medium,nucleus);
}
