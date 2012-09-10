#include "NucleonEMOperator.hpp"
using namespace std;

const FourVector<GammaStructure> NucleonEMOperator::gamma_mu=FourVector<GammaStructure>(GammaStructure(0.,0.,1.),GammaStructure(0.,0.,0.,1.),
				      GammaStructure(0.,0.,0.,0.,1.),GammaStructure(0.,0.,0.,0.,0.,1.));
const GammaStructure NucleonEMOperator::Id=GammaStructure(1.);


NucleonEMOperator::NucleonEMOperator(){
}

NucleonEMOperator::NucleonEMOperator(double q2, bool prot, int para):
proton(prot),parametrization(para),Q2(q2){
  tau=Q2/(4.*(proton? MASSP*MASSP:MASSN*MASSN));
  GE_null=proton? 1.:1.91304;
  GM_null=proton? 2.79:-1.91304;
  setGM();
  setGE();
  setF1();
  setF2();
}

NucleonEMOperator::~NucleonEMOperator(){
}

double NucleonEMOperator::getGE(){
  return GE;
}

double NucleonEMOperator::getGM(){
  return GM;
}

double NucleonEMOperator::getF1(){
  return F1;
}

double NucleonEMOperator::getF2(){
  return F2;
}


FourVector<GammaStructure> NucleonEMOperator::getCC1(FourVector<double> pi, FourVector<double> pf){
  return getGM()*gamma_mu-getF2()/(2.*(proton?MASSP:MASSN))*(pi+pf)*Id;
}
FourVector<GammaStructure> NucleonEMOperator::getCC2(FourVector<double> q){
  return gamma_mu*getF1()+getF2()/(4.*(proton?MASSP:MASSN))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
  
}
FourVector<GammaStructure> NucleonEMOperator::getCC3(FourVector<double> q, FourVector<double> pi, FourVector<double> pf){
  return getF1()/(2.*(proton?MASSP:MASSN))*(pi+pf)*Id + getGM()/(4.*(proton?MASSP:MASSN))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q)); 
}

FourVector<GammaStructure> NucleonEMOperator::getCC(int current, FourVector<double> q, 
						    FourVector<double> pi, FourVector<double> pf){
  switch(current){
    case(1):
      return getCC1(pi,pf);
      break;
    case(2):
      return getCC2(q);
      break;
    case(3):
      return getCC3(q,pi,pf);
      break;
    default:
      cerr << "Current operator not supported " << current << endl;
      exit(1);
  }
}



void NucleonEMOperator::setGE(){
  if(parametrization==0){
    if(!proton) GE=GE_null*pow(1+Q2/.71*1.E06,-2)*tau*0.942/(1+4.61*tau);
    else{
      if(Q2<6.*1.E06) GE=GE_null/(1. + 3.253 * Q2*1.E-06 + 1.422 * pow(Q2*1.E-06,2) + 0.08582 * pow(Q2*1.E-06,3) 
	    + 0.3318 * pow(Q2*1.E-06,4) +  (-0.09371)*pow(Q2*1.E-06,5) + 
	    0.01076*pow(Q2*1.E-06,6));
      else{
	NucleonEMOperator form6(6.,proton,parametrization);
	GE=GM/form6.getGM()*form6.getGE();
      }
    }
    return;
  }
  else if(parametrization==1){
    GE=proton? Get_Gdipole(Q2)*GE_null : GE_null*Get_Gdipole(Q2)*tau/(1+5.6*tau);
    return;
  }
}

void NucleonEMOperator::setGM(){
  if(parametrization==0){
     GM=proton? GM_null/(1. + 3.104 *Q2*1.E-06 + 1.428 * pow(Q2*1.E-06,2) + 0.1112 * pow(Q2*1.E-06,3)
	       +(-0.006981)*pow(Q2*1.E-06,4) + 0.0003705*pow(Q2*1.E-06,5) + 
	       (-0.000007063)*pow(Q2*1.E-06,6)):
	GM_null/(1. + 3.043 *Q2*1.E-06 + 0.8548 * pow(Q2*1.E-06,2) +
		0.6806 * pow(Q2*1.E-06,3)
		+(-0.1287)*pow(Q2*1.E-06,4) + 0.008912*pow(Q2*1.E-06,5));
     return;
  }
  else if(parametrization==1){
    GM=proton? GM_null*Get_Gdipole(Q2) : GM_null*Get_Gdipole(Q2);
    return;
  }
}


double NucleonEMOperator::Get_Gdipole(double Q2)  
{
  return pow(1+Q2/.71*1.E06,-2);
}


void NucleonEMOperator::setF1(){
  F1=(tau*GM+GE)/(1.+tau);
}
void NucleonEMOperator::setF2(){
  F2=(GM-GE)/(1.+tau);
}

