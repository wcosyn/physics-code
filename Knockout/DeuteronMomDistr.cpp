#include "DeuteronMomDistr.hpp"

#define MASS_N (MASSP+MASSN)*0.5

#include <cmath>
#include <Utilfunctions.hpp>

DeuteronMomDistr::DeuteronMomDistr(string name, double mass, int offshells, double sigmain, double betain, double epsilonin,
  double betaoffin, double lambdain):
sigma(sigmain/10.*INVHBARC*INVHBARC),
beta(betain*1.E-06),
epsilon(epsilonin),
betaoff(betaoffin*1.E-06),
lambda(lambdain*1.E+06),
massi(mass),
massr(massi==MASSP? MASSN:MASSP),
przprime(-9999.),
offshellset(offshells){
  wfref = TDeuteron::Wavefunction::CreateWavefunction(name);
  for(int i=0;i<=1000;i++){
    wf.AddUp(i,wfref->GetUp(i));
    wf.AddWp(i,wfref->GetWp(i));
  }
  
    
}

//for simple case
DeuteronMomDistr::DeuteronMomDistr(string name):
sigma(0.),
beta(0.),
epsilon(0.),
betaoff(0.),
lambda(0.),
massi(MASS_N),
massr(MASS_N),
przprime(-9999.),
offshellset(3){
  wfref = TDeuteron::Wavefunction::CreateWavefunction(name);
  for(int i=0;i<=2000.;i++){
    wf.AddUp(i,wfref->GetUp(i));
    wf.AddWp(i,wfref->GetWp(i));
  }
  
    
}


DeuteronMomDistr::~DeuteronMomDistr(){
  delete wfref;
}

//input in [mb],[GeV-2],[]
void DeuteronMomDistr::setScatter(double sigmain, double betain, double epsin){
 sigma=sigmain*0.1*INVHBARC*INVHBARC;
 beta=betain*1.E-06;
 epsilon=epsin;
 return; 
}


double DeuteronMomDistr::getMomDistrpw(TKinematics2to2 &kin, double phi) const{
  double pwtotal=0.;
  phi=PI/5.;
  double sintheta=sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      complex<double> wave=wf.DeuteronPState(M, -1, spinr, TVector3(kin.GetPklab()*sintheta*cos(phi),
								     kin.GetPklab()*sintheta*sin(phi),
								     kin.GetPklab()*kin.GetCosthklab()));
								    
      pwtotal+=norm(wave);
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  pwtotal*=2./3.;
  return pwtotal*MASSD/(2.*(MASSD-sqrt(kin.GetPklab()*kin.GetPklab()+kin.GetMesonMass()*kin.GetMesonMass())));
}


double DeuteronMomDistr::getMomDistrpw(TVector3 &pvec) const{
  double pwtotal=0.;
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      complex<double> wave=wf.DeuteronPState(M, -1, spinr, pvec);
      pwtotal+=norm(wave);
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  pwtotal*=2./3.;
  return pwtotal*MASSD/(2.*(MASSD-sqrt(pvec.Mag2()+massr*massr)));
}


complex<double> DeuteronMomDistr::getMomDistrpwcoh(TVector3 &pvec,TVector3 &pvec2,int M, int M2) const{
  complex<double> total=0.;
  for(int spinr=-1;spinr<=1;spinr+=2){
    for(int spins=-1;spins<=1;spins+=2){
      
    complex<double> wave= wf.DeuteronPState(M, spins, spinr, pvec)*conj(wf.DeuteronPState(M2,spins,spinr,pvec));
    //cout << M << " " << M2 << " " << spinr << " " << spins << " " << wave << endl;
    total+=wave;
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  return total*MASSD/(2.*sqrt((MASSD-sqrt(pvec.Mag2()+massr*massr))*(MASSD-sqrt(pvec2.Mag2()+massr*massr))));
}


double DeuteronMomDistr::getMomDistrfsi(TKinematics2to2 &kin, double phi){
  double fsitotal=0.;
  double sintheta=sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  double Er=sqrt(kin.GetPklab()*kin.GetPklab()+kin.GetMesonMass()*kin.GetMesonMass());
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      complex<double> wave=wf.DeuteronPState(M, -1, spinr, TVector3(kin.GetPklab()*sintheta*cos(phi),
								    kin.GetPklab()*sintheta*sin(phi),
								    kin.GetPklab()*kin.GetCosthklab()));
      complex<double> result;
      double qestimate=0.,thestimate=0.;
      rombergerN(this,&DeuteronMomDistr::totdens_qt,0.,1.E03,1,&result,PREC,3,10,&qestimate, 
				  &kin,M,spinr, Er,phi,&thestimate);
      wave+=I/(32.*PI*PI*kin.GetKlab()*sqrt(Er))*result/sqrt(MASSD/(2.*(MASSD-Er)));

      fsitotal+=norm(wave);
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  fsitotal*=2./3.;
  return fsitotal*MASSD/(2.*(MASSD-Er));
}

void DeuteronMomDistr::totdens_qt(const double qt, complex<double>* result, va_list ap){
  
  TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
  int M = va_arg(ap,int);
  int spinr = va_arg(ap,int);
  double Er = va_arg(ap,double);
  double phi = va_arg(ap,double);
  double *pthestimate = va_arg(ap,double*);
  

  rombergerN(this,&DeuteronMomDistr::totdens_qphi,0.,2*PI,1,result,PREC,3,10,pthestimate,pkin,qt,M,spinr,Er,phi);
  *result*=qt;
  
}

void DeuteronMomDistr::totdens_qphi(const double qphi, complex<double>* result, va_list ap){
  
  TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
  double qt = va_arg(ap,double);
  int M = va_arg(ap,int);
  int spinr = va_arg(ap,int);
  double Er = va_arg(ap,double);
  double phi = va_arg(ap,double);
 
  double prt=pkin->GetPklab()*sqrt(1.-pkin->GetCosthklab()*pkin->GetCosthklab()); //pr'_t
  double pt2=prt*prt+qt*qt-2.*prt*qt*cos(qphi-phi); //perp mom of initial spectator

  Wxprime2=-1.;  //invariant mass X' before rescattering
  get_przprime(pt2,Er,pkin);
  if(Wxprime2<0.) {*result= 0.; return;}
  double pprime = sqrt(pt2+przprime*przprime);
  double thetaprime=acos(przprime/pprime);
  double phiprime=atan2(prt*sin(phi)-qt*sin(qphi),prt*cos(phi)-qt*cos(qphi));
  double Erprime=sqrt(massr*massr+pprime*pprime);
  double t=(Er-Erprime)*(Er-Erprime)-(pkin->GetPklab()*pkin->GetCosthklab()-przprime)*(pkin->GetPklab()*pkin->GetCosthklab()-przprime)-qt*qt;
  double costhetaprime,sinthetaprime,cosphiprime,sinphiprime; 
  sincos(thetaprime,&sinthetaprime,&costhetaprime);
  sincos(phiprime,&sinphiprime,&cosphiprime);
  double chi=sqrt(pkin->GetS()*pkin->GetS()-2.*pkin->GetS()*(Wxprime2+massr*massr)+pow(massr*massr-Wxprime2,2.));
  double offshellness=0.;
  if(offshellset==0){
    double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*przprime;
    offshellness=exp(beta*massdiff);
  }
  if(offshellset==1) {
    double onshellm=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()+massi*massi+2.*pkin->GetWlab()*massi;
    offshellness=pow(lambda*lambda-onshellm,2.)/(pow(Wxprime2-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
  }
  if(offshellset==2) offshellness=exp((betaoff-beta)*t/2.);
  if(offshellset==3) offshellness=0.;
  if(offshellset==4) offshellness=1.;
  TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
  *result = chi*scatter(t)*(wf.DeuteronPState(M, -1, spinr, vecprime)
 	 +(offshellness*wfref->DeuteronPStateOff(M, -1, spinr, vecprime)))
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));
  
}

void DeuteronMomDistr::get_przprime(double pt2, double Er, TKinematics2to2 *pkin){
  przprime=0.;
  for(int i=0;i<50;i++){
    double f_Erprime=sqrt(massr*massr+pt2+przprime*przprime);  //guess for on-shell energy of initial spectator
    //guess for invariant mass initial X'
    double f_Wxprime2=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()
		      +pow(MASSD-f_Erprime,2.)-pt2-przprime*przprime
		      +2.*pkin->GetWlab()*(MASSD-f_Erprime)+2.*pkin->GetKlab()*przprime;
		      
    double f_przprime=pkin->GetPklab()*pkin->GetCosthklab()-(pkin->GetWlab()+MASSD)/pkin->GetKlab()*(Er-f_Erprime);  //new value for przprime 
    //add mass diff term if necessary
    if((pkin->GetHyperonMass()*pkin->GetHyperonMass())>f_Wxprime2)
      f_przprime-=((pkin->GetHyperonMass()*pkin->GetHyperonMass())-f_Wxprime2)/(2.*pkin->GetKlab());
    //check convergence
    if((abs((przprime-f_przprime)/f_przprime)<1e-03)) {Wxprime2= f_Wxprime2; przprime = f_przprime; return;}
    //start again
    przprime=f_przprime;
    Wxprime2=f_Wxprime2;  
  }
  //no real convergence
  return;
}

complex<double> DeuteronMomDistr::scatter(double t){
  return sigma*(I+epsilon)*exp(beta*t/2.); 
}

//all in MeV!
double DeuteronMomDistr::getMomDistrfsi(TVector3 &pvec, double nu, double qvec, double s, double massother){
  double fsitotal=0.;
  double Er=sqrt(pvec.Mag2()+massr*massr);
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      complex<double> wave=wf.DeuteronPState(M, -1, spinr, pvec);
      complex<double> result;
      double qestimate=0.,thestimate=0.;
      rombergerN(this,&DeuteronMomDistr::totdens_qt_simple,0.,1.E03,1,&result,PREC,3,7,&qestimate, 
				  &pvec, M,spinr, Er,nu,qvec,s,massother, &thestimate);
      wave+=I/(32.*PI*PI*qvec*sqrt(Er))*result/sqrt(MASSD/(2.*(MASSD-Er)));

      fsitotal+=norm(wave);
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  fsitotal*=2./3.;
  return fsitotal*MASSD/(2.*(MASSD-Er));
}


void DeuteronMomDistr::totdens_qt_simple(const double qt, complex<double>* result, va_list ap){
  
  TVector3 *p_pvec = va_arg(ap,TVector3*);
  int M = va_arg(ap,int);
  int spinr = va_arg(ap,int);
  double Er = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double s = va_arg(ap,double);
  double massother = va_arg(ap,double);
  double *pthestimate = va_arg(ap,double*);
  

  rombergerN(this,&DeuteronMomDistr::totdens_qphi_simple,0.,2*PI,1,result,PREC,3,7,
	      pthestimate,p_pvec,qt,M,spinr,Er,nu,qvec,massother);
  *result*=qt;
  
}

void DeuteronMomDistr::totdens_qphi_simple(const double qphi, complex<double>* result, va_list ap){
  
  TVector3 *p_pvec = va_arg(ap,TVector3*);
  double qt = va_arg(ap,double);
  int M = va_arg(ap,int);
  int spinr = va_arg(ap,int);
  double Er = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double s = va_arg(ap,double);
  double massother = va_arg(ap,double);
  
  double prt=p_pvec->Mag()*sqrt(1.-p_pvec->CosTheta()*p_pvec->CosTheta()); //pr'_t
  double pt2=prt*prt+qt*qt-2.*prt*qt*cos(qphi-p_pvec->Phi()); //perp mom of initial spectator

  przprime = p_pvec->Z()-(nu+MASSD)/qvec*(Er-sqrt(massr*massr+pt2));
  double pprime = sqrt(pt2+przprime*przprime);
  double thetaprime=acos(przprime/pprime);
  double phiprime=atan2(prt*sin(p_pvec->Phi())-qt*sin(qphi),prt*cos(p_pvec->Phi())-qt*cos(qphi));
  double Erprime=sqrt(massr*massr+pprime*pprime);
  if(Erprime>MASSD){
    *result=0.;
    return;
  }
  double t=(Er-Erprime)*(Er-Erprime)-(p_pvec->Z()-przprime)*(p_pvec->Z()-przprime)-qt*qt;
  double costhetaprime,sinthetaprime,cosphiprime,sinphiprime; 
  sincos(thetaprime,&sinthetaprime,&costhetaprime);
  sincos(phiprime,&sinphiprime,&cosphiprime);
  
  double chi=sqrt(s*s-2.*s*(massother*massother+massr*massr)+pow(massr*massr-massother*massother,2.));
  TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
  *result = chi*scatter(t)*(wf.DeuteronPState(M, -1, spinr, vecprime))
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));
  
}

