#include "DeuteronMomDistr.hpp"

#define MASS_N (MASSP+MASSN)*0.5

#include <cmath>
#include <Utilfunctions.hpp>
using namespace std;

DeuteronMomDistr::DeuteronMomDistr(std::string name, double mass, int offshells, double sigmain, double betain, double epsilonin,
  double betaoffin, double lambdain, int loop_limit):
sigma(sigmain/10.*INVHBARC*INVHBARC),  //to MeV^-2
beta(betain*1.E-06), // to MeV^-2
epsilon(epsilonin),
betaoff(betaoffin*1.E-06),// to MeV^-2
lambda(lambdain*1.E+06), // to MeV^2
massi(mass),
massr(massi==MASSP? MASSN:MASSP),
przprime(-9999.),
offshellset(offshells),
looplimit(loop_limit){
  wfref = TDeuteron::Wavefunction::CreateWavefunction(name);
  for(int i=0;i<=1000;i++){
    wf.AddUp(i,wfref->GetUp(i));
    wf.AddWp(i,wfref->GetWp(i));
  }
  
    
}

//for simple case
DeuteronMomDistr::DeuteronMomDistr(std::string name):
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
 sigma=sigmain*0.1*INVHBARC*INVHBARC;// to MeV^-2
 beta=betain*1.E-06;// to MeV^-2
 epsilon=epsin;
 return; 
}

double DeuteronMomDistr::getLCMomDistrpw(TKinematics2to2 &kin) const{
  //kaon translates to spectator nucleon
  double pwtotal=0.;
  double alpha=2.*(kin.GetEklab()-kin.GetPklab()*kin.GetCosthklab())/MASSD; //lightcone alpha_s=(E-p_z)/M_n
  double sintheta=sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  double pt=kin.GetPklab()*sintheta;
  double k=sqrt((M_NUCL*M_NUCL+pt*pt)/(alpha*(2.-alpha))-M_NUCL*M_NUCL); //lightcone momentum rescaling
  double k_z=sqrt(k*k-pt*pt);
//   for(int M=-2;M<=2;M+=2){
//     for(int spinr=-1;spinr<=1;spinr+=2){
//       complex<double> wave=wf.DeuteronPState(M, -1, spinr, TVector3(pt,
// 								     0.,
// 								     k_z));
// 								    
//       pwtotal+=norm(wave);
//     }
//   }
//   //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
//   pwtotal*=2./3.;
  //relativistic normalization
//   cout << k << " " << kin.GetPklab() << " " << sqrt(M_NUCL*M_NUCL+k*k)/(2.-alpha)/kin.GetEklab() << " " << alpha << endl;
  return (pow(wf.GetUp(k),2.)+pow(wf.GetWp(k),2.))/(4.*PI)*sqrt(massr*massr+k*k)/(2.-alpha)/kin.GetEklab(); 
}

double DeuteronMomDistr::getLCMomDistrpw(TVector3 p) const{
  //kaon translates to spectator nucleon
  double pwtotal=0.;
  double E=sqrt(p.Mag2()+massr*massr);
  double alpha=2.*(E-p[2])/MASSD; //lightcone alpha_s=(E-p_z)/M_n
  double pt2=p[0]*p[0]+p[1]*p[1];
  double k=sqrt((M_NUCL*M_NUCL+pt2)/(alpha*(2.-alpha))-M_NUCL*M_NUCL); //lightcone momentum rescaling
  double k_z=(1-alpha)*sqrt(M_NUCL*M_NUCL+k*k);
//   for(int M=-2;M<=2;M+=2){
//     for(int spinr=-1;spinr<=1;spinr+=2){
//       complex<double> wave=wf.DeuteronPState(M, -1, spinr, TVector3(pt,
// 								     0.,
// 								     k_z));
// 								    
//       pwtotal+=norm(wave);
//     }
//   }
//   //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
//   pwtotal*=2./3.;
  //relativistic normalization
//   cout << k << " " << kin.GetPklab() << " " << sqrt(M_NUCL*M_NUCL+k*k)/(2.-alpha)/kin.GetEklab() << " " << alpha << endl;
//   cout << p.Mag() << " " << k << " " << k_z << " " << alpha << " " << (pow(wf.GetUp(k),2.)+pow(wf.GetWp(k),2.))/(4.*PI)*sqrt(massr*massr+k*k)/(2.-alpha)/E 
//   << " " << (pow(wf.GetUp(p.Mag()),2.)+pow(wf.GetWp(p.Mag()),2.))/(4.*PI)*MASSD/(2.*(MASSD-E)) << " " << MASSD/(2.*(MASSD-E)) << " " << sqrt(massr*massr+k*k)/(2.-alpha)/E << endl; 
  return (pow(wf.GetUp(k),2.)+pow(wf.GetWp(k),2.))/(4.*PI)*sqrt(massr*massr+k*k)/(2.-alpha)/E; 
}

double DeuteronMomDistr::getMomDistrpw(TKinematics2to2 &kin) const{
  //kaon translates to spectator nucleon
  double pwtotal=0.;
//   double sintheta=sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
//   for(int M=-2;M<=2;M+=2){
//     for(int spinr=-1;spinr<=1;spinr+=2){
//       complex<double> wave=wf.DeuteronPState(M, -1, spinr, TVector3(kin.GetPklab()*sintheta,
// 								     0.,
// 								     kin.GetPklab()*kin.GetCosthklab()));
// 								    
//       pwtotal+=norm(wave);
//     }
//   }
//   //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
//   pwtotal*=2./3.;
  //relativistic normalization
  return (pow(wf.GetUp(kin.GetPklab()),2.)+pow(wf.GetWp(kin.GetPklab()),2.))/(4.*PI)
	    *MASSD/(2.*(MASSD-kin.GetEklab()));
}

double DeuteronMomDistr::getAzzDistrpw(TKinematics2to2 &kin) const{
  //kaon translates to spectator nucleon
  double pwtotal=0.;
  double sintheta=sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      complex<double> wave=wf.DeuteronPState(M, -1, spinr, TVector3(kin.GetPklab()*sintheta,
								     0.,
								     kin.GetPklab()*kin.GetCosthklab()));
								    
      pwtotal+=(M==0?-2.:1.)*norm(wave);
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  pwtotal*=2./3.;
  //relativistic normalization
  return pwtotal*MASSD/(2.*(MASSD-sqrt(kin.GetPklab()*kin.GetPklab()+kin.GetMesonMass()*kin.GetMesonMass())));
}


double DeuteronMomDistr::getMomDistrpw(TVector3 pvec) const{
//   double pwtotal=0.;
//   for(int M=-2;M<=2;M+=2){
//     for(int spinr=-1;spinr<=1;spinr+=2){
//       complex<double> wave=wf.DeuteronPState(M, -1, spinr, pvec);
//       pwtotal+=norm(wave);
//     }
//   }
//   //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
//   pwtotal*=2./3.;
  return (pow(wf.GetUp(pvec.Mag()),2.)+pow(wf.GetWp(pvec.Mag()),2.))/(4.*PI)*MASSD/(2.*(MASSD-sqrt(pvec.Mag2()+massr*massr)));
}


complex<double> DeuteronMomDistr::getMomDistrpwcoh(TVector3 pvec,TVector3 pvec2,int M, int M2) const{
  complex<double> total=0.;
  for(int spinr=-1;spinr<=1;spinr+=2){
    for(int spins=-1;spins<=1;spins+=2){
      
    complex<double> wave= wf.DeuteronPState(M, spins, spinr, pvec)*conj(wf.DeuteronPState(M2,spins,spinr,pvec));
    total+=wave;
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  return total*MASSD/(2.*sqrt((MASSD-sqrt(pvec.Mag2()+massr*massr))*(MASSD-sqrt(pvec2.Mag2()+massr*massr))));
}


double DeuteronMomDistr::getMomDistrfsi(TKinematics2to2 &kin, double phi){
  //kaon translates to spectator nucleon
  double fsitotal=0.;
  double sintheta=sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  double Er=sqrt(kin.GetPklab()*kin.GetPklab()+kin.GetMesonMass()*kin.GetMesonMass());
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      //plane-wave part
      complex<double> wave=wf.DeuteronPState(M, -1, spinr, TVector3(kin.GetPklab()*sintheta*cos(phi),
								    kin.GetPklab()*sintheta*sin(phi),
								    kin.GetPklab()*kin.GetCosthklab()));
      complex<double> result;
      double qestimate=0.,thestimate=0.;

      numint::array<double,2> lower = {{0.,0.}};
      numint::array<double,2> upper = {{1.E03,2.*PI}};
      DeuteronMomDistr::Ftor_FSI F;
      F.momdistr=this;
      F.kin=&kin;
      F.M = M;
      F.spinr=spinr;
      F.Er = Er;
      F.phi=phi;
      numint::mdfunction<numint::vector_z,2> mdf;
      mdf.func = &Ftor_FSI::exec;
      mdf.param = &F;
      numint::vector_z ret(1,0.);
      F.f=DeuteronMomDistr::FSI_int;
      int res=90;
      unsigned count=0;
    //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E03,ret,count,0);
      result=ret[0];
//       rombergerN(this,&DeuteronMomDistr::totdens_qt,0.,1.E03,1,&result,PREC,3,10,&qestimate, 
// 				  &kin,M,spinr, Er,phi,&thestimate);
      wave+=I_UNIT/(32.*PI*PI*kin.GetKlab()*sqrt(Er))*result/sqrt(MASSD/(2.*(MASSD-Er)));

      fsitotal+=norm(wave);
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  fsitotal*=2./3.;
  return fsitotal*MASSD/(2.*(MASSD-Er));
}

double DeuteronMomDistr::getAzzDistrfsi(TKinematics2to2 &kin, double phi){
  //kaon translates to spectator nucleon
  double fsitotal=0.;
  double sintheta=sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  double Er=sqrt(kin.GetPklab()*kin.GetPklab()+kin.GetMesonMass()*kin.GetMesonMass());
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      //plane-wave part
      complex<double> wave=wf.DeuteronPState(M, -1, spinr, TVector3(kin.GetPklab()*sintheta*cos(phi),
								    kin.GetPklab()*sintheta*sin(phi),
								    kin.GetPklab()*kin.GetCosthklab()));
      complex<double> result;
      double qestimate=0.,thestimate=0.;

      numint::array<double,2> lower = {{0.,0.}};
      numint::array<double,2> upper = {{1.E03,2.*PI}};
      DeuteronMomDistr::Ftor_FSI F;
      F.momdistr=this;
      F.kin=&kin;
      F.M = M;
      F.spinr=spinr;
      F.Er = Er;
      F.phi=phi;
      numint::mdfunction<numint::vector_z,2> mdf;
      mdf.func = &Ftor_FSI::exec;
      mdf.param = &F;
      numint::vector_z ret(1,0.);
      F.f=DeuteronMomDistr::FSI_int;
      int res=90;
      unsigned count=0;
    //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E03,ret,count,0);
      result=ret[0];
//       rombergerN(this,&DeuteronMomDistr::totdens_qt,0.,1.E03,1,&result,PREC,3,10,&qestimate, 
// 				  &kin,M,spinr, Er,phi,&thestimate);
      wave+=I_UNIT/(32.*PI*PI*kin.GetKlab()*sqrt(Er))*result/sqrt(MASSD/(2.*(MASSD-Er)));

      fsitotal+=(M==0?-2.:1.)*norm(wave);
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  fsitotal*=2./3.;
  return fsitotal*MASSD/(2.*(MASSD-Er));
}

void DeuteronMomDistr::FSI_int(numint::vector_z & result, double qt, double qphi, DeuteronMomDistr &momdistr,
			       TKinematics2to2 &kin, int M, int spinr, double Er, double phi){
  result=numint::vector_z(1,0.);
  if(qt<1.E-03) { result[0]=0.; return;}
  
  
  double prt=kin.GetPklab()*sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab()); //pr'_t
  double pt2=prt*prt+qt*qt-2.*prt*qt*cos(qphi-phi); //perp mom of initial spectator

  momdistr.Wxprime2=-1.;  //invariant mass X' before rescattering
  momdistr.get_przprime(pt2,Er,kin);
  if(momdistr.Wxprime2<0.) {result[0]= 0.; return;}
  double pprime = sqrt(pt2+momdistr.przprime*momdistr.przprime);
  if(pprime>1.E03) { result[0]=0.; return;}
  double thetaprime=acos(momdistr.przprime/pprime);
  double phiprime=atan2(prt*sin(phi)-qt*sin(qphi),prt*cos(phi)-qt*cos(qphi));
  double Erprime=sqrt(momdistr.massr*momdistr.massr+pprime*pprime);
  double t=(Er-Erprime)*(Er-Erprime)-(kin.GetPklab()*kin.GetCosthklab()-momdistr.przprime)
    *(kin.GetPklab()*kin.GetCosthklab()-momdistr.przprime)-qt*qt;
  double costhetaprime,sinthetaprime,cosphiprime,sinphiprime; 
  sincos(thetaprime,&sinthetaprime,&costhetaprime);
  sincos(phiprime,&sinphiprime,&cosphiprime);
  double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(momdistr.Wxprime2+momdistr.massr*momdistr.massr)
  +pow(momdistr.massr*momdistr.massr-momdistr.Wxprime2,2.));
  double offshellness=0.;
  if(momdistr.offshellset==0){
    double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-momdistr.massi*momdistr.massi
      +2.*kin.GetWlab()*(MASSD-Erprime-momdistr.massi)+2.*kin.GetKlab()*momdistr.przprime;
    offshellness=exp(momdistr.beta*massdiff);
  }
  if(momdistr.offshellset==1) {
    double onshellm=kin.GetWlab()*kin.GetWlab()-kin.GetKlab()*kin.GetKlab()
      +momdistr.massi*momdistr.massi+2.*kin.GetWlab()*momdistr.massi;
    offshellness=pow(momdistr.lambda*momdistr.lambda-onshellm,2.)/(pow(momdistr.Wxprime2-onshellm,2.)
    +pow(momdistr.lambda*momdistr.lambda-onshellm,2.));
  }
  if(momdistr.offshellset==2) offshellness=exp((momdistr.betaoff-momdistr.beta)*t/2.);
  if(momdistr.offshellset==3) offshellness=0.;
  if(momdistr.offshellset==4) offshellness=1.;
  TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
  result[0] = qt*chi*momdistr.scatter(t)*(momdistr.wf.DeuteronPState(M, -1, spinr, vecprime)
 	 +(offshellness*momdistr.wfref->DeuteronPStateOff(M, -1, spinr, vecprime)))
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));
  
  

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
  get_przprime(pt2,Er,*pkin);
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

void DeuteronMomDistr::get_przprime(double pt2, double Er, TKinematics2to2 & kin){
  przprime=0.;
  for(int i=0;i<looplimit;i++){
    double f_Erprime=sqrt(massr*massr+pt2+przprime*przprime);  //guess for on-shell energy of initial spectator
    //guess for invariant mass initial X'
    double f_Wxprime2=kin.GetWlab()*kin.GetWlab()-kin.GetKlab()*kin.GetKlab()
		      +pow(MASSD-f_Erprime,2.)-pt2-przprime*przprime
		      +2.*kin.GetWlab()*(MASSD-f_Erprime)+2.*kin.GetKlab()*przprime;
		      
    double f_przprime=kin.GetPklab()*kin.GetCosthklab()-(kin.GetWlab()+MASSD)/kin.GetKlab()*(Er-f_Erprime);  //new value for przprime 
    //add mass diff term if necessary
    if((kin.GetHyperonMass()*kin.GetHyperonMass())>f_Wxprime2)
      f_przprime-=((kin.GetHyperonMass()*kin.GetHyperonMass())-f_Wxprime2)/(2.*kin.GetKlab());
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
  return sigma*(I_UNIT+epsilon)*exp(beta*t/2.); 
}

//all in MeV!
//quasi-elastic fsi
double DeuteronMomDistr::getMomDistrfsi(TVector3 pvec, double nu, double qvec, double s, double massother){
  double fsitotal=0.;
  double Er=sqrt(pvec.Mag2()+massr*massr);
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      complex<double> wave=wf.DeuteronPState(M, -1, spinr, pvec);
      complex<double> result;
      double qestimate=0.,thestimate=0.;
      rombergerN(this,&DeuteronMomDistr::totdens_qt_simple,0.,1.E03,1,&result,PREC,3,6,&qestimate, 
				  &pvec, M,spinr, Er,nu,qvec,s,massother, &thestimate);
      wave+=I_UNIT/(32.*PI*PI*qvec*sqrt(Er))*result/sqrt(MASSD/(2.*(MASSD-Er)));

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
  

  rombergerN(this,&DeuteronMomDistr::totdens_qphi_simple,0.,2*PI,1,result,PREC,3,5,
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
  double thetaprime= acos(przprime/pprime);
  if(std::isnan(thetaprime)) thetaprime=SIGN(przprime);
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
  if(std::isnan(vecprime.Mag())) cout << pt2 << " " << przprime << " " << pprime << " " << thetaprime << " " << phiprime << " " << Erprime << endl;
  *result = chi*scatter(t)*(wf.DeuteronPState(M, -1, spinr, vecprime))
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));
  
}

