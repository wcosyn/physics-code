#include "DeuteronMomDistr.hpp"

#define MASS_N (MASSP+MASSN)*0.5

#include <cmath>
#include <Utilfunctions.hpp>
using namespace std;

DeuteronMomDistr::DeuteronMomDistr(std::string name, double mass, int offshells, double sigmain, double betain, double epsilonin,
  double betaoffin, double lambdain, int loop_limit):
wfname(name),
sigma(sigmain/10.*INVHBARC*INVHBARC),  //to MeV^-2
beta(betain*1.E-06), // to MeV^-2
epsilon(epsilonin),
betaoff(betaoffin*1.E-06),// to MeV^-2
lambda(lambdain*1.E+06), // to MeV^2
massi(mass),
massr(massi==MASSP? MASSN:MASSP),
przprime(-9999.),
offshellset(offshells),
looplimit(loop_limit)
{
  wfref = TDeuteron::Wavefunction::CreateWavefunction(name);
  for(int i=0;i<=1000;i++){
    wf.AddUp(i,wfref->GetUp(i));
    wf.AddWp(i,wfref->GetWp(i));
  }
//   double res=0.;
//   int m=-2;
//   for(int i=1;i<1000;i++){
//     double p=i;
//     for(int j=0;j<=60;j++){
//       double costh=-1+2./60.*j;
//       for(int k=0;k<=60;k++){
//         double phi = 2.*PI/60.*k;
//         double sinphi,cosphi;
//         sincos(phi,&sinphi,&cosphi);
//         TVector3 pmom=TVector3(p*sqrt(1-costh*costh)*cosphi,p*sqrt(1-costh*costh)*sinphi,p*costh);
//         double temp=0.;
//         temp+=norm(wf.DeuteronPState(m,-1,-1,pmom));
//         temp+=norm(wf.DeuteronPState(m,-1,1,pmom));
//         temp+=norm(wf.DeuteronPState(m,1,-1,pmom));
//         temp+=norm(wf.DeuteronPState(m,1,1,pmom));
//         res+=temp*p*p;
//       }
//     }
//   }
//   cout << res*2./60.*2.*PI/60. << endl;
//   exit(1);
//   
}

//for simple case
DeuteronMomDistr::DeuteronMomDistr(std::string name):
wfname(name),
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

DeuteronMomDistr::DeuteronMomDistr(){
  
}

DeuteronMomDistr::DeuteronMomDistr(const DeuteronMomDistr& rhs){
  wfname=rhs.wfname;
  sigma=rhs.sigma;
  beta=rhs.beta;
  epsilon=rhs.epsilon;
  betaoff=rhs.betaoff;
  lambda=rhs.lambda;
  massi=rhs.massi;
  massr=rhs.massr;
  przprime=rhs.przprime;
  offshellset=rhs.przprime;
  wfref = TDeuteron::Wavefunction::CreateWavefunction(wfname);
  wf=rhs.wf;
}

DeuteronMomDistr& DeuteronMomDistr::operator=(const DeuteronMomDistr& rhs){
  if(this!=&rhs) { // avoid self-assignment
    wfname=rhs.wfname;
    sigma=rhs.sigma;
    beta=rhs.beta;
    epsilon=rhs.epsilon;
    betaoff=rhs.betaoff;
    lambda=rhs.lambda;
    massi=rhs.massi;
    massr=rhs.massr;
    przprime=rhs.przprime;
    offshellset=rhs.przprime;
    wfref = TDeuteron::Wavefunction::CreateWavefunction(wfname);
    wf=rhs.wf;
  }
  return *this;
  
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
  return (pow(wf.GetUp(k),2.)+pow(wf.GetWp(k),2.))/(4.*PI)*sqrt(massr*massr+k*k)/pow(2.-alpha,2.)/kin.GetEklab(); 
}
// 
double DeuteronMomDistr::getLCMomDistrpw(TVector3 p) const{
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
  return (pow(wf.GetUp(k),2.)+pow(wf.GetWp(k),2.))/(4.*PI)*sqrt(massr*massr+k*k)/pow(2.-alpha,2.)/E; 
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

double DeuteronMomDistr::getMomDistrpwLC(LightConeKin2to2 &kin) const{
  //kaon translates to spectator nucleon
  double pwtotal=0.;
  return (pow(wf.GetUp(kin.getK()),2.)+pow(wf.GetWp(kin.getK()),2.))/(4.*PI)
	    *kin.getEk()/kin.getAlpha_s()/kin.getAlpha_i();
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
  fsitotal*=2./3.;  //factor 2 because of symmetry in spin_i
  return fsitotal*MASSD/(2.*(MASSD-Er));
}


double DeuteronMomDistr::getMomDistr_y_fsi(TKinematics2to2 &kin, double phi){
  //kaon translates to spectator nucleon
  double fsitotal=0.;
  double sintheta=sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  double Er=sqrt(kin.GetPklab()*kin.GetPklab()+kin.GetMesonMass()*kin.GetMesonMass());
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      //plane-wave part
      complex<double> wavedown=wf.DeuteronPState(M, -1, spinr, TVector3(kin.GetPklab()*sintheta*cos(phi),
								    kin.GetPklab()*sintheta*sin(phi),
								    kin.GetPklab()*kin.GetCosthklab()));
      complex<double> waveup=wf.DeuteronPState(M, 1, spinr, TVector3(kin.GetPklab()*sintheta*cos(phi),
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
      numint::vector_z ret(2,0.);
      F.f=DeuteronMomDistr::FSI_y_int;
      int res=90;
      unsigned count=0;
    //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E03,ret,count,0);
      result=ret[0];
//       rombergerN(this,&DeuteronMomDistr::totdens_qt,0.,1.E03,1,&result,PREC,3,10,&qestimate, 
// 				  &kin,M,spinr, Er,phi,&thestimate);
      wavedown+=I_UNIT/(32.*PI*PI*kin.GetKlab()*sqrt(Er))*ret[0]/sqrt(MASSD/(2.*(MASSD-Er)));
      waveup+=I_UNIT/(32.*PI*PI*kin.GetKlab()*sqrt(Er))*ret[0]/sqrt(MASSD/(2.*(MASSD-Er)));

      fsitotal+=2.*imag(conj(waveup)*wavedown);
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  fsitotal*=1./3.; 
  return fsitotal*MASSD/(2.*(MASSD-Er));
}


double DeuteronMomDistr::getMomDistrfsiLC(LightConeKin2to2 &kin){
  //kaon translates to spectator nucleon
  double fsitotal=0.;
  for(int M=-2;M<=2;M+=2){
    for(int spinr=-1;spinr<=1;spinr+=2){
      //plane-wave part
      complex<double> wave=wf.DeuteronPState(M, -1, spinr, kin.getKvec())*sqrt(kin.getEk());
      complex<double> result;

      numint::array<double,2> lower = {{0.,0.}};
      numint::array<double,2> upper = {{1.E03,2.*PI}};
      DeuteronMomDistr::Ftor_FSILC F;
      F.momdistr=this;
      F.kin=&kin;
      F.M = M;
      F.spinr=spinr;
      numint::mdfunction<numint::vector_z,2> mdf;
      mdf.func = &Ftor_FSILC::exec;
      mdf.param = &F;
      numint::vector_z ret(1,0.);
      F.f=DeuteronMomDistr::FSILC_int;
      int res=90;
      unsigned count=0;
    //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E03,ret,count,0);
      result=ret[0];
//       rombergerN(this,&DeuteronMomDistr::totdens_qt,0.,1.E03,1,&result,PREC,3,10,&qestimate, 
// 				  &kin,M,spinr, Er,phi,&thestimate);
      wave+=I_UNIT*kin.getAlpha_i()/(16*PI*PI*kin.getPA_plus()*(kin.getA_mu()[0]+kin.getQ_mu()[0]
      -kin.getA_mu()[3]-kin.getQ_mu()[3])/sqrt(2.))*result;
      fsitotal+=norm(wave);
    }
  }
  //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
  fsitotal*=2./3.; //factor 2 because of symmetry in spin_i
  return fsitotal/kin.getAlpha_i()/kin.getAlpha_s();
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
      wave-=I_UNIT/(32.*PI*PI*kin.GetKlab()*sqrt(Er))*result/sqrt(MASSD/(2.*(MASSD-Er)));  //minus sign!!!

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
//   cout << "VNA " << pprime << " " << kin.GetPklab() << endl;
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
 	 +(offshellness*(abs(offshellness)<1.E-05?1.:momdistr.wfref->DeuteronPStateOff(M, -1, spinr, vecprime))))
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));
  
  

}

void DeuteronMomDistr::FSI_y_int(numint::vector_z & result, double qt, double qphi, DeuteronMomDistr &momdistr,
			       TKinematics2to2 &kin, int M, int spinr, double Er, double phi){
  result=numint::vector_z(2,0.);
  if(qt<1.E-03) { result[0]=result[1]=0.; return;}
  
  
  double prt=kin.GetPklab()*sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab()); //pr'_t
  double pt2=prt*prt+qt*qt-2.*prt*qt*cos(qphi-phi); //perp mom of initial spectator

  momdistr.Wxprime2=-1.;  //invariant mass X' before rescattering
  momdistr.get_przprime(pt2,Er,kin);
  if(momdistr.Wxprime2<0.) {result[0]=result[1]=0.; return;}
  double pprime = sqrt(pt2+momdistr.przprime*momdistr.przprime);
//   cout << "VNA " << pprime << " " << kin.GetPklab() << endl;
  if(pprime>1.E03) { result[0]=result[1]=0.; return;}
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
 	 +(offshellness*(abs(offshellness)<1.E-05?1.:momdistr.wfref->DeuteronPStateOff(M, -1, spinr, vecprime))))
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));
  result[1] = qt*chi*momdistr.scatter(t)*(momdistr.wf.DeuteronPState(M, 1, spinr, vecprime)
 	 +(offshellness*(abs(offshellness)<1.E-05?1.:momdistr.wfref->DeuteronPStateOff(M, 1, spinr, vecprime))))
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));

}





//integration over perp LC rescaled momentum
// void DeuteronMomDistr::FSILC_int(numint::vector_z & result, double qt, double qphi, DeuteronMomDistr &momdistr,
// 			       LightConeKin2to2 &kin, int M, int spinr){
//   result=numint::vector_z(1,0.);
//   if(qt<1.E-03) { result[0]=0.; return;}
//   
//   
//   TVector3 k_perp=kin.getK_perp(); //k_t
//   TVector3 k_perp_prime=k_perp-TVector3(-qt*cos(qphi),-qt*sin(qphi),0.); //LC rescaled perp mom of initial spectator
//   
//   momdistr.Wxprime2=-1.;  //invariant mass X' before rescattering
//   momdistr.get_k_z_prime(k_perp_prime,kin);
//   if(momdistr.Wxprime2<0.) {result[0]= 0.; return;}
//   double kprime=sqrt(k_perp_prime.Mag2()+momdistr.k_z_prime*momdistr.k_z_prime);
// //   cout << "LC " << kprime << " " << kin.getK() << endl;
//   if(kprime>1.E03) { result[0]=0.; return;}
//   double Ekprime=sqrt(momdistr.massr*momdistr.massr+kprime*kprime);
//   double alpha_s_prime=1.+momdistr.k_z_prime/Ekprime; //WATCH SIGN!!
//   TVector3 psprime_perp=k_perp_prime+alpha_s_prime/2.*kin.getPA_perp();
//   double prprime_plus=alpha_s_prime/2.*kin.getPA_plus();
//   double prprime_min=(psprime_perp.Mag2()+momdistr.massr*momdistr.massr)/2./prprime_plus;
//   double przprime=(prprime_plus-prprime_min)/sqrt(2.);
//   double Erprime=sqrt(momdistr.massr*momdistr.massr+psprime_perp.Mag2()+przprime*przprime);
//   double thetaprime=acos(momdistr.k_z_prime/kprime);
//   FourVector<double> qmu=FourVector<double>(Erprime-kin.getPs_mu()[0],psprime_perp.X()-kin.getPs_perp().X(),
// 					    psprime_perp.Y()-kin.getPs_perp().Y(),przprime-kin.getPs_mu()[3]);
//   double t=qmu*qmu;
// //   double t=-qt*qt;
//   double chi=sqrt(kin.getS()*kin.getS()-2.*kin.getS()*(momdistr.Wxprime2+momdistr.massr*momdistr.massr)
//   +pow(momdistr.massr*momdistr.massr-momdistr.Wxprime2,2.));
//   double offshellness=0.;
//   if(momdistr.offshellset==2) offshellness=exp((momdistr.betaoff-momdistr.beta)*t/2.);
//   if(momdistr.offshellset==3) offshellness=0.;
//   if(momdistr.offshellset==4) offshellness=1.;
//   TVector3 vecprime(k_perp_prime.X(),k_perp_prime.Y(),momdistr.k_z_prime);
//   result[0] = qt*chi*momdistr.scatter(t)*(momdistr.wf.DeuteronPState(M, -1, spinr, vecprime)
//  	 -(offshellness*momdistr.wfref->DeuteronPStateOff(M, -1, spinr, vecprime)))
// 	*sqrt(Ekprime);
// }


//integration over perp LC rescaled momentum
void DeuteronMomDistr::FSILC_int(numint::vector_z & result, double qt, double qphi, DeuteronMomDistr &momdistr,
			       LightConeKin2to2 &kin, int M, int spinr){
  result=numint::vector_z(1,0.);
  if(qt<1.E-03) { result[0]=0.; return;}
  
  
  TVector3 k_perp=kin.getK_perp(); //k_t
  TVector3 k_perp_prime=k_perp-TVector3(-qt*cos(qphi),-qt*sin(qphi),0.); //LC rescaled perp mom of initial spectator
  
  momdistr.Wxprime2=-1.;  //invariant mass X' before rescattering
  momdistr.get_k_z_prime(k_perp_prime,kin);
  if(momdistr.Wxprime2<0.) {result[0]= 0.; return;}
  double kprime=sqrt(k_perp_prime.Mag2()+momdistr.k_z_prime*momdistr.k_z_prime);
//   cout << "LC " << kprime << " " << kin.getK() << endl;
  if(kprime>1.E03) { result[0]=0.; return;}
  double Ekprime=sqrt(momdistr.massr*momdistr.massr+kprime*kprime);
  double alpha_s_prime=1.-momdistr.k_z_prime/Ekprime; //WATCH SIGN!!..determined k_z_prime is that of non spectator nucleon!
  TVector3 psprime_perp=k_perp_prime+alpha_s_prime/2.*kin.getPA_perp();
  double prprime_plus=alpha_s_prime/2.*kin.getPA_plus();
  double prprime_min=(psprime_perp.Mag2()+momdistr.massr*momdistr.massr)/2./prprime_plus;
  double przprime=(prprime_plus-prprime_min)/sqrt(2.);
  double Erprime=sqrt(momdistr.massr*momdistr.massr+psprime_perp.Mag2()+przprime*przprime);
  double thetaprime=acos(-momdistr.k_z_prime/kprime);
  FourVector<double> qmu=FourVector<double>(Erprime-kin.getPs_mu()[0],psprime_perp.X()-kin.getPs_perp().X(),
					    psprime_perp.Y()-kin.getPs_perp().Y(),przprime-kin.getPs_mu()[3]);
  double t=qmu*qmu;
//   double t=-qt*qt;
  double chi=sqrt(kin.getS()*kin.getS()-2.*kin.getS()*(momdistr.Wxprime2+momdistr.massr*momdistr.massr)
  +pow(momdistr.massr*momdistr.massr-momdistr.Wxprime2,2.));
  double offshellness=0.;
  if(momdistr.offshellset==2) offshellness=exp((momdistr.betaoff-momdistr.beta)*t/2.);
  if(momdistr.offshellset==3) offshellness=0.;
  if(momdistr.offshellset==4) offshellness=1.;
  TVector3 vecprime(k_perp_prime.X(),k_perp_prime.Y(),-momdistr.k_z_prime);
  result[0] = qt*chi*momdistr.scatter(t)*(momdistr.wf.DeuteronPState(M, -1, spinr, vecprime)
 	 -(offshellness*momdistr.wfref->DeuteronPStateOff(M, -1, spinr, vecprime)))
	*sqrt(Ekprime);
}




// void DeuteronMomDistr::totdens_qt(const double qt, complex<double>* result, va_list ap){
//   
//   TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
//   int M = va_arg(ap,int);
//   int spinr = va_arg(ap,int);
//   double Er = va_arg(ap,double);
//   double phi = va_arg(ap,double);
//   double *pthestimate = va_arg(ap,double*);
//   
// 
//   rombergerN(this,&DeuteronMomDistr::totdens_qphi,0.,2*PI,1,result,PREC,3,10,pthestimate,pkin,qt,M,spinr,Er,phi);
//   *result*=qt;
//   
// }
// 
// void DeuteronMomDistr::totdens_qphi(const double qphi, complex<double>* result, va_list ap){
//   
//   TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
//   double qt = va_arg(ap,double);
//   int M = va_arg(ap,int);
//   int spinr = va_arg(ap,int);
//   double Er = va_arg(ap,double);
//   double phi = va_arg(ap,double);
//  
//   double prt=pkin->GetPklab()*sqrt(1.-pkin->GetCosthklab()*pkin->GetCosthklab()); //pr'_t
//   double pt2=prt*prt+qt*qt-2.*prt*qt*cos(qphi-phi); //perp mom of initial spectator
// 
//   Wxprime2=-1.;  //invariant mass X' before rescattering
//   get_przprime(pt2,Er,*pkin);
//   if(Wxprime2<0.) {*result= 0.; return;}
//   double pprime = sqrt(pt2+przprime*przprime);
//   double thetaprime=acos(przprime/pprime);
//   double phiprime=atan2(prt*sin(phi)-qt*sin(qphi),prt*cos(phi)-qt*cos(qphi));
//   double Erprime=sqrt(massr*massr+pprime*pprime);
//   double t=(Er-Erprime)*(Er-Erprime)-(pkin->GetPklab()*pkin->GetCosthklab()-przprime)*(pkin->GetPklab()*pkin->GetCosthklab()-przprime)-qt*qt;
//   double costhetaprime,sinthetaprime,cosphiprime,sinphiprime; 
//   sincos(thetaprime,&sinthetaprime,&costhetaprime);
//   sincos(phiprime,&sinphiprime,&cosphiprime);
//   double chi=sqrt(pkin->GetS()*pkin->GetS()-2.*pkin->GetS()*(Wxprime2+massr*massr)+pow(massr*massr-Wxprime2,2.));
//   double offshellness=0.;
//   if(offshellset==0){
//     double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*przprime;
//     offshellness=exp(beta*massdiff);
//   }
//   if(offshellset==1) {
//     double onshellm=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()+massi*massi+2.*pkin->GetWlab()*massi;
//     offshellness=pow(lambda*lambda-onshellm,2.)/(pow(Wxprime2-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
//   }
//   if(offshellset==2) offshellness=exp((betaoff-beta)*t/2.);
//   if(offshellset==3) offshellness=0.;
//   if(offshellset==4) offshellness=1.;
//   TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
//   *result = chi*scatter(t)*(wf.DeuteronPState(M, -1, spinr, vecprime)
//  	 +(offshellness*wfref->DeuteronPStateOff(M, -1, spinr, vecprime)))
// 	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));
//   
// }

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
//   cout << "VNA " << Wxprime2 << endl;
  return;
}

// void DeuteronMomDistr::get_k_z_prime(TVector3 &k_perp_prime, LightConeKin2to2 & kin){
//   k_z_prime=0.;
//   przprime=0.;
//   double factor=kin.getPA_plus()*(kin.getA_mu()[0]+kin.getQ_mu()[0]
//       -kin.getA_mu()[3]-kin.getQ_mu()[3])/sqrt(2.);
//   for(int i=0;i<looplimit;i++){
//     TVector3 k_prime=TVector3(k_perp_prime.X(),k_perp_prime.Y(),k_z_prime);
//     double Ekprime=sqrt(massr*massr+k_prime.Mag2());
//     double alpha_s_prime=1.+k_z_prime/Ekprime;
//     TVector3 psprime_perp=k_perp_prime+alpha_s_prime/2.*kin.getPA_perp();
//     double prprime_plus=alpha_s_prime/2.*kin.getPA_plus();
//     double prprime_min=(psprime_perp.Mag2()+massr*massr)/2./prprime_plus;
//     double przprime=(prprime_plus-prprime_min)/sqrt(2.);
//     double f_Erprime=(prprime_plus+prprime_min)/sqrt(2.);  //guess for on-shell energy of initial spectator
//     FourVector<double> ps_prime_mu=FourVector<double>(f_Erprime,psprime_perp.X(),psprime_perp.Y(),przprime);
//     //guess for invariant mass initial X'
//     double f_Wxprime2=(kin.getQ_mu()+kin.getA_mu()-ps_prime_mu)*(kin.getQ_mu()+kin.getA_mu()-ps_prime_mu);
// 		      
//     double f_k_z_prime=Ekprime/kin.getEk()*kin.getK_z()/*+Ekprime/factor*
// 	((sqrt(2.)*kin.getPA_plus()+(kin.getA_mu()[0]+kin.getQ_mu()[0]
//       +kin.getA_mu()[3]+kin.getQ_mu()[3]))*
// 	  (kin.getPs_mu()[0]-f_Erprime-kin.getPs_mu()[3]+przprime)
// 	  -2.*(kin.getPs_perp()-psprime_perp)*(kin.getPA_perp()+kin.getQ_perp()))*/;  //new value for kzprime 
//     double f_k_z_prime_ratio=	kin.getK_z()/kin.getEk()/*+1./factor*
// 	((sqrt(2.)*kin.getPA_plus()+(kin.getA_mu()[0]+kin.getQ_mu()[0]
//       +kin.getA_mu()[3]+kin.getQ_mu()[3]))*
// 	  (kin.getPs_mu()[0]-f_Erprime-kin.getPs_mu()[3]+przprime)
// 	  -2.*(kin.getPs_perp()-psprime_perp)*(kin.getPA_perp()+kin.getQ_perp()))*/;
//     //add mass diff term if necessary
// //     cout << i << " " << kin.getK_z() << " " << f_k_z_prime << " ";
//     if((kin.getMassX()*kin.getMassX())>f_Wxprime2){
//       f_k_z_prime+=(kin.getMassX()*kin.getMassX()-f_Wxprime2)*Ekprime/factor;
//       f_k_z_prime_ratio+=(kin.getMassX()*kin.getMassX()-f_Wxprime2)/factor;
//     }
// //     cout << f_k_z_prime << endl;
//     
//     //check convergence
//     double f_k_z_prime2=sqrt((massr*massr+k_perp_prime.Mag2())/(1./f_k_z_prime_ratio/f_k_z_prime_ratio-1.))*SIGN(f_k_z_prime_ratio);
//     if(abs(f_k_z_prime_ratio)>1){ Wxprime2=-1.; k_z_prime=1.E04; return;}
//     
//     if((abs((k_z_prime-f_k_z_prime)/f_k_z_prime)<1e-03)) {Wxprime2= f_Wxprime2; k_z_prime = f_k_z_prime; return;}
//     //start again
//     k_z_prime=f_k_z_prime;
//     Wxprime2=f_Wxprime2;  
//   }
// //   no real convergence
// //   cout << endl;
// //  k_z_prime=f_k_z_prime;//sqrt((massr*massr+k_perp_prime.Mag2())/(1./pow(alphaprime-1,2.)-1.))*SIGN(alphaprime-1);
// //  Wxprime2=f_Wxprime2;
// //  cout << f_Wxprime2 << " " << f_Wxprime22 << " " << alphaprime << " " << k_z_prime << endl;
//   return;
// }

void DeuteronMomDistr::get_k_z_prime(TVector3 &k_perp_prime, LightConeKin2to2 & kin){
  k_z_prime=0.; //k_z of INITIAL struck nucleon
  przprime=0.;
  double factor=kin.getPA_plus()*(kin.getQ_mu()[0]
      -kin.getQ_mu()[3])/sqrt(2.);
  for(int i=0;i<looplimit;i++){
    TVector3 k_prime=TVector3(k_perp_prime.X(),k_perp_prime.Y(),-k_z_prime); //SPECTATOR 3-mom (LC)
    double Ekprime=sqrt(massr*massr+k_prime.Mag2());
    double alpha_s_prime=1.-k_z_prime/Ekprime;  //WATCH SIGN
    TVector3 psprime_perp=k_perp_prime+alpha_s_prime/2.*kin.getPA_perp();
    double prprime_plus=alpha_s_prime/2.*kin.getPA_plus();
    double prprime_min=(psprime_perp.Mag2()+massr*massr)/2./prprime_plus;
    double przprime=(prprime_plus-prprime_min)/sqrt(2.);
    double f_Erprime=(prprime_plus+prprime_min)/sqrt(2.);  //guess for on-shell energy of initial spectator
    FourVector<double> ps_prime_mu=FourVector<double>(f_Erprime,psprime_perp.X(),psprime_perp.Y(),przprime);
    //guess for invariant mass initial X'
    double f_Wxprime2=(kin.getQ_mu()+kin.getA_mu()-ps_prime_mu)*(kin.getQ_mu()+kin.getA_mu()-ps_prime_mu);
		      
    double f_k_z_prime=-Ekprime/kin.getEk()*kin.getK_z()/*-Ekprime*(kin.getQ_mu()[0]+kin.getQ_mu()[3])
      *(kin.getPs_mu()[0]-f_Erprime-kin.getPs_mu()[3]+przprime)/factor-Ekprime/factor*2.*kin.getA_mu()*(ps_prime_mu-kin.getPs_mu())*/;  //new value for kzprime 
    double f_k_z_prime_ratio=	-kin.getK_z()/kin.getEk();
    //add mass diff term if necessary
//     cout << i << " " << -kin.getK_z() << " " << f_k_z_prime << " ";
    if((kin.getMassX()*kin.getMassX())>f_Wxprime2){
      f_k_z_prime+=-(kin.getMassX()*kin.getMassX()-f_Wxprime2)*Ekprime/factor;
      f_k_z_prime_ratio+=-(kin.getMassX()*kin.getMassX()-f_Wxprime2)/factor;
    }
//     cout << f_k_z_prime << " " << (kin.getMassX()*kin.getMassX()-f_Wxprime2)*Ekprime/factor << endl;
    
    //check convergence
    double f_k_z_prime2=sqrt((massr*massr+k_perp_prime.Mag2())/(1./f_k_z_prime_ratio/f_k_z_prime_ratio-1.))*SIGN(f_k_z_prime_ratio);
    if(abs(f_k_z_prime_ratio)>1){ Wxprime2=-1.; k_z_prime=1.E04; return;}
    
    if((abs((k_z_prime-f_k_z_prime)/f_k_z_prime)<1e-03)) {Wxprime2= f_Wxprime2; k_z_prime = f_k_z_prime; return;}
    //start again
    k_z_prime=f_k_z_prime;
    Wxprime2=f_Wxprime2;  
  }
//   no real convergence
//   cout << endl;
//  k_z_prime=f_k_z_prime;//sqrt((massr*massr+k_perp_prime.Mag2())/(1./pow(alphaprime-1,2.)-1.))*SIGN(alphaprime-1);
//  Wxprime2=f_Wxprime2;
//  cout << f_Wxprime2 << " " << f_Wxprime22 << " " << alphaprime << " " << k_z_prime << endl;
  return;
}


complex<double> DeuteronMomDistr::scatter(double t){
  return sigma*(I_UNIT+epsilon)*exp(beta*t/2.); 
}

//all in MeV!
//quasi-elastic fsi
// double DeuteronMomDistr::getMomDistrfsi(TVector3 pvec, double nu, double qvec, double s, double massother){
//   double fsitotal=0.;
//   double Er=sqrt(pvec.Mag2()+massr*massr);
//   for(int M=-2;M<=2;M+=2){
//     for(int spinr=-1;spinr<=1;spinr+=2){
//       complex<double> wave=wf.DeuteronPState(M, -1, spinr, pvec);
//       complex<double> result;
//       double qestimate=0.,thestimate=0.;
//       rombergerN(this,&DeuteronMomDistr::totdens_qt_simple,0.,1.E03,1,&result,PREC,3,6,&qestimate, 
// 				  &pvec, M,spinr, Er,nu,qvec,s,massother, &thestimate);
//       wave+=I_UNIT/(32.*PI*PI*qvec*sqrt(Er))*result/sqrt(MASSD/(2.*(MASSD-Er)));
// 
//       fsitotal+=norm(wave);
//     }
//   }
//   //cout << kin.GetMesonMass() << " " << kin.GetPklab() << " " << kin.GetCosthklab() << endl;
//   fsitotal*=2./3.;
//   return fsitotal*MASSD/(2.*(MASSD-Er));
// }
// 
// 
// void DeuteronMomDistr::totdens_qt_simple(const double qt, complex<double>* result, va_list ap){
//   
//   TVector3 *p_pvec = va_arg(ap,TVector3*);
//   int M = va_arg(ap,int);
//   int spinr = va_arg(ap,int);
//   double Er = va_arg(ap,double);
//   double nu = va_arg(ap,double);
//   double qvec = va_arg(ap,double);
//   double s = va_arg(ap,double);
//   double massother = va_arg(ap,double);
//   double *pthestimate = va_arg(ap,double*);
//   
// 
//   rombergerN(this,&DeuteronMomDistr::totdens_qphi_simple,0.,2*PI,1,result,PREC,3,5,
// 	      pthestimate,p_pvec,qt,M,spinr,Er,nu,qvec,massother);
//   *result*=qt;
//   
// }
// 
// void DeuteronMomDistr::totdens_qphi_simple(const double qphi, complex<double>* result, va_list ap){
//   
//   TVector3 *p_pvec = va_arg(ap,TVector3*);
//   double qt = va_arg(ap,double);
//   int M = va_arg(ap,int);
//   int spinr = va_arg(ap,int);
//   double Er = va_arg(ap,double);
//   double nu = va_arg(ap,double);
//   double qvec = va_arg(ap,double);
//   double s = va_arg(ap,double);
//   double massother = va_arg(ap,double);
//   
//   double prt=p_pvec->Mag()*sqrt(1.-p_pvec->CosTheta()*p_pvec->CosTheta()); //pr'_t
//   double pt2=prt*prt+qt*qt-2.*prt*qt*cos(qphi-p_pvec->Phi()); //perp mom of initial spectator
// 
//   przprime = p_pvec->Z()-(nu+MASSD)/qvec*(Er-sqrt(massr*massr+pt2));
//   double pprime = sqrt(pt2+przprime*przprime);
//   double thetaprime= acos(przprime/pprime);
//   if(std::isnan(thetaprime)) thetaprime=SIGN(przprime);
//   double phiprime=atan2(prt*sin(p_pvec->Phi())-qt*sin(qphi),prt*cos(p_pvec->Phi())-qt*cos(qphi));
//   double Erprime=sqrt(massr*massr+pprime*pprime);
//   if(Erprime>MASSD){
//     *result=0.;
//     return;
//   }
//   double t=(Er-Erprime)*(Er-Erprime)-(p_pvec->Z()-przprime)*(p_pvec->Z()-przprime)-qt*qt;
//   double costhetaprime,sinthetaprime,cosphiprime,sinphiprime; 
//   sincos(thetaprime,&sinthetaprime,&costhetaprime);
//   sincos(phiprime,&sinphiprime,&cosphiprime);
//   double chi=sqrt(s*s-2.*s*(massother*massother+massr*massr)+pow(massr*massr-massother*massother,2.));
//   TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
//   if(std::isnan(vecprime.Mag())) cout << pt2 << " " << przprime << " " << pprime << " " << thetaprime << " " << phiprime << " " << Erprime << endl;
//   *result = chi*scatter(t)*(wf.DeuteronPState(M, -1, spinr, vecprime))
// 	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));
//   
// }
// 
