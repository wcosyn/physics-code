#include "InclusiveCross.hpp"
#include <Utilfunctions.hpp>
#include <gsl/gsl_poly.h>

using namespace std;

InclusiveCross::InclusiveCross(bool proton, string strucname, string wavename, TElectronKinematics &elec, int sym, int offshell, 
		double sigmain,double betain, double epsilonin, double betaoffin, double lambdain)
:sigma(sigmain/10.*INVHBARC*INVHBARC),
beta(betain*1.E-06),
epsilon(epsilonin),
betaoff(betaoffin*1.E-06),
lambda(lambdain*1.E+06),
massi(proton? MASSP:MASSN),
massr(proton? MASSN:MASSP),
symm(sym),
prz(-9999.),
offshellset(offshell),
electron(elec),
structure(proton,strucname){

  wf = TDeuteron::Wavefunction::CreateWavefunction(wavename);
  
}


InclusiveCross::~InclusiveCross(){
 delete wf; 
}


double InclusiveCross::calc_F2Dinc(double Q2,double x){
  
  double result;
  double prestimate=0.,cosestimate=0.;
//   rombergerN(this,&InclusiveCross::int_pr,0.,1.e03,1,&result,PREC,3,7,&prestimate, x, Q2, &cosestimate);
//   return 2.*PI*2.*massi/MASSD*result;
  numint::array<double,2> lower = {{0.,-1.}};
  numint::array<double,2> upper = {{1.E03,1.}};
  InclusiveCross::Ftor_planewave F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  numint::mdfunction<numint::vector_d,2> mdf;
  mdf.func = &Ftor_planewave::exec;
  mdf.param = &F;
  numint::vector_d ret(1,0.);
  F.f=InclusiveCross::planewave_int;
  int res=90;
  unsigned count=0;
//     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
  res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,2E04,ret,count,0);
//     cout << res << " " << count << endl;
  return 2.*PI*2.*massi/MASSD*ret[0];
  
}
void InclusiveCross::planewave_int(numint::vector_d & result, double prnorm, double costheta,
				   InclusiveCross& cross, double Q2, double x){
  
  result=numint::vector_d(1,0.);
  if(prnorm<1.E-03) { result[0]=0.; return;}
  double Er=sqrt(cross.massr*cross.massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*cross.massi*x),2.));
  //Wx^2> massi^2
  double lowerlimit= -(-cross.massi*cross.massi-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*cross.massi*x))/(2.*qvec*prnorm);
  if(lowerlimit>costheta) {result[0]=0.;return;}
  if(1.<lowerlimit){ result[0]=0.; return;}
  double nu=Q2/(2.*cross.massi*x);
  double Wx=sqrt(MASSD*MASSD-Q2+cross.massr*cross.massr+2.*MASSD*(nu-Er)-2.*nu*Er+2.*qvec*prnorm*costheta);
  TKinematics2to2 kin("","",MASSD,cross.massr,Wx,"qsquared:wlab:pklab",Q2,nu,prnorm);
  double structfactor=cross.structure.getInclStructure(kin,MASSD-Er);
  if(abs(structfactor<1E-09)||isnan(structfactor)) {result[0]=0.; return; }
  //no phi integration, symmetry is already used in getInclStructure()
  result[0]=prnorm*prnorm*(pow(cross.wf->GetUp(prnorm),2.)+pow(cross.wf->GetWp(prnorm),2.))/(4.*PI)*(MASSD/(2.*(MASSD-Er)))*structfactor;
  return;
  

}

void InclusiveCross::int_pr(double prnorm, double *result, va_list ap){
  
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double *pcosestimate = va_arg(ap,double*);

  if(prnorm<1.E-03) { *result=0.; return;}
  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double lo=-1;
  double hi=1.;
  //Wx^2> massi^2
  double lowerlimit= -(-massi*massi-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*massi*x))/(2.*qvec*prnorm);
  if(lowerlimit>lo) {lo=lowerlimit;}
  if(hi<lo){ *result=0.; return;}
  rombergerN(this, &InclusiveCross::int_costheta_incl,lo+1e-04,hi-1.e-04,1,result,PREC,3,7,pcosestimate,prnorm, x,Q2);
  if(prnorm<1.e-04) *result=0.; 
  else *result*=prnorm*prnorm*(pow(wf->GetUp(prnorm),2.)+pow(wf->GetWp(prnorm),2.))/(4.*PI)*(MASSD/(2.*(MASSD-Er)));
  return;
}


void InclusiveCross::int_costheta_incl(double costheta, double *result, va_list ap){
  
  double prnorm = va_arg(ap,double);
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);

  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double nu=Q2/(2.*massi*x);
  double Wx=sqrt(MASSD*MASSD-Q2+massr*massr+2.*MASSD*(nu-Er)-2.*nu*Er+2.*qvec*prnorm*costheta);
  TKinematics2to2 kin("","",MASSD,massr,Wx,"qsquared:wlab:pklab",Q2,nu,prnorm);
  double structfactor=structure.getInclStructure(kin,MASSD-Er);
  if(structfactor==0.) {*result=0.; return; }
  //no phi integration, symmetry is already used in getInclStructure()
  *result=structfactor;
  return;
  
}



void InclusiveCross::calc_F2DincFSI(double &fsi1, double &fsi2, double Q2,double x){
  
//   double results[2];
//   double prestimate=0.,cosestimate=0.;
//   rombergerN(this,&InclusiveCross::int_pr_fsi,0.,1.e03,2,results,PREC,3,7,&prestimate, x, Q2, &cosestimate);
//   //fsi1: integration spec momentum of x1 first
//   //fsi2: integration spec momentum of x2 first
//   fsi1=results[0]*2.*PI*2.*massi/MASSD;
//   fsi2=results[1]*2.*PI*2.*massi/MASSD;
 //cout << fsi1 << " " << fsi2 << endl;
  
  numint::array<double,4> lower = {{0.,-1.,0.,0.}};
  numint::array<double,4> upper = {{1.E03,1.,1.E03,2.*PI}};
  InclusiveCross::Ftor_FSI F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  numint::mdfunction<numint::vector_d,4> mdf;
  mdf.func = &Ftor_FSI::exec;
  mdf.param = &F;
  numint::vector_d ret(2,0.);
  F.f=InclusiveCross::FSI_int;
  int res=90;
  unsigned count=0;
//     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
  res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,2E04,ret,count,0);
//   cout << res << " " << count << endl;
  fsi1= 2.*PI*2.*massi/MASSD*ret[0];
  fsi2= 2.*PI*2.*massi/MASSD*ret[1];
  return;
}

void InclusiveCross::FSI_int(numint::vector_d & result, double prnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x){
			
  result=numint::vector_d(2,0.);
  if(prnorm<1.E-03||qt<1.E-03) { result[0]=result[1]=0.; return;}
  double Er=sqrt(cross.massr*cross.massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*cross.massi*x),2.));
  //Wx^2> massi^2
  double lowerlimit= -(-cross.massi*cross.massi-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*cross.massi*x))/(2.*qvec*prnorm);
  if(lowerlimit>costheta) {result[0]=result[1]=0.;return;}
  if(1.<lowerlimit){ result[0]=result[1]=0.; return;}
  double nu=Q2/(2.*cross.massi*x);
  cross.Wxprime2=MASSD*MASSD-Q2+cross.massr*cross.massr+2.*MASSD*(nu-Er)-2.*nu*Er+2.*qvec*prnorm*costheta;
  TKinematics2to2 kin("","",MASSD,cross.massr,sqrt(cross.Wxprime2),"qsquared:wlab:pklab",Q2,nu,prnorm);
  double structfactor=cross.structure.getInclStructure(kin,MASSD-Er);
  if(abs(structfactor<1E-09)||isnan(structfactor)) {result[0]=result[1]=0.; return; }
  //no phi integration, symmetry is already used in getInclStructure()
  double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(cross.Wxprime2+cross.massr*cross.massr)+pow(cross.massr*cross.massr-cross.Wxprime2,2.));
  complex<double> wave[6];
  for(int i=0;i<6;i++){
    wave[i] = cross.wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, TVector3(prnorm*sqrt(1.-costheta*costheta),0.,
								   prnorm*costheta)); 
  }
  result[0]=result[1]=prnorm*prnorm*structfactor/(kin.GetKlab()*32.*PI*PI*3.)
		*sqrt(MASSD/(2.*(MASSD-Er)*Er))*qt;
  double prt=kin.GetPklab()*sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  double sinqphi,cosqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double pt2=prt*prt+qt*qt-2.*prt*qt*cosqphi;  
  cross.prz=kin.GetPklab()*kin.GetCosthklab();
  double cosphiprime = (prt-qt*cosqphi)/sqrt(pt2);
  double sinphiprime = -qt*sinqphi/sqrt(pt2);
  if(abs(pt2)<1.E-03){
    pt2=0.;
    cosphiprime=1.;
    sinphiprime=0.;
  }
  cross.otherWx2=-1.;
  cross.get_prz2(pt2,&kin);
  if(cross.otherWx2<0.) {result[0]=result[1]=0.;return;}
  double pprime = sqrt(pt2+cross.przprime*cross.przprime);
  double costhetaprime=cross.przprime/pprime;
  double sinthetaprime=sqrt(pt2)/pprime;
  if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
  double Erprime=sqrt(cross.massr*cross.massr+pprime*pprime);
  double xprime=kin.GetQsquared()/(2.*((MASSD-Erprime)*kin.GetWlab()+cross.przprime*kin.GetKlab()));
  if((Erprime>MASSD)||xprime>1) {result[0]=result[1]=0.; return;}
  //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
  //t = (ps-ps')^2
  double t=(Er-Erprime)*(Er-Erprime)-(cross.prz-cross.przprime)*(cross.prz-cross.przprime)-qt*qt;
  cross.sigma=(25.3*1.E-06*1.8E06+53*(sqrt(cross.Wxprime2>5.76E06?5.76E06:cross.Wxprime2)-MASSP)*1.E-03)
    /(1.E-05*1.8E06)*INVHBARC*INVHBARC;
  complex<double> wave2[6];
  TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,cross.przprime);
  for(int i=0;i<6;i++){
    wave2[i] = cross.wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime);
  }
  complex<double> temp=0.;
  for(int i=0;i<6;i++){
    temp+=wave[i]*conj(wave2[i]);
  }
  result[0]*= imag(cross.scatter(t)*2.*temp //factor 2 due to symmetry
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)))*chi;
  double chi2=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(cross.otherWx2+cross.massr*cross.massr)
    +pow(cross.massr*cross.massr-cross.otherWx2,2.));
  result[1]*=imag(cross.scatter(t)*2.*conj(temp) //factor 2 due to symmetry
	  *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)))*chi2;
  return;
  

}

void InclusiveCross::calc_F2DincFSI_off(double &fsi1, double &fsi2, double Q2,double x){
  
  double Wmax=massi;
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double nu=Q2/(2.*massi*x);
  for(int i=0;i<1000;i++){
    double W=sqrt(-Q2+MASSD*MASSD+2.*MASSD*nu+massr*massr-2.*sqrt(massr*massr+i*i)*(nu+MASSD)+2.*qvec*i);
    if(W>Wmax) Wmax=W;
  }
//   cout << Wmax <<endl;
  
  numint::array<double,4> lower = {{0.,massi,0.,0.}};
  numint::array<double,4> upper = {{1.E03,Wmax,1.E03,2.*PI}};
  InclusiveCross::Ftor_FSI F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  numint::mdfunction<numint::vector_d,4> mdf;
  mdf.func = &Ftor_FSI::exec;
  mdf.param = &F;
  numint::vector_d ret(2,0.);
  F.f=InclusiveCross::FSI_int_off;
  int res=90;
  unsigned count=0;
//     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
  res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,2E04,ret,count,0);
  fsi1= 2.*PI*2.*massi/MASSD*ret[0];
  fsi2= 2.*PI*2.*massi/MASSD*ret[1];
//   cout << res << " " << count << " " << fsi1 << " " << fsi2 << endl;
  return;
}

void InclusiveCross::FSI_int_off(numint::vector_d & result, double prt, double W, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x){
			
  result=numint::vector_d(2,0.);
  if((qt<1.E-03)||(prt<1.E-03)) { result[0]=result[1]=0.; return;}
  
  double qvec=sqrt(Q2+pow(Q2/(2.*cross.massi*x),2.));
  double nu=Q2/(2.*cross.massi*x);
  cross.Wxprime2=W*W;
  cross.get_prz_first(prt*prt,Q2,qvec,nu);
  double prnorm=sqrt(prt*prt+cross.prz*cross.prz);
  double Er=sqrt(cross.massr*cross.massr+prnorm*prnorm);
  double xprime = Q2/(2.*((MASSD-Er)*nu+cross.prz*qvec));
  if((Er>MASSD)||(xprime>1)){result[1]=result[0]=0.;return;}
  double costheta=cross.prz/prnorm;
  TKinematics2to2 kin("","",MASSD,cross.massr,W,"qsquared:wlab:pklab",Q2,nu,prnorm);
  
  double structfactor=cross.structure.getInclStructure(kin,MASSD-Er);
  if(abs(structfactor<1E-09)||isnan(structfactor)) {result[0]=result[1]=0.; return; }
  //no phi integration, symmetry is already used in getInclStructure()
  double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(cross.Wxprime2+cross.massr*cross.massr)+pow(cross.massr*cross.massr-cross.Wxprime2,2.));
  complex<double> wave[6];
  for(int i=0;i<6;i++){
    wave[i] = cross.wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, TVector3(prnorm*sqrt(1.-costheta*costheta),0.,
								  prnorm*costheta)); 
  }
  result[0]=result[1]=prt*W*structfactor/(kin.GetKlab()*kin.GetKlab()*32.*PI*PI*3.)
		*sqrt(MASSD/(2.*(MASSD-Er)*Er))*qt;
  double sinqphi,cosqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double pt2=prt*prt+qt*qt-2.*prt*qt*cosqphi;  
  double cosphiprime = (prt-qt*cosqphi)/sqrt(pt2);
  double sinphiprime = -qt*sinqphi/sqrt(pt2);
  if(abs(pt2)<1.E-03){
    pt2=0.;
    cosphiprime=1.;
    sinphiprime=0.;
  }
  cross.otherWx2=-1.;
  cross.get_prz2(pt2,&kin);
  double pprime = sqrt(pt2+cross.przprime*cross.przprime);
  double costhetaprime=cross.przprime/pprime;
  double sinthetaprime=sqrt(pt2)/pprime;
  if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
  double Erprime=sqrt(cross.massr*cross.massr+pprime*pprime);
  double Ex=MASSD+nu-Er;
  double Exprime=MASSD+nu-Erprime;
  cross.otherWx2=cross.Wxprime2;
  if(cross.otherWx2<0.) {result[0]=result[1]=0.;return;}
  xprime=kin.GetQsquared()/(2.*((MASSD-Erprime)*kin.GetWlab()+cross.przprime*kin.GetKlab()));;
  if((Erprime>MASSD)||xprime>1) {result[0]=result[1]=0.; return;}
  //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
  //t = (ps-ps')^2
  double t=(Er-Erprime)*(Er-Erprime)-(cross.prz-cross.przprime)*(cross.prz-cross.przprime)-qt*qt;
//   cout << t << " " << t2 << " " << sqrt(cross.Wxprime2) << " " << sqrt(cross.otherWx2) << " " << Ex << " " << Exprime << endl;
  
  cross.sigma=(25.3*1.E-06*kin.GetQsquared()+53*(sqrt(cross.Wxprime2>5.76E06?5.76E06:cross.Wxprime2)-MASSP)*1.E-03)
    /(1.E-05*kin.GetQsquared())*INVHBARC*INVHBARC;
  complex<double> wave2[6];
  TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,cross.przprime);
  for(int i=0;i<6;i++){
    wave2[i] = -cross.wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime);
  }
  complex<double> temp=0.;
  for(int i=0;i<6;i++){
    temp+=wave[i]*conj(wave2[i]);
  }
  double offshellness=0.;
  if(cross.offshellset==0){
    double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-cross.massi*cross.massi
      +2.*kin.GetWlab()*(MASSD-Erprime-cross.massi)+2.*kin.GetKlab()*cross.prz;
    offshellness=exp(cross.beta*massdiff);
  }
  if(cross.offshellset==1) {
    double onshellm=-kin.GetQsquared()+cross.massi*cross.massi+2.*kin.GetWlab()*cross.massi;
    offshellness=pow(cross.lambda*cross.lambda-onshellm,2.)/(pow(cross.Wxprime2-onshellm,2.)+pow(cross.lambda*cross.lambda-onshellm,2.));
  }
  if(cross.offshellset==2){
    cross.betaoff=16.8*1.E-06;
    cross.sigma=(25.3*1.E-06*kin.GetQsquared()+53*(sqrt(cross.Wxprime2>5.76E06?5.76E06:cross.Wxprime2)-MASSP)*1.E-03)/(1.E-05*kin.GetQsquared())*INVHBARC*INVHBARC;
    offshellness=0.;//exp((cross.betaoff-cross.beta)*t/2.);
  }
  if(cross.offshellset==3){
    cross.betaoff=16.8*1.E-06;
    cross.sigma=(25.3*1.E-06*kin.GetQsquared()+53*(sqrt(cross.Wxprime2>5.76E06?5.76E06:cross.Wxprime2)-MASSP)*1.E-03)/(1.E-05*kin.GetQsquared())*INVHBARC*INVHBARC;
    offshellness=exp((cross.betaoff-cross.beta)*t/2.);
  }
  if(cross.offshellset==4){
    cross.sigma=(25.3*1.E-06*kin.GetQsquared()+53*(sqrt(cross.Wxprime2>5.76E06?5.76E06:cross.Wxprime2)-MASSP)*1.E-03)/(1.E-05*kin.GetQsquared())*INVHBARC*INVHBARC;
    offshellness=1.;
  }
  
  result[0]*= imag(cross.scatter(t)*2.*temp //factor 2 due to symmetry
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)))*chi*offshellness;
  double chi2=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(cross.otherWx2+cross.massr*cross.massr)
    +pow(cross.massr*cross.massr-cross.otherWx2,2.));
  result[1]*=imag(cross.scatter(t)*2.*conj(temp) //factor 2 due to symmetry
      *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)))*chi2*offshellness;

}




void InclusiveCross::int_pr_fsi(double prnorm, double *results, va_list ap){
  
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double *pcosestimate = va_arg(ap,double*);
  
  if(prnorm<1.E-03) { results[0]=0.; results[1]=0.; return;}
  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double lo=-1;
  double hi=1.;
  //Wx^2> massi^2
  double lowerlimit= -(-massi*massi-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*massi*x))/(2.*qvec*prnorm);
  if(lowerlimit>lo) {lo=lowerlimit;}
  if(hi<lo){ results[0]=0.; results[1]=0.; return;}
  rombergerN(this, &InclusiveCross::int_costheta_incl_fsi,lo+1e-04,hi-1.e-04,2,results,PREC,3,7,pcosestimate,prnorm, x,Q2);

  results[0]*=prnorm*prnorm;
  results[1]*=prnorm*prnorm;
  return;
}


void InclusiveCross::int_costheta_incl_fsi(double costheta, double *results, va_list ap){
  
  double prnorm = va_arg(ap,double);
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);

  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double nu=Q2/(2.*massi*x);
  Wxprime2=MASSD*MASSD-Q2+massr*massr+2.*MASSD*(nu-Er)-2.*nu*Er+2.*qvec*prnorm*costheta;
  TKinematics2to2 kin("","",MASSD,massr,sqrt(Wxprime2),"qsquared:wlab:pklab",Q2,nu,prnorm);
  double structfactor=structure.getInclStructure(kin,MASSD-Er);
  if(abs(structfactor<1E-09)||isnan(structfactor)) {results[0]=0.; results[1]=0.; return; }
  double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(Wxprime2+massr*massr)+pow(massr*massr-Wxprime2,2.));
  
  complex<double> wave[6];
  complex<double> waveoff[6];
  for(int i=0;i<6;i++){
    wave[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, TVector3(prnorm*sqrt(1.-costheta*costheta),0.,
								   prnorm*costheta)); 
    waveoff[i] = wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, TVector3(prnorm*sqrt(1.-costheta*costheta),0.,
								   prnorm*costheta)); 
  }
  
  double qresults[2];
  double qestimate=0.,qphiestimate=0.;
  rombergerN(this,&InclusiveCross::int_qt,0.,1.E03,2,qresults,PREC,3,6,&qestimate,&kin,wave,waveoff,&qphiestimate);
  
  results[0]= structfactor/(kin.GetKlab()*32.*PI*PI*3.)
		*sqrt(MASSD/(2.*(MASSD-Er)*Er));
  results[1]= results[0];
  results[0]*=qresults[0]*chi;
  results[1]*=qresults[1];
  return;
  
}

void InclusiveCross::int_qt(double qt, double *results, va_list ap){

  TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
  complex<double>* wave = va_arg(ap,complex<double>*);
  complex<double>* waveoff = va_arg(ap,complex<double>*);
  double *pqphiestimate = va_arg(ap,double*);
  
  if(qt<1.E-03) {results[0]=results[1]=0.; return;}
  rombergerN(this, &InclusiveCross::int_qphi,0.,2*PI,2,results,PREC,3,6,pqphiestimate,qt,pkin,wave,waveoff);
  
  results[0]*=qt;
  results[1]*=qt;
  return;  
}


void InclusiveCross::int_qphi(double qphi, double *results, va_list ap){

  double qt = va_arg(ap,double);
  TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
  complex<double>* wave = va_arg(ap,complex<double>*);
  complex<double>* waveoff = va_arg(ap,complex<double>*);
  
  double prt=pkin->GetPklab()*sqrt(1.-pkin->GetCosthklab()*pkin->GetCosthklab());
  double sinqphi,cosqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double pt2=prt*prt+qt*qt-2.*prt*qt*cosqphi;  
  double Er=sqrt(massr*massr+pkin->GetPklab()*pkin->GetPklab());
  prz=pkin->GetPklab()*pkin->GetCosthklab();
  double cosphiprime = (prt-qt*cosqphi)/sqrt(pt2);
  double sinphiprime = -qt*sinqphi/sqrt(pt2);
  if(abs(pt2)<1.E-03){
    pt2=0.;
    cosphiprime=1.;
    sinphiprime=0.;
  }

  //spectator integration over X1
  otherWx2=-1.;
  get_prz2(pt2,pkin);
  if(otherWx2<0.) {results[0]=results[1]=0.;return;}
  else{
    double pprime = sqrt(pt2+przprime*przprime);
    double costhetaprime=przprime/pprime;
    double sinthetaprime=sqrt(pt2)/pprime;
    if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
    double Erprime=sqrt(massr*massr+pprime*pprime);
    double x=pkin->GetQsquared()/(2.*(MASSD-Erprime)*pkin->GetWlab());
    if((Erprime>MASSD)||x>1) {results[0]=results[1]=0.; return;}
    else{
      //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
      //t = (ps-ps')^2
      double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
      double offshellness=0.;
      if(offshellset==0){
	double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*prz;
	offshellness=exp(beta*massdiff);
      }
      if(offshellset==1) {
	double onshellm=-pkin->GetQsquared()+massi*massi+2.*pkin->GetWlab()*massi;
	offshellness=pow(lambda*lambda-onshellm,2.)/(pow(Wxprime2-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
      }
      if(offshellset==2){
	betaoff=16.8*1.E-06;
	sigma=(25.3*1.E-06*pkin->GetQsquared()+53*(sqrt(otherWx2>5.76E06?5.76E06:otherWx2)-MASSP)*1.E-03)/(1.E-05*pkin->GetQsquared())*INVHBARC*INVHBARC;
	offshellness=0.;//exp((betaoff-beta)*t/2.);
      }
      if(offshellset==3){
	sigma=(25.3*1.E-06*pkin->GetQsquared()+53*(sqrt(Wxprime2>5.76E06?5.76E06:Wxprime2)-MASSP)*1.E-03)/(1.E-05*pkin->GetQsquared())*INVHBARC*INVHBARC;
	offshellness=0.;
      }
      if(offshellset==4){
	sigma=(25.3*1.E-06*pkin->GetQsquared()+53*(sqrt(Wxprime2>5.76E06?5.76E06:Wxprime2)-MASSP)*1.E-03)/(1.E-05*pkin->GetQsquared())*INVHBARC*INVHBARC;
	offshellness=1.;
      }
      complex<double> wave2[6];
      complex<double> wave2off[6];
      TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,przprime);
      for(int i=0;i<6;i++){
	wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime);
	wave2off[i] = wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 			
      }
      complex<double> temp=0.;
      for(int i=0;i<6;i++){
	temp+=wave[i]*conj(wave2[i])+offshellness*waveoff[i]*conj(wave2off[i]);
      }
      //cout << przprime << " " << pkin->GetQsquared()/(2.*(MASSD-Erprime)*pkin->GetWlab()) << " " << abs(temp) << endl;
      results[0]= imag(scatter(t)*2.*temp //factor 2 due to symmetry
	    *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));
    double chi=sqrt(pkin->GetS()*pkin->GetS()-2.*pkin->GetS()*(otherWx2+massr*massr)+pow(massr*massr-otherWx2,2.));
//       cout << pprime << " " << temp << " " << results[0] << " " << results[0]*chi << " ";
    results[1]=imag(scatter(t)*2.*conj(temp) //factor 2 due to symmetry
	    *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)))*chi;
// 	    cout << results[1] << " ";
	    return;
    }
  }
  
    //spectator integration over X1
//   otherWx2=-1.;
//   get_prz2(pt2,Er,pkin,0);   
//   if(otherWx2<0.) {results[1]=0.; return;}
//   double chi=sqrt(pkin->GetS()*pkin->GetS()-2.*pkin->GetS()*(otherWx2+massr*massr)+pow(massr*massr-otherWx2,2.));
//   double pprime = sqrt(pt2+przprime*przprime);
//   double costhetaprime=przprime/pprime;
//   double sinthetaprime=sqrt(pt2)/pprime;
//   if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
//   double Erprime=sqrt(massr*massr+pprime*pprime);
//   if(Erprime>MASSD) {results[1]=0.; return;}
//   //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
//   double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
//   double offshellness=0.;
//   if(offshellset==0){
//     double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*przprime;
//     offshellness=exp(beta*massdiff);
//   }
//   if(offshellset==1) {
//     double onshellm=-pkin->GetQsquared()+massi*massi+2.*pkin->GetWlab()*massi;
//     offshellness=pow(lambda*lambda-onshellm,2.)/(pow(Wxprime2-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
//   }
//   if(offshellset==2){
//     betaoff=16.8*1.E-06;
//     sigma=(25.3*1.E-06*pkin->GetQsquared()+53*(sqrt(Wxprime2>5.76E06?5.76E06:Wxprime2)-MASSP)*1.E-03)/(1.E-05*pkin->GetQsquared())*INVHBARC*INVHBARC;
//    offshellness=0.;//exp((betaoff-beta)*t/2.);
//   }
//   if(offshellset==3){
//     sigma=(25.3*1.E-06*pkin->GetQsquared()+53*(sqrt(Wxprime2>5.76E06?5.76E06:Wxprime2)-MASSP)*1.E-03)/(1.E-05*pkin->GetQsquared())*INVHBARC*INVHBARC;
//     offshellness=0.;
//   }
//   if(offshellset==4){
//     sigma=(25.3*1.E-06*pkin->GetQsquared()+53*(sqrt(Wxprime2>5.76E06?5.76E06:Wxprime2)-MASSP)*1.E-03)/(1.E-05*pkin->GetQsquared())*INVHBARC*INVHBARC;
//     offshellness=1.;
//   }
//   complex<double> wave2[6];
//   complex<double> wave2off[6];
//   TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,przprime);
//   for(int i=0;i<6;i++){
//     wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime);
//     wave2off[i] = wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime);
// //      cout << wave2[i] << " " << wave2off[i] << " " << vecprime[0] << " " << vecprime[1] << " " << vecprime[2] <<  endl;
//   }	    
//   complex<double> temp=0.,temp1=0.,temp2=0.;
//   for(int i=0;i<6;i++) {temp+=conj(wave[i])*wave2[i]+offshellness*conj(waveoff[i])*wave2off[i];
//   temp1+=conj(wave[i])*wave2[i];
//   temp2+=conj(waveoff[i])*wave2off[i];
// // 	cout << qphi << " " << przprime << " " << (i/4)*2-2 << " " << ((i/2)%2)*2-1 << " " << (i%2)*2-1<< " " << wave[i]*conj(wave2[i]) << " " << waveoff[i]*conj(wave2off[i]) << endl;
//   }
// //   cout << abs(temp1) << " " << abs(temp2) << endl;
//   results[1]= imag(scatter(t)*2.*temp  //factor 2 due to symmetry
// 	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)))*chi;
// //   cout << pprime << " " << temp << " " << results[1] << " " << results[1]/chi << endl;
// 	cout << results[1] << endl;
  return;  
}


void InclusiveCross::get_prz(double pt2, double Er, TKinematics2to2 *pkin, int first){
  przprime=0;
  for(int i=0;i<50;i++){
    double f_Erprime=sqrt(massr*massr+pt2+przprime*przprime);  //guess for on-shell energy of initial spectator
    //guess for invariant mass initial X' Wx'2=(q+p_D-p_s')^2
    double f_Wxprime2=-pkin->GetQsquared()
		      +pow(MASSD-f_Erprime,2.)-pt2-przprime*przprime
		      +2.*pkin->GetWlab()*(MASSD-f_Erprime)+2.*pkin->GetKlab()*przprime;
		      
    double f_prz=prz-(pkin->GetWlab()+MASSD)/pkin->GetKlab()*(Er-f_Erprime);  //new value for prz 
    //add mass diff term if necessary
    if(symm==0){
      if(first/*&&(Wxprime2<f_Wxprime2)*/) f_prz-=(Wxprime2-f_Wxprime2)/(2.*pkin->GetKlab());
      if((!first)/*&&(Wxprime2>f_Wxprime2)*/) f_prz-=(Wxprime2-f_Wxprime2)/(2.*pkin->GetKlab());
    }
    if(symm==-1){
      if(first&&(Wxprime2<f_Wxprime2)) f_prz-=(Wxprime2-f_Wxprime2)/(2.*pkin->GetKlab());
      if((!first)&&(Wxprime2>f_Wxprime2)) f_prz-=(Wxprime2-f_Wxprime2)/(2.*pkin->GetKlab());
    }
    //check convergence
    if((abs((przprime-f_prz)/f_prz)<1e-03)) {otherWx2= f_Wxprime2; przprime = f_prz;/* cout << prz << " " << przprime << endl;*/ return;}
    //start again
    przprime=f_prz;
    otherWx2=f_Wxprime2;  
  }
  //no real convergence
//   cout << "bleeeep " << prz << " " << przprime << endl;
  return;
}


void InclusiveCross::get_prz2(double pt2, TKinematics2to2 *pkin){
  if(symm==1){
    przprime=0.;
    double x=0.;
    for(int i=0;i<1;i++){
      double f_Erprime=sqrt(massr*massr+pt2+przprime*przprime);  //guess for on-shell energy of initial spectator
//       x=pkin->GetQsquared()/(2.*((MASSD-f_Erprime)*pkin->GetWlab()+przprime*pkin->GetKlab()));
      x=pkin->GetQsquared()/(2.*(MASSN*pkin->GetWlab()));
    //guess for invariant mass initial X' Wx'2=(q+p_D-p_s')^2
//       double f_Wxprime2=-pkin->GetQsquared()
// 			+pow(MASSD-f_Erprime,2.)-pt2-przprime*przprime
// 			+2.*pkin->GetWlab()*(MASSD-f_Erprime)+2.*pkin->GetKlab()*przprime;
			
      double f_prz=-massi*(1-x);
      f_Erprime=sqrt(massr*massr+pt2+f_prz*f_prz);
      double f_Wxprime2=-pkin->GetQsquared()
			+pow(MASSD-f_Erprime,2.)-pt2-f_prz*f_prz
			+2.*pkin->GetWlab()*(MASSD-f_Erprime)+2.*pkin->GetKlab()*f_prz;
			
      if((abs((przprime-f_prz)/f_prz)<1e-03)) {otherWx2= f_Wxprime2; przprime = f_prz; /*cout << i << " " << prz << " " << przprime << " " << x << endl;*/ return;}
      //start again
      przprime=f_prz;
      otherWx2=f_Wxprime2;  
    }
    //no real convergence
  //   cout << "bleeeep " << prz << " " << przprime << " " << x << endl;
  //   cout << first << " " << x << " " << przprime << " " << otherWx2 << endl;
  }
  if(symm==0){
    gsl_poly_complex_workspace * space =  gsl_poly_complex_workspace_alloc(5);
    double z[8];
    double a[5];
    double C=massi*pkin->GetQsquared()/2.;
    double B=massr*massr+pt2;
    a[0]=C*C-2.*C*massi*MASSD*pkin->GetWlab()+massi*massi*(MASSD*MASSD-B)*pkin->GetWlab()*pkin->GetWlab();
    a[1]=2.*C*(pkin->GetKlab()*massi-pkin->GetWlab()*MASSD)-2.*massi*massi*MASSD*pkin->GetKlab()*pkin->GetWlab()
    +2.*massi*pkin->GetWlab()*pkin->GetWlab()*(MASSD*MASSD-B);
    a[2]=2.*C*pkin->GetKlab()+(MASSD*MASSD-B)*pkin->GetWlab()*pkin->GetWlab()-4.*massi*MASSD*pkin->GetKlab()*pkin->GetWlab()+massi*massi*pkin->GetQsquared();
    a[3]=2.*massi*pkin->GetQsquared()-2.*MASSD*massi*pkin->GetKlab()*pkin->GetWlab();
    a[4]=pkin->GetQsquared();
    int s = gsl_poly_complex_solve(a, 5, space, z);
    cout << endl << endl << s << " " << z[0] << " " << z[1] << " " << z[2] << " " << z[3] << " " << z[4] << " " << z[5] << " " << z[6] << " " <<z[7] << endl;
    for(int i=0;i<4;i++){
      if((z[2*i+1]==0)&&(abs(z[2*i])<1.E03)) {
      double Er = sqrt(massr*massr+pt2+z[2*i]*z[2*i]);
      double x=pkin->GetQsquared()/(2.*((MASSD-Er)*pkin->GetWlab()-z[2*i]*pkin->GetKlab()));
  	cout << Er << " " << x << " " << massi*(x-1) << " " << z[2*i] << endl;
	    cout << i << " " << a[4]*pow(z[2*i],4.)+a[3]*pow(z[2*i],3.)+a[2]*pow(z[2*i],2.)+a[1]*z[2*i]+a[0] << endl;

      }
    }
    
  }
  return;
}

void InclusiveCross::get_prz_first(double pt2, double Q2, double qvec, double nu){
  if(symm==1){
    prz=0.;
    double x=0.;
    for(int i=0;i<1;i++){
      double f_Er=sqrt(massr*massr+pt2+prz*prz);  //guess for on-shell energy of initial spectator
//       x=Q2/(2.*((MASSD-f_Er)*nu+prz*qvec));
      x=Q2/(2.*((MASSN*nu)));
    //guess for invariant mass initial X' Wx'2=(q+p_D-p_s')^2
      double f_Wxprime2=-Q2
			+pow(MASSD-f_Er,2.)-pt2-prz*prz
			+2.*nu*(MASSD-f_Er)+2.*qvec*prz;
			
      double f_prz=-massi*(1-x);  //new value for prz 
      f_Er=sqrt(massr*massr+pt2+f_prz*f_prz);
      f_Wxprime2=-Q2
			+pow(MASSD-f_Er,2.)-pt2-f_prz*f_prz
			+2.*nu*(MASSD-f_Er)+2.*qvec*f_prz;
      
       if((abs((prz-f_prz)/f_prz)<1e-03)) {przprime = f_prz; /*cout << i << " " << prz << " " << przprime << " " << x << endl;*/ return;}
      //start again
      prz=f_prz;
    }
    //no real convergence
  //   cout << "bleeeep " << prz << " " << przprime << " " << x << endl;
  //   cout << first << " " << x << " " << przprime << " " << otherWx2 << endl;

  }
}

void InclusiveCross::get_prz_first2(double pt2, double Q2, double qvec, double nu, double W1){
  if(symm==1){
    prz=0.;
    for(int i=0;i<50;i++){
      double f_Er=sqrt(massr*massr+pt2+prz*prz);  //guess for on-shell energy of initial spectator
      double f_prz=-(-Q2+MASSD*MASSD+2.*nu*MASSD+massr*massr-W1*W1)/(2.*qvec)+f_Er/qvec*(nu+MASSD);  //new value for prz 
       if((abs((prz-f_prz)/f_prz)<1e-03)) {przprime = f_prz; /*cout << i << " " << prz <<  endl;*/ return;}
      //start again
      prz=f_prz;
    }
    //no real convergence
  //   cout << "bleeeep " << prz << " " << przprime << " " << x << endl;
  //   cout << first << " " << x << " " << przprime << " " << otherWx2 << endl;

  }
}


complex<double> InclusiveCross::scatter(double t){
  return sigma*(I_UNIT+epsilon)*exp(beta*t/2.); 
}

/*
void InclusiveCross::calc_F2DincFSI2(double &fsi1, double &fsi2, double Q2,double x){
  
  double results[2];
  double prestimate=0.,qtestimate=0.,qtphiestimate=0.;
  rombergerN(this,&InclusiveCross::int_pperp_fsi,0.,1.e03,2,results,PREC,3,7,&prestimate, x, Q2, &qtestimate,&qtphiestimate);
  fsi1= results[0]*2.*PI*2.*massi/MASSD;
  fsi2= results[1]*2.*PI*2.*massi/MASSD;
 cout << fsi1 << " " << fsi2 << endl;
  
}


void InclusiveCross::int_pperp_fsi(double pperp, double *results, va_list ap){
  
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double *pqtestimate = va_arg(ap,double*);
  double *pqtphiestimate = va_arg(ap,double*);

  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double nu=Q2/(2.*massi*x);
  rombergerN(this,&InclusiveCross::int_qt_bis,0.,1.E03,2,results,PREC,3,6,pqtestimate,pperp,qvec,nu,pqtphiestimate);
  results[0]*= pperp/(qvec*32.*PI*PI*3.);
  results[1]*= pperp/(qvec*32.*PI*PI*3.);
// cout << results[0] << " " << results[1] << endl;
    
  return;
  
}

void InclusiveCross::int_qt_bis(double qt, double *results, va_list ap){

  double pperp = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double *pqtphiestimate = va_arg(ap,double*);
  
  if(qt<1.E-03) {results[0]=results[1]=0.; return;}
  rombergerN(this, &InclusiveCross::int_qphi_bis,0.,2*PI,2,results,PREC,3,6,pqtphiestimate,qt,pperp,qvec,nu);
  
  results[0]*=qt;
  results[1]*=qt;
  return;  
}


void InclusiveCross::int_qphi_bis(double qphi, double *results, va_list ap){

  double qt = va_arg(ap,double);
  double pperp = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double nu = va_arg(ap,double);
  
  double pperp2other=pperp*pperp+qt*qt-2.*pperp*qt*cos(qphi);  
  //double Er=sqrt(massr*massr+pkin->GetPklab()*pkin->GetPklab());
  //double prz=pkin->GetPklab()*pkin->GetCosthklab();
  double sinqphi,cosqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double cosphiprime = (pperp-qt*cosqphi)/sqrt(pperp2other);
  double sinphiprime = -qt*sinqphi/sqrt(pperp2other);
  if(abs(pperp2other)<1.E-03){
    pperp2other=0.;
    cosphiprime=1.;
    sinphiprime=0.;
  }
    
  otherWx2=-1.;
  Wxprime2=-1;
  get_przs(pperp*pperp,pperp2other,qvec,nu,0);
//   cout << prz << " " << przprime << " " << Wxprime2 << " " << otherWx2 << endl;
  if(otherWx2<0.||Wxprime2<0.) results[0]=0.;
  else{
    double pprime = sqrt(pperp2other+przprime*przprime);
    double costhetaprime=przprime/pprime;
    double sinthetaprime=sqrt(pperp2other)/pprime;
    if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=1.;}
    double Erprime=sqrt(massr*massr+pprime*pprime);
    double pp = sqrt(pperp*pperp+prz*prz);
    double costhetapp=prz/pp;
    double sinthetapp=sqrt(pperp*pperp)/pp;
    if(pp<1.E-03) {costhetapp=1.;sinthetapp=1.;}
    double Erpp=sqrt(massr*massr+pp*pp);    
    if(Erprime>MASSD||Erpp>MASSD) {results[0]=0.;}    
    else{
      TKinematics2to2 kin("","",MASSD,massr,sqrt(Wxprime2),"qsquared:wlab:pklab",qvec*qvec-nu*nu,nu,pp);
      double structfactor=structure.getInclStructure(kin,MASSD-Erpp);
      if(abs(structfactor<1E-09)||isnan(structfactor)) {results[0]=0.; }
      else{
	double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(Wxprime2+massr*massr)+pow(massr*massr-Wxprime2,2.));        
	double t=(Erpp-Erprime)*(Erpp-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
	double offshellness=0.;
	if(offshellset==0){
	  double massdiff=(MASSD-Erpp)*(MASSD-Erpp)-pp*pp-massi*massi+2.*nu*(MASSD-Erpp-massi)+2.*qvec*prz;
	  offshellness=exp(beta*massdiff);
	}
	if(offshellset==1) {
	  double onshellm=nu*nu-qvec*qvec+massi*massi+2.*nu*massi;
	  offshellness=pow(lambda*lambda-onshellm,2.)/(pow(Wxprime2-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
	}
	if(offshellset==2){
	betaoff=16.8*1.E-06;
	sigma=(25.3+53*(sqrt(Wxprime2)-MASSP)*1.E-03)/(1.E-05*kin.GetQsquared())*INVHBARC*INVHBARC;
	offshellness=exp((betaoff-beta)*t/2.);
	}
	if(offshellset==3) offshellness=0.;
	if(offshellset==4) offshellness=1.;
	complex<double> wave1[12],wave2[12];
	TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
	TVector3 vecpp(pp*sinthetapp,0.,pp*costhetapp);
	for(int i=0;i<12;i++){
	  wave1[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecpp)
			  +offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecpp); 
	  wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
			  -offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
	}
	complex<double> temp=0.;
	for(int i=0;i<12;i++) temp+=wave1[i]*conj(wave2[i]);
	results[0]= imag(scatter(t)*temp
	      *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime))*sqrt(MASSD/(2.*(MASSD-Erpp)*Erpp)))*chi*structfactor;
	
      }
    }
  }
  
  otherWx2=-1.;
  Wxprime2=-1;
  get_przs(pperp*pperp,pperp2other,qvec,nu,0);
  if(otherWx2<0.||Wxprime2<0.) {results[1]=0.;return;}
  double pprime = sqrt(pperp2other+przprime*przprime);
  double costhetaprime=przprime/pprime;
  double sinthetaprime=sqrt(pperp2other)/pprime;
  if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=1.;}
  double Erprime=sqrt(massr*massr+pprime*pprime);
  double pp = sqrt(pperp*pperp+prz*prz);
  double costhetapp=prz/pp;
  double sinthetapp=sqrt(pperp*pperp)/pp;
  if(pp<1.E-03) {costhetapp=1.;sinthetapp=1.;}
  double Erpp=sqrt(massr*massr+pp*pp);    
  if(Erprime>MASSD||Erpp>MASSD) {results[1]=0.; return;}    

  TKinematics2to2 kin("","",MASSD,massr,sqrt(otherWx2),"qsquared:wlab:pklab",qvec*qvec-nu*nu,nu,pprime);
  double structfactor=structure.getInclStructure(kin,MASSD-Erprime);
  if(abs(structfactor<1E-09)||isnan(structfactor)) {results[1]=0.; return;}
  double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(otherWx2+massr*massr)+pow(massr*massr-otherWx2,2.));        
  double t=(Erpp-Erprime)*(Erpp-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
  double offshellness=0.;
  if(offshellset==0){
    double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*nu*(MASSD-Erprime-massi)+2.*qvec*przprime;
    offshellness=exp(beta*massdiff);
  }
  if(offshellset==1) {
    double onshellm=nu*nu-qvec*qvec+massi*massi+2.*nu*massi;
    offshellness=pow(lambda*lambda-onshellm,2.)/(pow(Wxprime2-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
  }
  if(offshellset==2){
  betaoff=16.8*1.E-06;
  sigma=(25.3+53*(sqrt(otherWx2)-MASSP)*1.E-03)/(1.E-05*kin.GetQsquared())*INVHBARC*INVHBARC;
  offshellness=exp((betaoff-beta)*t/2.);
  }
  if(offshellset==3) offshellness=0.;
  if(offshellset==4) offshellness=1.;
  complex<double> wave1[12],wave2[12];
  TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
  TVector3 vecpp(pp*sinthetapp,0.,pp*costhetapp);
  for(int i=0;i<12;i++){
    wave1[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecpp)
		    -offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecpp); 
    wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
		    +offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
  }
  complex<double> temp=0.;
  for(int i=0;i<12;i++) temp+=wave2[i]*conj(wave1[i]);
  results[1]= imag(scatter(t)*temp
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime))*sqrt(MASSD/(2.*(MASSD-Erpp)*Erpp)))*chi*structfactor;
  
  return;  
}


void InclusiveCross::get_przs(double pperp2, double pperp2other, double qvec, double nu, int first){
  prz=0.;
  przprime=0.;
  for(int i=0;i<50;i++){
    double f_Erprime=sqrt(massr*massr+pperp2+prz*prz);  //guess for on-shell energy of initial spectator
    double f_Erprimeother=sqrt(massr*massr+pperp2other+przprime*przprime);  //guess for on-shell energy of initial spectator
    //guess for invariant mass initial X'
    double f_Wxprime2=nu*nu-qvec*qvec
		      +pow(MASSD-f_Erprime,2.)-pperp2-prz*prz
		      +2.*nu*(MASSD-f_Erprime)+2.*qvec*prz;
    double f_Wxprime2other=nu*nu-qvec*qvec
		      +pow(MASSD-f_Erprimeother,2.)-pperp2other-przprime*przprime
		      +2.*nu*(MASSD-f_Erprimeother)+2.*qvec*przprime;
		      
    double f_prz=przprime-(nu+MASSD)/qvec*(f_Erprimeother-f_Erprime);  //new value for prz 
    double f_przprime=prz-(nu+MASSD)/qvec*(f_Erprime-f_Erprimeother);  //new value for prz 
    //add mass diff term if necessary
    if(!symm){
      if(first&&(f_Wxprime2other<f_Wxprime2)){  //other one is first
	f_prz-=(f_Wxprime2other-f_Wxprime2)/(2.*qvec);
	f_przprime-=(f_Wxprime2-f_Wxprime2other)/(2.*qvec);
      }
      if((!first)&&(f_Wxprime2other>f_Wxprime2)) {  //regular one is first
	f_prz-=(f_Wxprime2other-f_Wxprime2)/(2.*qvec);
	f_przprime-=(f_Wxprime2-f_Wxprime2other)/(2.*qvec);
      }
    }
    //check convergence
    if((abs((prz-f_prz)/f_prz)<1e-03)&&(abs((przprime-f_przprime)/f_przprime)<1e-03)) {
      Wxprime2= f_Wxprime2; 
      prz = f_prz; 
      otherWx2= f_Wxprime2other; 
      przprime = f_przprime; 
      return;      
    }
    //start again
    Wxprime2= f_Wxprime2; 
    prz = f_prz; 
    otherWx2= f_Wxprime2other; 
    przprime = f_przprime; 
  }
  //no real convergence
  return;
}

*/



