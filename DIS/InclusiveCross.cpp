#include "InclusiveCross.hpp"
#include <Utilfunctions.hpp>

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
structure(elec,proton,strucname){

  wf = TDeuteron::Wavefunction::CreateWavefunction(wavename);
  
}


InclusiveCross::~InclusiveCross(){
 delete wf; 
}


double InclusiveCross::calc_F2Dinc(double Q2,double x){
  
  double result;
  double prestimate=0.,cosestimate=0.;
  rombergerN(this,&InclusiveCross::int_pr,0.,1.e03,1,&result,PREC,3,7,&prestimate, x, Q2, &cosestimate);
  return 2.*PI*2.*massi/MASSD*result;
  
  
}


void InclusiveCross::int_pr(double prnorm, double *result, va_list ap){
  
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double *pcosestimate = va_arg(ap,double*);

  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double lo=-1;
  double hi=1.;
  double lowerlimit= -(-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*massi*x))/(2.*qvec*prnorm);
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
  //cout << *result << " " << structfactor << " " << Wx << endl;
  return;
  
}



void InclusiveCross::calc_F2DincFSI(double &fsi1, double &fsi2, double Q2,double x){
  
  double results[2];
  double prestimate=0.,cosestimate=0.;
  rombergerN(this,&InclusiveCross::int_pr_fsi,0.,1.e03,2,results,PREC,3,7,&prestimate, x, Q2, &cosestimate);
  fsi1=results[0]*2.*PI*2.*massi/MASSD;
  fsi2=results[1]*2.*PI*2.*massi/MASSD;
 cout << fsi1 << " " << fsi2 << endl;
  
  
}


void InclusiveCross::int_pr_fsi(double prnorm, double *results, va_list ap){
  
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double *pcosestimate = va_arg(ap,double*);

  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double lo=-1;
  double hi=1.;
  double lowerlimit= -(-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*massi*x))/(2.*qvec*prnorm);
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
  
  complex<double> wave[12];
  for(int i=0;i<12;i++) wave[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, TVector3(prnorm*sqrt(1.-costheta*costheta),0.,
								   prnorm*costheta)); 
  double qresults[2];
  double qestimate=0.,qphiestimate=0.;
  rombergerN(this,&InclusiveCross::int_qt,0.,1.E03,2,qresults,PREC,3,6,&qestimate,&kin,wave,&qphiestimate);
  
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
  double *pqphiestimate = va_arg(ap,double*);
  
  if(qt<1.E-03) {results[0]=results[1]=0.; return;}
  rombergerN(this, &InclusiveCross::int_qphi,0.,2*PI,2,results,PREC,3,6,pqphiestimate,qt,pkin,wave);
  
  results[0]*=qt;
  results[1]*=qt;
  return;  
}


void InclusiveCross::int_qphi(double qphi, double *results, va_list ap){

  double qt = va_arg(ap,double);
  TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
  complex<double>* wave = va_arg(ap,complex<double>*);
  
  double prt=pkin->GetPklab()*sqrt(1.-pkin->GetCosthklab()*pkin->GetCosthklab());
  double pt2=prt*prt+qt*qt-2.*prt*qt*cos(qphi);  
  double Er=sqrt(massr*massr+pkin->GetPklab()*pkin->GetPklab());
  prz=pkin->GetPklab()*pkin->GetCosthklab();
  double sinqphi,cosqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double cosphiprime = (prt-qt*cosqphi)/sqrt(pt2);
  double sinphiprime = -qt*sinqphi/sqrt(pt2);
  if(abs(pt2)<1.E-03){
    pt2=0.;
    cosphiprime=1.;
    sinphiprime=0.;
  }
    
  otherWx2=-1.;
  get_prz(pt2,Er,pkin,1);
  if(otherWx2<0.) results[0]=0.;
  else{
    double pprime = sqrt(pt2+przprime*przprime);
    double costhetaprime=przprime/pprime;
    double sinthetaprime=sqrt(pt2)/pprime;
    if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
    double Erprime=sqrt(massr*massr+pprime*pprime);
    if(Erprime>MASSD) {results[0]=0.;}
    else{
      //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
      double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
      double offshellness=0.;
      if(offshellset==0){
	double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*prz;
	offshellness=exp(beta*massdiff);
      }
      if(offshellset==1) {
	double onshellm=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()+massi*massi+2.*pkin->GetWlab()*massi;
	offshellness=pow(lambda*lambda-onshellm,2.)/(pow(Wxprime2-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
      }
      if(offshellset==2){
	betaoff=16.8*1.E-06;
	sigma=(25.3*1.E-06*pkin->GetQsquared()+53*(sqrt(otherWx2)-MASSP)*1.E-03)/(1.E-05*pkin->GetQsquared())*INVHBARC*INVHBARC;
	offshellness=exp((betaoff-beta)*t/2.);
      }
      if(offshellset==3) offshellness=0.;
      if(offshellset==4) offshellness=1.;
      complex<double> wave2[12];
      TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
      for(int i=0;i<12;i++) wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
			-offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
      complex<double> temp=0.;
      for(int i=0;i<12;i++) temp+=wave[i]*conj(wave2[i]);
      results[0]= imag(scatter(t)*temp
	    *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));

    }
  }
  
  otherWx2=-1.;
  get_prz(pt2,Er,pkin,0);
//   
  if(otherWx2<0.) {results[1]=0.; return;}
  double chi=sqrt(pkin->GetS()*pkin->GetS()-2.*pkin->GetS()*(otherWx2+massr*massr)+pow(massr*massr-otherWx2,2.));
  double pprime = sqrt(pt2+przprime*przprime);
  double costhetaprime=przprime/pprime;
  double sinthetaprime=sqrt(pt2)/pprime;
  if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
  double Erprime=sqrt(massr*massr+pprime*pprime);
  if(Erprime>MASSD) {results[1]=0.; return;}
  //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
  double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
  double offshellness=0.;
  if(offshellset==0){
    double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*przprime;
    offshellness=exp(beta*massdiff);
  }
  if(offshellset==1) {
    double onshellm=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()+massi*massi+2.*pkin->GetWlab()*massi;
    offshellness=pow(lambda*lambda-onshellm,2.)/(pow(Wxprime2-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
  }
  if(offshellset==2){
    betaoff=16.8*1.E-06;
    sigma=(25.3*1.E-06*pkin->GetQsquared()+53*(sqrt(Wxprime2)-MASSP)*1.E-03)/(1.E-05*pkin->GetQsquared())*INVHBARC*INVHBARC;
   offshellness=exp((betaoff-beta)*t/2.);
  }
  if(offshellset==3) offshellness=0.;
  if(offshellset==4) offshellness=1.;
  complex<double> wave2[12];
  TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
  for(int i=0;i<12;i++) wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
		    +offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
  complex<double> temp=0.;
  for(int i=0;i<12;i++) temp+=conj(wave[i])*wave2[i];
  results[1]= imag(scatter(t)*temp
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)))*chi;

  return;  
}


void InclusiveCross::get_prz(double pt2, double Er, TKinematics2to2 *pkin, int first){
  przprime=0.;
  for(int i=0;i<50;i++){
    double f_Erprime=sqrt(massr*massr+pt2+przprime*przprime);  //guess for on-shell energy of initial spectator
    //guess for invariant mass initial X'
    double f_Wxprime2=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()
		      +pow(MASSD-f_Erprime,2.)-pt2-przprime*przprime
		      +2.*pkin->GetWlab()*(MASSD-f_Erprime)+2.*pkin->GetKlab()*przprime;
		      
    double f_prz=prz-(pkin->GetWlab()+MASSD)/pkin->GetKlab()*(Er-f_Erprime);  //new value for prz 
    //add mass diff term if necessary
    if(!symm){
      if(first/*&&(Wxprime2<f_Wxprime2)*/) f_prz-=(Wxprime2-f_Wxprime2)/(2.*pkin->GetKlab());
      if((!first)/*&&(Wxprime2>f_Wxprime2)*/) f_prz-=(Wxprime2-f_Wxprime2)/(2.*pkin->GetKlab());
    }
    //check convergence
    if((abs((przprime-f_prz)/f_prz)<1e-03)) {otherWx2= f_Wxprime2; przprime = f_prz; return;}
    //start again
    przprime=f_prz;
    otherWx2=f_Wxprime2;  
  }
  //no real convergence
  return;
}

complex<double> InclusiveCross::scatter(double t){
  return sigma*(I+epsilon)*exp(beta*t/2.); 
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


InclusiveCrossRes::InclusiveCrossRes(bool proton, string strucname, string wavename, TElectronKinematics &elec, int sym, int offshell)
:massi(proton? MASSP:MASSN),
massr(proton? MASSN:MASSP),
symm(sym),
prz(-9999.),
offshellset(offshell),
electron(elec),
structure(elec,proton,strucname){

  wf = TDeuteron::Wavefunction::CreateWavefunction(wavename);
  
}


InclusiveCrossRes::~InclusiveCrossRes(){
 delete wf; 
}


void InclusiveCrossRes::addResonance(Resonance &res){
  
  resonance_vec.push_back(res);
}


double InclusiveCrossRes::calc_F2Dinc(double Q2,double x){
  
  double result;
  double prestimate=0.,cosestimate=0.;
  rombergerN(this,&InclusiveCrossRes::int_pr,0.,1.e03,1,&result,PREC,3,7,&prestimate, x, Q2, &cosestimate);
  return 2.*PI*2.*massi/MASSD*result;
  
  
}


void InclusiveCrossRes::int_pr(double prnorm, double *result, va_list ap){
  
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double *pcosestimate = va_arg(ap,double*);

  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double lo=-1;
  double hi=1.;
  double lowerlimit= -(-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*massi*x))/(2.*qvec*prnorm);
  if(lowerlimit>lo) {lo=lowerlimit;}
  if(hi<lo){ *result=0.; return;}
  rombergerN(this, &InclusiveCrossRes::int_costheta_incl,lo+1e-04,hi-1.e-04,1,result,PREC,3,7,pcosestimate,prnorm, x,Q2);
  if(prnorm<1.e-04) *result=0.; 
  else *result*=prnorm*prnorm*(pow(wf->GetUp(prnorm),2.)+pow(wf->GetWp(prnorm),2.))/(4.*PI)*(MASSD/(2.*(MASSD-Er)));
  return;
}


void InclusiveCrossRes::int_costheta_incl(double costheta, double *result, va_list ap){
  
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
  //cout << *result << " " << structfactor << " " << Wx << endl;
  return;
  
}



void InclusiveCrossRes::calc_F2DincFSI(double &fsi1, double &fsi2, double Q2,double x){
  
  fsi1=fsi2=0.;
  double s = Q2*(MASSD/(massi*x)-1.)+MASSD*MASSD;
  size_t max=getResonance_vec().size();
  for(size_t res1=0; res1<getResonance_vec().size(); res1++) {
    if(pow(getResonance_vec()[res1].getMass()+massr,2.)>s){
      max = res1;
      break;
    }
  }
  for(size_t res1=0; res1<max; res1++){
    for(size_t res2=0; res2<max; res2++){
      double results[2];
      double prestimate=0.,cosestimate=0.;
      rombergerN(this,&InclusiveCrossRes::int_pr_fsi,0.,1.e03,2,results,PREC,3,7,&prestimate, x, Q2, res1, res2, &cosestimate);
      double chi2 = sqrt(s*s-2.*s*(getResonance_vec()[res1].getMass2()+massr*massr)+pow(massr*massr-getResonance_vec()[res1].getMass2(),2.));
      //cout << res1 << " " << res2 << " " << results[0]*chi2 << " " << results[1]*chi2 << endl;
      //cout << res1 << " " << res2 << " " << s*1.E-06 << " " << chi2*1.E-06 << endl;
      fsi1+=results[0]*2.*PI*2.*massi/MASSD*chi2;
      fsi2+=results[1]*2.*PI*2.*massi/MASSD*chi2;                
    }
  }
 //cout << fsi1 << " " << fsi2 << endl;
  
  
}


void InclusiveCrossRes::int_pr_fsi(double prnorm, double *results, va_list ap){
  
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  size_t res1 = va_arg(ap, size_t);
  size_t res2 = va_arg(ap, size_t);
  double *pcosestimate = va_arg(ap,double*);

  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double lo=-1;
  double hi=1.;
  double lowerlimit= -(-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*massi*x))/(2.*qvec*prnorm);
  if(lowerlimit>lo) {lo=lowerlimit;}
  if(hi<lo){ results[0]=0.; results[1]=0.; return;}
  rombergerN(this, &InclusiveCrossRes::int_costheta_incl_fsi,lo+1e-04,hi-1.e-04,2,results,PREC,3,7,pcosestimate,prnorm, x,Q2,res1,res2);

  results[0]*=prnorm*prnorm;
  results[1]*=prnorm*prnorm;
  return;
}


void InclusiveCrossRes::int_costheta_incl_fsi(double costheta, double *results, va_list ap){
  
  double prnorm = va_arg(ap,double);
  double x = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  size_t res1 = va_arg(ap, size_t);
  size_t res2 = va_arg(ap, size_t);

  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double nu=Q2/(2.*massi*x);
  Wxprime2=MASSD*MASSD-Q2+massr*massr+2.*MASSD*(nu-Er)-2.*nu*Er+2.*qvec*prnorm*costheta;
  //TKinematics2to2 kin("","",MASSD,massr,sqrt(getResonance_vec()[res1].getMass2()),"qsquared:wlab:pklab",Q2,nu,prnorm);
  TKinematics2to2 kin("","",MASSD,massr,sqrt(Wxprime2),"qsquared:wlab:pklab",Q2,nu,prnorm);
  double structfactor=structure.getInclStructure(kin,MASSD-Er);
  if(abs(structfactor<1E-09)||isnan(structfactor)) {results[0]=0.; results[1]=0.; return; }
//   double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*
// 	      (getResonance_vec()[res1].getMass2()+massr*massr)+pow(massr*massr-getResonance_vec()[res1].getMass2(),2.));
  complex<double> wave[12];
  for(int i=0;i<12;i++) wave[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, TVector3(prnorm*sqrt(1.-costheta*costheta),0.,
								   prnorm*costheta)); 
  double qresults[2];
  double qestimate=0.,qphiestimate=0.;
  rombergerN(this,&InclusiveCrossRes::int_qt,0.,1.E03,2,qresults,PREC,3,6,&qestimate,&kin,wave,res1,res2,&qphiestimate);
  
  results[0]= structfactor/(kin.GetKlab()*32.*PI*PI*3.)
		*sqrt(MASSD/(2.*(MASSD-Er)*Er))/**chi*/;
  results[1]= results[0];
  results[0]*=qresults[0];
  results[1]*=qresults[1];
  return;
  
}

void InclusiveCrossRes::int_qt(double qt, double *results, va_list ap){

  TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
  complex<double>* wave = va_arg(ap,complex<double>*);
  size_t res1 = va_arg(ap, size_t);
  size_t res2 = va_arg(ap, size_t);
  double *pqphiestimate = va_arg(ap,double*);
  
  if(qt<1.E-03) {results[0]=results[1]=0.; return;}
  rombergerN(this, &InclusiveCrossRes::int_qphi,0.,2*PI,2,results,PREC,3,6,pqphiestimate,qt,pkin,wave,res1,res2);
  
  results[0]*=qt;
  results[1]*=qt;
  return;  
}


void InclusiveCrossRes::int_qphi(double qphi, double *results, va_list ap){

//   double qt = va_arg(ap,double);
//   TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
//   complex<double>* wave = va_arg(ap,complex<double>*);
//   size_t res1 = va_arg(ap, size_t);
//   size_t res2 = va_arg(ap, size_t);
//   
//   double prt=pkin->GetPklab()*sqrt(1.-pkin->GetCosthklab()*pkin->GetCosthklab());
//   double pt2=prt*prt+qt*qt-2.*prt*qt*cos(qphi);  
//   double Er=sqrt(massr*massr+pkin->GetPklab()*pkin->GetPklab());
//   prz=pkin->GetPklab()*pkin->GetCosthklab();
//   double sinqphi,cosqphi;
//   sincos(qphi,&sinqphi,&cosqphi);
//   double cosphiprime = (prt-qt*cosqphi)/sqrt(pt2);
//   double sinphiprime = -qt*sinqphi/sqrt(pt2);
//   if(abs(pt2)<1.E-03){
//     pt2=0.;
//     cosphiprime=1.;
//     sinphiprime=0.;
//   }
//     
//   get_prz(pt2,Er,pkin,1,res1,res2);
//   double pprime = sqrt(pt2+przprime*przprime);
//   double costhetaprime=przprime/pprime;
//   double sinthetaprime=sqrt(pt2)/pprime;
//   if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=1.;}
//   double Erprime=sqrt(massr*massr+pprime*pprime);
//   if(Erprime>MASSD) {results[0]=0.;}
//   else{
//     //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
//     double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
//     double offshellness=0.;
//     if(offshellset==0){
//       double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*prz;
//       offshellness=exp((getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.*massdiff);
//     }
//     if(offshellset==1) {
//       double onshellm=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()+massi*massi+2.*pkin->GetWlab()*massi;
//       double lambda=(getResonance_vec()[res1].getLambda()+getResonance_vec()[res2].getLambda());
//       offshellness=pow(lambda*lambda-onshellm,2.)/(pow(getResonance_vec()[res1].getMass2()-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
//     }
//     if(offshellset==2){
// //       betaoff=16.8*1.E-06;
//       //sigma=(25.2+(53*(sqrt(otherWx2)-MASSP))/(1.E-03*pkin->GetQsquared()))/10*INVHBARC*INVHBARC;
//       offshellness=exp(((getResonance_vec()[res1].getBetaoff()+getResonance_vec()[res2].getBetaoff())
// 		      -(getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta()))*t/4.);
//     }
//     if(offshellset==3) offshellness=0.;
//     if(offshellset==4) offshellness=1.;
//     complex<double> wave2[12];
//     TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
//     for(int i=0;i<12;i++) wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
// 		      -offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
//     complex<double> temp=0.;
//     for(int i=0;i<12;i++) temp+=wave[i]*conj(wave2[i]);
//     results[0]= imag(scatter(t,pkin->GetQsquared(),res1,res2)*temp
// 	  *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));
// 
//   }
// 
//   
//   //otherWx2=-1.;
//   get_prz(pt2,Er,pkin,0,res1,res2);
// //   
//   pprime = sqrt(pt2+przprime*przprime);
//   costhetaprime=przprime/pprime;
//   sinthetaprime=sqrt(pt2)/pprime;
//   if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=1.;}
//   Erprime=sqrt(massr*massr+pprime*pprime);
//   if(Erprime>MASSD) {results[1]=0.; return;}
//   //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
//   double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
//   double offshellness=0.;
//   if(offshellset==0){
//     double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*przprime;
//     offshellness=exp((getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.*massdiff);
//   }
//   if(offshellset==1) {
//     double onshellm=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()+massi*massi+2.*pkin->GetWlab()*massi;
//     double lambda=(getResonance_vec()[res1].getLambda()+getResonance_vec()[res2].getLambda());
//     offshellness=pow(lambda*lambda-onshellm,2.)/(pow(getResonance_vec()[res1].getMass2()-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
//   }
//   if(offshellset==2){
// //     betaoff=16.8*1.E-06;
// //     sigma=4.*INVHBARC*INVHBARC;//(25.3+53*(sqrt(Wxprime2)-MASSP)*1.E-03)/(1.E-05*pkin->GetQsquared())*INVHBARC*INVHBARC;
//     offshellness=exp(((getResonance_vec()[res1].getBetaoff()+getResonance_vec()[res2].getBetaoff())
// 		      -(getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta()))*t/4.);
//   }
//   if(offshellset==3) offshellness=0.;
//   if(offshellset==4) offshellness=1.;
//   complex<double> wave2[12];
//   TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
//   for(int i=0;i<12;i++) wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
// 		    +offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
//   complex<double> temp=0.;
//   for(int i=0;i<12;i++) temp+=conj(wave[i])*wave2[i];
//   results[1]= imag(scatter(t,pkin->GetQsquared(),res1,res2)*temp
// 	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));

  double qt = va_arg(ap,double);
  TKinematics2to2 *pkin = va_arg(ap, TKinematics2to2*);
  complex<double>* wave = va_arg(ap,complex<double>*);
  size_t res1 = va_arg(ap, size_t);
  size_t res2 = va_arg(ap, size_t);
  
  double prt=pkin->GetPklab()*sqrt(1.-pkin->GetCosthklab()*pkin->GetCosthklab());
  double pt2=prt*prt+qt*qt-2.*prt*qt*cos(qphi);  
  double Er=sqrt(massr*massr+pkin->GetPklab()*pkin->GetPklab());
  prz=pkin->GetPklab()*pkin->GetCosthklab();
  double sinqphi,cosqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double cosphiprime = (prt-qt*cosqphi)/sqrt(pt2);
  double sinphiprime = -qt*sinqphi/sqrt(pt2);
  if(abs(pt2)<1.E-03){
    pt2=0.;
    cosphiprime=1.;
    sinphiprime=0.;
  }
    
  otherWx2=-1.;
  get_prz(pt2,Er,pkin,1,res1,res2);
  if(otherWx2<0.) results[0]=0.;
  else{
    double pprime = sqrt(pt2+przprime*przprime);
    double costhetaprime=przprime/pprime;
    double sinthetaprime=sqrt(pt2)/pprime;
    if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
    double Erprime=sqrt(massr*massr+pprime*pprime);
    if(Erprime>MASSD) {results[0]=0.;}
    else{
      double Ex = sqrt(getResonance_vec()[res1].getMass2()+(pkin->GetKlab()-prz)*(pkin->GetKlab()-prz)+prt*prt);
      double Exprime = sqrt(getResonance_vec()[res2].getMass2()+(pkin->GetKlab()-przprime)*(pkin->GetKlab()-przprime)+pt2);
      //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
      //double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
      double t = (Ex-Exprime)*(Ex-Exprime)-(prz-przprime)*(prz-przprime)-qt*qt;
      double offshellness=0.;
      if(offshellset==0){
	double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*prz;
	offshellness=exp((getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.*massdiff);
      }
      if(offshellset==1) {
	double onshellm=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()+massi*massi+2.*pkin->GetWlab()*massi;
	double lambda=(getResonance_vec()[res1].getLambda()+getResonance_vec()[res2].getLambda());
	offshellness=pow(lambda*lambda-onshellm,2.)/(pow(getResonance_vec()[res1].getMass2()-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
      }
      if(offshellset==2){
	offshellness=exp(((getResonance_vec()[res1].getBetaoff()+getResonance_vec()[res2].getBetaoff())
			  -(getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta()))*t/4.);
      }
      if(offshellset==3) offshellness=0.;
      if(offshellset==4) offshellness=1.;
      complex<double> wave2[12];
      TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
      for(int i=0;i<12;i++) wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
			-offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
      complex<double> temp=0.;
      for(int i=0;i<12;i++) temp+=wave[i]*conj(wave2[i]);
      results[0]= imag(scatter(t,pkin->GetQsquared(),res1,res2)*temp
	    *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));

    }
  }
  
  otherWx2=-1.;
  get_prz(pt2,Er,pkin,0,res1,res2);
//   
  if(otherWx2<0.) {results[1]=0.; return;}
  double pprime = sqrt(pt2+przprime*przprime);
  double costhetaprime=przprime/pprime;
  double sinthetaprime=sqrt(pt2)/pprime;
  if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
  double Erprime=sqrt(massr*massr+pprime*pprime);
  if(Erprime>MASSD) {results[1]=0.; return;}
  double Ex = sqrt(getResonance_vec()[res2].getMass2()+(pkin->GetKlab()-prz)*(pkin->GetKlab()-prz)+prt*prt);
  double Exprime = sqrt(getResonance_vec()[res1].getMass2()+(pkin->GetKlab()-przprime)*(pkin->GetKlab()-przprime)+pt2);
  //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
  //double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
  double t = (Ex-Exprime)*(Ex-Exprime)-(prz-przprime)*(prz-przprime)-qt*qt;
  double offshellness=0.;
  if(offshellset==0){
    double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin->GetWlab()*(MASSD-Erprime-massi)+2.*pkin->GetKlab()*przprime;
    offshellness=exp((getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.*massdiff);
  }
  if(offshellset==1) {
    double onshellm=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()+massi*massi+2.*pkin->GetWlab()*massi;
    double lambda=(getResonance_vec()[res1].getLambda()+getResonance_vec()[res2].getLambda());
    offshellness=pow(lambda*lambda-onshellm,2.)/(pow(getResonance_vec()[res1].getMass2()-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
  }
  if(offshellset==2){
    offshellness=exp(((getResonance_vec()[res1].getBetaoff()+getResonance_vec()[res2].getBetaoff())
		      -(getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta()))*t/4.);
  }
  if(offshellset==3) offshellness=0.;
  if(offshellset==4) offshellness=1.;
  complex<double> wave2[12];
  TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,pprime*costhetaprime);
  for(int i=0;i<12;i++) wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
		    +offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
  complex<double> temp=0.;
  for(int i=0;i<12;i++) temp+=conj(wave[i])*wave2[i];
  results[1]= imag(scatter(t,pkin->GetQsquared(),res1,res2)*temp
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));

  return;  
}


void InclusiveCrossRes::get_prz(double pt2, double Er, TKinematics2to2 *pkin, int first, size_t res1, size_t res2){
//   przprime=0.;
//   for(int i=0;i<50;i++){
//     double f_Erprime=sqrt(massr*massr+pt2+przprime*przprime);  //guess for on-shell energy of initial spectator
//     //guess for invariant mass initial X'
// //     double f_Wxprime2=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()
// // 		      +pow(MASSD-f_Erprime,2.)-pt2-przprime*przprime
// // 		      +2.*pkin->GetWlab()*(MASSD-f_Erprime)+2.*pkin->GetKlab()*przprime;
// 		      
//     double f_prz=prz-(pkin->GetWlab()+MASSD)/pkin->GetKlab()*(Er-f_Erprime);  //new value for prz 
//     //add mass diff term if necessary
//     if(!symm){
//       if(first/*&&(getResonance_vec()[res1].getMass2()<getResonance_vec()[res1].getMass2())*/) 
// 	f_prz-=(getResonance_vec()[res1].getMass2()-getResonance_vec()[res2].getMass2())/(2.*pkin->GetKlab());
//       if((!first)/*&&(getResonance_vec()[res1].getMass2()>getResonance_vec()[res1].getMass2())*/) 
// 	f_prz-=(getResonance_vec()[res2].getMass2()-getResonance_vec()[res1].getMass2())/(2.*pkin->GetKlab());
//     }
//     //check convergence
//     if((abs((przprime-f_prz)/f_prz)<1e-03)) {przprime = f_prz; return;}
//     //start again
//     przprime=f_prz;
//     //otherWx2=f_Wxprime2;  
//   }
//   //no real convergence


  przprime=0.;
  for(int i=0;i<50;i++){
    double f_Erprime=sqrt(massr*massr+pt2+przprime*przprime);  //guess for on-shell energy of initial spectator
    //guess for invariant mass initial X'
    double f_Wxprime2=pkin->GetWlab()*pkin->GetWlab()-pkin->GetKlab()*pkin->GetKlab()
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
    if((abs((przprime-f_prz)/f_prz)<1e-03)) {otherWx2= f_Wxprime2; przprime = f_prz; return;}
    //start again
    przprime=f_prz;
    otherWx2=f_Wxprime2;  
  }


  return;
}

complex<double> InclusiveCrossRes::scatter(double t, double Q2, size_t res1, size_t res2){
  double sigma = (getResonance_vec()[res1].getSigma(Q2)+getResonance_vec()[res2].getSigma(Q2))/2.;
  double beta = (getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.;
  double epsilon = (getResonance_vec()[res1].getEpsilon()+getResonance_vec()[res2].getEpsilon())/2.;
  return getResonance_vec()[res1].getCoeff()*conj(getResonance_vec()[res2].getCoeff())*sigma*(I+epsilon)*exp(beta*t/2.); 
}




Resonance::Resonance(double m, complex<double> coefficient, double s0, double slope,
		     double betain,double epsilonin,double betaoffin, double lambdain): 
mass(m),
mass2(m*m),
sigma0(s0),
sigmaslope(slope),
beta(betain*1.E-06),
epsilon(epsilonin),
betaoff(betaoffin*1.E-06),
lambda(lambdain*1.E+06),
coeff(coefficient){
  
}

Resonance::Resonance(const Resonance &copy): 
mass(copy.getMass()),
mass2(copy.getMass2()),
sigma0(copy.getSigma0()),
sigmaslope(copy.getSigmaslope()),
beta(copy.getBeta()),
epsilon(copy.getEpsilon()),
betaoff(copy.getBetaoff()),
lambda(copy.getLambda()),
coeff(copy.getCoeff()){
  
}

Resonance& Resonance::operator=(const Resonance& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
     (*this).mass=rhs.getMass();
     (*this).mass2=rhs.getMass2();
     (*this).sigma0=rhs.getSigma0();
     (*this).sigmaslope=rhs.getSigmaslope();
     (*this).beta=rhs.getBeta();
     (*this).epsilon=rhs.getEpsilon();
     (*this).betaoff=rhs.getBetaoff();
     (*this).lambda=rhs.getLambda();
     (*this).coeff=rhs.getCoeff();

     
     
  }

  return *this;
}



Resonance::~Resonance(){}

double Resonance::getSigma(double Q2) const{
    return (getSigma0()+getSigmaslope()*(getMass()-MASSP)/(1.E-03*Q2))/10*INVHBARC*INVHBARC;

}