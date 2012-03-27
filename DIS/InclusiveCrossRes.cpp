#include "InclusiveCrossRes.hpp"
#include "Resonance.hpp"
#include <Utilfunctions.hpp>


InclusiveCrossRes::InclusiveCrossRes(bool proton, string strucname, string wavename, 
				     TElectronKinematics &elec, int sym, int offshell, bool fixpropp, int t_c)
:massi(proton? MASSP:MASSN),
massr(proton? MASSN:MASSP),
symm(sym),
fixprop(fixpropp),
t_choice(t_c),
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
  rombergerN(this,&InclusiveCrossRes::int_pr,0.,1.e03,1,&result,PREC*1.E-03,3,10,&prestimate, x, Q2, &cosestimate);
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
  //Wx^2> massi^2
  double lowerlimit= -(-massi*massi-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*massi*x))/(2.*qvec*prnorm);
  if(lowerlimit>lo) {lo=lowerlimit;}
  if(hi<lo){ *result=0.; return;}
  rombergerN(this, &InclusiveCrossRes::int_costheta_incl,lo+1e-04,hi-1.e-04,1,result,PREC*1.E-03,3,10,pcosestimate,prnorm, x,Q2);
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
      rombergerN(this,&InclusiveCrossRes::int_pr_fsi,0.,1.e03,2,results,PREC,3,10,&prestimate, x, Q2, res1, res2, &cosestimate);
      double chi2 = sqrt(s*s-2.*s*(getResonance_vec()[res1].getMass2()+massr*massr)+pow(massr*massr-getResonance_vec()[res1].getMass2(),2.));
      cout << res1 << " " << res2 << " " << results[0]*chi2 << " " << results[1]*chi2 << endl;
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

  if(prnorm<1.E-03) { results[0]=0.; results[1]=0.; return;}
  double Er=sqrt(massr*massr+prnorm*prnorm);
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double lo=-1;
  double hi=1.;
  //Wx^2> massi^2
  double lowerlimit= -(-massi*massi-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*Q2/(2.*massi*x))/(2.*qvec*prnorm);
  if(lowerlimit>lo) {lo=lowerlimit;}
  if(hi<lo){ results[0]=0.; results[1]=0.; return;}
  rombergerN(this, &InclusiveCrossRes::int_costheta_incl_fsi,lo+1e-04,hi-1.e-04,2,results,PREC,3,8,pcosestimate,prnorm, x,Q2,res1,res2);

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
  TKinematics2to2 *kin1=NULL,*kin2=NULL;
  if(fixprop) {
    kin1 = new TKinematics2to2("","",MASSD,massr,sqrt(getResonance_vec()[res1].getMass2()),"qsquared:wlab:pklab",Q2,nu,prnorm);
    kin2 = new TKinematics2to2("","",MASSD,massr,sqrt(getResonance_vec()[res2].getMass2()),"qsquared:wlab:pklab",Q2,nu,prnorm);
  }
  else kin1= new TKinematics2to2("","",MASSD,massr,sqrt(Wxprime2),"qsquared:wlab:pklab",Q2,nu,prnorm);
  double structfactor1=structure.getInclStructure(*kin1,MASSD-Er);
  double structfactor2=structure.getInclStructure(*kin2,MASSD-Er);
  if(abs(structfactor1<1E-09)||isnan(structfactor1)) {results[0]=0.; results[1]=0.; return; }
  if(abs(structfactor2<1E-09)||isnan(structfactor2)) {results[0]=0.; results[1]=0.; return; }
//   double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*
// 	      (getResonance_vec()[res1].getMass2()+massr*massr)+pow(massr*massr-getResonance_vec()[res1].getMass2(),2.));
  complex<double> wave[12];
  for(int i=0;i<12;i++) wave[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, TVector3(prnorm*sqrt(1.-costheta*costheta),0.,
								   prnorm*costheta)); 
  double qresults[2];
  double qestimate=0.,qphiestimate=0.;
  rombergerN(this,&InclusiveCrossRes::int_qt,0.,1.E03,2,qresults,PREC,3,7,&qestimate,kin1,kin2,wave,res1,res2,&qphiestimate);
  
  results[0]= 1./(32.*PI*PI*3.)
		*sqrt(MASSD/(2.*(MASSD-Er)*Er))/**chi*/;
  results[1]= results[0];
  results[0]*=structfactor1*qresults[0]/kin1->GetKlab();
  results[1]*=structfactor2*qresults[1]/kin2->GetKlab();
  delete kin1;
  if(kin2) delete kin2;
  return;
  
}

void InclusiveCrossRes::int_qt(double qt, double *results, va_list ap){

  TKinematics2to2 *pkin1 = va_arg(ap, TKinematics2to2*);
  TKinematics2to2 *pkin2 = va_arg(ap, TKinematics2to2*);
  complex<double>* wave = va_arg(ap,complex<double>*);
  size_t res1 = va_arg(ap, size_t);
  size_t res2 = va_arg(ap, size_t);
  double *pqphiestimate = va_arg(ap,double*);
  
  if(qt<1.E-03) {results[0]=results[1]=0.; return;}
  rombergerN(this, &InclusiveCrossRes::int_qphi,0.,2*PI,2,results,PREC,3,7,pqphiestimate,qt,pkin1,pkin2,wave,res1,res2);
  
  results[0]*=qt;
  results[1]*=qt;
  return;  
}


void InclusiveCrossRes::int_qphi(double qphi, double *results, va_list ap){

  double qt = va_arg(ap,double);
  TKinematics2to2 *pkin1 = va_arg(ap, TKinematics2to2*);
  TKinematics2to2 *pkin2 = va_arg(ap, TKinematics2to2*);
  complex<double>* wave = va_arg(ap,complex<double>*);
  size_t res1 = va_arg(ap, size_t);
  size_t res2 = va_arg(ap, size_t);
  
  double prt=pkin1->GetPklab()*sqrt(1.-pkin1->GetCosthklab()*pkin1->GetCosthklab());
  double sinqphi,cosqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double pt2=prt*prt+qt*qt-2.*prt*qt*cosqphi;  
  double Er=sqrt(massr*massr+pkin1->GetPklab()*pkin1->GetPklab());
  prz=pkin1->GetPklab()*pkin1->GetCosthklab();
  double cosphiprime = (prt-qt*cosqphi)/sqrt(pt2);
  double sinphiprime = -qt*sinqphi/sqrt(pt2);
  if(abs(pt2)<1.E-03){
    pt2=0.;
    cosphiprime=1.;
    sinphiprime=0.;
  }
 
  //fixed masses in the propagator
  if(fixprop){
    
      
    get_prz(pt2,Er,pkin1,1,res1,res2);
    double pprime = sqrt(pt2+przprime*przprime);
    double costhetaprime=przprime/pprime;
    double sinthetaprime=sqrt(pt2)/pprime;
    if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
    double Erprime=sqrt(massr*massr+pprime*pprime);
    if(Erprime>MASSD) {results[0]=0.;}  //too high spectator momentum
    else{
      double t=0.;
      if(t_choice==0) {
	t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*sqrt(pt2)*cosphiprime;
	double t2=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
	cout << t << " " << t2 << endl;
      }
      
      if(t_choice==1){
	double Ex = sqrt(getResonance_vec()[res1].getMass2()+(pkin1->GetKlab()-prz)*(pkin1->GetKlab()-prz)+prt*prt);
	double Exprime = sqrt(getResonance_vec()[res2].getMass2()+(pkin1->GetKlab()-przprime)*(pkin1->GetKlab()-przprime)+pt2);
	t = (Ex-Exprime)*(Ex-Exprime)-(prz-przprime)*(prz-przprime)-qt*qt;
      }

      if(t_choice==2) t = getResonance_vec()[res1].getMass2() + getResonance_vec()[res2].getMass2()
			  -2.*(-pkin1->GetQsquared()+2.*MASSD*pkin1->GetWlab()+MASSD*MASSD
			  -(MASSD+pkin1->GetWlab())*(Er+Erprime)+pkin1->GetKlab()*(prz+przprime)
			  +Er*Erprime-prt*sqrt(pt2)*cosphiprime-prz*przprime);
      
      double offshellness=0.;
      if(offshellset==0){
	double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin1->GetWlab()*(MASSD-Erprime-massi)+2.*pkin1->GetKlab()*prz;
	offshellness=exp((getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.*massdiff);
      }
      if(offshellset==1) {
	double onshellm=-pkin1->GetQsquared()+massi*massi+2.*pkin1->GetWlab()*massi;
	double lambda=(getResonance_vec()[res1].getLambda()+getResonance_vec()[res2].getLambda())/2.;
	offshellness=pow(lambda*lambda-onshellm,2.)/(pow(getResonance_vec()[res1].getMass2()-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
      }
      if(offshellset==2){
	offshellness=exp(((getResonance_vec()[res1].getBetaoff()+getResonance_vec()[res2].getBetaoff())
			-(getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta()))*t/4.);
      }
      if(offshellset==3) offshellness=0.;
      if(offshellset==4) offshellness=1.;
      complex<double> wave2[12];
      TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,przprime);
      for(int i=0;i<12;i++) wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
			-offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
      complex<double> temp=0.;
      for(int i=0;i<12;i++) temp+=wave[i]*conj(wave2[i]);
      results[0]= imag(scatter(t,pkin1->GetQsquared(),res1,res2)*temp
	    *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));

    }

    
    //otherWx2=-1.;
    get_prz(pt2,Er,pkin2,0,res1,res2);
  //   
    pprime = sqrt(pt2+przprime*przprime);
    costhetaprime=przprime/pprime;
    sinthetaprime=sqrt(pt2)/pprime;
    if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
    Erprime=sqrt(massr*massr+pprime*pprime);
    if(Erprime>MASSD) {results[1]=0.; return;}
    double t=0.;
    if(t_choice==0) {
      t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*sqrt(pt2)*cosphiprime;
      double t2=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
      cout << t << " " << t2 << endl;
    }
    
    if(t_choice==1){
      double Ex = sqrt(getResonance_vec()[res1].getMass2()+(pkin1->GetKlab()-prz)*(pkin1->GetKlab()-prz)+prt*prt);
      double Exprime = sqrt(getResonance_vec()[res2].getMass2()+(pkin1->GetKlab()-przprime)*(pkin1->GetKlab()-przprime)+pt2);
      t = (Ex-Exprime)*(Ex-Exprime)-(prz-przprime)*(prz-przprime)-qt*qt;
    }

    if(t_choice==2) t = getResonance_vec()[res1].getMass2() + getResonance_vec()[res2].getMass2()
			-2.*(-pkin1->GetQsquared()+2.*MASSD*pkin1->GetWlab()+MASSD*MASSD
			-(MASSD+pkin1->GetWlab())*(Er+Erprime)+pkin1->GetKlab()*(prz+przprime)
			+Er*Erprime-prt*sqrt(pt2)*cosphiprime-prz*przprime);
      
    double offshellness=0.;
    if(offshellset==0){
      double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin2->GetWlab()*(MASSD-Erprime-massi)+2.*pkin2->GetKlab()*przprime;
      offshellness=exp((getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.*massdiff);
    }
    if(offshellset==1) {
      double onshellm=-pkin2->GetQsquared()+massi*massi+2.*pkin2->GetWlab()*massi;
      double lambda=(getResonance_vec()[res1].getLambda()+getResonance_vec()[res2].getLambda())/2.;
      offshellness=pow(lambda*lambda-onshellm,2.)/(pow(getResonance_vec()[res1].getMass2()-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
    }
    if(offshellset==2){
      offshellness=exp(((getResonance_vec()[res1].getBetaoff()+getResonance_vec()[res2].getBetaoff())
			-(getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta()))*t/4.);
    }
    if(offshellset==3) offshellness=0.;
    if(offshellset==4) offshellness=1.;
    complex<double> wave2[12];
    TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,przprime);
    for(int i=0;i<12;i++) wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
		      +offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
    complex<double> temp=0.;
    for(int i=0;i<12;i++) temp+=conj(wave[i])*wave2[i];
    results[1]= imag(scatter(t,pkin2->GetQsquared(),res1,res2)*temp
	  *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));
  }
    
    
  else{
    
    otherWx2=-1.;
    get_prz(pt2,Er,pkin1,1,res1,res2);
    if(otherWx2<0.) results[0]=0.;
    else{
      double pprime = sqrt(pt2+przprime*przprime);
      double costhetaprime=przprime/pprime;
      double sinthetaprime=sqrt(pt2)/pprime;
      if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
      double Erprime=sqrt(massr*massr+pprime*pprime);
      if(Erprime>MASSD) {results[0]=0.;}
      else{
	double t=0.;
	if(t_choice==0) {
	  t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*sqrt(pt2)*cosphiprime;
	  double t2=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
	  cout << t << " " << t2 << endl;
	}
	
	if(t_choice==1){
	  double Ex = sqrt(getResonance_vec()[res1].getMass2()+(pkin1->GetKlab()-prz)*(pkin1->GetKlab()-prz)+prt*prt);
	  double Exprime = sqrt(getResonance_vec()[res2].getMass2()+(pkin1->GetKlab()-przprime)*(pkin1->GetKlab()-przprime)+pt2);
	  t = (Ex-Exprime)*(Ex-Exprime)-(prz-przprime)*(prz-przprime)-qt*qt;
	}

	if(t_choice==2) t = getResonance_vec()[res1].getMass2() + getResonance_vec()[res2].getMass2()
			    -2.*(-pkin1->GetQsquared()+2.*MASSD*pkin1->GetWlab()+MASSD*MASSD
			    -(MASSD+pkin1->GetWlab())*(Er+Erprime)+pkin1->GetKlab()*(prz+przprime)
			    +Er*Erprime-prt*sqrt(pt2)*cosphiprime-prz*przprime);
      	double offshellness=0.;
	if(offshellset==0){
	  double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin1->GetWlab()*(MASSD-Erprime-massi)+2.*pkin1->GetKlab()*prz;
	  offshellness=exp((getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.*massdiff);
	}
	if(offshellset==1) {
	  double onshellm=-pkin1->GetQsquared()+massi*massi+2.*pkin1->GetWlab()*massi;
	  double lambda=(getResonance_vec()[res1].getLambda()+getResonance_vec()[res2].getLambda())/2.;
	  offshellness=pow(lambda*lambda-onshellm,2.)/(pow(getResonance_vec()[res1].getMass2()-onshellm,2.)+pow(lambda*lambda-onshellm,2.));
	}
	if(offshellset==2){
	  offshellness=exp(((getResonance_vec()[res1].getBetaoff()+getResonance_vec()[res2].getBetaoff())
			    -(getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta()))*t/4.);
	}
	if(offshellset==3) offshellness=0.;
	if(offshellset==4) offshellness=1.;
	complex<double> wave2[12];
	TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,przprime);
	for(int i=0;i<12;i++) wave2[i] = wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime)
			  -offshellness*wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime); 
	complex<double> temp=0.;
	for(int i=0;i<12;i++) temp+=wave[i]*conj(wave2[i]);
	results[0]= imag(scatter(t,pkin1->GetQsquared(),res1,res2)*temp
	      *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));

      }
    }
    
    otherWx2=-1.;
    get_prz(pt2,Er,pkin1,0,res1,res2);
  //   
    if(otherWx2<0.) {results[1]=0.; return;}
    double pprime = sqrt(pt2+przprime*przprime);
    double costhetaprime=przprime/pprime;
    double sinthetaprime=sqrt(pt2)/pprime;
    if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
    double Erprime=sqrt(massr*massr+pprime*pprime);
    if(Erprime>MASSD) {results[1]=0.; return;}
    double t=0.;
    if(t_choice==0) {
      t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*sqrt(pt2)*cosphiprime;
      double t2=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
      cout << t << " " << t2 << endl;
    }
    
    if(t_choice==1){
      double Ex = sqrt(getResonance_vec()[res1].getMass2()+(pkin1->GetKlab()-prz)*(pkin1->GetKlab()-prz)+prt*prt);
      double Exprime = sqrt(getResonance_vec()[res2].getMass2()+(pkin1->GetKlab()-przprime)*(pkin1->GetKlab()-przprime)+pt2);
      t = (Ex-Exprime)*(Ex-Exprime)-(prz-przprime)*(prz-przprime)-qt*qt;
    }

    if(t_choice==2) t = getResonance_vec()[res1].getMass2() + getResonance_vec()[res2].getMass2()
			-2.*(-pkin1->GetQsquared()+2.*MASSD*pkin1->GetWlab()+MASSD*MASSD
			-(MASSD+pkin1->GetWlab())*(Er+Erprime)+pkin1->GetKlab()*(prz+przprime)
			+Er*Erprime-prt*sqrt(pt2)*cosphiprime-prz*przprime);
    

    double offshellness=0.;
    if(offshellset==0){
      double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*pkin1->GetWlab()*(MASSD-Erprime-massi)+2.*pkin1->GetKlab()*przprime;
      offshellness=exp((getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.*massdiff);
    }
    if(offshellset==1) {
      double onshellm=-pkin1->GetQsquared()+massi*massi+2.*pkin1->GetWlab()*massi;
      double lambda=(getResonance_vec()[res1].getLambda()+getResonance_vec()[res2].getLambda())/2.;
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
    results[1]= imag(scatter(t,pkin1->GetQsquared(),res1,res2)*temp
	  *sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)));
  }
    
  return;  
}


void InclusiveCrossRes::get_prz(double pt2, double Er, TKinematics2to2 *pkin, int first, size_t res1, size_t res2){

  if(fixprop){
    przprime=0.;
    for(int i=0;i<50;i++){
      double f_Erprime=sqrt(massr*massr+pt2+przprime*przprime);  //guess for on-shell energy of initial spectator
      double f_prz=prz-(pkin->GetWlab()+MASSD)/pkin->GetKlab()*(Er-f_Erprime);  //new value for prz 
      //add mass diff term if necessary
      if(symm==0){
	if(first/*&&(getResonance_vec()[res1].getMass2()<getResonance_vec()[res1].getMass2())*/) 
	  f_prz-=(getResonance_vec()[res1].getMass2()-getResonance_vec()[res2].getMass2())/(2.*pkin->GetKlab());
	if((!first)/*&&(getResonance_vec()[res1].getMass2()>getResonance_vec()[res1].getMass2())*/) 
	  f_prz-=(getResonance_vec()[res2].getMass2()-getResonance_vec()[res1].getMass2())/(2.*pkin->GetKlab());
      }
      if(symm==-1){
	if(first&&(getResonance_vec()[res1].getMass2()<getResonance_vec()[res2].getMass2())) 
	  f_prz-=(getResonance_vec()[res1].getMass2()-getResonance_vec()[res2].getMass2())/(2.*pkin->GetKlab());
	if((!first)&&(getResonance_vec()[res1].getMass2()<getResonance_vec()[res2].getMass2())) 
	  f_prz-=(getResonance_vec()[res2].getMass2()-getResonance_vec()[res1].getMass2())/(2.*pkin->GetKlab());
      }
      //check convergence
      if((abs((przprime-f_prz)/f_prz)<1e-03)) {przprime = f_prz; return;}
      //start again
      przprime=f_prz;
      //otherWx2=f_Wxprime2;  
    }
    //no real convergence
  }
  else{
      
    przprime=0.;
    for(int i=0;i<50;i++){
      double f_Erprime=sqrt(massr*massr+pt2+przprime*przprime);  //guess for on-shell energy of initial spectator
      //guess for invariant mass initial X'
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
      if((abs((przprime-f_prz)/f_prz)<1e-03)) {otherWx2= f_Wxprime2; przprime = f_prz; return;}
      //start again
      przprime=f_prz;
      otherWx2=f_Wxprime2;  
    }
  }

  return;
}

complex<double> InclusiveCrossRes::scatter(double t, double Q2, size_t res1, size_t res2){
  double sigma = (getResonance_vec()[res1].getSigma(Q2)+getResonance_vec()[res2].getSigma(Q2))/2.;
  double beta = (getResonance_vec()[res1].getBeta()+getResonance_vec()[res2].getBeta())/2.;
  double epsilon = (getResonance_vec()[res1].getEpsilon()+getResonance_vec()[res2].getEpsilon())/2.;
  return getResonance_vec()[res1].getCoeff()*conj(getResonance_vec()[res2].getCoeff())*sigma*(I+epsilon)*exp(beta*t/2.); 
}



