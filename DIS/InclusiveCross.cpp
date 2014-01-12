#include "InclusiveCross.hpp"
#include <Utilfunctions.hpp>
// #include <gsl/gsl_poly.h>

using namespace std;

InclusiveCross::InclusiveCross(bool proton, string strucname, string wavename, TElectronKinematics &elec, 
			       std::vector<double> & res,int offshell, 
				double sigmain,double betain, double epsilonin, double betaoffin, double lambdain)
:sigma(sigmain/10.*INVHBARC*INVHBARC),
beta(betain*1.E-06),
epsilon(epsilonin),
betaoff(betaoffin*1.E-06),
lambda(lambdain*1.E+06),
massi(proton? MASSP:MASSN),
massr(proton? MASSN:MASSP),
resonances(res),
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
  if((abs(structfactor)<1E-09)||isnan(structfactor)||isinf(structfactor)) {result[0]=0.; return; }
  //no phi integration, symmetry is already used in getInclStructure()
  result[0]=prnorm*prnorm*(pow(cross.wf->GetUp(prnorm),2.)+pow(cross.wf->GetWp(prnorm),2.))/(4.*PI)*(MASSD/(2.*(MASSD-Er)))*structfactor;
//   /*if(cross.massr==MASSN) */cout << prnorm << " " << costheta << " " << result[0] << " " << result[0]/prnorm/prnorm << endl;
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
  fsi1=fsi2=0.;
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
  for(size_t it=0; it<resonances.size(); it++){
    F.it=it;
  int res=90;
    unsigned count=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
    res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,6E04,ret,count,0);
  //   cout << res << " " << count << endl;
    fsi1+= 2.*PI*2.*massi/MASSD*ret[0];
    fsi2+= 2.*PI*2.*massi/MASSD*ret[1];
  }
  return;
}

void InclusiveCross::FSI_int(numint::vector_d & result, double prnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, size_t itt, size_t it2=0){
			
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
  if((abs(structfactor)<1E-09)||isnan(structfactor)||isinf(structfactor)) {result[0]=result[1]=0.; return; }
  //no phi integration, symmetry is already used in getInclStructure()
  double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(cross.Wxprime2+cross.massr*cross.massr)
	      +pow(cross.massr*cross.massr-cross.Wxprime2,2.));
  complex<double> wave[6];
  for(int i=0;i<6;i++){
    wave[i] = cross.wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, TVector3(prnorm*sqrt(1.-costheta*costheta),0.,
								   prnorm*costheta)); 
  }
  double prt=kin.GetPklab()*sqrt(1.-kin.GetCosthklab()*kin.GetCosthklab());
  double sinqphi,cosqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double pt2=prt*prt+qt*qt-2.*prt*qt*cosqphi;  
  double prz=kin.GetPklab()*kin.GetCosthklab();
  double cosphiprime = (prt-qt*cosqphi)/sqrt(pt2);
  double sinphiprime = -qt*sinqphi/sqrt(pt2);
  if(abs(pt2)<1.E-03){
    pt2=0.;
    cosphiprime=1.;
    sinphiprime=0.;
  }

  vector<double> res_result;
  res_result=vector<double>(2,0.);
  cross.otherWx2=cross.resonances[itt]*cross.resonances[itt];
  double przprime = cross.get_prz_res(pt2,cross.otherWx2,kin);
  if(abs(przprime+1.E09)<1.E-02){result[0]=result[1]=0.; return;}
  double pprime = sqrt(pt2+przprime*przprime);
  double costhetaprime=przprime/pprime;
  double sinthetaprime=sqrt(pt2)/pprime;
  if(pprime<1.E-03) {costhetaprime=1.;sinthetaprime=0.;}
  TKinematics2to2 kin2("","",MASSD,cross.massr,sqrt(cross.otherWx2),"qsquared:wlab:pklab",Q2,nu,pprime);
  double Erprime=sqrt(cross.massr*cross.massr+pprime*pprime);
  double xprime=kin.GetQsquared()/(2.*((MASSD-Erprime)*kin.GetWlab()+przprime*kin.GetKlab())
      -cross.massi*cross.massi+pow(MASSD-Erprime,2.)-pprime*pprime);
  if((Erprime>MASSD)/*||xprime>1*/) {result[0]=result[1]=0.; return;}
  double structfactor2=cross.structure.getInclStructure(kin2,MASSD-Erprime);
  if((abs(structfactor2)<1E-09)||isnan(structfactor2)||isinf(structfactor2)) {result[0]=result[1]=0.; return; }
  res_result[0]=res_result[1]=prnorm*prnorm*sqrt(structfactor*structfactor2)/(abs(kin.GetKlab()-przprime/Erprime*(MASSD+nu))*32.*PI*PI*3.)
		*MASSD/sqrt(4.*(MASSD-Er)*Er*(MASSD-Erprime)*Erprime)*qt;
  //double t=2.*massr-2.*Er*Erprime+2.*prz*prz+2.*prt*pt*cos(phiprime-phi);
  //t = (ps-ps')^2
  double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
  complex<double> wave2[6];
  TVector3 vecprime(pprime*sinthetaprime*cosphiprime,pprime*sinthetaprime*sinphiprime,przprime);
  for(int i=0;i<6;i++){
    wave2[i] = cross.wf->DeuteronPState((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, vecprime);
  }
  complex<double> temp=0.;
  for(int i=0;i<6;i++){
    temp+=wave[i]*conj(wave2[i]);
  }
  cross.sigma=InclusiveCross::sigmaparam(cross.Wxprime2,kin.GetQsquared());
  res_result[0]*= imag(cross.scatter(t)*2.*temp)*chi; //factor 2 due to symmetry	
  cross.sigma=InclusiveCross::sigmaparam(cross.otherWx2,kin.GetQsquared());
  double chi2=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(cross.otherWx2+cross.massr*cross.massr)
	  +pow(cross.massr*cross.massr-cross.otherWx2,2.));
  res_result[1]*=imag(cross.scatter(t)*2.*conj(temp))*chi2; //factor 2 due to symmetry
	  
  result[0]+=res_result[0]; result[1]+= res_result[1];  
//   if(cross.massr==MASSN) cout << prnorm << " " << costheta << " " << pprime << " " << costhetaprime << " " << t << " " << qt << " " << 
//     sqrt(MASSD*MASSD+cross.massr*cross.massr-Q2+2.*MASSD*(nu-Erprime)-2.*nu*Erprime+qvec*przprime) << " " << result[0] << " " << result[0]/prnorm/prnorm << endl;
  
  return;
  

}

void InclusiveCross::calc_F2DincFSI_off(double &fsi1, double &fsi2, double Q2,double x){
  
  fsi1=fsi2=0.;
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
  for(size_t it=0; it<resonances.size(); it++){
    F.it=it;
    /*for(size_t it2=0; it2<resonances.size(); it2++)*/{
      F.it2=it;
      int res=90;
      unsigned count=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,6E04,ret,count,0);
    //   cout << res << " " << count << endl;
      fsi1+= 2.*PI*2.*massi/MASSD*ret[0];
      fsi2+= 2.*PI*2.*massi/MASSD*ret[1];
    }
  }
  return;
}


void InclusiveCross::FSI_int_off(numint::vector_d & result, double prt, double W, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x,size_t itt, size_t itt2){
			
  result=numint::vector_d(2,0.);
  if((qt<1.E-03)||(prt<1.E-03)) { result[0]=result[1]=0.; return;}
  
  //photon kinematics
  double qvec=sqrt(Q2+pow(Q2/(2.*cross.massi*x),2.));
  double nu=Q2/(2.*cross.massi*x);

  double sinqphi,cosqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double pt=sqrt(prt*prt+qt*qt-2.*prt*qt*cosqphi);  
  double cosphiprime = (prt-qt*cosqphi)/pt;
  double sinphiprime = -qt*sinqphi/pt;
  if(pt<1.E-03){
    pt=0.;
    cosphiprime=1.;
    sinphiprime=0.;
  }
  
  //1st vertex structure function part.  On-shell with W!!!
  double prz=cross.get_prz_res(prt*prt,W*W,Q2,nu,qvec); //on-shell value for structure functions!! 1st vertex
  if(abs(prz+1.E09)<1.E-02){result[0]=result[1]=0.; return;}
  double prnorm=sqrt(prt*prt+prz*prz);
  double Er=sqrt(cross.massr*cross.massr+prnorm*prnorm);
  double xprime = Q2/(2.*((MASSD-Er)*nu+prz*qvec) -cross.massi*cross.massi+pow(MASSD-Er,2.)-prnorm*prnorm);
  if(xprime>1){result[1]=result[0]=0.;return;}
  TKinematics2to2 kin("","",MASSD,cross.massr,W,"qsquared:wlab:pklab",Q2,nu,prnorm);
  double structfactor=cross.structure.getInclStructure_off(kin,W*W,MASSD-Er);
  if((abs(structfactor)<1E-09)||isnan(structfactor)||isinf(structfactor)) {result[0]=result[1]=0.; return; }
  //no phi integration, symmetry is already used in getInclStructure()

  //2nd vertex structure function part.  On-shell with W!!!
  double przprime=cross.get_prz_res(pt*pt,W*W,Q2,nu,qvec);//on-shell value for structure functions!! 2nd vertex
  if(abs(przprime+1.E09)<1.E-02){result[0]=result[1]=0.; return;}
  double pprime = sqrt(pt*pt+przprime*przprime);
  double Erprime=sqrt(cross.massr*cross.massr+pprime*pprime);
  xprime=Q2/(2.*((MASSD-Erprime)*nu+przprime*qvec)-cross.massi*cross.massi+pow(MASSD-Erprime,2.)-pprime*pprime);
  if(xprime>1) {result[0]=result[1]=0.; return;}
  TKinematics2to2 kin2("","",MASSD,cross.massr,W,"qsquared:wlab:pklab",Q2,nu,pprime);
  double structfactor2=cross.structure.getInclStructure_off(kin2,W*W,MASSD-Erprime);
  if((abs(structfactor2)<1E-09)||isnan(structfactor2)||isinf(structfactor2)) {result[0]=result[1]=0.; return; }
  
  
  
  //off-shell wave function first vertex (on-shell tilde)
  cross.Wxprime2=cross.resonances[itt]*cross.resonances[itt];
  double prztilde=cross.get_prz_res(prt*prt,cross.Wxprime2,Q2,nu,qvec); //"on-shell" prop value (needed in wf evaluation)
  if(abs(prztilde+1.E09)<1.E-02){result[0]=result[1]=0.; return;}
  double prnormtilde=sqrt(prt*prt+prztilde*prztilde);
  double Ertilde=sqrt(cross.massr*cross.massr+prnormtilde*prnormtilde);
  if((Ertilde>MASSD)){result[1]=result[0]=0.;return;}
  double costhetatilde=prztilde/prnormtilde;
  
  complex<double> wave[6];
  for(int i=0;i<6;i++){
    wave[i] = cross.wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, 
		TVector3(prnormtilde*sqrt(1.-costhetatilde*costhetatilde),0.,prnormtilde*costhetatilde)); 
  }

  //off-shell wave function second vertex (on-shell tilde)
  cross.otherWx2=cross.resonances[itt2]*cross.resonances[itt2];
  double prztildeprime=cross.get_prz_res(pt*pt,cross.otherWx2,Q2,nu,qvec);
  if(abs(prztildeprime+1.E09)<1.E-02){result[0]=result[1]=0.; return;}
  double ptildeprime = sqrt(pt*pt+prztildeprime*prztildeprime);
  double costhetatildeprime=prztildeprime/ptildeprime;
  double sinthetatildeprime=pt/ptildeprime;
  if(ptildeprime<1.E-03) {costhetatildeprime=1.;sinthetatildeprime=0.;}
  double Ertildeprime=sqrt(cross.massr*cross.massr+ptildeprime*ptildeprime);
  if((Ertildeprime>MASSD)) {result[0]=result[1]=0.; return;}
  
  complex<double> wave2[6];
  for(int i=0;i<6;i++){
    wave2[i] = -cross.wf->DeuteronPStateOff((i/4)*2-2, ((i/2)%2)*2-1, (i%2)*2-1, 
      TVector3(ptildeprime*sinthetatildeprime*cosphiprime,ptildeprime*sinthetatildeprime*sinphiprime,prztildeprime));
  }
  complex<double> temp=0.;
  for(int i=0;i<6;i++){
    temp+=wave[i]*conj(wave2[i]);
  }

  result[0]=result[1]=prt*W*sqrt(structfactor*structfactor2)/(abs(qvec-prztilde/Ertilde*(MASSD+nu))
			*abs(qvec-prztildeprime/Ertildeprime*(MASSD+nu))*32.*PI*PI*3.)
		      *MASSD/sqrt(4.*(MASSD-Ertilde)*Ertilde*(MASSD-Ertildeprime)*Ertildeprime)*qt;
		      
  //t = (ps-ps')^2
  double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
  
  double offshellness=0.;
  if(cross.offshellset==0){
    double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-cross.massi*cross.massi
      +2.*kin.GetWlab()*(MASSD-Erprime-cross.massi)+2.*kin.GetKlab()*prztildeprime;
    offshellness=exp(cross.beta*massdiff);
  }
  if(cross.offshellset==1) {
    double onshellm=-kin.GetQsquared()+cross.massi*cross.massi+2.*kin.GetWlab()*cross.massi;
    offshellness=pow(cross.lambda*cross.lambda-onshellm,2.)/(pow(cross.Wxprime2-onshellm,2.)+pow(cross.lambda*cross.lambda-onshellm,2.));
  }
  if(cross.offshellset==2){
    cross.betaoff=16.8*1.E-06;
    offshellness=0.;//exp((cross.betaoff-cross.beta)*t/2.);
  }
  if(cross.offshellset==3){
    cross.betaoff=16.8*1.E-06;
    offshellness=exp(-cross.beta*abs(W*W-cross.Wxprime2)/2.);
  }
  if(cross.offshellset==4){
    offshellness=1.;
  }
  
  double chi=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(cross.Wxprime2+cross.massr*cross.massr)
		+pow(cross.massr*cross.massr-cross.Wxprime2,2.));
  cross.sigma=InclusiveCross::sigmaparam(cross.Wxprime2,kin.GetQsquared());
  result[0]*= offshellness*imag(cross.scatter(t)*2.*temp)*chi; //factor 2 due to symmetry
	
// 	if(isnan(result[0])) cout << offshellness << " " << chi << " " << cross.scatter(t) << " " << sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime)) << endl;
//   if(isnan(result[0])) cout << kin.GetS() << " " << kin2.GetS() << " " << pow(cross.massr+sqrt(cross.Wxprime2),2.) <<  " " << kin.GetS()*kin.GetS()-2.*kin.GetS()*(cross.Wxprime2+cross.massr*cross.massr)
//     +pow(cross.massr*cross.massr-cross.Wxprime2,2.) << endl;
	cross.sigma=InclusiveCross::sigmaparam(cross.otherWx2,kin.GetQsquared());
  double chi2=sqrt(kin.GetS()*kin.GetS()-2.*kin.GetS()*(cross.otherWx2+cross.massr*cross.massr)
    +pow(cross.massr*cross.massr-cross.otherWx2,2.));
  result[1]*=offshellness*imag(cross.scatter(t)*2.*conj(temp))*chi2; //factor 2 due to symmetry

}



void InclusiveCross::calc_F2DincFSI_Gauss(double &fsi1, double &fsi2, double Q2,double x, 
					   std::vector<double> & centrals, std::vector<double> & widths){
  
  fsi1=fsi2=0.;
  InclusiveCross::Ftor_FSI_distr F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.centrals=centrals;
  F.widths=widths;
  
  numint::mdfunction<numint::vector_d,5> mdf;
  mdf.func = &Ftor_FSI_distr::exec;
  mdf.param = &F;
  numint::vector_d ret(2,0.);
  F.f=InclusiveCross::FSI_int_Gauss;
  for(size_t it=0; it<resonances.size(); it++){
    numint::array<double,5> lower = {{centrals[it]-3.*widths[it],0.,-1.,0.,0.}};
    numint::array<double,5> upper = {{centrals[it]+3.*widths[it],1.E03,1.,1.E03,2.*PI}};
    F.it=it;
    F.it2=0;
  int res=90;
  unsigned count=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
    res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E07,ret,count,0);
  //   cout << res << " " << count << endl;
    fsi1+= 2.*PI*2.*massi/MASSD*ret[0];
    fsi2+= 2.*PI*2.*massi/MASSD*ret[1];
  }
  return;
}

void InclusiveCross::FSI_int_Gauss(numint::vector_d & result, double mass, double prnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, size_t itt, size_t it2,
		      std::vector<double> & centrals, std::vector<double> & widths){

  cross.setResonance(itt,mass);
  FSI_int(result,prnorm,costheta,qt,qphi,cross,Q2,x,itt,it2);
  double coef=exp(-(mass-centrals[itt])*(mass-centrals[itt])/(2.*widths[itt]*widths[itt]))/(widths[itt]*sqrt(2.*PI));
  result[0]*=coef;
  result[1]*=coef;
  return;
}


void InclusiveCross::calc_F2DincFSI_Gauss_off(double &fsi1, double &fsi2, double Q2,double x, 
						std::vector<double> & centrals, std::vector<double> & widths){
  fsi1=fsi2=0.;
  double Wmax=massi;
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double nu=Q2/(2.*massi*x);
  for(int i=0;i<1000;i++){
    double W=sqrt(-Q2+MASSD*MASSD+2.*MASSD*nu+massr*massr-2.*sqrt(massr*massr+i*i)*(nu+MASSD)+2.*qvec*i);
    if(W>Wmax) Wmax=W;
  }
//   cout << Wmax <<endl;
  
  InclusiveCross::Ftor_FSI_distr F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.centrals=centrals;
  F.widths=widths;
  
  numint::mdfunction<numint::vector_d,5> mdf;
  mdf.func = &Ftor_FSI_distr::exec;
  mdf.param = &F;
  numint::vector_d ret(2,0.);
  F.f=InclusiveCross::FSI_int_Gauss_off;
  for(size_t it=0; it<resonances.size(); it++){
    numint::array<double,5> lower = {{centrals[it]-3.*widths[it],0.,massi,0.,0.}};
    numint::array<double,5> upper = {{centrals[it]+3.*widths[it],1.E03,Wmax,1.E03,2.*PI}};
    F.it=it;
    /*for(size_t it2=0; it2<resonances.size(); it2++)*/{
      F.it2=it;
      int res=90;
      unsigned count=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E07,ret,count,0);
    //   cout << res << " " << count << endl;
      fsi1+= 2.*PI*2.*massi/MASSD*ret[0];
      fsi2+= 2.*PI*2.*massi/MASSD*ret[1];
    }
  }
  return;
}


void InclusiveCross::FSI_int_Gauss_off(numint::vector_d & result, double mass, double prt, double W, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x,size_t itt, size_t itt2,
		      std::vector<double> & centrals, std::vector<double> & widths){

  cross.setResonance(itt,mass);
  FSI_int_off(result,prt,W,qt,qphi,cross,Q2,x,itt,itt2);
  double coef=exp(-(mass-centrals[itt])*(mass-centrals[itt])/(2.*widths[itt]*widths[itt]))/(widths[itt]*sqrt(2.*PI));
  result[0]*=coef;
  result[1]*=coef;
  return;
  
}

void InclusiveCross::calc_F2DincFSI_uniform(double &fsi1, double &fsi2, double Q2,double x, 
					    std::vector<double> & centrals, std::vector<double> & widths){
  
  fsi1=fsi2=0.;
  InclusiveCross::Ftor_FSI_distr F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.centrals=centrals;
  F.widths=widths;
  
  numint::mdfunction<numint::vector_d,5> mdf;
  mdf.func = &Ftor_FSI_distr::exec;
  mdf.param = &F;
  numint::vector_d ret(2,0.);
  F.f=InclusiveCross::FSI_int_uniform;
  for(size_t it=0; it<resonances.size(); it++){
    numint::array<double,5> lower = {{centrals[it]-widths[it],0.,-1.,0.,0.}};
    numint::array<double,5> upper = {{centrals[it]+widths[it],1.E03,1.,1.E03,2.*PI}};
    F.it=it;
    F.it2=0;
    int res=90;
    unsigned count=0;
    //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
    res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E07,ret,count,0);
  //   cout << res << " " << count << endl;
    fsi1+= 2.*PI*2.*massi/MASSD*ret[0];
    fsi2+= 2.*PI*2.*massi/MASSD*ret[1];
  }
  return;
}

void InclusiveCross::FSI_int_uniform(numint::vector_d & result, double mass, double prnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, size_t itt, size_t it2,
		      std::vector<double> & centrals, std::vector<double> & widths){

  cross.setResonance(itt,mass);
  FSI_int(result,prnorm,costheta,qt,qphi,cross,Q2,x,itt,it2);
  double coef=1./(2.*widths[itt]);
  result[0]*=coef;
  result[1]*=coef;
  return;
}


void InclusiveCross::calc_F2DincFSI_uniform_off(double &fsi1, double &fsi2, double Q2,double x, 
					      std::vector<double> & centrals, std::vector<double> & widths){
  
  fsi1=fsi2=0.;
  double Wmax=massi;
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  double nu=Q2/(2.*massi*x);
  for(int i=0;i<1000;i++){
    double W=sqrt(-Q2+MASSD*MASSD+2.*MASSD*nu+massr*massr-2.*sqrt(massr*massr+i*i)*(nu+MASSD)+2.*qvec*i);
    if(W>Wmax) Wmax=W;
  }
//   cout << Wmax <<endl;
  
  InclusiveCross::Ftor_FSI_distr F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.centrals=centrals;
  F.widths=widths;
  
  numint::mdfunction<numint::vector_d,5> mdf;
  mdf.func = &Ftor_FSI_distr::exec;
  mdf.param = &F;
  numint::vector_d ret(2,0.);
  F.f=InclusiveCross::FSI_int_uniform_off;
  for(size_t it=0; it<resonances.size(); it++){
    numint::array<double,5> lower = {{centrals[it]-widths[it],0.,massi,0.,0.}};
    numint::array<double,5> upper = {{centrals[it]+widths[it],1.E03,Wmax,1.E03,2.*PI}};
    F.it=it;
    /*for(size_t it2=0; it2<resonances.size(); it2++)*/{
      F.it2=it;
      int res=90;
      unsigned count=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E07,ret,count,0);
    //   cout << res << " " << count << endl;
      fsi1+= 2.*PI*2.*massi/MASSD*ret[0];
      fsi2+= 2.*PI*2.*massi/MASSD*ret[1];
    }
  }
  return;
}


void InclusiveCross::FSI_int_uniform_off(numint::vector_d & result, double mass, double prt, double W, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x,size_t itt, size_t itt2,
		      std::vector<double> & centrals, std::vector<double> & widths){

  cross.setResonance(itt,mass);
  FSI_int_off(result,prt,W,qt,qphi,cross,Q2,x,itt,itt2);
  double coef=1./(2.*widths[itt]);
  result[0]*=coef;
  result[1]*=coef;
  return;
  
}

double InclusiveCross::sigmaparam(double W_sq, double Q2){
  /*cout << sqrt(W_sq) << " " << 65/10.*1.8E06/Q2*INVHBARC*INVHBARC << " " << (25.3*1.E-06*2.3+53*(sqrt(W_sq>5.76E06?5.76E06:W_sq)-MASSP)*1.E-03)
	/(1.E-05*Q2)*INVHBARC*INVHBARC << " " << endl;
  */if(abs(W_sq-1.232*1.232E06)<2.5E5) return 65/10.*1.8E06/Q2*INVHBARC*INVHBARC;
  return (25.3*1.E-06*2.3+53*(sqrt(W_sq>5.76E06?5.76E06:W_sq)-MASSP)*1.E-03)
	/(1.E-05*Q2)*INVHBARC*INVHBARC;
}



double InclusiveCross::get_prz_res(double pt2, double W_sq, TKinematics2to2 & kin){
//   double Er=sqrt(massr*massr+pt2);
//     przprime = (W_sq-massr*massr-MASSD*MASSD+kin.GetQsquared()-2.*MASSD*kin.GetWlab()
//     + 2.*Er*(MASSD+kin.GetWlab()))/(2.*kin.GetKlab());         
  double A=(W_sq-massr*massr-MASSD*MASSD+kin.GetQsquared()-2.*MASSD*kin.GetWlab())/(-2.*(MASSD+kin.GetWlab()));
  double B= kin.GetKlab()/(MASSD+kin.GetWlab());
  double aa=B*B-1.;
  double bb=2.*A*B;
  double cc=A*A-massr*massr-pt2;
  double discr=bb*bb-4.*aa*cc;
  int check=0;
  double result=-1.E09;
  if(abs(discr)<1.E-09) {result = -bb/(2.*aa);}
  if(discr>0.){
    double p1 = (-bb+sqrt(discr))/(2.*aa);
    double Er=sqrt(massr*massr+pt2+p1*p1);
    if(SIGN(A+B*p1)==1){check++; result=p1;}
    double p2 = (-bb-sqrt(discr))/(2.*aa);
    double Er2=sqrt(massr*massr+pt2+p2*p2);
    if(SIGN(A+B*p2)==1){check++; result=p2;}
//     cout << check << " " << A+B*p1 << " " << A+B*p2 << " " << Er << " " << Er2 << " " << p1 << " " << p2 << endl;
    if(check==2) result=p1; //other one is too big to give physical results....
  }
  return result;
}
double InclusiveCross::get_prz_res(double pt2, double W_sq, double Q2, double nu, double qvec){
  double A=(W_sq-massr*massr-MASSD*MASSD+Q2-2.*MASSD*nu)/(-2.*(MASSD+nu));
  double B= qvec/(MASSD+nu);
  double aa=B*B-1.;
  double bb=2.*A*B;
  double cc=A*A-massr*massr-pt2;
  double discr=bb*bb-4.*aa*cc;
  int check=0;
  double result=-1.E09;
  if(abs(discr)<1.E-09) {result = -bb/(2.*aa);}
  if(discr>0.){
    double p1 = (-bb+sqrt(discr))/(2.*aa);
    double Er=sqrt(massr*massr+pt2+p1*p1);
    if(SIGN(A+B*p1)==1){check++; result=p1;}
    double p2 = (-bb-sqrt(discr))/(2.*aa);
    double Er2=sqrt(massr*massr+pt2+p2*p2);
    if(SIGN(A+B*p2)==1){check++; result=p2;}
//     cout << check << " " << A+B*p1 << " " << A+B*p2 << " " << Er << " " << Er2 << " " << p1 << " " << p2 << endl;
    if(check==2) result=p1; //other one is too big to give physical results....
  }
  return result;
}




complex<double> InclusiveCross::scatter(double t){
  return sigma*(I_UNIT+epsilon)*exp(beta*t/2.); 
}


