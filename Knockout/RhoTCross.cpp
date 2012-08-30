#include "RhoTCross.hpp"
#include <Utilfunctions.hpp>


RhoTCross::RhoTCross(const int nucleus, const double p_max, const string dir, const bool no_cuts,
  const bool user_set, const double user_sigma, const double precision
):
homedir(dir),
pmax(p_max),
nocuts(no_cuts),
userset(user_set),
usersigma(user_sigma),
nucleusthick(nucleus,dir),
pdistgrid(NULL),
pfsigrid(NULL),
prec(precision){

  pfsigrid = new GlauberDecayGridThick*[nucleusthick.getTotalLevels()];
  pdistgrid = new DistMomDistrGrid*[nucleusthick.getTotalLevels()];
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
    pfsigrid[i] = new GlauberDecayGridThick(60,20,5,&nucleusthick,getPrec(),homedir);
    pdistgrid[i] = new DistMomDistrGrid(i, pmax, 30,20,5,pfsigrid[i],getPrec(),homedir);
  }
}

RhoTCross::~RhoTCross(){
  
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
    delete pdistgrid[i];
    delete pfsigrid[i];
  }
  delete [] pdistgrid;
  delete [] pfsigrid;
}

//input in GeV!!!!
//cross section is obtained in units GeV.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoTCross::getCrosst(double *results, const double Ebeam, const double Q2, const double nu, const double t){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
  double pmestimate=0.,cthestimate=0.,phiestimate=0.;
  rombergerN(this,&RhoTCross::intPmt,0.,pmax*1.E-03,NROFRES,
	      results,getPrec(),3,7,&pmestimate,Q2,nu,qvec,t,&cthestimate, &phiestimate);
  for(int i=0;i<NROFRES;i++) results[i]*= ALPHA*(Ebeam-nu)/(2.*pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon))*pow(INVHBARC*1.E03,3.)/nucleusthick.getA();
  
}

void RhoTCross::intPmt(const double pm, double *results, va_list ap){
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double t = va_arg(ap,double);
  double *pcthestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  
  if(pm==0.){
    for(int i=0;i<NROFRES;i++) results[i]=0.;
    return;
  }

  rombergerN(this,&RhoTCross::intCosThetat,-1.,1.,NROFRES,
	    results,getPrec(),3,10,pcthestimate, pm, Q2,nu,qvec,t,pphiestimate);
  for(int i=0;i<NROFRES;i++){
    results[i]*=pm*pm;
  }
 
}

void RhoTCross::intCosThetat(const double costheta, double *results, va_list ap){
  double pm = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double t = va_arg(ap,double);
  double *pphiestimate = va_arg(ap,double*);
  
  double sintheta = sqrt(1.-costheta*costheta);
  rombergerN(this,&RhoTCross::intPhit,0.,2.*PI,NROFRES,
	    results,getPrec(),3,10,pphiestimate,pm, costheta, sintheta, Q2,nu,qvec,t);
 
}

void RhoTCross::intPhit(const double phi, double *results, va_list ap){
  double pm = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double t = va_arg(ap,double);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  
  for(int i=0;i<NROFRES;i++) results[i]=0.;
  
  //ALL IN GeV to avoid some numeric overflow almost zero shit!
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
    double massi = nucleusthick.getMassA()-nucleusthick.getMassA_min_1(i)-nucleusthick.getExcitation()[i];
    massi*=1.E-03;
    double mN = i<nucleusthick.getPLevels()? MASSP:MASSN;
    mN*=1.E-03;
    //determine kinematics
    double A = (t+Q2-MASSRHO*MASSRHO*1.E-06)*0.5;
    double C0 = nu+massi;
    double Cz = qvec+pm*costheta;
    double Cy = pm*sintheta*sinphi;
    double Cx = pm*sintheta*cosphi;
    double s = C0*C0-Cz*Cz-Cx*Cx-Cy*Cy; //mandelstam 
    double D = (s-mN*mN+MASSRHO*MASSRHO*1.E-06)/2.+C0*A/nu;
    double E = Cz-C0*qvec/nu;
    double a = E*E-Q2*Cx*Cx/(nu*nu);
    double b = 2.*(D*E+A*qvec*Cx*Cx/(nu*nu));
    double c = D*D-A*A*Cx*Cx/(nu*nu)+Cx*Cx*MASSRHO*MASSRHO*1.E-06;
    double discr = b*b-4.*a*c;
    if(discr>-1.E-09){
      //zero discriminant
      if(abs(discr)<1.E-09){
	double pzrho=-b/(2.*a);
	double Erho = (pzrho*qvec-A)/nu;
	//check cuts
	if(nocuts||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	  if(!isnan(pxrho)){
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double costhetarho = pzrho/prho;
	    double intresults[NROFRES];
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);	  
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t)*prho;
	  }
	}
      }
      else{
	discr = sqrt(discr);
	double pzrho = (-b+discr)/(2.*a);
	double Erho = (pzrho*qvec-A)/nu;
	double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(nocuts||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	  if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double costhetarho = pzrho/prho;
	    
	    double intresults[NROFRES];
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t)*prho;
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	Erho = (pzrho*qvec-A)/nu;
	t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(nocuts||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	  if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double costhetarho = pzrho/prho;
	    double intresults[NROFRES];	  
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t)*prho;
	  }
	}
      }    
    }
  }
}

//input in GeV!!!!
//cross section is obtained in units GeV.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoTCross::getCrossz(double *results, const double Ebeam,  const double Q2, const double nu, const double z){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
  double Erho=nu*z;
  double prho=sqrt(Erho*Erho-MASSRHO*MASSRHO*1.E-06);
  double pmestimate=0.,cthestimate=0.,phiestimate=0.;
  rombergerN(this,&RhoTCross::intPmz,0.,pmax*1.E-03,NROFRES,
	      results,getPrec(),3,10,&pmestimate,Q2,nu,qvec,Erho,prho,&cthestimate, &phiestimate);
  
  for(int i=0;i<NROFRES;i++) results[i]*= ALPHA*(Ebeam-nu)*prho/(2.*pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon))*pow(INVHBARC*1.E03,3.)/nucleusthick.getA();
  
}

void RhoTCross::intPmz(const double pm, double *results, va_list ap){
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double Erho = va_arg(ap,double);
  double prho = va_arg(ap,double);
  double *pcthestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  if(pm==0.){
    for(int i=0;i<NROFRES;i++) results[i]=0.;
    return;
  }
  rombergerN(this,&RhoTCross::intCosThetaz,-1.,1.,NROFRES,
	    results,getPrec(),3,8,pcthestimate, pm, Q2,nu,qvec,Erho,prho,pphiestimate);
  for(int i=0;i<NROFRES;i++){
    results[i]*=pm*pm;
  }
}

void RhoTCross::intCosThetaz(const double costheta, double *results, va_list ap){
  double pm = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double Erho = va_arg(ap,double);
  double prho = va_arg(ap,double);
  double *pphiestimate = va_arg(ap,double*);
  
  double sintheta = sqrt(1.-costheta*costheta);
  rombergerN(this,&RhoTCross::intPhiz,0.,2.*PI,NROFRES,
	    results,getPrec(),3,8,pphiestimate,pm, costheta, sintheta, Q2,nu,qvec,Erho,prho);
 
}

void RhoTCross::intPhiz(const double phi, double *results, va_list ap){
  double pm = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double Erho = va_arg(ap,double);
  double prho = va_arg(ap,double);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  for(int i=0;i<NROFRES;i++) results[i]=0.;  
  
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
    //offshell initial energy
    double massi = nucleusthick.getMassA()-nucleusthick.getMassA_min_1(i)-nucleusthick.getExcitation()[i];
    massi*=1.E-03;
    double mN = i<nucleusthick.getPLevels()? MASSP:MASSN;
    mN*=1.E-03;
    double C0 = nu+massi;
    double Cz = qvec+pm*costheta;
    double Cy = pm*sintheta*sinphi;
    double Cx = pm*sintheta*cosphi;
    double s = C0*C0-Cz*Cz-Cx*Cx-Cy*Cy; //mandelstam 
    double En = C0-Erho;
    double pn = sqrt(En*En-massi*massi);
    
    double A=(Cx*Cx+Cy*Cy+Cz*Cz-pn*pn+prho*prho)*0.5;
    double a = Cz*Cz+Cx*Cx;
    double b = -2.*A*Cz;
    double c = A*A-Cx*Cx*prho*prho;
    double discr = b*b-4.*a*c;
    if(discr>-1.E-09){
      if(abs(discr)<1.E-09) {
	double pzrho = -b/(2.*a);
	double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(prho>pzrho&&(nocuts||(t<-0.1&&t>-0.4))){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);      
	  double costhetarho = pzrho/prho;
	  double intresults[NROFRES];
	  getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	  double front=1.;
	  //correct for not completely full shells
	  if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	  if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	  for(int dd=0;dd<NROFRES;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t);

	}
      }
      else{
	discr = sqrt(discr);
	double pzrho = (-b+discr)/(2.*a);
	double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(pzrho<prho&&(nocuts||(t<-0.1&&t>-0.4))){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);
	  if(SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho)){
	    double costhetarho = pzrho/prho;
	    double intresults[NROFRES];
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);

	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t);
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(pzrho<prho&&(nocuts||(t<-0.1&&t>-0.4))){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);	  
	  if(SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho)){
	    double costhetarho = pzrho/prho;
	    double intresults[NROFRES];
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t);

	  }
	}
      }
    }
  }
}

//all energies in MeV please
void RhoTCross::getMomdistr(double *results, double prho, double thetarho, double Q2, int shell, 
			    double pm, double pmcostheta, double pmphi){
  FastParticle rho(4, 0, prho,thetarho,0.,Q2,145.,homedir);
  if(userset) rho.setScatter(usersigma,6.,-0.2);
  pfsigrid[shell]->clearParticles();
  pfsigrid[shell]->addParticle(rho);
  pfsigrid[shell]->updateGrids();
  pfsigrid[shell]->clearKnockout();
  //update the grid if necessary
  pdistgrid[shell]->updateGrids(pfsigrid[shell],shell);
  //plane-wave
  results[NROFRES-1] = pdistgrid[shell]->getRhopwGridFull_interp(pm);
  //fsi
  for(int i=0;i<NROFRES-1;i++) results[i]= pdistgrid[shell]->getRhoGridFull_interp3(i, pm, pmcostheta, pmphi);
    
}
  
  
double RhoTCross::getfrontfactor(double nu, double qvec, double Erho, double prho, double pzrho, double pxrho,
				 double s, double Q2, double mN, double t){
  //X is residual nucleus system
  double EX = nu+nucleusthick.getMassA()*1.E-03-Erho;
  double massX = sqrt(EX*EX-pxrho*pxrho-(qvec-pzrho)*(qvec-pzrho));
  //elementary rho production cross section parametrized as exponential e^{beta*t}
  return massX*(s*s-2.*s*(mN*mN-Q2)+pow(mN*mN+Q2,2.))/(mN*mN*abs(EX+Erho*(1.-qvec*pzrho/prho/prho)))*exp(t*6.);
  
}

  