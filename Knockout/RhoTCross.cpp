#include "RhoTCross.hpp"
#include <Utilfunctions.hpp>
#include <TRotation.h>

using namespace std;

RhoTCross::RhoTCross(const int nucleus, const double p_max, const string dir, const bool no_cuts,
  const bool user_set, const double user_sigma, const double precision, const int integr, const int max_Eval,
  const double p_dil, const int particletype
):
homedir(dir),
pmax(p_max),
nocuts(no_cuts),
userset(user_set),
usersigma(user_sigma),
nucleusthick(nucleus,dir),
pdistgrid(NULL),
pfsigrid(NULL),
prec(precision),
integrator(integr),
abserror(1.E-012),
maxEval(max_Eval),
pdil(p_dil),
particle_type(particletype){

  pfsigrid = new GlauberDecayGridThick*[nucleusthick.getTotalLevels()];
//   pdistgrid = new DistMomDistrGrid*[nucleusthick.getTotalLevels()];
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
    pfsigrid[i] = new GlauberDecayGridThick(60,20,5,&nucleusthick,prec,2,homedir);
//     pdistgrid[i] = new DistMomDistrGrid(i, pmax, 30,20,5,pfsigrid[i],1.E-03,2,2E04,0.,homedir);
  }
  nrofcross=pfsigrid[0]->getNumber_of_grids()+1;
}

RhoTCross::~RhoTCross(){
  
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
//     delete pdistgrid[i];
    delete pfsigrid[i];
  }
//   delete [] pdistgrid;
  delete [] pfsigrid;
}

//input in GeV!!!!
//cross section is obtained in units GeV.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoTCross::getCrosst(double *results, const double Ebeam, const double Q2, const double nu, const double t){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
//   double A = (t+Q2-MASSRHO*MASSRHO*1.E-06)*0.5;
//   double MA=getNucleusthick().getMassA()*1.E-03;
//   double C0=nu+MA;
//   double pzrho=(MA*MA+qvec*qvec-C0*C0-MASSRHO*MASSRHO*1.E-06-2.*C0*A/nu)/(2.*qvec*(1.-C0/nu));
//   double Erho = (pzrho*qvec-A)/nu;
//   double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06); //pxrho defines x-axis (so only positive sqrt is considered)
//   double prho = sqrt(Erho*Erho-MASSRHO*MASSRHO*1.E-06); //pxrho defines x-axis (so only positive sqrt is considered)
//   if(isnan(pxrho)) pxrho=0.;
//   double pAz=qvec-pzrho;
//   double z=Erho/nu;
//   if(z>1.||z<0.9){for(int i=0;i<nrofcross;i++) results[i]=0.;}
//   //else for(int i=0;i<nrofcross;i++){ results[i]=1./(2.*qvec*abs(C0-Erho+Erho*(1.-qvec*pzrho/prho/prho)));}
//   else { for(int i=0;i<nrofcross;i++) results[i]=1./(2.*MA);}
//   return;
  
  
  int res=90;
  unsigned count=0;
  if(integrator==0){
    double pmestimate=0.,cthestimate=0.,phiestimate=0.;
    rombergerN(this,&RhoTCross::intPmt,0.,pmax*1.E-03,nrofcross,
		results,getPrec(),3,7,&pmestimate,Q2,nu,qvec,t,&cthestimate, &phiestimate);
  }
  else if(integrator==1||integrator==2){
    numint::vector_d ret(nrofcross,0.);
    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{pmax*1.E-03,1.,2.*PI}};
    
    RhoTCross::Ftor_rho F;
    F.cross = this;
    F.Q2=Q2;
    F.nu=nu;
    F.qvec=qvec;
    F.t=t;
    F.Erho=0.;
    F.prho=0.;
    numint::mdfunction<numint::vector_d,3> mdf;
    mdf.func = &Ftor_rho::exec;
    mdf.param = &F;
    F.f=klaas_rho_t;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,abserror,1.E-03,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,abserror,1.E-03,maxEval,ret,count,0);
    for(int i=0;i<nrofcross;++i) results[i]=ret[i];
    if(ret[2]*1.E-04>abserror) {abserror=ret[2]*1.E-04; }
  }
  for(int i=0;i<nrofcross;i++) results[i]*= ALPHA*(Ebeam-nu)/(2.*pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon))*pow(INVHBARC*1.E03,3.)/nucleusthick.getA()/(2.*qvec);
//   cout << nu << " " << t << " " << results[0] << " " << res << " " << count << " " << abserror << " " << endl;
}

void RhoTCross::intPmt(const double pm, double *results, va_list ap){
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double t = va_arg(ap,double);
  double *pcthestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  
  if(pm==0.){
    for(int i=0;i<nrofcross;i++) results[i]=0.;
    return;
  }

  rombergerN(this,&RhoTCross::intCosThetat,-1.,1.,nrofcross,
	    results,getPrec(),3,10,pcthestimate, pm, Q2,nu,qvec,t,pphiestimate);
  for(int i=0;i<nrofcross;i++){
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
  rombergerN(this,&RhoTCross::intPhit,0.,2.*PI,nrofcross,
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
  
  for(int i=0;i<nrofcross;i++) results[i]=0.;
  
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
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06); //pxrho defines x-axis (so only positive sqrt is considered)
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double intresults[nrofcross];
	    getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);	  
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	    for(int dd=0;dd<nrofcross;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,1);
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
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);//pxrho defines x-axis (so only positive sqrt is considered)
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    
	    double intresults[nrofcross];
	    getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	    for(int dd=0;dd<nrofcross;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,1);
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	Erho = (pzrho*qvec-A)/nu;
	t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(nocuts||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);//pxrho defines x-axis (so only positive sqrt is considered)
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double intresults[nrofcross];	  
	    getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	    for(int dd=0;dd<nrofcross;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,1);
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
  
//   double MA=getNucleusthick().getMassA()*1.E-03;
//   double EA=-Erho+MA+nu;
//   double pA=sqrt(EA*EA-MA*MA);
//   if(isnan(pA)) pA=0.;
//   double pzrho = (prho*prho+qvec*qvec-pA*pA)/(2.*qvec);
//   double pxrho=sqrt(prho*prho-pzrho*pzrho);
//   if(isnan(pxrho)) pxrho=0.;
//   double pxA=-pxrho;
//   double pzA=qvec-pzrho;
//   double t = -Q2+MASSRHO*MASSRHO*1.E-06-2.*Erho*nu+2.*qvec*pzrho;
//   if(t>-0.1||t<-0.4){for(int i=0;i<nrofcross;i++) results[i]=0.;}
//   else for(int i=0;i<nrofcross;i++){ results[i]=nu;}
//   return;
  
  int res=90;
  unsigned count=0;
  if(integrator==0){
    double pmestimate=0.,cthestimate=0.,phiestimate=0.;
    rombergerN(this,&RhoTCross::intPmz,0.,pmax*1.E-03,nrofcross,
		results,getPrec(),3,10,&pmestimate,Q2,nu,qvec,Erho,prho,&cthestimate, &phiestimate);
  }
  else if(integrator==1||integrator==2){
    numint::vector_d ret(nrofcross,0.);
    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{pmax*1.E-03,1.,2.*PI}};
    
    RhoTCross::Ftor_rho F;
    F.cross = this;
    F.Q2=Q2;
    F.nu=nu;
    F.qvec=qvec;
    F.t=0.;
    F.Erho=Erho;
    F.prho=prho;
    numint::mdfunction<numint::vector_d,3> mdf;
    mdf.func = &Ftor_rho::exec;
    mdf.param = &F;
    F.f=klaas_rho_z;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,abserror,1.E-03,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,abserror,1.E-03,maxEval,ret,count,0);
    for(int i=0;i<nrofcross;++i) results[i]=ret[i];
    if(ret[2]*1.E-04>abserror)  {abserror=ret[2]*1.E-04; }
  }
  
  for(int i=0;i<nrofcross;i++) results[i]*= ALPHA*(Ebeam-nu)/(2.*pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon))*pow(INVHBARC*1.E03,3.)/nucleusthick.getA()*nu/qvec;
//   cout << nu << " " << z << " " << results[0] << " " << res << " " << count << " " << abserror << endl;
  
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
    for(int i=0;i<nrofcross;i++) results[i]=0.;
    return;
  }
  rombergerN(this,&RhoTCross::intCosThetaz,-1.,1.,nrofcross,
	    results,getPrec(),3,8,pcthestimate, pm, Q2,nu,qvec,Erho,prho,pphiestimate);
  for(int i=0;i<nrofcross;i++){
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
  rombergerN(this,&RhoTCross::intPhiz,0.,2.*PI,nrofcross,
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
  for(int i=0;i<nrofcross;i++) results[i]=0.;  
  
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
    double pn = sqrt(En*En-mN*mN);
    
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
	if(nocuts||(t<-0.1&&t>-0.4)){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho); //pxrho defines x-axis (so only positive sqrt is considered)     
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	    double intresults[nrofcross];
	    getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	    for(int dd=0;dd<nrofcross;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,0);
	  }
	}
      }
      else{
	discr = sqrt(discr);
	double pzrho = (-b+discr)/(2.*a);
	double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(nocuts||(t<-0.1&&t>-0.4)){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);//pxrho defines x-axis (so only positive sqrt is considered)  
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	    double intresults[nrofcross];
	    getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);

	    for(int dd=0;dd<nrofcross;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,0);
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(nocuts||(t<-0.1&&t>-0.4)){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);//pxrho defines x-axis (so only positive sqrt is considered)  	  
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	    double intresults[nrofcross];
	    getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=1.;
	    //correct for not completely full shells
	    if(i==nucleusthick.getPLevels()-1) front=double((nucleusthick.getFinalMProton()+1)*2-nucleusthick.getOnlyOneProton())/(nucleusthick.getJ_array()[i]+1);
	    if(i==nucleusthick.getTotalLevels()-1) front=double((nucleusthick.getFinalMNeutron()+1)*2-nucleusthick.getOnlyOneNeutron())/(nucleusthick.getJ_array()[i]+1);
	    for(int dd=0;dd<nrofcross;dd++) results[dd]+=front*intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,0);
	  }
	}
      }
    }
  }
}

//all energies in MeV please
void RhoTCross::getMomdistr(double *results, double prho, double thetarho, double Q2, int shell, 
			    double pm, double pmcostheta, double pmphi){
  //if userset there's only one possible glaubergrid, so we don't have to check every time!
  if(!(userset&&pfsigrid[shell]->getFilledallgrid())){
    FastParticle rho(particle_type, 0, prho,0.,0.,Q2,145.,homedir);
    if(userset) rho.setScatter(usersigma,6.,-0.2,pdil);
    pfsigrid[shell]->clearParticles();
    pfsigrid[shell]->addParticle(rho);
    pfsigrid[shell]->updateGrids();
    pfsigrid[shell]->clearKnockout();
  }
  //update the grid if necessary
  TRotation rot;
  rot.Rotate(-thetarho,TVector3(0.,1.,0.));
  double key=double(shell)/getNucleusthick().getTotalLevels()+int(round(thetarho*10.));
  map<double,DistMomDistrGrid>::iterator it=distgridmap.find(key);
//   cout << key << " " << thetarho << " " << shell  << endl;
  if(it==distgridmap.end()){
    distgridmap[key]=DistMomDistrGrid(shell, pmax, 30,20,5,pfsigrid[shell],1.E-03,2,5E04,0.,homedir);
//     distgridmap.insert(pair<double,DistMomDistrGrid>(key,DistMomDistrGrid(shell, pmax, 30,20,5,pfsigrid[shell],1.E-03,2,2E04,0.,homedir)));
    distgridmap[key].updateGrids(pfsigrid[shell],shell,rot);
    distgridmap[key].printRho_grid(0);
    cout << endl << endl;
    it=distgridmap.find(key);
  }
  else distgridmap[key].updateGrids(pfsigrid[shell],shell,rot);

  results[nrofcross-1] = distgridmap[key].getRhopwGridFull_interp(pm);
  //fsi
  for(int i=0;i<nrofcross-1;i++) results[i] = distgridmap[key].getRhoGridFull_interp3(i, pm, pmcostheta, pmphi);
//    if(pm>200.&&results[0]>50.) {distgridmap[key].printRho_grid(0); cout  << endl << endl << endl;}
  
  
//   pdistgrid[shell]->updateGrids(pfsigrid[shell],shell,rot);
//   //plane-wave
//   results[nrofcross-1] = pdistgrid[shell]->getRhopwGridFull_interp(pm);
//   //fsi
//   for(int i=0;i<nrofcross-1;i++) results[i]= pdistgrid[shell]->getRhoGridFull_interp3(i, pm, pmcostheta, pmphi);
}
  
  
double RhoTCross::getfrontfactor(double nu, double qvec, double Erho, double prho, double pzrho, double pxrho,
				 double s, double Q2, double mN, double t, bool torz){
  //X is residual nucleus system
  double EX = nu+nucleusthick.getMassA()*1.E-03-Erho;
  double massX = sqrt(EX*EX-pxrho*pxrho-(qvec-pzrho)*(qvec-pzrho));
  //elementary rho production cross section parametrized as exponential e^{beta*t}
  return (torz? 1.:massX)*(s*s-2.*s*(mN*mN-Q2)+pow(mN*mN+Q2,2.))/(mN*mN)*exp(t*6.);
  
}

void RhoTCross::klaas_rho_t(numint::vector_d & results, double pm, double costheta, double phi, RhoTCross & cross, double Q2, double nu, double qvec,
    double t, double dummy, double dummy2){
  results=numint::vector_d(cross.getNrofcross(),0.);
  if(pm==0.){
    return;
  }
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  double sintheta = sqrt(1.-costheta*costheta);
  
  
  //ALL IN GeV to avoid some numeric overflow almost zero shit!
  for(int i=0;i<cross.getNucleusthick().getTotalLevels();i++){
    double massi = cross.getNucleusthick().getMassA()-cross.getNucleusthick().getMassA_min_1(i)-cross.getNucleusthick().getExcitation()[i];
    massi*=1.E-03;
    double mN = i<cross.getNucleusthick().getPLevels()? MASSP:MASSN;
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
	if(cross.getNocuts()||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(isnan(pxrho)) pxrho=0.; //underflow
	  if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);	  
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,1);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
      }
      else{
	discr = sqrt(discr);
	double pzrho = (-b+discr)/(2.*a);
	double Erho = (pzrho*qvec-A)/nu;
	//check cuts
	if(cross.getNocuts()||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,1);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	Erho = (pzrho*qvec-A)/nu;
	//check cuts
	if(cross.getNocuts()||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double intresults[cross.getNrofcross()];	  
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,1);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
      }    
    }
  }
  for(int i=0;i<cross.getNrofcross();i++){
    results[i]*=pm*pm;
  }
}

void RhoTCross::klaas_rho_z(numint::vector_d & results, double pm, double costheta, double phi, RhoTCross & cross, double Q2, double nu, double qvec,
    double dummy, double Erho, double prho){
  results=numint::vector_d(cross.getNrofcross(),0.);
  if(pm==0.){
    return;
  }
  double sintheta = sqrt(1.-costheta*costheta);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  
  for(int i=0;i<cross.getNucleusthick().getTotalLevels();i++){
    //offshell initial energy
    double massi = cross.getNucleusthick().getMassA()-cross.getNucleusthick().getMassA_min_1(i)-cross.getNucleusthick().getExcitation()[i];
    massi*=1.E-03;
    double mN = i<cross.getNucleusthick().getPLevels()? MASSP:MASSN;
    mN*=1.E-03;
    double C0 = nu+massi;
    double Cz = qvec+pm*costheta;
    double Cy = pm*sintheta*sinphi;
    double Cx = pm*sintheta*cosphi;
    double s = C0*C0-Cz*Cz-Cx*Cx-Cy*Cy; //mandelstam 
    double En = C0-Erho;
    double pn = sqrt(En*En-mN*mN);
    
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
	if(cross.getNocuts()||(t<-0.1&&t>-0.4)){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);//pxrho defines x-axis (so only positive sqrt is considered)           
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,0);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
      }
      else{
	discr = sqrt(discr);
	double pzrho = (-b+discr)/(2.*a);
	double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(cross.getNocuts()||(t<-0.1&&t>-0.4)){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 

	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,0);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(cross.getNocuts()||(t<-0.1&&t>-0.4)){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	  
	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,0);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
      }
    }
  }
  
  for(int i=0;i<cross.getNrofcross();i++){
    results[i]*=pm*pm;
  }
}
  