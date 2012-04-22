#include "RhoTCross.hpp"
#include <Utilfunctions.hpp>


RhoTCross::RhoTCross(const int nucleus, const double p_max, const string dir):
homedir(dir),
pmax(p_max),
nucleusthick(nucleus,dir),
pdistgrid(NULL),
pfsigrid(NULL),
prho(NULL){

  
  
}

RhoTCross::~RhoTCross(){
  
}

//input in GeV!!!!
void RhoTCross::getCrosst(double *results, const double Ebeam, const double Q2, const double nu, const double t){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
  double pmestimate=0.,cthestimate=0.,phiestimate=0.;
  rombergerN(this,&RhoTCross::intPmt,0.,pmax*1.-03,NROFRES,
	      results,PREC,3,7,&pmestimate,Q2,nu,qvec,t,&cthestimate, &phiestimate);
  for(int i=0;i<NROFRES;i++) results[i]*= ALPHA*(Ebeam-nu)/(2.*pow(2.*PI,3.)*Ebeam*Q2*(1.-epsilon));
  
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
	    results,PREC,3,7,pcthestimate, pm, Q2,nu,qvec,t,pphiestimate);
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
	    results,PREC,3,5,pphiestimate,pm, costheta, sintheta, Q2,nu,qvec,t);
 
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
      if(abs(discr)<1.E-09){
	double pzrho=-b/(2.*a);
	double Erho = (pzrho*qvec-A)/nu;
	double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	if(Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	  if(!isnan(pxrho)){
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double costhetarho = pzrho/prho;
  	  double pnz = Cz-pzrho;
  	  double pnx = Cx-pxrho;
  	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
  	  double En = sqrt(mN*mN+pn*pn);
	  double EX = nu+nucleusthick.getMassA()*1.E-03-Erho;
	  double massX = sqrt(EX*EX-pxrho*pxrho-(qvec-pzrho)*(qvec-pzrho));
//   	  if(Erho/nu>0.9) cout << "zero " << pm << " " << costheta << " " << phi << " " << i << " " <<  
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << pn<< " " << t << 
//   	    " " << massX - nucleusthick.getMassA()*1.E-03 << endl;
	    double intresults[NROFRES];
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);	  
	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t)*prho;
	  }
	}
      }
      else{
	discr = sqrt(discr);
	double pzrho = (-b+discr)/(2.*a);
	double Erho = (pzrho*qvec-A)/nu;
	double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	if(Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	  if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double costhetarho = pzrho/prho;
	    double pnz = Cz-pzrho;
	    double pnx = Cx-pxrho;
	    double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	    double En = sqrt(mN*mN+pn*pn);
	    double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	    double EX = nu+nucleusthick.getMassA()*1.E-03-Erho;
	    double massX = sqrt(EX*EX-pxrho*pxrho-(qvec-pzrho)*(qvec-pzrho));
//   	  if(Erho/nu>0.9) cout << "bla1 " << pm << " " << costheta << " " << phi << " " << i << " " <<  
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << pn<< " " << t <<  
//   	    " " << massX - nucleusthick.getMassA()*1.E-03 << endl;
	    
	    double intresults[NROFRES];
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t)*prho;
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	Erho = (pzrho*qvec-A)/nu;
	t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	if(Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	  if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double costhetarho = pzrho/prho;
	    double pnz = Cz-pzrho;
	    double pnx = Cx-pxrho;
	    double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	    double En = sqrt(mN*mN+pn*pn);
	    double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	    double EX = nu+nucleusthick.getMassA()*1.E-03-Erho;
	    double massX = sqrt(EX*EX-pxrho*pxrho-(qvec-pzrho)*(qvec-pzrho));
//   	  if(Erho/nu>0.9) cout << "bla2 " << pm << " " << costheta << " " << phi << " " << i << " " <<  
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << pn<< " " << t <<  
//   	    " " << massX - nucleusthick.getMassA()*1.E-03 << endl;
	    double intresults[NROFRES];	  
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t)*prho;
	  }
	}
      }    
    }
  }
}

//input in GeV!!!!
void RhoTCross::getCrossz(double *results, const double Ebeam,  const double Q2, const double nu, const double z){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
  double Erho=nu*z;
  double prho=sqrt(Erho*Erho-MASSRHO*MASSRHO*1.E-06);
  double pmestimate=0.,cthestimate=0.,phiestimate=0.;
  rombergerN(this,&RhoTCross::intPmz,0.,pmax*1.E-03,NROFRES,
	      results,PREC,3,7,&pmestimate,Q2,nu,qvec,Erho,prho,&cthestimate, &phiestimate);
  
  for(int i=0;i<NROFRES;i++) results[i]*= ALPHA*(Ebeam-nu)*prho/(2.*pow(2.*PI,3.)*Ebeam*Q2*(1.-epsilon));
  
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
	    results,PREC,3,7,pcthestimate, pm, Q2,nu,qvec,Erho,prho,pphiestimate);
  cout << pm << " " ;
  for(int i=0;i<NROFRES;i++){
    results[i]*=pm*pm;
    cout << results[i] << " ";
  }
 cout << endl;
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
	    results,PREC,3,5,pphiestimate,pm, costheta, sintheta, Q2,nu,qvec,Erho,prho);
 
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
	if(prho>pzrho&&t<-0.1&&t>-0.4){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);      
	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
	  double pnn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	  double EX = nu+nucleusthick.getMassA()*1.E-03-Erho;
	  double massX = sqrt(EX*EX-pxrho*pxrho-(qvec-pzrho)*(qvec-pzrho));
// 	  if(t<-0.1&&t>-0.4) cout << "zero " << pm << " " << costheta << " " << phi << " " << i << " " <<  
// 	    Erho << " " << costhetarho << " " << Erho/nu << " " << t << " " << pnn<< " " << pn  <<
// 	    " " << massX - nucleusthick.getMassA()*1.E-03 << endl;
	  double intresults[NROFRES];
	  getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	  for(int dd=0;dd<NROFRES;dd++) results[dd]+=intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t);
	}
      }
      else{
	discr = sqrt(discr);
	double pzrho = (-b+discr)/(2.*a);
	double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	if(pzrho<prho&&t<-0.1&&t>-0.4){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);
	  if(SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho)){
	    double costhetarho = pzrho/prho;
	    double pnz = Cz-pzrho;
	    double pnx = Cx-pxrho;
	    double pnn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	    double EX = nu+nucleusthick.getMassA()*1.E-03-Erho;
	    double massX = sqrt(EX*EX-pxrho*pxrho-(qvec-pzrho)*(qvec-pzrho));
// 	    if(t<-0.1&&t>-0.4) cout << "bla1 " << pm << " " << costheta << " " << phi << " " << i << " " <<  
// 	      Erho << " " << costhetarho << " " << Erho/nu << " " << t << " " << pnn<< " " << pn  << 
// 	      " " << massX - nucleusthick.getMassA()*1.E-03 << endl;
	    double intresults[NROFRES];
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t);
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	if(pzrho<prho&&t<-0.1&&t>-0.4){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);	  
	  if(SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho)){
	    double costhetarho = pzrho/prho;
	    double pnz = Cz-pzrho;
	    double pnx = Cx-pxrho;
	    double pnn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	    double EX = nu+nucleusthick.getMassA()*1.E-03-Erho;
	    double massX = sqrt(EX*EX-pxrho*pxrho-(qvec-pzrho)*(qvec-pzrho));
// 	    if(t<-0.1&&t>-0.4) cout << "bla2 " << pm << " " << costheta << " " << phi << " " << i << " " <<  
// 	      Erho << " " << costhetarho << " " << Erho/nu << " " << t << " " << pnn<< " " << pn <<
// 	      " " << massX - nucleusthick.getMassA()*1.E-03 << endl;
	    double intresults[NROFRES];
	    getMomdistr(intresults,prho*1.E03,acos(costhetarho),Q2,i,pm*1.E03,costheta,phi);
	    for(int dd=0;dd<NROFRES;dd++) results[dd]+=intresults[dd]*getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t);
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
  pfsigrid = new GlauberDecayGridThick(60,20,5,&nucleusthick,homedir);
  pfsigrid->addParticle(rho);
  pfsigrid->fillGrids();
  pfsigrid->clearKnockout();
  pdistgrid = new DistMomDistrGrid(shell, pmax, 1,1,1,pfsigrid,homedir);
  for(int i=0;i<NROFRES-1;i++) results[i]+= pdistgrid->getRhoGridFull_interp3(i, pm, pmcostheta, pmphi);
  results[NROFRES-1] += pdistgrid->getRhopwGridFull_interp(pm);
  
  delete pdistgrid;
  delete pfsigrid;
  
}
  
  
double RhoTCross::getfrontfactor(double nu, double qvec, double Erho, double prho, double pzrho, double pxrho,
				 double s, double Q2, double mN, double t){
  
  double EX = nu+nucleusthick.getMassA()*1.E-03-Erho;
  double massX = sqrt(EX*EX-pxrho*pxrho-(qvec-pzrho)*(qvec-pzrho));
  
  return massX*(s*s-2.*s*(mN*mN-Q2)+pow(mN*mN+Q2,2.))/(mN*mN*abs(EX+Erho*(1.-qvec*pzrho/prho/prho)))*exp(t*6.);
  
}

  