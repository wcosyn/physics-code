#include "RhoTCross.hpp"
#include <Utilfunctions.hpp>


RhoTCross::RhoTCross(const int nucleus, const double p_max, const string dir):
homedir(dir),
pmax(p_max),
nucleusthick(nucleus,dir),
pdistgrid(NULL),
pfsigrid(NULL),
pkin(NULL),
prho(NULL){

  
  
}

RhoTCross::~RhoTCross(){
  
}


void RhoTCross::getCrosst(double *results, const double Q2, const double nu, const double t){

  double qvec = sqrt(Q2+nu*nu);
  double pmestimate=0.,cthestimate=0.,phiestimate=0.;
  rombergerN(this,&RhoTCross::intPmt,0.,pmax,NROFRES,
	      results,PREC,3,7,&pmestimate,Q2,nu,qvec,t,&cthestimate, &phiestimate);
  
  
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
	    results,PREC,3,7,pphiestimate,pm, costheta, sintheta, Q2,nu,qvec,t);
 
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
    double D = (C0*C0-Cz*Cz-Cx*Cx-Cy*Cy-mN*mN+MASSRHO*MASSRHO*1.E-06)/2.+C0*A/nu;
    double E = Cz-C0*qvec/nu;
    double a = E*E-Q2*Cx*Cx/(nu*nu);
    double b = 2.*(D*E+A*qvec*Cx*Cx/(nu*nu));
    double c = D*D-A*A*Cx*Cx/(nu*nu)+Cx*Cx*MASSRHO*MASSRHO*1.E-06;
    double discr = b*b-4.*a*c;
    if(discr<0.){
    }
    else{
      if(abs(discr)<1.E-09){
	double pzrho=-b/(2.*a);
	double Erho = (pzrho*qvec-A)/nu;
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(isnan(pxrho)){
	}
	else{
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	  double En = sqrt(mN*mN+pn*pn);
	  double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	  if(Erho/nu>0.9) cout << "zero " << pm << " " << costheta << " " << phi << " " << i << " " <<  Erho << " " << costhetarho << " " << Erho/nu << " " << pn<< " " << t << " " << Erho+En << " " << C0 << endl;
	}
      }
      else{
	int kak=0;
	discr = sqrt(discr);
	if(isnan(discr)) cout << "tralalal " << b*b-4.*a*c << endl;
	double pzrho = (-b+discr)/(2.*a);
	double Erho = (pzrho*qvec-A)/nu;
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	  kak++;
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	  double En = sqrt(mN*mN+pn*pn);
	  double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	  if(Erho/nu>0.9) cout << "bla1 " << pm << " " << costheta << " " << phi << " " << i << " " <<  Erho << " " << costhetarho << " " << Erho/nu << " " << pn<< " " << t << " " << Erho+En << " " << C0 << endl;
	}
	pzrho = (-b-discr)/(2.*a);
	Erho = (pzrho*qvec-A)/nu;
	pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	  kak++;
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	  double En = sqrt(mN*mN+pn*pn);
	  double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	  if(Erho/nu>0.9) cout << "bla2 " << pm << " " << costheta << " " << phi << " " << i << " " <<  Erho << " " << costhetarho << " " << Erho/nu << " " << pn<< " " << t << " " << Erho+En << " " << C0 << endl;
	}
      }    
    }
  }
  for(int i=0;i<NROFRES;i++) results[i]=1.;
}


void RhoTCross::getCrossz(double *results, const double Q2, const double nu, const double z){

  double qvec = sqrt(Q2+nu*nu);
  double Erho=nu*z;
  double prho=sqrt(Erho*Erho-MASSRHO*MASSRHO*1.E-06);
  double pmestimate=0.,cthestimate=0.,phiestimate=0.;
  rombergerN(this,&RhoTCross::intPmz,0.,pmax,NROFRES,
	      results,PREC,3,7,&pmestimate,Q2,nu,qvec,Erho,prho,&cthestimate, &phiestimate);
  
  
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
	    results,PREC,3,7,pphiestimate,pm, costheta, sintheta, Q2,nu,qvec,Erho,prho);
 
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
  for(int i=0;i<NROFRES;i++) results[i]=1.;  
  
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
    double massi = nucleusthick.getMassA()-nucleusthick.getMassA_min_1(i)-nucleusthick.getExcitation()[i];
    massi*=1.E-03;
    double mN = i<nucleusthick.getPLevels()? MASSP:MASSN;
    mN*=1.E-03;
    double C0 = nu+massi;
    double Cz = qvec+pm*costheta;
    double Cy = pm*sintheta*sinphi;
    double Cx = pm*sintheta*cosphi;
    double En = C0-Erho;
    double pn = sqrt(En*En-massi*massi);
    
    double A=(Cx*Cx+Cy*Cy+Cz*Cz-pn*pn+prho*prho)*0.5;
    double a = Cz*Cz+Cx*Cx;
    double b = -2.*A*Cz;
    double c = A*A-Cx*Cx*prho*prho;
    double discr = b*b-4.*a*c;
    if(abs(discr)<1.E-09) {
      double pzrho = -b/(2.*a);
      if(prho>pzrho){
	double pxrho = sqrt(prho*prho-pzrho*pzrho);      
	double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	double costhetarho = pzrho/prho;
	double pnz = Cz-pzrho;
	double pnx = Cx-pxrho;
	double pnn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	if(t<-0.1&&t>-0.4) cout << "zero " << pm << " " << costheta << " " << phi << " " << i << " " <<  Erho << " " << costhetarho << " " << Erho/nu << " " << t << " " << pnn<< " " << pn  << " " << A-Cz*pzrho << " " << pxrho*Cx <<  endl;
      }
    }
    else{
      discr = sqrt(discr);
      if(isnan(discr)) {
	for(int i=0;i<NROFRES;i++) results[i]=1.;  
	return;
      }
      double pzrho = (-b+discr)/(2.*a);
      double pxrho = sqrt(prho*prho-pzrho*pzrho);
      if(pzrho<prho){
  //       if(isnan(pxrho)) cout << "tttttt " << pzrho << " " << prho << endl;
	if(SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho)){
	  double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
	  double pnn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	  if(t<-0.1&&t>-0.4) cout << "bla1 " << pm << " " << costheta << " " << phi << " " << i << " " <<  Erho << " " << costhetarho << " " << Erho/nu << " " << t << " " << pnn<< " " << pn  << " " << A-Cz*pzrho << " " << pxrho*Cx << endl;

	}
      }
      pzrho = (-b-discr)/(2.*a);
      pxrho = sqrt(prho*prho-pzrho*pzrho);
      if(pzrho<prho){
  // 	if(isnan(pxrho)) cout << "tttttt " << pzrho << " " << prho << endl;
	if(SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho)){
	  double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
	  double pnn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	   if(t<-0.1&&t>-0.4) cout << "bla2 " << pm << " " << costheta << " " << phi << " " << i << " " <<  Erho << " " << costhetarho << " " << Erho/nu << " " << t << " " << pnn<< " " << pn << " " << A-Cz*pzrho << " " << pxrho*Cx << endl;

	}
      }
    }        
  }
  
}

