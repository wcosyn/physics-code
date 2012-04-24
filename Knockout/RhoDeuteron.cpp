#include "RhoDeuteron.hpp"
#include <Utilfunctions.hpp>
#include <FastParticle.hpp>
#define MASS_N (MASSP+MASSN)*0.5E-03

RhoDeuteron::RhoDeuteron(const string name, const double p_max, const string dir):
pmax(p_max),
homedir(dir),
deuteron(name){

}

RhoDeuteron::~RhoDeuteron(){
}

//input in GeV!!!!
//cross section is obtained in units fm^3 GeV^4.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoDeuteron::getCrosst(double *result, const double Ebeam, const double Q2, const double nu, const double t){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  double pmestimate=0.,cthestimate=0.,phiestimate=0.;
  rombergerN(this,&RhoDeuteron::intPmt,0.,pmax*1.E-03,2,
	      result,PREC,3,7,&pmestimate,Q2,nu,qvec,t,&cthestimate, &phiestimate);
  for(int i=0;i<2;i++) result[i]*=ALPHA*(Ebeam-nu)/(pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon)*MASSD*1.E-03);
}

void RhoDeuteron::intPmt(const double pm, double *result, va_list ap){
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double t = va_arg(ap,double);
  double *pcthestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  
  if(pm==0.){
    for(int i=0;i<2;i++) result[i]=0.;
    return;
  }
  double EN = sqrt(MASS_N*MASS_N+pm*pm);
  rombergerN(this,&RhoDeuteron::intCosThetat,-1.,1.,2,
	    result,PREC,3,7,pcthestimate, pm, Q2,nu,qvec,t,EN,pphiestimate);
    for(int i=0;i<2;i++) result[i]*=pm*pm;
  
}

void RhoDeuteron::intCosThetat(const double costheta, double *result, va_list ap){
  double pm = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double t = va_arg(ap,double);
  double EN = va_arg(ap,double);
  double *pphiestimate = va_arg(ap,double*);
  
  double sintheta = sqrt(1.-costheta*costheta);
  rombergerN(this,&RhoDeuteron::intPhit,0.,2.*PI,2,
	    result,PREC,3,7,pphiestimate,pm, costheta, sintheta, Q2,nu,qvec,t,EN);
 
}

void RhoDeuteron::intPhit(const double phi, double *result, va_list ap){
  double pm = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double t = va_arg(ap,double);
  double EN = va_arg(ap,double);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TVector3 pmvec = TVector3(pm*sintheta*cosphi*1.E03,pm*sintheta*sinphi*1.E03,pm*costheta*1.E03);  
  
  for(int i=0;i<2;i++) result[i]=0.;
  
  //ALL IN GeV to avoid some numeric overflow almost zero shit!
  double A = (t+Q2-MASSRHO*MASSRHO*1.E-06)*0.5;
  double C0 = nu+MASSD*1.E-03-EN;
  double Cz = qvec-pm*costheta;
  double Cy = -pm*sintheta*sinphi;
  double Cx = -pm*sintheta*cosphi;
  double s = C0*C0-Cz*Cz-Cx*Cx-Cy*Cy; //mandelstam rhoN system
  double D = (s-MASS_N*MASS_N+MASSRHO*MASSRHO*1.E-06)/2.+C0*A/nu;
  double E = Cz-C0*qvec/nu;
  double a = E*E-Q2*Cx*Cx/(nu*nu);
  double b = 2.*(D*E+A*qvec*Cx*Cx/(nu*nu));
  double c = D*D-A*A*Cx*Cx/(nu*nu)+Cx*Cx*MASSRHO*MASSRHO*1.E-06;
  double discr = b*b-4.*a*c;
  if(discr>-1.E-09){
    if(abs(discr)<1.E-09){
      double pzrho=-b/(2.*a);
      double Erho = (pzrho*qvec-A)/nu;
      double tt = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(Erho/nu>0.9){
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(!isnan(pxrho)){
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	  double costhetarho = pzrho/prho;
	double pnz = Cz-pzrho;
	double pnx = Cx-pxrho;
	double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	double En = sqrt(MASS_N*MASS_N+pn*pn);	
//   	  /*if(Erho/nu>0.9)*/ cout << "zero " << pm << " " << costheta << " " << phi << " " << 
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
//   	    " " << C0 << " " << En+Erho << endl;
// 	  *result +=getMomdistr(prho*1.E03,acos(costhetarho),Q2,pm*1.E03,costheta,phi);	  
// 		    *getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En)*prho;
	double intresults[2];
	getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pzrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En)*prho;
	}
      }
    }
    else{
      discr = sqrt(discr);
      double pzrho = (-b+discr)/(2.*a);
      double Erho = (pzrho*qvec-A)/nu;
      double tt = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(Erho/nu>0.9){
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	  double En = sqrt(MASS_N*MASS_N+pn*pn);
	  
//  	  /*if(Erho/nu>0.9)*/ cout << "2nd " << pm << " " << costheta << " " << phi << " " << 
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
//   	    " " << C0 << " " << En+Erho << endl;
	  
	  double intresults[2];
	  getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pzrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En)*prho;
	}
      }
      pzrho = (-b-discr)/(2.*a);
      Erho = (pzrho*qvec-A)/nu;
      tt = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(Erho/nu>0.9){
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	  double En = sqrt(MASS_N*MASS_N+pn*pn);
//  	  /*if(Erho/nu>0.9)*/ cout << "3rd " << pm << " " << costheta << " " << phi << " " <<  
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
//   	    " " << C0 << " " << En+Erho << endl;
	  
	  double intresults[2];
	  getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pzrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En)*prho;
	}
      }
    }    
  }

}

//input in GeV!!!!
//cross section is obtained in units fm^3 GeV^4.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoDeuteron::getCrossz(double *result, const double Ebeam,  const double Q2, const double nu, const double z){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
  double Erho=nu*z;
  double prho=sqrt(Erho*Erho-MASSRHO*MASSRHO*1.E-06);
  double pmestimate=0.,cthestimate=0.,phiestimate=0.;
  rombergerN(this,&RhoDeuteron::intPmz,0.,pmax*1.E-03,2,
	      result,PREC,2,2,&pmestimate,Q2,nu,qvec,Erho,prho,&cthestimate, &phiestimate);
  
  for(int i=0;i<2;i++) result[i]*=ALPHA*(Ebeam-nu)*prho/(2.*pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon));
  
}

void RhoDeuteron::intPmz(const double pm, double *result, va_list ap){
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double Erho = va_arg(ap,double);
  double prho = va_arg(ap,double);
  double *pcthestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  if(pm==0.){
    for(int i=0;i<2;i++) result[i]=0.;
    return;
  }
  double EN = sqrt(MASS_N*MASS_N+pm*pm);
  rombergerN(this,&RhoDeuteron::intCosThetaz,-1.,1.,2,
	    result,PREC,2,2,pcthestimate, pm, Q2,nu,qvec,Erho,prho,EN,pphiestimate);
  for(int i=0;i<2;i++){
    result[i]*=pm*pm;
  }
  cout << pm << " " << result[0] << " " << result[1] << endl;

}

void RhoDeuteron::intCosThetaz(const double costheta, double *result, va_list ap){
  double pm = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double Erho = va_arg(ap,double);
  double prho = va_arg(ap,double);
  double EN = va_arg(ap,double);
  double *pphiestimate = va_arg(ap,double*);
  
  double sintheta = sqrt(1.-costheta*costheta);
  rombergerN(this,&RhoDeuteron::intPhiz,0.,2.*PI,2,
	    result,PREC,2,2,pphiestimate,pm, costheta, sintheta, Q2,nu,qvec,Erho,prho,EN);
 
}

void RhoDeuteron::intPhiz(const double phi, double *result, va_list ap){
  double pm = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double Erho = va_arg(ap,double);
  double prho = va_arg(ap,double);
  double EN = va_arg(ap,double);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TVector3 pmvec = TVector3(pm*sintheta*cosphi*1.E03,pm*sintheta*sinphi*1.E03,pm*costheta*1.E03);
  for(int i=0;i<2;i++) result[i]=0.;  
  
  double C0 = nu+MASSD*1.E-03-EN;
  double Cz = qvec-pm*costheta;
  double Cy = -pm*sintheta*sinphi;
  double Cx = -pm*sintheta*cosphi;
  double s = C0*C0-Cz*Cz-Cx*Cx-Cy*Cy; //mandelstam rhoN
  double En = C0-Erho;
  double pn = sqrt(En*En-MASS_N*MASS_N);
  
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
	double Enn = sqrt(MASS_N*MASS_N+pnx*pnx+pnz*pnz+Cy*Cy);
// 	/*if(t<-0.1&&t>-0.4)*/ cout << "zero " << pm << " " << costheta << " " << phi << " " <<
// 	  Erho << " " << costhetarho << " " << Erho/nu << " " << t << " " << En<< " " << Enn  << endl;
// 	  *result +=getMomdistr(prho*1.E03,acos(costhetarho),Q2,pm*1.E03,costheta,phi);	  
// 		    *getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En);
	double intresults[2];
	getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pzrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En);

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
	  double Enn = sqrt(MASS_N*MASS_N+pnx*pnx+pnz*pnz+Cy*Cy);
// 	  /*if(t<-0.1&&t>-0.4)*/ cout << "2nd " << pm << " " << costheta << " " << phi << " " <<  
// 	    Erho << " " << costhetarho << " " << Erho/nu << " " << t << " " << En<< " " << Enn  << endl;
// 	  *result +=getMomdistr(prho*1.E03,acos(costhetarho),Q2,pm*1.E03,costheta,phi);	  
// 		    *getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En);
	double intresults[2];
	getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pzrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En);
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
	  double Enn = sqrt(MASS_N*MASS_N+pnx*pnx+pnz*pnz+Cy*Cy);
// 	 /* if(t<-0.1&&t>-0.4)*/ cout << "2nd " << pm << " " << costheta << " " << phi << " " << 
// 	    Erho << " " << costhetarho << " " << Erho/nu << " " << t << " " << En<< " " << Enn  << endl;
// 	  *result +=getMomdistr(prho*1.E03,acos(costhetarho),Q2,pm*1.E03,costheta,phi);	  
// 		    *getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En);

	  double intresults[2];
	  getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pzrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En);
	}
      }
    }
  }
//   cout << pm << " " << discr << " end ";
//   for(int dd=0;dd<2;dd++) cout << result[dd] << " ";
//   cout << endl;
}

//all energies in MeV please
void RhoDeuteron::getMomdistr(double *results, double Erho, double prho, double pzrho, 
			      double pxrho, TVector3 &pmvec, double nu, double qvec){
  
  double s = MASS_N*MASS_N*1.E06+MASSRHO*MASSRHO+2.*Erho*sqrt(MASS_N*MASS_N*1.E06+pmvec.Mag2())-2.*(pzrho*pmvec.Z()+pxrho*pmvec.X());
  double sigmap, beta2p, epsp, sigman, beta2n, epsn;
  FastParticle::setPionGlauberData(prho,sigmap,beta2p,epsp,sigman,beta2n,epsn,homedir);
  deuteron.setScatter(0.5*10.*(sigmap+sigman),0.5*(beta2n+beta2p)*INVHBARC*INVHBARC*1.E06,0.5*(epsn+epsp));//careful with units!
  /*results[1] = */results[0] = deuteron.getMomDistrpw(pmvec)*pow(HBARC,3.);
  results[1] = deuteron.getMomDistrfsi(pmvec, nu, qvec, s, MASSRHO)*pow(HBARC,3.);
  //cout << pmvec.Mag() << " " << nu << " " << qvec << " " << s << " " << results[0] << " " << results[1] << endl;
}  
  
  
double RhoDeuteron::getfrontfactor(double qvec, double Erho, double prho, double pzrho, double s, double Q2, double t, double En){
  
  
  return (s*s-2.*s*(MASS_N*MASS_N-Q2)+pow(MASS_N*MASS_N+Q2,2.))/(abs(En+Erho*(1.-qvec*pzrho/prho/prho)))*exp(t*6.);
  
}

  