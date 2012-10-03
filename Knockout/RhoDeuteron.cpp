#include "RhoDeuteron.hpp"
#include <Utilfunctions.hpp>
#include <FastParticle.hpp>
#define MASS_N (MASSP+MASSN)*0.5E-03

using namespace std;


RhoDeuteron::RhoDeuteron(const string name, const double p_max,const bool no_cuts, const int integr,
	      const double precision, const int max_Eval, const bool FSI):
pmax(p_max),
nocuts(no_cuts),
deuteron(name),
prec(precision),
integrator(integr),
abserror(1.E-012),
maxEval(max_Eval),
fsi(FSI){
}

RhoDeuteron::~RhoDeuteron(){
}

//input in GeV!!!!
//cross section is obtained in units GeV.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoDeuteron::getCrosst(double *result, const double Ebeam, const double Q2, const double nu, const double t){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  int res=90;
  unsigned count=0;
  if(integrator==0){
    double pmestimate=0.,cthestimate=0.,phiestimate=0.;
    rombergerN(this,&RhoDeuteron::intPmt,0.,pmax*1.E-03,2,
		result,PREC,3,7,&pmestimate,Q2,nu,qvec,t,&cthestimate, &phiestimate);
  }
  else if(integrator==1||integrator==2){
    numint::vector_d ret(2,0.);
    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{pmax*1.E-03,1.,PI}};
    
    RhoDeuteron::Ftor_rhoD F;
    F.deuteron = this;
    F.Q2=Q2;
    F.nu=nu;
    F.qvec=qvec;
    F.t=t;
    F.Erho=0.;
    F.prho=0.;
    numint::mdfunction<numint::vector_d,3> mdf;
    mdf.func = &Ftor_rhoD::exec;
    mdf.param = &F;
    F.f=klaas_rhoD_t;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,abserror,1.E-03,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,abserror,1.E-03,maxEval,ret,count,0);
    for(int i=0;i<2;++i) result[i]=ret[i];
    if(ret[1]*1.E-04>abserror) {abserror=ret[1]*1.E-04; }
  }
    
    
  for(int i=0;i<2;i++) result[i]*=ALPHA*(Ebeam-nu)/(pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon)*MASSD*1.E-03)*1.E09/(2.*qvec);
  cout << nu << " " << t << " " << result[0] << " " << result[1] << " " << res << " " << count << " " << abserror << " " << endl;
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
    for(int i=0;i<2;i++) result[i]*=2.*pm*pm;
  
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
  rombergerN(this,&RhoDeuteron::intPhit,0.,PI,2,
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
  double Eoff = MASSD*1.E-03-EN;
  
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
      if(nocuts||Erho/nu>0.9){
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(!isnan(pxrho)){
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
// 	  double costhetarho = pzrho/prho;
	double pnz = Cz-pzrho;
	double pnx = Cx-pxrho;
// 	double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	double En = sqrt(MASS_N*MASS_N+pnx*pnx+pnz*pnz+Cy*Cy);	
//   	   cout << "zero " << pm << " " << costheta << " " << phi << " " << 
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
//   	    " " << C0 << " " << En+Erho << endl;
// 	  *result +=getMomdistr(prho*1.E03,acos(costhetarho),Q2,pm*1.E03,costheta,phi);	  
// 		    *getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En)*prho;
	double intresults[2];
	getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,1);
	}
      }
    }
    else{
      discr = sqrt(discr);
      double pzrho = (-b+discr)/(2.*a);
      double Erho = (pzrho*qvec-A)/nu;
      double tt = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(nocuts||Erho/nu>0.9){
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
// 	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
// 	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	double En = sqrt(MASS_N*MASS_N+pnx*pnx+pnz*pnz+Cy*Cy);	
	  
//  	  cout << "2nd " << pm << " " << costheta << " " << phi << " " << 
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
//   	    " " << C0 << " " << En+Erho << endl;
	  
	  double intresults[2];
	  getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,1);
	}
      }
      pzrho = (-b-discr)/(2.*a);
      Erho = (pzrho*qvec-A)/nu;
      tt = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(nocuts||Erho/nu>0.9){
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(SIGN(D+E*pzrho)==SIGN(-Cx*pxrho)){
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
// 	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
// 	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	double En = sqrt(MASS_N*MASS_N+pnx*pnx+pnz*pnz+Cy*Cy);	
//  	   cout << "3rd " << pm << " " << costheta << " " << phi << " " <<  
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
//   	    " " << C0 << " " << En+Erho << endl;
	  
	  double intresults[2];
	  getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,1);
	}
      }
    }    
  }

}

//input in GeV!!!!
//cross section is obtained in units GeV.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoDeuteron::getCrossz(double *result, const double Ebeam,  const double Q2, const double nu, const double z){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  double Erho=nu*z;
  double prho=sqrt(Erho*Erho-MASSRHO*MASSRHO*1.E-06);
  int res=90;
  unsigned count=0;
  if(integrator==0){
    double pmestimate=0.,cthestimate=0.,phiestimate=0.;
    rombergerN(this,&RhoDeuteron::intPmz,0.,pmax*1.E-03,2,
		result,PREC,3,7,&pmestimate,Q2,nu,qvec,Erho,prho,&cthestimate, &phiestimate);
   }
  else if(integrator==1||integrator==2){
    numint::vector_d ret(2,0.);
    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{pmax*1.E-03,1.,PI}};
    
    RhoDeuteron::Ftor_rhoD F;
    F.deuteron = this;
    F.Q2=Q2;
    F.nu=nu;
    F.qvec=qvec;
    F.t=0.;
    F.Erho=Erho;
    F.prho=prho;
    numint::mdfunction<numint::vector_d,3> mdf;
    mdf.func = &Ftor_rhoD::exec;
    mdf.param = &F;
    F.f=klaas_rhoD_z;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,abserror,1.E-03,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,abserror,1.E-03,maxEval,ret,count,0);
    for(int i=0;i<2;++i) result[i]=ret[i];
    if(ret[1]*1.E-04>abserror)  {abserror=ret[1]*1.E-04; }
  }
  
 
  for(int i=0;i<2;i++) result[i]*=ALPHA*(Ebeam-nu)/(pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon)*MASSD*1.E-03)*1.E09*nu/qvec;
   cout << nu << " " << z << " " << result[0] << " " << result[1] << " " << res << " " << count << " " << abserror << " " << endl;
 
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
	    result,PREC,3,7,pcthestimate, pm, Q2,nu,qvec,Erho,prho,EN,pphiestimate);
  for(int i=0;i<2;i++){
    result[i]*=2.*pm*pm;  //2. because of symmetry in phi
  }

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
  rombergerN(this,&RhoDeuteron::intPhiz,0.,PI,2,
	    result,PREC,3,7,pphiestimate,pm, costheta, sintheta, Q2,nu,qvec,Erho,prho,EN);

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
  double Eoff = MASSD*1.E-03-EN;
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
      if(prho>pzrho&&(nocuts||(t<-0.1&&t>-0.4))){
	double pxrho = sqrt(prho*prho-pzrho*pzrho);      
	double intresults[2];
	getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,0);

      }
    }
    else{
      discr = sqrt(discr);
      double pzrho = (-b+discr)/(2.*a);
      double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(pzrho<prho&&(nocuts||(t<-0.1&&t>-0.4))){
	double pxrho = sqrt(prho*prho-pzrho*pzrho);
	if(SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho)){
	double intresults[2];
	getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,0);
	}
      }
      pzrho = (-b-discr)/(2.*a);
      t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(pzrho<prho&&(nocuts||(t<-0.1&&t>-0.4))){
	double pxrho = sqrt(prho*prho-pzrho*pzrho);	  
	if(SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho)){
	  double intresults[2];
	  getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) result[dd]+=intresults[dd]*getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,0);
	}
      }
    }
  }
}

//all energies in MeV please
void RhoDeuteron::getMomdistr(double *results, double Erho, double prho, double pzrho, 
			      double pxrho, TVector3 &pmvec, double nu, double qvec){
  
  double s = MASS_N*MASS_N*1.E06+MASSRHO*MASSRHO+2.*Erho*sqrt(MASS_N*MASS_N*1.E06+pmvec.Mag2())-2.*(pzrho*pmvec.Z()+pxrho*pmvec.X());
  double sigmap, beta2p, epsp, sigman, beta2n, epsn;
  FastParticle::interpPionGlauberData(4,prho,sigmap,beta2p,epsp,sigman,beta2n,epsn);
  deuteron.setScatter(10.*sigmap,beta2p*INVHBARC*INVHBARC*1.E06,epsp);//careful with units!
  results[0] = deuteron.getMomDistrpw(pmvec);
  results[1] = fsi? deuteron.getMomDistrfsi(pmvec, nu, qvec, s, MASSRHO) : results[0];
}  
  
  
double RhoDeuteron::getfrontfactor(double qvec, double Erho, double prho, double pzrho, double s, double Q2, double t, 
				   double En, double Eoff, bool torz){
  
  
  return (s*s-2.*s*(MASS_N*MASS_N-Q2)+pow(MASS_N*MASS_N+Q2,2.))*exp(t*6.)/(torz? Eoff:1.);
  
}

//  input in GeV!!!!
//cross section is obtained in units GeV.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoDeuteron::getCrosst_coh(double *result, const double Ebeam, const double Q2, const double nu, const double t){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
  //ALL IN GeV to avoid some numeric overflow almost zero shit!
  double A = (t+Q2-MASSRHO*MASSRHO*1.E-06)*0.5;
  double C0 = nu+MASSD*1.E-03;
  double Cz = qvec;
  double s = C0*C0-Cz*Cz; //mandelstam rhoD system
  double D = (s-MASSD*MASSD*1.E-06+MASSRHO*MASSRHO*1.E-06)/2.+C0*A/nu;
  double E = Cz-C0*qvec/nu;
  double pzrho=-D/E;
  double Erho = (pzrho*qvec-A)/nu;
//   double tt = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
  if(Erho/nu>0.9){
    double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
    if(isnan(pxrho)) pxrho=0.;
    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
    double costhetarho = pzrho/prho;
    double pDz = Cz-pzrho;
    double pDx = -pxrho;
    double ss = -Q2+MASS_N*MASS_N+2.*MASS_N*nu;
    double ED = sqrt(MASSD*MASSD*1.E-06+pDx*pDx+pDz*pDz);	
//   	  /*if(Erho/nu>0.9)*/ cout << "zero " <<
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
//   	    " " << C0 << " " << ED+Erho << endl;
      TVector3 pDvec = TVector3(pDx*1.E-03,0.,pDx*1.E-03);
//       *result = 1./(2.*MASSD*1.E-03);
//       return;
      *result = calcCross_coh(pDvec)
		*ALPHA*(Ebeam-nu)*prho*ss*exp(6*t)/(pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon)*MASSD*1.E-03*abs(ED+Erho*(1-qvec*pzrho/prho/prho)));
//       cout << nu << " " << t << " " << *result << endl;
    
  }
  else *result=0.;
}

void RhoDeuteron::getCrossz_coh(double *result, const double Ebeam,  const double Q2, const double nu, const double z){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
  double Erho=nu*z;
  double prho=sqrt(Erho*Erho-MASSRHO*MASSRHO*1.E-06);
  double ED = nu+MASSD*1.E-03-Erho;
  double pD = sqrt(ED*ED-MASSD*MASSD*1.E-06);
  double pzrho = (qvec*qvec+prho*prho-pD*pD)/(2.*qvec);
  double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
  if(t<-0.1&&t>-0.4)
  {
    double pxrho = sqrt(prho*prho-pzrho*pzrho);
    if(isnan(pxrho)) pxrho=0.;
    double ss = -Q2+MASS_N*MASS_N+2.*MASS_N*nu;
    TVector3 pDvec = TVector3(-pxrho*1.E03,0.,(qvec-pzrho)*1.E03);
//     cout << t << " " << pD*pD << " " << pxrho*pxrho+(qvec-pzrho)*(qvec-pzrho) << " " << pzrho/prho << endl;
//     *result = nu;
//     return;
    *result = calcCross_coh(pDvec)
    *ALPHA*(Ebeam-nu)*prho*ss*exp(6*t)/(pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon)*MASSD*1.E-03*abs(ED+Erho*(1-qvec*pzrho/prho/prho)));    
//     cout << nu << " " << z << " " << *result << endl;
  }
  else *result=0.;
}

double RhoDeuteron::calcCross_coh(TVector3 &pDvec){
  double result=0.;
  double pmestimate=0.,cthestimate=0.,phiestimate=0.;
  //because of symmetry between sd and -sd for both incoming and outgoing deuteron we can simplify the summation
  complex<double> intresult;
  rombergerN(this,&RhoDeuteron::intPmcoh,0.,pmax,1,
	  &intresult,PREC,3,7,&pmestimate,&pDvec,-2,-2,&cthestimate, &phiestimate);
  result+=2.*norm(intresult);
  rombergerN(this,&RhoDeuteron::intPmcoh,0.,pmax,1,
	  &intresult,PREC,3,7,&pmestimate,&pDvec,-2,2,&cthestimate, &phiestimate);
  result+=2.*norm(intresult);
  rombergerN(this,&RhoDeuteron::intPmcoh,0.,pmax,1,
	  &intresult,PREC,3,7,&pmestimate,&pDvec,-2,0,&cthestimate, &phiestimate);
  result+=4.*norm(intresult);
  rombergerN(this,&RhoDeuteron::intPmcoh,0.,pmax,1,
	  &intresult,PREC,3,7,&pmestimate,&pDvec,0,0,&cthestimate, &phiestimate);
  result+=norm(intresult);
  return result*4./3.*pow(HBARC*1.E-03,3.);
}


void RhoDeuteron::intPmcoh(const double pm, complex<double> *result, va_list ap){
  TVector3 *p_pDvec = va_arg(ap,TVector3*);
  int M = va_arg(ap,int);
  int M2 = va_arg(ap,int);
  double *pcthestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  
  if(pm==0.){
    for(int i=0;i<2;i++) result[i]=0.;
    return;
  }
  rombergerN(this,&RhoDeuteron::intCosThetacoh,-1.,1.,1,
	    result,PREC,3,7,pcthestimate, pm, p_pDvec,M,M2,pphiestimate);
  *result*=pm*pm;
  
}

void RhoDeuteron::intCosThetacoh(const double costheta, complex<double> *result, va_list ap){
  double pm = va_arg(ap,double);
  TVector3 *p_pDvec = va_arg(ap,TVector3*);
  int M = va_arg(ap,int);
  int M2 = va_arg(ap,int);
  double *pphiestimate = va_arg(ap,double*);
  
  double sintheta = sqrt(1.-costheta*costheta);
  rombergerN(this,&RhoDeuteron::intPhicoh,0.,2.*PI,1,
	    result,PREC,3,7,pphiestimate,pm, costheta, sintheta, p_pDvec,M,M2);
 
}

void RhoDeuteron::intPhicoh(const double phi, complex<double> *result, va_list ap){
  double pm = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  TVector3 *p_pDvec = va_arg(ap,TVector3*);
  int M = va_arg(ap,int);
  int M2 = va_arg(ap,int);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TVector3 ps = TVector3(pm*sintheta*cosphi,pm*sintheta*sinphi,pm*costheta);
  TVector3 ps2 = TVector3((ps-*p_pDvec)*0.5);
  
  *result=deuteron.getMomDistrpwcoh(ps,ps2,M,M2);
}  

void RhoDeuteron::klaas_rhoD_t(numint::vector_d & results, double pm, double costheta, double phi, RhoDeuteron & deuteron, double Q2, double nu, double qvec,
    double t, double dummy, double dummy2){
  results=numint::vector_d(2,0.);
  if(pm==0.){
    return;
  }
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  double sintheta = sqrt(1.-costheta*costheta);

  double EN = sqrt(MASS_N*MASS_N+pm*pm);
  TVector3 pmvec = TVector3(pm*sintheta*cosphi*1.E03,pm*sintheta*sinphi*1.E03,pm*costheta*1.E03);  
  double Eoff = MASSD*1.E-03-EN;
  
  double A = (t+Q2-MASSRHO*MASSRHO*1.E-06)*0.5;
  double C0 = nu+Eoff;
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
      //zero discriminant
    if(abs(discr)<1.E-09){
      double pzrho=-b/(2.*a);
      double Erho = (pzrho*qvec-A)/nu;
      
      double tt = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      
      if(deuteron.getNocuts()||Erho/nu>0.9){
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(isnan(pxrho)) pxrho=0.; //underflow
	if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
  // 	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
  // 	double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	  double En = sqrt(MASS_N*MASS_N+pnx*pnx+pnz*pnz+Cy*Cy);	
  //   	   cout << "zero " << pm << " " << costheta << " " << phi << " " << 
  //   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
  //   	    " " << C0 << " " << En+Erho << endl;
  // 	  *result +=deuteron.getMomdistr(prho*1.E03,acos(costhetarho),Q2,pm*1.E03,costheta,phi);	  
  // 		    *deuteron.getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En)*prho;
	  double front=deuteron.getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,1);
	  double intresults[2];
	  deuteron.getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) results[dd]+=intresults[dd]*front;	  
	}
      }
    }
    else{
      discr = sqrt(discr);
      double pzrho = (-b+discr)/(2.*a);
      double Erho = (pzrho*qvec-A)/nu;
      double tt = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
      if(deuteron.getNocuts()||Erho/nu>0.9){
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(isnan(pxrho)) pxrho=0.;//underflow
	if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
// 	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
// 	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	double En = sqrt(MASS_N*MASS_N+pnx*pnx+pnz*pnz+Cy*Cy);	
	  
//  	  cout << "2nd " << pm << " " << costheta << " " << phi << " " << 
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
//   	    " " << C0 << " " << En+Erho << endl;
	  
	  double front=deuteron.getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,1);
	  double intresults[2];
	  deuteron.getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) results[dd]+=intresults[dd]*front;	  
	}
      }
      pzrho = (-b-discr)/(2.*a);
      Erho = (pzrho*qvec-A)/nu;
      tt = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(deuteron.getNocuts()||Erho/nu>0.9){
	double pxrho = sqrt(Erho*Erho-pzrho*pzrho-MASSRHO*MASSRHO*1.E-06);
	if(isnan(pxrho)) pxrho=0.;//underflow
	if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	  double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
// 	  double costhetarho = pzrho/prho;
	  double pnz = Cz-pzrho;
	  double pnx = Cx-pxrho;
// 	  double pn = sqrt(pnx*pnx+pnz*pnz+Cy*Cy);
	double En = sqrt(MASS_N*MASS_N+pnx*pnx+pnz*pnz+Cy*Cy);	
//  	   cout << "3rd " << pm << " " << costheta << " " << phi << " " <<  
//   	    Erho << " " << costhetarho << " " << Erho/nu << " " << tt << 
//   	    " " << C0 << " " << En+Erho << endl;
	  
	  double front=deuteron.getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,1);
	  double intresults[2];
	  deuteron.getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) results[dd]+=intresults[dd]*front;	  
	}
      }
    }    
  }
  
  for(int i=0;i<2;i++) results[i]*=2.*pm*pm;
}


void RhoDeuteron::klaas_rhoD_z(numint::vector_d & results, double pm, double costheta, double phi, RhoDeuteron & deuteron, double Q2, double nu, double qvec,
    double dummy, double Erho, double prho){
  results=numint::vector_d(2,0.);
  if(pm==0.){
    return;
  }
  double sintheta = sqrt(1.-costheta*costheta);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  
  double EN = sqrt(MASS_N*MASS_N+pm*pm);
  TVector3 pmvec = TVector3(pm*sintheta*cosphi*1.E03,pm*sintheta*sinphi*1.E03,pm*costheta*1.E03);  
  double Eoff = MASSD*1.E-03-EN;

  double C0 = nu+Eoff;
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
      if(deuteron.getNocuts()||(t<-0.1&&t>-0.4)){
	double pxrho = sqrt(prho*prho-pzrho*pzrho);//pxrho defines x-axis (so only positive sqrt is considered)           
	if(isnan(pxrho)) pxrho=0.;//underflow
	if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	  double front=deuteron.getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,0);
	  double intresults[2];
	  deuteron.getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) results[dd]+=intresults[dd]*front;	  
	}
      }
    }
    else{
      discr = sqrt(discr);
      double pzrho = (-b+discr)/(2.*a);
      double t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(deuteron.getNocuts()||(t<-0.1&&t>-0.4)){
	double pxrho = sqrt(prho*prho-pzrho*pzrho);
	if(isnan(pxrho)) pxrho=0.;//underflow
	if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	  double front=deuteron.getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,0);
	  double intresults[2];
	  deuteron.getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) results[dd]+=intresults[dd]*front;	  
	}
      }
      pzrho = (-b-discr)/(2.*a);
      t =  -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*pzrho;
      if(deuteron.getNocuts()||(t<-0.1&&t>-0.4)){
	double pxrho = sqrt(prho*prho-pzrho*pzrho);	  
	if(isnan(pxrho)) pxrho=0.;//underflow
	if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	  double front=deuteron.getfrontfactor(qvec,Erho,prho,pzrho,s,Q2,t,En,Eoff,0);
	  double intresults[2];
	  deuteron.getMomdistr(intresults,Erho*1.E03,prho*1.E03,pzrho*1.E03,pxrho*1.E03,pmvec,nu*1.E03,qvec*1.E03);	  
	  for(int dd=0;dd<2;dd++) results[dd]+=intresults[dd]*front;	  
	}
      }
    }
  }
  
  for(int i=0;i<2;i++) results[i]*=2.*pm*pm;
  
}