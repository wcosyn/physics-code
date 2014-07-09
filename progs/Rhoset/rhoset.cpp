#include <iostream>
#include <cstdlib>

using namespace std;

#include <RhoTCross.hpp>
#include <Utilfunctions.hpp>
#include <string>

void tzint(double tz, double* result, va_list ap);
void nuint(double nu, double* result, va_list ap);
void klaas_t(numint::vector_d & res, double nu, double t, RhoTCross & cross, double Q2, double Ebeam){
  res=numint::vector_d(cross.getNrofcross(),0.);
  double results[cross.getNrofcross()];
  cross.getCrosst(results,Ebeam,Q2,nu,t);
  for(int i=0;i<cross.getNrofcross();++i) res[i]=results[i];

//   double qvec=sqrt(Q2+nu*nu);
//   double massA = cross.getNucleusthick().getMassA();
//   massA*=1.E-03;
//   //determine kinematics
//   double A = (t+Q2-MASSRHO*MASSRHO*1.E-06)*0.5;
//   double B=(2.*massA*nu+MASSRHO*MASSRHO*1.E-06-Q2)/2.+A*(massA+nu)/nu;
//   double prz=B/((massA/nu)*qvec);
//   double Erho=(prz*qvec-A)/nu;
//   double prx=sqrt(Erho*Erho-prz*prz-MASSRHO*MASSRHO*1.E-06);
//   double EA=sqrt(massA*massA+prx*prx+(qvec-prz)*(qvec-prz));
//   double z=Erho/nu;
//   if((!std::isnan(prx))&&(z>0.9&&z<1.)){
// //     cout << EA+Erho << " " << massA+nu << endl;
//     for(int i=0;i<cross.getNrofcross();++i) res[i]=0.5*Erho/qvec;    
//   }
//   else for(int i=0;i<cross.getNrofcross();++i) res[i]=0.;
  
}
void klaas_z(numint::vector_d & res, double nu, double z, RhoTCross & cross, double Q2, double Ebeam){
  res=numint::vector_d(cross.getNrofcross(),0.);
  double results[cross.getNrofcross()];
  cross.getCrossz(results,Ebeam,Q2,nu,z);
  for(int i=0;i<cross.getNrofcross();++i) res[i]=results[i];

  
  //   double qvec=sqrt(Q2+nu*nu);
//   double massA = cross.getNucleusthick().getMassA();
//   massA*=1.E-03;
//   double Erho=z*nu;
//   double EA=nu+massA-Erho;
//   double prho=sqrt(Erho*Erho-MASSRHO*MASSRHO*1.E-06);
//   double pA=sqrt(EA*EA-massA*massA);
//   double prhoz=(prho*prho+qvec*qvec-pA*pA)/(2.*qvec);
//   double pAz=(-prho*prho+qvec*qvec+pA*pA)/(2.*qvec);
//   double prhox=sqrt(prho*prho-prhoz*prhoz);
//   double t = -Q2 + MASSRHO*MASSRHO*1.E-06-2*Erho*nu+2.*qvec*prhoz;
//   if((!std::isnan(prhox))&&(t<-0.1&&t>-0.4)) {
//     for(int i=0;i<cross.getNrofcross();++i) res[i]=nu*massA*Erho/qvec;
//     if(std::isnan(prhox)) cout << "Blaaa" << endl;
// //   cout << pA << " " << sqrt(prhox*prhox+pAz*pAz) << endl;
//   }
//   else for(int i=0;i<cross.getNrofcross();++i) res[i]=0.;
}

int main(int argc, char *argv[])
{

  int torz = atoi(argv[1]);
  int nucleus = atoi(argv[2]);
  double Q2 = atof(argv[3]);
  double nu_min = atof(argv[4]);
  double nu_max = atof(argv[5]);
  double p_dil = atof(argv[6]);
  int nocuts = atoi(argv[7]);
  double usersigma = atof(argv[8]);
  int integrator = atoi(argv[9]);
  int maxEval = atoi(argv[10]);
  double prec=atof(argv[11]);
  int ptype=atoi(argv[12]);
//   int maxEval2 = atoi(argv[11]);
  string homedir=argv[13];
  double Ebeam = 5.014;
  

  RhoTCross test = RhoTCross(nucleus,330,homedir,nocuts,1,usersigma,prec,2,maxEval,p_dil,ptype);  
  
  if(integrator==0){
    double result[test.getNrofcross()];
    double nuestimate=0.,tzestimate=0.;
    rombergerN(nuint,nu_min,nu_max,test.getNrofcross(),result,PREC,3,7,&nuestimate, 
				  torz,Q2,nocuts,Ebeam, &test, &tzestimate);
    cout << Q2 << " ";
    for (int l=0;l<test.getNrofcross();l++) cout << result[l] << " ";
    cout << endl;
  }
  else if(integrator==1||integrator==2){
    /*! struct that is used for integrators (clean ones)*/
    struct Ftor_tz {

      /*! integrandum function */
      static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
	Ftor_tz &p = * (Ftor_tz *) param;
	p.f(ret,x[0],x[1],*p.cross,p.Q2,p.Ebeam);
      }
      RhoTCross *cross;/*!< pointer to the grid where the integration is performed */
      double Q2;/*!< [GeV^2] Qsquared */
      double Ebeam; /*!< [GeV] rho momentum */
      /*! integrandum 
      * \param res results
      * \param pm first integration variable
      * \param costheta second integration variable
      * \param phi third integration variable
      * \param cross the RhoTCross instance
      */
      void (*f)(numint::vector_d & res, double nu, double other, RhoTCross & cross, double Q2, double Ebeam);
    };
    int res=90;
    unsigned count=0;
    numint::vector_d ret(test.getNrofcross(),0.);
    numint::array<double,2> lower = {{nu_min,torz?-0.4:0.9}};
    numint::array<double,2> upper = {{nu_max,torz?-0.1:1.}};
    
    Ftor_tz F;
    F.cross = &test;
    F.Q2=Q2;
    F.Ebeam=Ebeam;
    numint::mdfunction<numint::vector_d,2> mdf;
    mdf.func = &Ftor_tz::exec;
    mdf.param = &F;
    F.f=torz?klaas_t:klaas_z;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-12,1.E-03,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,1.E-12,1.E-03,20000,ret,count,0);
    cout << Q2 << " ";
    for (int l=0;l<test.getNrofcross();l++) cout << ret[l] << " ";
    cout << count << " " << res << endl;
  }
  
  
}


void nuint(double nu, double* result, va_list ap){
  
  int torz = va_arg(ap,int);
  double Q2 = va_arg(ap,double);
  int nocuts = va_arg(ap,int);
  double Ebeam = va_arg(ap,double);
  RhoTCross *ptest = va_arg(ap,RhoTCross*);
  double *ptzestimate = va_arg(ap,double*);

  if(torz) rombergerN(tzint,-0.4,-0.1,ptest->getNrofcross(),result,PREC,3,6,
	      ptzestimate,nu,torz,Q2,nocuts,Ebeam,ptest);
  else rombergerN(tzint,0.9,1.,ptest->getNrofcross(),result,PREC,3,6,
	      ptzestimate,nu,torz,Q2,nocuts,Ebeam,ptest);
  return;
}

void tzint(double tz, double* result, va_list ap){
  
  double nu = va_arg(ap,double); 
  int torz = va_arg(ap,int);
  double Q2 = va_arg(ap,double);
  int nocuts = va_arg(ap,int);
  double Ebeam = va_arg(ap,double);
  RhoTCross *ptest = va_arg(ap,RhoTCross*);  
  
  if(torz) ptest->getCrosst(result,Ebeam,Q2,nu,tz);
  else ptest->getCrossz(result,Ebeam,Q2,nu,tz);
  return;
}


