#include <iostream>
#include <cstdlib>

using namespace std;

#include <RhoDeuteron.hpp>
#include <Utilfunctions.hpp>
#include <string>

void tzint(double tz, double* result, va_list ap);
void nuint(double nu, double* result, va_list ap);
void klaas_t(numint::vector_d & res, double nu, double t, RhoDeuteron & deuteron, double Q2, double Ebeam){
  res=numint::vector_d(2,0.);
  double results[2];
  deuteron.getCrosst(results,Ebeam,Q2,nu,t);
  for(int i=0;i<2;++i) res[i]=results[i];
}
void klaas_z(numint::vector_d & res, double nu, double z, RhoDeuteron & deuteron, double Q2, double Ebeam){
  res=numint::vector_d(2,0.);
  double results[2];
  deuteron.getCrossz(results,Ebeam,Q2,nu,z);
  for(int i=0;i<2;++i) res[i]=results[i];
}

int main(int argc, char *argv[])
{

  int torz = atoi(argv[1]);
  double Q2 = atof(argv[2]);
  double nu_min = atof(argv[3]);
  double nu_max = atof(argv[4]);
  int nocuts = atoi(argv[5]);
  int integrator = atoi(argv[6]);
  int maxEval = atoi(argv[7]);
  bool fsi = atoi(argv[8]);
//   double prec=atof(argv[9]);
  //int maxEval2 = atoi(argv[10]);
  //string homedir=argv[10];
  double Ebeam = 5.014;
  

  RhoDeuteron test = RhoDeuteron("paris",400,nocuts,2,1.E-03,maxEval,fsi);  
  
  if(integrator==0){
    double result[2];
    double nuestimate=0.,tzestimate=0.;
    rombergerN(nuint,nu_min,nu_max,2,result,PREC,3,7,&nuestimate, 
				  torz,Q2,nocuts,Ebeam, &test, &tzestimate);
    cout << Q2 << " ";
    for (int l=0;l<2;l++) cout << result[l] << " ";
    cout << endl;
  }
  else if(integrator==1||integrator==2){
    /*! struct that is used for integrators (clean ones)*/
    struct Ftor_tz {

      /*! integrandum function */
      static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
	Ftor_tz &p = * (Ftor_tz *) param;
	p.f(ret,x[0],x[1],*p.deuteron,p.Q2,p.Ebeam);
      }
      RhoDeuteron *deuteron;/*!< pointer to the grid where the integration is performed */
      double Q2;/*!< [GeV^2] Qsquared */
      double Ebeam; /*!< [GeV] rho momentum */
      /*! integrandum 
      * \param res results
      * \param pm first integration variable
      * \param costheta second integration variable
      * \param phi third integration variable
      * \param deuteron the RhoDeuteron instance
      */
      void (*f)(numint::vector_d & res, double nu, double other, RhoDeuteron & deuteron, double Q2, double Ebeam);
    };
    int res=90;
    unsigned count=0;
    numint::vector_d ret(2,0.);
    numint::array<double,2> lower = {{nu_min,torz?-0.4:0.9}};
    numint::array<double,2> upper = {{nu_max,torz?-0.1:1.}};
    
    Ftor_tz F;
    F.deuteron = &test;
    F.Q2=Q2;
    F.Ebeam=Ebeam;
    numint::mdfunction<numint::vector_d,2> mdf;
    mdf.func = &Ftor_tz::exec;
    mdf.param = &F;
    F.f=torz?klaas_t:klaas_z;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-12,1.E-03,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,1.E-12,1.E-03,1e05,ret,count,0);
    cout << Q2 << " ";
    for (int l=0;l<2;l++) cout << ret[l] << " ";
    cout << count << " " << res << endl;
  }
  
  
}


void nuint(double nu, double* result, va_list ap){
  
  int torz = va_arg(ap,int);
  double Q2 = va_arg(ap,double);
  int nocuts = va_arg(ap,int);
  double Ebeam = va_arg(ap,double);
  RhoDeuteron *ptest = va_arg(ap,RhoDeuteron*);
  double *ptzestimate = va_arg(ap,double*);

  if(torz) rombergerN(tzint,-0.4,-0.1,2,result,PREC,3,6,
	      ptzestimate,nu,torz,Q2,nocuts,Ebeam,ptest);
  else rombergerN(tzint,0.9,1.,2,result,PREC,3,6,
	      ptzestimate,nu,torz,Q2,nocuts,Ebeam,ptest);
  return;
}

void tzint(double tz, double* result, va_list ap){
  
  double nu = va_arg(ap,double); 
  int torz = va_arg(ap,int);
  double Q2 = va_arg(ap,double);
  int nocuts = va_arg(ap,int);
  double Ebeam = va_arg(ap,double);
  RhoDeuteron *ptest = va_arg(ap,RhoDeuteron*);  
  
  if(torz) ptest->getCrosst(result,Ebeam,Q2,nu,tz);
  else ptest->getCrossz(result,Ebeam,Q2,nu,tz);
  return;
}


