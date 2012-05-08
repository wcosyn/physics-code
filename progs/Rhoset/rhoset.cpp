#include <iostream>
#include <cstdlib>

using namespace std;

#include <RhoTCross.hpp>
#include <Utilfunctions.hpp>


void tzint(double tz, double* result, va_list ap);
void nuint(double nu, double* result, va_list ap);

int main(int argc, char *argv[])
{

  int torz = atoi(argv[1]);
  int nucleus = atoi(argv[2]);
  double Q2 = atof(argv[3]);
  double nu_min = atof(argv[4]);
  double nu_max = atof(argv[5]);
  int nocuts = atoi(argv[6]);
  double usersigma = atof(argv[7]);
  double Ebeam = 5.014;
  

  RhoTCross test = RhoTCross(nucleus,400,argv[8],nocuts,1,usersigma);
  double result[NROFRES];
  double nuestimate=0.,tzestimate=0.;
  rombergerN(nuint,nu_min,nu_max,NROFRES,result,PREC,3,7,&nuestimate, 
				torz,Q2,nocuts,Ebeam, &test, &tzestimate);
  cout << Q2 << " ";
  for (int l=0;l<NROFRES;l++) cout << result[l] << " ";
  cout << endl;
  
}


void nuint(double nu, double* result, va_list ap){
  
  int torz = va_arg(ap,int);
  double Q2 = va_arg(ap,double);
  int nocuts = va_arg(ap,int);
  double Ebeam = va_arg(ap,double);
  RhoTCross *ptest = va_arg(ap,RhoTCross*);
  double *ptzestimate = va_arg(ap,double*);

  if(torz) rombergerN(tzint,-0.4,-0.1,NROFRES,result,PREC,3,6,
	      ptzestimate,nu,torz,Q2,nocuts,Ebeam,ptest);
  else rombergerN(tzint,0.9,1.,NROFRES,result,PREC,3,6,
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


