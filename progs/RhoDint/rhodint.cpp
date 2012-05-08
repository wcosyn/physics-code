#include <iostream>
#include <cstdlib>

using namespace std;

#include <RhoDeuteron.hpp>
#include <Utilfunctions.hpp>

void tzint(double tz, double* result, va_list ap);
void nuint(double nu, double* result, va_list ap);

int main(int argc, char *argv[])
{

  int torz = atoi(argv[1]);
  double Q2 = atof(argv[2]);
  double nu_min = atof(argv[3]);
  double nu_max = atof(argv[4]);
  int nocuts = atoi(argv[5]);
  double Ebeam = 5.014;
  

  RhoDeuteron test = RhoDeuteron("paris",400.,nocuts);
  double result[2];
  double nuestimate=0.,tzestimate=0.;
  rombergerN(nuint,nu_min,nu_max,2,result,PREC,3,7,&nuestimate, 
				torz,Q2,nocuts,Ebeam, &test, &tzestimate);
  cout << Q2 << " " << result[0] << endl;
  
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

