#include <iostream>
#include <cstdlib>

using namespace std;

#include <numint/array.hpp>
#include <numint/numint.hpp>
#include <numint/typedef.hpp>
#include <cmath>
#include <vector>
#include <adaptive/cubature.h>

int main(int argc, char *argv[])
{

  void bol(numint::vector_d &, double r, double ctheta, double phi, double c) ;
  void bol(unsigned n, const double *x, void *param, unsigned fdim, double *ret) ;
  
  struct Ftor {

    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor &p = * (Ftor *) param;
      p.f(ret,x[0],x[1],x[2],p.c);
    }
    double c;
    void (*f)(numint::vector_d &, double x, double y, double z, double c);
  };

  Ftor F;
  F.c = 1.;
  F.f = bol;

  numint::mdfunction<numint::vector_d,3> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

  unsigned neval = 0;
  int error;
  numint::array<double,3> lower = {{0.,-1.,0.}};
  numint::array<double,3> upper = {{1.,1.,2.*PI}};
  
  double SIGNIF = 1.E-05;
  vector<double> ret(2,0.);
  
//   error = numint::cube_romb(mdf,lower,upper,SIGNIF*1.E-03,SIGNIF,ret,neval,0);
//   cout << ret[0] << " " << ret[1] << " " << neval << " " << error << " " << 4./3.*PI << endl;
  
  double xmin[3]={0.,-1.,0.},xmax[3]={1.,1.,2.*PI}, c=1.;
  double res[2];
  double err[2];
  error = adapt_integrate(2,bol,&c,3,xmin,xmax,0,SIGNIF*1.E-03,SIGNIF,res,err);
  cout << res[0] << " " << res[1] << " " << err[0] << " " << err[1] << " " << error << endl;
  
}

void bol( numint::vector_d &ret, double r, double ctheta, double phi, double c)
{
  ret=numint::vector_d(2,0.);
  ret[0]=(c*r*r*cos(phi)*cos(phi)+pow(phi,0.5));
//   if(phi>PI) ret[1]=(c*r*cos(phi)/(4.+ctheta));
//   else ret[1]=0.;
//   if(phi>3.*PI/2.) ret[1]=(c*r*r);
//   else ret[1]=0.;
  
}

void bol(unsigned n, const double *x, void *param, unsigned fdim, double *ret){
  double c=*((double*)param);
  ret[0]=c*x[0]*x[0]*cos(x[2])*cos(x[2])+pow(x[2],0.5);
//   if(x[2]>PI/3.) ret[1]=(c*x[0]*cos(x[2])/(4.+x[1]));
//   else ret[1]=0.;
//   if(x[2]>sqrt(2.)/2.*PI) ret[1]=c*x[0]*x[0];
//   else ret[1]=0.;
}  
  