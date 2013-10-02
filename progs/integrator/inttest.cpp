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

  //integration function declaration
  void bol(numint::vector_d &, double r, double ctheta, double phi, double c) ;
  
  //structure needed for integration
  struct Ftor {
    //function that is executed by integration algorithm
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor &p = * (Ftor *) param; //refer to yourself to get the parameters
      p.f(ret,x[0],x[1],x[2],p.c); //use the integration function supplied through "f"
    }
    //parameters
    double c;
    //integration function prototype, first param is ref to vector with results, then integration vars, then parameter
    void (*f)(numint::vector_d &, double r, double theta, double phi, double c);
  };

  Ftor F;
  F.c = 1.; //initialize parameter
  F.f = bol; //initialize integration function

  numint::mdfunction<numint::vector_d,3> mdf; //multidim function (3d) that returns multiple values
  mdf.func = &Ftor::exec; //initialize integration function (static)
  mdf.param = &F; //initialize parameters pointer

  unsigned neval = 0; //# of function evaluations when convergence reached, return value via reference
  int error; //evaluation code
  numint::array<double,3> lower = {{0.,-1.,0.}}; //upper and lower bounds
  numint::array<double,3> upper = {{1.,1.,2.*PI}};
  
  double SIGNIF = 1.E-05; //convergence criterion
  vector<double> ret(2,0.); //return values
  
  error = numint::cube_romb(mdf,lower,upper,SIGNIF*1.E-03,SIGNIF,ret,neval,0);
  cout << ret[0] << " " << ret[1] << " " << neval << " " << error << " " << 4./3.*PI << endl;
  
  int maxEval=1E07; //max # of function evaluations
  error = numint::cube_adaptive(mdf,lower,upper,SIGNIF*1.E-03,SIGNIF,maxEval,ret,neval,0);
  cout << ret[0] << " " << ret[1] << " " << neval << " " << error << endl;
  
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

  