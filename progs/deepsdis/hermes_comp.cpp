#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <NuclStructure.hpp>
#include <NucleonStructure.hpp>

#include <TDeuteron.h>
#include <numint/numint.hpp>


void k_int(numint::vector_d & res, double knorm, double kcosth, TDeuteron::Wavefunction *wfref,
          double Q2, double x, double nu, double qvec, double gamma, NucleonStructure &strp);


int main(int argc, char *argv[])
{
  double Q2 = atof(argv[1])*1.E06; // parse from argv or something
  double x=atof(argv[2]);
  double Ein=27.6E03;
  double Azz=atof(argv[3])*0.01;
  double b1=atof(argv[4])*0.01;
  
//   string wf = argv[2];
  string nuclstruc = argv[5];
//   
//   TDeuteron::Wavefunction *wfref;
//   wfref = TDeuteron::Wavefunction::CreateWavefunction(wf);

  NucleonStructure strp(nuclstruc);
  
  
  
  struct Ftor_b1 {

    /*! integrandum function */
    static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
      Ftor_b1 &p = * (Ftor_b1 *) param;
      p.f(ret,x[0],x[1], p.wfref,p.Q2,p.x,p.nu,p.qvec,p.gamma,p.strp);
    }
    TDeuteron::Wavefunction *wfref;
    double Q2;
    double x;
    double nu;
    double qvec;
    double gamma;
    NucleonStructure strp;
    
    
    
    void (*f)(numint::vector_d & res, double knorm, double kcosth, TDeuteron::Wavefunction *wfref,
              double Q2, double x, double nu, double qvec, double gamma, NucleonStructure &strp);
  };
 
  double nu=Q2/(2.*MASSn*x);
  double Eout=Ein-nu;
  double qvec=sqrt(Q2+nu*nu);
  double y=nu/Ein;
  double gamma=sqrt(Q2)/nu;
  double eps=(1-y-gamma*gamma*y*y/4)/(1-y+y*y/2+gamma*gamma*y*y/4);
  double thetae=asin(sqrt(Q2/4/Ein/Eout))*2;
  double thetaq=acos((Ein*Ein+qvec*qvec-Eout*Eout)/2/Ein/qvec);
  double kap=1+gamma*gamma;
  double R=NucleonStructure::getr1998(x,Q2);
  
  double term1=(1.+3.*cos(2*thetaq))*(2.*(1-eps)*kap-gamma*gamma/3./kap*(1+R)+eps*(kap*kap+kap+1)*2*(1+R)/3/kap);
  double term2=3.*sin(2.*thetaq)*sqrt(2.*eps*(1+eps))*gamma*(kap+2)*(1+R)/3./kap;
  double term3=3.*(1-2.*cos(2.*thetaq))*eps*gamma*gamma*(1+R)/3/kap;

  double factor=1./8./(1+eps*R)*(term1+term2+term3);

  cout << thetae/PI*180<< " " << thetaq/PI*180<< " " << nu<< " " << qvec<< " " << y<< " " << eps<< " " << gamma*gamma<< " " << pow(2.*MASSn*x/sqrt(Q2),2.)<< " " << R << " "<< factor << endl;
  
  
  
}