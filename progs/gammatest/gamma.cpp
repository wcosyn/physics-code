#include <iostream>
#include <cstdlib>

using namespace std;


#include <iostream>
#include <cmath>
#include "constants.hpp"
#include "FourVector.h"
#include "GammaStructure.h"

int main(int argc, char *argv[]){


  const FourVector<GammaStructure> gamma_mu=FourVector<GammaStructure>(GammaStructure(0.,0.,1.),
											GammaStructure(0.,0.,0.,1.),
				      GammaStructure(0.,0.,0.,0.,1.),GammaStructure(0.,0.,0.,0.,0.,1.));


  const GammaStructure gamma_5(0.,1.);
  
  
  double Ein=4454.;
  double Q2=1.9E06;
  double x=1.2;
  double mp=938.272;
  double nu=Q2/(2.*mp*x);
  double qvec=sqrt(Q2+nu*nu);
  double Eout=Ein-nu;
  double thetae=2.*asin(sqrt(Q2/(4.*Ein*Eout)));
  double tanth=tan(thetae/2.);
  double kx=sqrt(Q2/qvec/qvec*Ein*Eout*pow(cos(thetae/2.),2.));
  double kz=Q2/(2.*qvec)+nu*Ein/qvec;
  double kprimez=-Q2/(2.*qvec)+nu*Eout/qvec;
  double denom=4.*Ein*Eout*pow(cos(thetae/2.),2.);
  
  FourVector<complex<double> > q(nu,0.,0.,qvec);
  FourVector<complex<double> > ein(Ein,kx,0.,kz);
  FourVector<complex<double> > eout(Eout,kx,0.,kprimez);
  
  FourVector<complex<double> > polVectorPlus(0.,
                                                    -1./sqrt(2.),
                                                    complex<double>(0.,-1./sqrt(2.)),
                                                 0.);
  FourVector<complex<double> > polVectorMin(0.,
                                                   1./sqrt(2.),
                                                   complex<double>(0.,-1./sqrt(2.)),
                                                   0.);
  FourVector<complex<double> > polVector0(qvec/sqrt(Q2),0.,0.,nu/sqrt(Q2));

  
  
  cout << Trace(((polVector0*gamma_mu)*(ein*gamma_mu)*(polVectorPlus*gamma_mu)*(eout*gamma_mu)).value())/denom/2. << " " << sqrt(Q2)/qvec*sqrt(Q2/qvec/qvec+tanth*tanth)/sqrt(2.)<< endl;

  cout << Trace(((polVector0*gamma_mu)*(ein*gamma_mu)*(polVectorPlus*gamma_mu)*(eout*gamma_mu)*gamma_5).value())/denom/2. << endl;
      
      
      
}