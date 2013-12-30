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
  double m=1.4455E03;
  double p=3.45346E03;
  double E=sqrt(m*m+p*p);
  
  FourVector<complex<double> > Z(E,0.,0.,p);
  
  FourVector<complex<double> > polVectorPlus(0.,
                                                    -1./sqrt(2.),
                                                    complex<double>(0.,-1./sqrt(2.)),
                                                 0.);
  FourVector<complex<double> > polVectorMin(0.,
                                                   1./sqrt(2.),
                                                   complex<double>(0.,-1./sqrt(2.)),
                                                   0.);
  FourVector<complex<double> > polVector0(qvec/sqrt(Q2),0.,0.,nu/sqrt(Q2));
  FourVector<complex<double> > polVectorZ(p/m,0.,0.,E/m);
  
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      cout << i << " " << j << " " << -1.*polVectorMin[j]*conj(polVectorMin[i])-conj(polVectorZ[i])*polVectorZ[j]-1.*polVectorPlus[j]*conj(polVectorPlus[i]) << " " << 
	  ((i==j?(i==0?1.:-1.):(0.))-Z[i]*Z[j]/(m*m)) << endl;
    }
  }
  cout << endl << endl;
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      cout << i << " " << j << " " << -1.*polVectorMin[j]*conj(polVectorMin[i])+conj(polVector0[i])*polVector0[j]-1.*polVectorPlus[j]*conj(polVectorPlus[i]) << " " << 
	  ((i==j?(i==0?1.:-1.):(0.))+q[i]*q[j]/Q2) << endl;
    }
  }
  cout << endl << endl;
  
  cout << Trace(((polVector0*gamma_mu)*(ein*gamma_mu)*(polVectorPlus*gamma_mu)*(eout*gamma_mu)).value())/denom/2. << " " << sqrt(Q2)/qvec*sqrt(Q2/qvec/qvec+tanth*tanth)/sqrt(2.)<< endl;

  cout << Trace(((polVector0*gamma_mu)*gamma_5*(ein*gamma_mu)*(polVectorPlus*gamma_mu)*(eout*gamma_mu)).value())/denom/2. << " " 
  << sqrt(Q2/2)/qvec*tanth << endl;
  cout << -Trace(((polVectorMin*gamma_mu)*gamma_5*(ein*gamma_mu)*(polVectorPlus*gamma_mu)*(eout*gamma_mu)).value())/denom/2. << " " 
  << tanth*sqrt(Q2/qvec/qvec+tanth*tanth) << endl;
      
  cout << Trace((gamma_mu[0]*gamma_mu[1]*gamma_mu[2]*gamma_mu[3]*gamma_5).value()) << endl;
      
}