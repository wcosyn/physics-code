//small program that checked kinematics of the hermes b1 experiment and what influence the angle between beam (=polarization axis) and virtual photon had

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


int main(int argc, char *argv[])
{
  double Q2 = atof(argv[1])*1.E06; // parse from argv or something
  double x=atof(argv[2]);
  double Ein=atof(argv[5]); //27.6E03 for hermes 11E03 for JLab12
  double Azz=atof(argv[3])*0.01;
  double b1=atof(argv[4])*0.01;
  
//   string wf = argv[2];
  
     
    
 
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
  
  double R_b=R;
  
  double term1=(1.+3.*cos(2*thetaq))*(2.*(1-eps)*kap-gamma*gamma/3./kap+eps*(kap*kap+kap+1)*2/3/kap);
  double term2=sin(2.*thetaq)*sqrt(2.*eps*(1+eps))*gamma*(kap+2)/kap;
  double term3=(1-2.*cos(2.*thetaq))*eps*gamma*gamma/kap;

  double factor=1./8./(1+eps*R)*(term1+term2+term3);

  cout << thetae/PI*180<< " " << thetaq/PI*180<< " " << eps<< " " << gamma*gamma<< " " << " " << R << " "<< factor << " " << term1/8. << " " << term2/8. << " " << term3/8. << " " << 1+eps*R << endl;
  
  
  
}