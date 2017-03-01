#include <iostream>
#include <cstdlib>

using namespace std;


#include <iostream>
#include <cmath>
#include "constants.hpp"
#include "FourVector.h"
#include "GammaStructure.h"
#include "NuclStructure.hpp"
#include "WeakQECross.hpp"
#include "TLeptonKinematics.h"
#include <MeanFieldNucleusThick.hpp>


int main(int argc, char *argv[]){


//   for(int i=1;i<100;i++){
//     double x=0.01*i;
//     NuclStructure F2CTEQ=NuclStructure(1,200.E06,x,0,"CTEQ");
//     NuclStructure F2SLAC=NuclStructure(1,200.E06,x,0,"SLAC");
// //     NuclStructure F2CB=NuclStructure(1,200.E06,x,0,"CB");
//     NuclStructure F2alekhin=NuclStructure(1,200.E06,x,0,"Alekhin");
//     cout << x << " " << F2CTEQ.getF2() << " " << F2SLAC.getF2() << 	" " << F2alekhin.getF2() << endl;
//     
//   }
//   exit(1);
  
  const FourVector<GammaStructure> gamma_mu=FourVector<GammaStructure>(GammaStructure(0.,0.,1.),
											GammaStructure(0.,0.,0.,1.),
				      GammaStructure(0.,0.,0.,0.,1.),GammaStructure(0.,0.,0.,0.,0.,1.));


  const GammaStructure gamma_5(0.,1.);
  
  
  double Ein=4454.;
  double Q2=0.9E06;
  double x=1.;
  double mp=MASSn;//938.272;
  double nu=Q2/(2.*mp*x);
  double qvec=sqrt(Q2+nu*nu);
  double Eout=Ein-nu;
  double mmu=105.;
  double pout=sqrt(Eout*Eout-mmu*mmu);
  double zeta=pout/Eout;
//   cout << zeta << " " << nu << " " << Eout << " " << pout << endl;
  double thetae=acos((1.-(Q2+mmu*mmu)/(2.*Ein*Eout))/zeta);
//   cout << Q2 << " " << 2.*Ein*Eout*(1-zeta*cos(thetae))-mmu*mmu <<  endl;;
  
  double tanth=tan(thetae/2.);
//   double kx=sqrt(Q2/qvec/qvec*Ein*Eout*pow(cos(thetae/2.),2.));
//   double kz=Q2/(2.*qvec)+nu*Ein/qvec;
//   double kprimez=-Q2/(2.*qvec)+nu*Eout/qvec;
  double denom=4.*Ein*Eout*pow(cos(thetae/2.),2.);
//   
  
  double kx=zeta*Ein*Eout*sin(thetae)/qvec;
  double kz=(Q2+mmu*mmu)/2/qvec+nu*Ein/qvec;
  double kprimez=(-Q2+mmu*mmu)/2/qvec+nu*Eout/qvec;
  
  double tanPRC=Q2/(pow(Ein+Eout,2.)-qvec*qvec);
  double delta=mmu/sqrt(Q2);
  double rho=Q2/qvec/qvec;
  double rho2=qvec/(Ein+Eout);
  
  double ratio=(pow(Ein+Eout,2.)-qvec*qvec)/2./Ein/Eout;
  
  double vlPRC=1-delta*delta*tanPRC;
  double vllPRC=nu*nu/qvec/qvec+(1+2.*nu/qvec/rho2+rho*delta*delta)*delta*delta*tanPRC;
  double vtPRC =tanPRC+rho/2.-delta*delta/rho2*(nu/qvec+rho*rho2*delta*delta/2.)*tanPRC;
  double vtprimePRC=1./rho2*(1.-nu*rho2/qvec*delta*delta)*tanPRC;
  
//   cout << vtprimePRC*ratio << " " << (Ein+Eout)/qvec*(1.-zeta*cos(thetae))-mmu*mmu/qvec/Eout << endl;
  double M_A=1.03E03;
  string homedir="/home/wim/Code/trunk/share/";
  MeanFieldNucleusThick Nucleus(1,homedir);
  TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(TLeptonKinematics::muon,Ein);
//   WeakQECross obs(lepton,&Nucleus,1.E-03,2,homedir,1, 1.03E03, 0., 0.);  
  double crossHanu=WeakQECross::getElWeakQECross(Q2,Ein,2,1,1,M_A);
  cout << crossHanu << endl;

  exit(1);
  
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
      
  cout << Trace((gamma_5*gamma_mu[0]*gamma_mu[1]*gamma_mu[2]*gamma_mu[3]).value()) << endl;
  cout << (I_UNIT*gamma_mu[0]*gamma_mu[1]*gamma_mu[2]*gamma_mu[3]).value() << endl;
  cout << gamma_5.value() << endl;
}
