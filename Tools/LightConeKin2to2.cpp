#include "LightConeKin2to2.hpp"

#include <cassert>
#include <constants.hpp>

using namespace std;


LightConeKin2to2::LightConeKin2to2(double massA_, double Q2_, double massC_, 
				   TVector3 &vecA, TVector3 &vecq, TVector3 &vecC, TVector3 &vecBeam):
				   massA(massA_),Q2(Q2_),massN(massC_){
  A_mu=FourVector<double>(sqrt(massA*massA+vecA.Mag2()),vecA.X(),vecA.Y(),vecA.Z());
  q_mu=FourVector<double>(sqrt(vecq.Mag2()-Q2),vecq.X(),vecq.Y(),vecq.Z());
  ps_mu=FourVector<double>(sqrt(massN*massN+vecC.Mag2()),vecC.X(),vecC.Y(),vecC.Z());
  pi_mu=A_mu-ps_mu;
  X_mu=(A_mu+q_mu-ps_mu);
  massX=sqrt(X_mu*X_mu);
  if(isnan(massX)){
    cerr << "invalid kinematics, X has no physical mass!" << endl;
    assert(1==0);
  }
  Beam_mu=FourVector<double>(vecBeam.Mag(),vecBeam.X(),vecBeam.Y(),vecBeam.Z());
  isCollinear=0;
  
  yA=(A_mu*q_mu)/(A_mu*Beam_mu);
  xA=Q2/(A_mu*q_mu);
  yN=(pi_mu*q_mu)/(pi_mu*Beam_mu);
  xN=Q2/(2.*(pi_mu*q_mu));
  pA_plus=(A_mu[0]+A_mu[3])/sqrt(2.);
  pA_perp=TVector3(vecA.X(),vecA.Y(),0.);
  ps_perp=TVector3(vecC.X(),vecC.Y(),0.);
  q_perp=TVector3(vecq.X(),vecq.Y(),0.);
  alpha_s=(ps_mu[0]+ps_mu[3])*sqrt(2.)/pA_plus;
  alpha_i=2.-alpha_s;
  epsilon=(1.-yA-pow(xA*yA*massA,2.)/(4.*Q2))/(1.-yA+yA*yA/2.+pow(xA*yA*massA,2.)/(4.*Q2));
  k_perp=ps_perp-alpha_s/2.*pA_perp;
  k=sqrt((massN*massN+k_perp.Mag2())/(alpha_s*alpha_i)-massN*massN);
  Ek=sqrt(k*k+massN*massN);
  k_z=(alpha_s-1)*Ek; //sign could also be taken +!!!! (check this later...)
  zs=(A_mu*ps_mu)/(A_mu*q_mu);
  kvec=TVector3(k_perp.X(),k_perp.Y(),k_z);
}

LightConeKin2to2::LightConeKin2to2(double massA_, double Q2_, double massC_, 
				   double pA, double vecq, TVector3 &vecC, TVector3 &vecBeam):
				   massA(massA_),Q2(Q2_),massN(massC_){
  A_mu=FourVector<double>(sqrt(massA*massA+pA*pA),0.,0.,pA);
  q_mu=FourVector<double>(sqrt(vecq*vecq-Q2),0.,0.,vecq);
  ps_mu=FourVector<double>(sqrt(massN*massN+vecC.Mag2()),vecC.X(),vecC.Y(),vecC.Z());
  pi_mu=A_mu-ps_mu;
  X_mu=(A_mu+q_mu-ps_mu);
  massX=sqrt(X_mu*X_mu);
  if(isnan(massX)){
    cerr << "invalid kinematics, X has no physical mass!" << endl;
    assert(1==0);
  }
  Beam_mu=FourVector<double>(vecBeam.Mag(),vecBeam.X(),vecBeam.Y(),vecBeam.Z());
  isCollinear=1;
  
  yA=(A_mu*q_mu)/(A_mu*Beam_mu);
  xA=Q2/(A_mu*q_mu);
  yN=(pi_mu*q_mu)/(pi_mu*Beam_mu);
  xN=Q2/(2.*(pi_mu*q_mu));
  pA_plus=(A_mu[0]+A_mu[3])/sqrt(2.);
  pA_perp=TVector3(0.,0.,0.);
  q_perp=TVector3(0.,0.,0.);
  ps_perp=TVector3(vecC.X(),vecC.Y(),0.);
  alpha_s=(ps_mu[0]+ps_mu[3])*sqrt(2.)/pA_plus;
  alpha_i=2.-alpha_s;
  epsilon=(1.-yA-pow(xA*yA*massA,2.)/(4.*Q2))/(1.-yA+yA*yA/2.+pow(xA*yA*massA,2.)/(4.*Q2));
  k_perp=ps_perp-alpha_s/2.*pA_perp;
  k=sqrt((massN*massN+k_perp.Mag2())/(alpha_s*alpha_i)-massN*massN);
  Ek=sqrt(k*k+massN*massN);
  k_z=(1-alpha_s)*Ek; //sign could also be taken +!!!! (check this later...)
  zs=(A_mu*ps_mu)/(A_mu*q_mu);
  kvec=TVector3(k_perp.X(),k_perp.Y(),k_z);
}

