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
  if(std::isnan(massX)){
    cerr << "invalid kinematics, X has no physical mass!" << X_mu*X_mu << endl;
    assert(1==0);
  }
  Beam_mu=FourVector<double>(vecBeam.Mag(),vecBeam.X(),vecBeam.Y(),vecBeam.Z());
  isCollinear=0;
  if(vecA.Perp2()+vecq.Perp2()==0.) isCollinear=1.;  //check if this is collinear kinematics
  
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
  if(std::isnan(massX)){
    cerr << "invalid kinematics, X has no physical mass! " << X_mu*X_mu << endl;
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
  k_z=(alpha_s-1)*Ek; //sign could also be taken +!!!! (check this later...)
  zs=(A_mu*ps_mu)/(A_mu*q_mu);
  kvec=TVector3(k_perp.X(),k_perp.Y(),k_z);
}

LightConeKin2to2::LightConeKin2to2(LightConeKin2to2 &rhs){
  massA=rhs.massA;
  Q2=rhs.Q2;
  massN=rhs.massN;
  A_mu=rhs.A_mu;
  q_mu=rhs.q_mu;
  ps_mu=rhs.ps_mu;
  pi_mu=rhs.pi_mu;
  X_mu=rhs.X_mu;
  massX=rhs.massX;
  Beam_mu=rhs.Beam_mu;
  isCollinear=rhs.isCollinear;
  yA=rhs.yA;
  xA=rhs.xA;
  yN=rhs.yN;
  xN=rhs.xN;
  pA_plus=rhs.pA_plus;
  pA_perp=rhs.pA_perp;
  q_perp=rhs.q_perp;
  ps_perp=rhs.ps_perp;
  alpha_s=rhs.alpha_s;
  alpha_i=rhs.alpha_i;
  epsilon=rhs.epsilon;
  k_perp=rhs.k_perp;
  k=rhs.k;
  Ek=rhs.Ek;
  k_z=rhs.k_z;
  zs=rhs.zs;
  kvec=rhs.kvec;
}

LightConeKin2to2& LightConeKin2to2::operator=(LightConeKin2to2 &rhs){
 if(this!=&rhs) { // avoid self-assignment
    massA=rhs.massA;
    Q2=rhs.Q2;
    massN=rhs.massN;
    A_mu=rhs.A_mu;
    q_mu=rhs.q_mu;
    ps_mu=rhs.ps_mu;
    pi_mu=rhs.pi_mu;
    X_mu=rhs.X_mu;
    massX=rhs.massX;
    Beam_mu=rhs.Beam_mu;
    isCollinear=rhs.isCollinear;
    yA=rhs.yA;
    xA=rhs.xA;
    yN=rhs.yN;
    xN=rhs.xN;
    pA_plus=rhs.pA_plus;
    pA_perp=rhs.pA_perp;
    q_perp=rhs.q_perp;
    ps_perp=rhs.ps_perp;
    alpha_s=rhs.alpha_s;
    alpha_i=rhs.alpha_i;
    epsilon=rhs.epsilon;
    k_perp=rhs.k_perp;
    k=rhs.k;
    Ek=rhs.Ek;
    k_z=rhs.k_z;
    zs=rhs.zs;
    kvec=rhs.kvec;
  }
  return *this;
  
}
