#include "Resonance.hpp"
#include <constants.hpp>

using namespace std;

Resonance::Resonance(double m, complex<double> coefficient, double s0, double slope,
		     double betain,double epsilonin,double betaoffin, double lambdain): 
mass(m),
mass2(m*m),
sigma0(s0),
sigmaslope(slope),
beta(betain*1.E-06), //units conversion from GeV^-2 to MeV^-2
epsilon(epsilonin),
betaoff(betaoffin*1.E-06),
lambda(lambdain*1.E+06),
coeff(coefficient){
  
}

Resonance::Resonance(const Resonance &copy): 
mass(copy.getMass()),
mass2(copy.getMass2()),
sigma0(copy.getSigma0()),
sigmaslope(copy.getSigmaslope()),
beta(copy.getBeta()),
epsilon(copy.getEpsilon()),
betaoff(copy.getBetaoff()),
lambda(copy.getLambda()),
coeff(copy.getCoeff()){
  
}

Resonance& Resonance::operator=(const Resonance& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
     (*this).mass=rhs.getMass();
     (*this).mass2=rhs.getMass2();
     (*this).sigma0=rhs.getSigma0();
     (*this).sigmaslope=rhs.getSigmaslope();
     (*this).beta=rhs.getBeta();
     (*this).epsilon=rhs.getEpsilon();
     (*this).betaoff=rhs.getBetaoff();
     (*this).lambda=rhs.getLambda();
     (*this).coeff=rhs.getCoeff();

     
     
  }

  return *this;
}



Resonance::~Resonance(){}

double Resonance::getSigma(double Q2) const{
    return (getSigma0()+getSigmaslope()*((getMass()>2.4E03?2.4E03:getMass())-MASSP)/(1.E-03*Q2))/10*INVHBARC*INVHBARC;

}



