#ifndef RESONANCE_HPP
#define RESONANCE_HPP

#include <complex>
using namespace std;

class Resonance{
public:
  Resonance(double m, complex<double> coeff, double sigma0=40., double sigmaslope=0.,
	    double betain=6.,double epsilonin=-0.5,double betaoffin=8., double lambdain=1.2);
  Resonance(const Resonance &copy);
  Resonance& operator=(const Resonance&);
  ~Resonance();
  double getSigma(double Q2) const;
  double getBeta(double Q2=0.) const{return beta;}
  double getBetaoff() const{return betaoff;}
  double getSigma0() const{return sigma0;}
  double getSigmaslope() const{return sigmaslope;}
  double getEpsilon() const{return epsilon;}
  double getLambda() const{return lambda;}
  double getMass2() const{return mass2;}
  double getMass() const{return mass;}
  complex<double> getCoeff() const{return coeff;}
  
  void setCoeff(complex<double> coeff_in){coeff=coeff_in;}
  
private:
  double mass;
  double mass2;
  double sigma0;
  double sigmaslope;
  double beta;
  double betaoff;
  double epsilon;
  double lambda;
  complex<double> coeff;  
  
};
#endif