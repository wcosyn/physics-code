/** \file
 * This is a simple test file to try out different constructors and member functions
 * of some classes.
 * \author Camille Colle
 * \date 07/01/2013
 */


#include <FourVector.h>
#include <TSpinor.h>
#include "MeanFieldNucleus.hpp"
#include "TMFSpinor.hpp"
#include "MeanFieldNucleusThick.hpp"
#include <iostream>

using std::string;
using std::cout;
using std::endl;
using std::complex;

int main(int argc, char* argv[])
{
  FourVector<double> f(1.,0.,0.,0.); // simple four vector
  double mass=1.; //mass of particle, if on shell this could be derived from the four vector, f[0]^2-f[1]^2-f[2]^2-f[3]^2
  TSpinor s(f,mass,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity); //plane wave dirac spinor
  
  // from here on this resembles very much ../glauber/testglauber.cpp
  string homedir = argv[1];
  complex<double> wave[4];
  MeanFieldNucleusThick CarbonThick(1,homedir); // make a MF Carbon nucleus (thick = include nuclear densities)
  int shellindex = 0;
  int m= -1;
  double r= 2.;
  double costheta= 0.5;
  double phi= 0.;
  CarbonThick.getWaveFunction(wave,shellindex,m,r,costheta,phi); // get wavefunction of nucleon with spec. quantum numbers
  cout << wave[0] << " " << wave[1] << " " << wave[2] << " " << wave[3] << endl;
  cout << TMFSpinor(CarbonThick,shellindex,m,r,costheta,phi) << endl; // mean field spinor of nucleus
  
  return 0;
}
