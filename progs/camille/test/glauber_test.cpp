/** \file
 * This is a simple test file to try out different constructors and member functions
 * of some classes.
 * \author Camille Colle
 * \date 07/01/2013
 */


#include <FourVector.h>
#include "MeanFieldNucleus.hpp"
#include "MeanFieldNucleusThick.hpp"
#include <iostream>
#include <FastParticle.hpp>
#include <GlauberGridThick.hpp>

using std::string;
using std::cout;
using std::endl;
using std::complex;

int main(int argc, char* argv[])
{

  string homedir = argv[1];
  complex<double> wave[4];
  MeanFieldNucleusThick CarbonThick(1,homedir); // make a MF Carbon nucleus (thick = include nuclear densities)



  int type=0; // proton
  int inc=0; // is particle a beam particle
  double momentum=20e3; //Momentum of particle [MeV] along z-axis
  double ptheta = 0.; // momentum theta angle
  double pphi = 0.; // momentum phi angle
  double hard_scale = 8. ; // hard_scale ?? -> seems only to matter in CT-calculations
  double Gamma = 0. ; // decay width of particle, stable proton -> Gamma = 0
  FastParticle proton(type,inc,momentum,ptheta,pphi,hard_scale,Gamma,homedir);
  proton.printParticle();
  
  
  int r_grid = 80 ; // rgrid
  int cth_grid = 18; // cos theta grid
  int phi_grid = 5; // phi grid
  double prec = 1.;
  int integrator = 1; // 0 Win fubini, 1 Klaas romberg?, 2 MIT adaptive
  GlauberGridThick gridthick(r_grid,cth_grid,phi_grid,&CarbonThick,prec,integrator,homedir);
  gridthick.printParticles();
  gridthick.addParticle(proton);
  gridthick.printParticles();
  gridthick.updateGrids();
  gridthick.printFsi_src_ct_grid();
}
