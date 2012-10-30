/*! \mainpage Glauber ISI/FSI RMSGA C++ Project
 * \author Wim Cosyn
 * \date 16/08/2011
 * \brief This code implements classes for the Glauber RMSGA formalism, including CT and SRC effects.
 * 
 * \details 
 * - It contains classes for a mean field nucleus, a mean field nucleus with densities (used in thickness calculations). <BR>
 * - A class for correlated FSI calculations containing all the needed functions and a grid for the gamma functions. <BR>
 * - Four abstract classes (one for a general FSI grid, one that adds a CT grid, one for a general thickness FSI grid and one that adds a thickness CT grid). <BR>
 * - Two glauber classes : one without thickness, one with (also adds SRC to the ISI/FSI). <BR>
 * - A special class for the glauber grid of one particle, exploiting the symmetry along the z-axis (no phi dependence).  <BR>
*/
#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <FsiCorrelator.hpp>
#include <FastParticle.hpp>
#include <constants.hpp>
#include <GlauberGrid.hpp>
#include <OneGlauberGrid.hpp>
#include <GlauberDecayGridThick.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <Cross.hpp>
#include <DistMomDistrGrid.hpp>
#include <TRotation.h>

//run ./onenucl [Q2 [MeV^2]] [omega] [missing momentum]
int main(int argc, char *argv[])
{
  double prec=atof(argv[1]);
  int integrator=atoi(argv[2]);
  int maxEval=atoi(argv[3]);

  string homedir="/home/wim/Code/share";
  MeanFieldNucleusThick CarbonThick(1,homedir);
  double thetarho=0.;
  double phirho=0.;
  FastParticle rho(4, 0, 1821.82,0.,0.,1.01,145.,homedir);
  rho.setScatter(25.,6.,-0.2);
  GlauberDecayGridThick grid(60,20,5,&CarbonThick,1.E-05,2,homedir);
  DistMomDistrGrid rhogrid(1, 400., 30,20,5,&grid,prec,integrator,maxEval,abs(thetarho),homedir);
  TVector3 rhovec(sin(thetarho)*cos(phirho),sin(thetarho)*sin(phirho),cos(thetarho));
  rhovec*=1500.;
  TRotation rot;
  TVector3 axis(-sin(thetarho)*sin(phirho),sin(thetarho)*cos(phirho),0.);
// //   axis.Print();
  if(thetarho>0.) rot.Rotate(thetarho,-axis);
  else rot.Rotate(-thetarho,axis);
// //   rhovec.Transform(rot).Print();
  
  grid.clearParticles();
  grid.addParticle(rho);
  grid.updateGrids();
//   grid.printFsi_grid();
//   grid.printFsi_decay_grid();
  grid.clearKnockout();
    rhogrid.updateGrids(&grid,1,rot);
    rhogrid.printRho_grid(0);
//   rhogrid.printRho_grid(1);
//   rhogrid.printRho_grid(2);
//   rhogrid.printRho_grid(4);
  
  
}


