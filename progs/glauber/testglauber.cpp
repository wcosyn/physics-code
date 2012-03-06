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

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <FsiCorrelator.hpp>
#include <FastParticle.hpp>
#include <constants.hpp>
#include <GlauberGrid.hpp>
#include <OneGlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include <TMFSpinor.hpp>

int main(int argc, char *argv[])
{
  
  string homedir="/home/wim/Code/share";

  MeanFieldNucleusThick CarbonThick(1,homedir);
  complex<double> wave[4];
  CarbonThick.getWaveFunction(wave,0, -1, 2.,0.5,1.);
  cout << wave[0] << " " << wave[1] << " " << wave[2] << " " << wave[3] << endl;
  cout << TMFSpinor(CarbonThick,0,-1,2.,0.5,1.) << endl;
  CarbonThick.getWaveFunction(wave,1, 1, 1.,0.5,1.);
  cout << wave[0] << " " << wave[1] << " " << wave[2] << " " << wave[3] << endl;
  cout << TMFSpinor(CarbonThick,1,1,1.,0.5,1.) << endl;
 // exit(1);
  //FsiCorrelator CarbonCorr(&CarbonThick,90,90,homedir);
  //FastParticle proton(0, 1, 1234,0.,0.,3.,homedir);
  FastParticle proton2(0, 0, 4878,0.,0.,8.,homedir);
  proton2.printParticle();
  //FastParticle pion(2, 0, 1234,-1.,0.,3.,homedir);
//   OneGlauberGrid grid(60,18,&CarbonThick,homedir);
//   grid.addParticle(&proton2);
//   //grid.addParticle(&pion);
// //  grid.printParticles();
//   grid.fillGrids();
//   grid.addKnockout(0,0);
//    grid.printFsi_grid();
   GlauberGridThick gridthick(60,18,5,&CarbonThick,homedir);
  gridthick.addParticle(&proton2);
//   gridthick.printParticles();
  gridthick.fillGrids();
   GlauberGrid grid(60,18,5,&CarbonThick,homedir);
  grid.addParticle(&proton2);
//   gridthick.printParticles();
  grid.fillGrids();

 /*  MeanFieldNucleusThick CupperThick(6,homedir);
  OneGlauberGrid grid3(60,18,&CupperThick,homedir);
  FastParticle pion(2, 0, 2570,0.,0.,3.,homedir);
  grid3.addParticle(&pion);
  grid3.printParticles();
  grid3.fillGrids();
*/
  
/*  GlauberGrid grid2(60,18,5,&CarbonThick,homedir);
  FastParticle proton0(0, 0, 1453,0.,0.,2.21,homedir);
  FastParticle proton1(0, 0, 1453,-113.99*DEGRTORAD,0.,2.21,homedir);
  grid2.addParticle(&proton0);
  grid2.addParticle(&proton1);
  grid2.printParticles();
  grid2.fillGrids();*/
  
  
  return 0;
}