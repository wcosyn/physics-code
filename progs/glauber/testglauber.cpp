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
#include <GlauberDecayGridThick.hpp>
#include <TMFSpinor.hpp>

int main(int argc, char *argv[])
{
  
  string homedir=argv[1];

//   MeanFieldNucleusThick PbThick(4,homedir);
//   FastParticle proton2(0, 0, 2000.,0.,0.,3.,0.,homedir);
//   GlauberGridThick onegrid(60,18,5,&PbThick,1.E-03,homedir);
//   onegrid.addParticle(proton2);
//   onegrid.updateGrids();
//   onegrid.printFsi_grid();
//   proton2.setScreening(1.03);
//   onegrid.clearParticles();
//   onegrid.addParticle(proton2);
//   onegrid.updateGrids();
//   onegrid.printFsi_grid();
  
  
  
  MeanFieldNucleusThick CarbonThick(1,homedir);
  complex<double> wave[4];
  CarbonThick.getWaveFunction(wave,0, -1, 2.,0.5,1.);
  cout << wave[0] << " " << wave[1] << " " << wave[2] << " " << wave[3] << endl;
  cout << TMFSpinor(CarbonThick,0,-1,2.,0.5,1.) << endl;
  CarbonThick.getWaveFunction(wave,1, 1, 1.,0.5,1.);
  cout << wave[0] << " " << wave[1] << " " << wave[2] << " " << wave[3] << endl;
  cout << TMFSpinor(CarbonThick,1,1,1.,0.5,1.) << endl;
  //FastParticle proton(0, 1, 1234,0.,0.,3.,0.,homedir);
  FastParticle proton2(0, 0, 4878,0.,0.,8.,0.,homedir);
  proton2.printParticle();
  GlauberGridThick gridthick(60,18,5,&CarbonThick,homedir);
  gridthick.addParticle(proton2);
  gridthick.printParticles();
  gridthick.updateGrids();
  gridthick.printFsi_src_ct_grid();
  GlauberGrid grid(60,18,5,&CarbonThick,homedir);
  grid.addParticle(proton2);
  grid.printParticles();
  grid.updateGrids();
  grid.printFsi_grid();
  OneGlauberGrid onegrid(60,18,&CarbonThick,homedir);
  onegrid.addParticle(proton2);
  onegrid.printParticles();
  onegrid.updateGrids();
  onegrid.printFsi_ct_grid();

  
  
  
  //   
//   FastParticle rho(4, 0, 2356,0.,0.,5.,145.,homedir);
//   GlauberDecayGridThick gridthick(3,3,2,&CarbonThick,homedir);  
//   gridthick.addParticle(rho);
//   gridthick.printParticles();
//   gridthick.fillGrids();
//   gridthick.printFsi_grid();
//   gridthick.printFsi_src_grid();
//   gridthick.printFsi_ct_grid();
//   gridthick.printFsi_src_ct_grid();
//   gridthick.printFsi_decay_grid();
//   gridthick.printFsi_src_decay_grid();
//   gridthick.printFsi_ct_decay_grid();
//   gridthick.printFsi_src_ct_decay_grid();
//   cout << gridthick.getNumber_of_grids() << endl;
  
//   FastParticle pion(2, 0, 2356,0.,0.,5.,0.,homedir);
//   OneGlauberGrid grid(60,18,&CarbonThick,homedir);  
//   grid.addParticle(pion);
//   grid.printParticles();
//   grid.fillGrids();
//   grid.printFsi_grid();
//   grid.printFsi_ct_grid();
  
  
//   FastParticle rho(4, 0, 2356,0.,0.,5.,145.,homedir);
//   FastParticle rho2(4, 0, 2056,0.,0.,5.,145.,homedir);
//   FastParticle rho3(4, 0, 2356,0.,0.,3.,145.,homedir);
//   FastParticle rho4(4, 0, 2056,0.,0.,3.,145.,homedir);
//   GlauberDecayGridThick gridthick(3,3,2,&CarbonThick,homedir);  
//   gridthick.addParticle(rho);
//   gridthick.updateGrids();
// //   gridthick.print_grid(0);
// //   gridthick.print_grid(1);
//   gridthick.clearParticles();
//   gridthick.addParticle(rho2);
//   gridthick.updateGrids();
//   gridthick.print_grid(0);
//   gridthick.print_grid(1);
//   gridthick.printFsi_ct_grid();
//   gridthick.clearParticles();
//   gridthick.addParticle(rho4);
//   gridthick.updateGrids();
//   gridthick.printFsi_ct_grid();
  
  return 0;
}