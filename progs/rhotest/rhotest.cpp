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
#include <DistMomDistrGrid.hpp>
#include <TMFSpinor.hpp>

int main(int argc, char *argv[])
{
  
  string homedir="/home/wim/Code/share";

  MeanFieldNucleusThick CarbonThick(3,homedir);
  FastParticle rho(4, 0, 2220,0.,0.,1.1,145.,homedir);
  GlauberDecayGridThick gridthick(60,20,5,&CarbonThick,homedir);  
  gridthick.addParticle(rho);
  //gridthick.printParticles();
  gridthick.updateGrids();
//   gridthick.print_grid(0);
//   gridthick.print_grid(1);
  //   cout << gridthick.getNumber_of_grids() << endl;
   DistMomDistrGrid Momgrid(9, 400., 30,20,5,&gridthick,homedir);
   Momgrid.printRhopw_grid();
   Momgrid.fillGrids();
   cout << endl << endl;
   Momgrid.printRho_grid(0);
   
//   MeanFieldNucleusThick CarbonThick(1,homedir);
//   FastParticle rho(4, 0, 2710,0.,0.,1.4,145.,homedir);
//   GlauberDecayGridThick gridthick(10,10,3,&CarbonThick,homedir);  
//   gridthick.addParticle(rho);
//   //gridthick.printParticles();
//   gridthick.updateGrids();
// //   gridthick.print_grid(0);
// //   gridthick.print_grid(1);
//   //   cout << gridthick.getNumber_of_grids() << endl;
//    DistMomDistrGrid Momgrid(0, 400., 2,1,1,&gridthick,homedir);
//    Momgrid.printRhopw_grid();
//    Momgrid.fillGrids();
//    cout << endl << endl;

//    MeanFieldNucleusThick CarbonThick(1,homedir);
//   FastParticle rho(4, 0, 2010,0.,0.,1.4,145.,homedir);
//   FastParticle rho2(4, 0, 2210,0.,0.,1.4,145.,homedir);
//   FastParticle rho3(4, 0, 2010,0.,0.,2.4,145.,homedir);
//   GlauberDecayGridThick gridthick(5,5,3,&CarbonThick,homedir);  
//   gridthick.addParticle(rho);
//   //gridthick.printParticles();
//   gridthick.updateGrids();
// //   gridthick.print_grid(0);
// //   gridthick.print_grid(1);
//   //   cout << gridthick.getNumber_of_grids() << endl;
//    DistMomDistrGrid Momgrid(0, 400., 2,1,1,&gridthick,homedir);
//    Momgrid.fillGrids();
//    Momgrid.printRhopw_grid();
//    Momgrid.printRho_grid(0);
//    Momgrid.printRho_grid(4);
// 
//    gridthick.clearParticles();
//    gridthick.addParticle(rho3);
//    gridthick.updateGrids();
// //   gridthick.print_grid(0);
// //   gridthick.print_grid(1);
//    Momgrid.updateGrids(&gridthick,1);
//    Momgrid.printRhopw_grid();
//    Momgrid.printRho_grid(0);
//    Momgrid.printRho_grid(4);  
   
   return 0;
}