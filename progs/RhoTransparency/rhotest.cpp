//program used to test Rho glauber grids and distorted momentum distributions used in the cross section calculations

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
#include <TRotation.h>

int main(int argc, char *argv[])
{
  
  string homedir=HOMEDIR;
  cout << homedir << endl;
// 
//   MeanFieldNucleusThick CarbonThick(3,homedir);
//   FastParticle rho(4, 0, 2220,0.,0.,1.1,145.,homedir);
//   GlauberDecayGridThick gridthick(60,20,5,&CarbonThick,homedir);  
//   gridthick.addParticle(rho);
//   //gridthick.printParticles();
//   gridthick.updateGrids();
// //   gridthick.print_grid(0);
// //   gridthick.print_grid(1);
//   //   cout << gridthick.getNumber_of_grids() << endl;
//    DistMomDistrGrid Momgrid(9, 400., 30,20,5,&gridthick,homedir);
//    Momgrid.printRhopw_grid();
//    Momgrid.fillGrids();
//    cout << endl << endl;
//    Momgrid.printRho_grid(0);
   
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
   
  MeanFieldNucleusThick FeThick(3,homedir);
  FastParticle rho(4, 0, 2340,0.,0.,1.0,145.,1.,1.,homedir);
  GlauberDecayGridThick gridthick(60,20,5,&FeThick,2,1.E-04,homedir);  
  gridthick.addParticle(rho);
  //gridthick.printParticles();
  gridthick.updateGrids();
//   gridthick.print_grid(0);
//   gridthick.print_grid(1);
  //   cout << gridthick.getNumber_of_grids() << endl;
   DistMomDistrGrid** pdistgrid = new DistMomDistrGrid*[FeThick.getTotalLevels()];
   for(int i=0;i<FeThick.getTotalLevels();i++){
     pdistgrid[i] = new DistMomDistrGrid(i, 400., 30,20,5,&gridthick,1.E-03,2,1E04,0.,homedir);
   }
   for(int i=0;i<=400;i+=10){
     double p=i;
     double tot =0.,totpw=0.;
//      cout << p << " "; 
     for(int shell=0;shell<FeThick.getTotalLevels();shell++){
      TRotation rot;
      rot.Rotate(0.,TVector3(0.,1.,0.));
         pdistgrid[shell]->updateGrids(&gridthick,shell,rot);
	 double pw = pdistgrid[shell]->getRhopwGridFull_interp(p);
	 double fsi = pdistgrid[shell]->getRhoGridFull_interp3(0, p, 1., 0.);
	 tot+=fsi;
	 totpw+=pw;
// 	 cout << fsi << " " << pw << " " << fsi/pw << " ";
     }
     cout << tot/totpw << endl;
   }
  
   return 0;
}