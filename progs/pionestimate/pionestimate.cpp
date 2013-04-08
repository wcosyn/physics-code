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
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <Cross.hpp>
#include <TRotation.h>


//run ./onenucl [Q2 [MeV^2]] [omega] [missing momentum]
int main(int argc, char *argv[])
{
  int nucleus=atoi(argv[1]);
  double prec=(nucleus==1? 1.E-05: 1.E-06);//atof(argv[7]);
  int integrator=2;//atoi(argv[8]);
  int thick=1;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  bool ppion=atoi(argv[2]);
  
  string homedir="/home/wim/Code/share";

  MeanFieldNucleusThick nucl(nucleus,homedir);
  DistMomDistrGrid *pgrid[nucl.getTotalLevels()];
  double prho=3E03;
  double ppi=sqrt((MASSRHO*MASSRHO+prho*prho)/4.-MASSPI*MASSPI);
  double ppiz=prho/2.
  double thetapi=2.*acos(ppiz/ppi);
  cout << thetapi*RADTODEGR << endl;
  
  
  if(ppion){
    GlauberGridThick grid = GlauberGridThick(120,36,5,&nucl,prec,integrator,homedir),
    FastParticle pion1(2, 0, ppi,0.,0.,3.,0.,homedir);
    FastParticle pion2(3, 0, ppi,0.,0.,3.,0.,homedir);
    grid.clearParticles();
    grid.addParticle(pion1);
    grid.addParticle(pion2);
    grid.updateGrids();
    grid.clearKnockout();
    for(int i=0;i<nucl.getTotalLevels();i++){
      pgrid[i] = new DistMomDistrGrid(i, 300., 30,20,5,&grid,1.E-03,2,5E04,0.,homedir);
      TRotation rot;
      rot.Rotate(0.,TVector3(0.,1.,0.));
      pgrid.updateGrids(grid,i,rot);
      pgrid.printRho_grid(0);
    }
  }
  else{
    OneGlauberGrid grid = OneGlauberGrid(120,36,&nucl,prec,integrator,homedir)
    FastParticle rhopion(7, 0, prho,0.,0.,3.,0.,homedir);
    grid.clearParticles();
    grid.addParticle(rhopion);
    grid.updateGrids();
    grid.clearKnockout();
    for(int i=0;i<nucl.getTotalLevels();i++){
      pgrid[i] = new DistMomDistrGrid(i, 300., 30,20,5,&grid,1.E-03,2,5E04,0.,homedir);
      TRotation rot;
      rot.Rotate(0.,TVector3(0.,1.,0.));
      pgrid.updateGrids(grid,i,rot);
      pgrid.printRho_grid(0);
    }
  }
  
