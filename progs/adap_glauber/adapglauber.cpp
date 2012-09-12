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
  
  int nucleus=atof(argv[1]);
  double prec=atof(argv[2]);
  int integr=atoi(argv[3]);
  string homedir=argv[4];

  MeanFieldNucleusThick CarbonThick(nucleus,homedir);
  FastParticle proton2(0, 0, 1503.85,0.,0.,1.696,0.,homedir);
//   proton2.printParticle();
  GlauberGridThick grid(60,18,5,&CarbonThick,prec,integr,homedir);
//   OneGlauberGrid grid(60,18,&CarbonThick,prec,integr,homedir);
  grid.addParticle(proton2);
//   gridthick.printParticles();
  grid.updateGrids();
  grid.printFsi_grid();
  grid.printFsi_ct_grid();
  
}