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
  FastParticle proton2(0, 0, 4878,0.,0.,8.,0.,homedir);
  proton2.printParticle();
  GlauberGridThick gridthick(10,9,5,&CarbonThick,prec,integr,homedir);
  gridthick.addParticle(proton2);
  gridthick.printParticles();
  gridthick.updateGrids();
  gridthick.printFsi_src_ct_grid();
  
}