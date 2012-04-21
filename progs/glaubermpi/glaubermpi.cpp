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

#define NROFRES 9

using namespace std;

#include <TMPI.h>
#include <stdio.h>
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
  TMPI mpi(&argc,&argv);
  int n=8;
  MeanFieldNucleusThick CarbonThick(1,homedir);
  double **results= new double*[n];
  for(int i=0;i<n;i++) results[i] = new double[NROFRES];
  for(int i=TMPI::Rank(); i<n; i += TMPI::NumberOfProcesses() ) {
    cout << TMPI::Rank() << endl;
    for(int j=0;j<NROFRES;j++) {
      FastParticle rho(4, 0, 2000+i*100.+j*20.,0.,0.,5.,145.,homedir);
      OneGlauberGrid grid = OneGlauberGrid(60,18,&CarbonThick,homedir);
      grid.addParticle(rho);
      grid.fillGrids();
      results[i][j]=abs(grid.getFsiGridFull_interp3(2.,0.,0.));
      grid.clearParticles();
    }
  }
  TMPI::GatherResults(n,results);

  TMPI::SilenceSlaves();
  for(int i=0; i<n; ++i){
    for (int j=0;j<NROFRES;j++) cout << results[i][j] << " ";
    cout << endl;
  }
  TMPI::SilenceSlaves(false);
  
  for(int i=0;i<n;i++) delete [] results[i];
  delete [] results;
  //return 0;
}