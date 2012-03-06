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
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <Cross.hpp>


//run ./onenucl [Q2 [MeV^2]] [omega] [missing momentum]
int main(int argc, char *argv[])
{
  
  string homedir="/home/wim/Code";

  MeanFieldNucleusThick Carbon(1,homedir);
  TKinematics2to2 kin("","",Carbon.getMassA(),Carbon.getMassA_min_proton(),MASSP,"qsquared:wlab:pklab",atof(argv[1]),atof(argv[2]),atof(argv[3]));
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(4627.);
  Cross obs(*elec,&Carbon,homedir);
  double cross=obs.getDiffCross(kin, 0, 0, 0, 1, 0.);
  double free=obs.getElCross(kin,0.)*HBARC*HBARC;
  cout << free/HBARC/HBARC << endl;
  cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " <<  cross << " " << free << " " << cross/free << endl;
  //cout << obs.getElCross(kin,0.) << " " << obs.getElCross(kin,PI) <<endl;
//   Model test(&Carbon,0,homedir);
//   cout << test.getMatrixEl(kin,-1,0,0,-1,1,1) << endl;
//   
  delete elec;
  return 0;
}


