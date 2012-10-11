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
  double Ein=atof(argv[1]);
  double Q2=atof(argv[2]);
  double omega=atof(argv[3]);
  int shell=atoi(argv[4]);
  double pm=atof(argv[5]);
  bool screening=0;//atoi(argv[4]);
  double scr=1.;//atof(argv[5]);
  int nucleus=1;//atoi(argv[6]);
  double prec=1.E-05;//atof(argv[7]);
  int integrator=2;//atoi(argv[8]);
  int thick=1;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  
  string homedir="/home/wim/Code/share";

  MeanFieldNucleusThick Carbon(nucleus,homedir);
  TKinematics2to2 kin("","",Carbon.getMassA(),Carbon.getMassA_min_proton(),MASSP,"qsquared:wlab:pklab",Q2,omega,pm);
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
  Cross obs(*elec,&Carbon,prec,integrator,homedir,screening,scr);
  double free=obs.getElCross(kin,2,PI)*HBARC*HBARC;
  vector<double> cross;

//   obs.getAllDiffCross(cross,kin,2,1,thick,0.,maxEval,1);
  cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " "  
      /*<< cross[0] << " " << cross[1] << " " << cross[thick?4:2] << " " << free << " " << cross[0]/free << " " << cross[1]/free << endl*/;

  vector<double> observ;
  obs.getAllObs(observ,kin,2,shell,thick,PI,maxEval,1);
  //cout << observ[0] << " " << observ[8*(thick?4:2)] << " " << observ[3] << " " << observ[8*(thick?4:2)+3] << endl;
  for(int i=0;i<40;i++) cout << observ[i] << " ";
  cout << free*1.E-09 << endl;
  delete elec;
  return 0;
}


