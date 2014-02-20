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
  double Q2=atof(argv[1]);
  double omega=atof(argv[2]);
  double pm=atof(argv[3]);
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
//   TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(4627.);  //Monaghan Data
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(3245.); //Dutta Data
  Cross obs(*elec,&Carbon,prec,integrator,homedir,screening,scr);
  double free=obs.getElCross(kin,2,0.,2E03)*HBARC*HBARC;
  vector<double> cross;

  
// comparison with pascal's code (input_compare_onenucl.dat input file) 
//  ../../bin/onenucl 1696000.0 851.5 190.0 with dutta beam inq
  double crossp=obs.getDiffCross(kin, 2, 0, 0, 0, 1, 1, 0.,2E04,1);
  double crosss=obs.getDiffCross(kin, 2, 0, 0, 0, 1,0,0.,2E04,1);
  cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << 
  acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << endl;
  cout << sqrt(Q2+omega*omega) << " " << crossp << " " << crosss << endl;

  MeanFieldNucleusThick Fe(3,homedir);
  TKinematics2to2 kin2("","",Fe.getMassA(),Fe.getMassA_min_proton(),MASSP,"qsquared:wlab:pklab",Q2,omega,pm);
  Cross obs2(*elec,&Fe,prec,integrator,homedir,screening,scr);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 0, 0.,2E04,1);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 1, 0.,2E04,1);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 2, 0.,2E04,1);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 3, 0.,2E04,1);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 4, 0.,2E04,1);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 5, 0.,2E04,1);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 6, 0.,2E04,1);
  
//   exit(1);
  

 /* obs.getAllDiffCross(cross,kin,2,1,thick,0.,maxEval,1);
  cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " <<  
      cross[0] << " " << cross[1] << " " << cross[4] << " " << free << " " << cross[0]/free << " " << cross[1]/free << endl;
  
      
  free=obs.getElCross(kin,2,PI,2E03)*HBARC*HBARC;
  obs.getAllDiffCross(cross,kin,2,1,thick,PI,maxEval,1);
  cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*RADTODEGR << " " << -kin.GetPklab() << " " <<  
      cross[0] << " " << cross[1] << " " << cross[4] << " " << free << " " << cross[0]/free << " " << cross[1]/free << endl;
 */     
      
  vector<double> observ;
  obs.getAllObs_tnl(observ,kin,2,1,thick,0,1.,2E04,1);
//   cout << cross[4] << " " << observ[8*4] << endl;
  for(int i=0;i<8;i++) cout << observ[i] << endl;
  for(int i=0;i<8;i++) cout << observ[32+i] << endl;
  //cout << observ[0] << " " << observ[8*(thick?4:2)] << " " << observ[3] << " " << observ[8*(thick?4:2)+3] << endl;
//   for(int i=0;i<5;i++) {
//     for(int j=0;j<8;j++) cout << observ[i*8+j] << " ";
//     cout << endl;
//   }
      //   double crossp=obs.getDiffCross(kin, 2, 0, 0, 0, 0, 1, 0.);
// // //   double crosss=obs.getDiffCross(kin, 0, 0, 0, 0, 0.);
//   double crosspsrc=obs.getDiffCross(kin, 2, 0, 1, 0, 0, 1, 0.);
// // //   double crosspct=obs.getDiffCross(kin, 0, 1, 0, 1, 0.);
// // //   double crosspsrcct=obs.getDiffCross(kin, 1, 1, 0, 1, 0.);
// // //   double crossppw=obs.getDiffCross(kin, 0, 0, 1, 1, 0.);
// // //   cout << free/HBARC/HBARC << endl;
//   cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " <<  
//       crossp << " " << crosspsrc << " " << free << " " << crossp/free << " " << crosspsrc/free << endl;
// 
//       //   cout << crossp << " " << crosss << " " << crosspsrc << " " << crosspct << " " << crosspsrcct << " " << crossppw << endl;
//   //cout << obs.getElCross(kin,0.) << " " << obs.getElCross(kin,PI) <<endl;
// //   Model test(&Carbon,0,homedir);
// //   cout << test.getMatrixEl(kin,-1,0,0,-1,1,1) << endl;
// //   
  delete elec;
  return 0;
}


