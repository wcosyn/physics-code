#include <iostream>
#include <cstdlib>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <FsiCorrelator.hpp>
#include <FastParticle.hpp>
#include <constants.hpp>
#include <GlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include <TKinematics2to3.h>
#include <TElectronKinematics.h>
#include <DoubleNCross.hpp>


int main(int argc, char *argv[])
{
  
  string homedir="/home/wim/Code";
  int shellindex1=0;
  int shellindex2=2;
  bool corr=atoi(argv[1]);
  cout << corr << endl;
  MeanFieldNucleusThick Carbon(1,homedir);
  TKinematics2to3WithLabAngles kin("","",7,TKinematics2to3::kN,"qsquared:wlab:pn:thn:thk:phik",1631389.67942,580,1000,15*DEGRTORAD,31.17*DEGRTORAD,PI,
		      Carbon.getMassA(), MASSP, MASSN, Carbon.getMassA_min_pn()+Carbon.getExcitation()[shellindex1]+Carbon.getExcitation()[shellindex2]);
  
  cout << kin.GetNumberOfSolutions() << " " << kin.GetMd() << " " << kin.GetMn() << " " << kin.GetMk() << " " << kin.GetMy() << endl;
  cout << kin.GetNumberOfSolutions() << " " << kin.GetPn() << " " << kin.GetPk() << " " << kin.GetPy() << endl;
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(4627.);
  DoubleNCross obs(*elec,&Carbon,homedir);
  double cross=obs.getDiffCross(kin,0,0,0,corr,shellindex1,shellindex2,0.);
  delete elec;
  return 0;
}


