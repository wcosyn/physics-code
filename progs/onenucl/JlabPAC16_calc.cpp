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

int main(int argc, char *argv[])
{

  double Ebeam=atof(argv[1]);
  double Q2=atof(argv[2]);
  double omega=atof(argv[3]);
  double pm=atof(argv[4]);
  double phi=atof(argv[5]);
  string sharedir=argv[6];
  phi=acos(phi);
  bool screening=0;//atoi(argv[4]);
  double scr=1.;//atof(argv[5]);
  //int nucleus=1;//atoi(argv[6]);
  double prec=1.E-05;//atof(argv[7]);
  int integrator=2;//atoi(argv[8]);
  int thick=1;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  
//   string homedir="/home/wim/Code/share";
//   cout << "bla" << endl;
  MeanFieldNucleusThick Fe54(MeanFieldNucleus::Fe54,sharedir);
  MeanFieldNucleusThick Ca40(MeanFieldNucleus::Ca40,sharedir);
  MeanFieldNucleusThick Ca48(MeanFieldNucleus::Ca48,sharedir);

  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ebeam);
  cout << Ebeam << " " << Q2 << " " << omega << " " << pm << " " << phi << " ";
  vector<Cross> obs;
  obs.push_back(Cross(*elec,&Fe54,prec,integrator,sharedir,screening,scr));
  obs.push_back(Cross(*elec,&Ca40,prec,integrator,sharedir,screening,scr));
  obs.push_back(Cross(*elec,&Ca48,prec,integrator,sharedir,screening,scr));
  cout << endl;
  for(unsigned int i=0;i<obs.size();i++){
    double cross=0.;
    double pw=0.;
    for(int plevel=0;plevel<obs[i].getPnucl()->getPLevels();plevel++){
      cout << i << " " << plevel << endl;
      TKinematics2to2 kin("","",obs[i].getPnucl()->getMassA(),obs[i].getPnucl()->getMassA_min_proton()+(obs[i].getPnucl()->getExcitation())[plevel]
                          ,MASSP,"qsquared:wlab:pklab",Q2,omega,pm);
//       cout << kin.GetPYlab() << " " << acos(kin.GetCosthYlab()) << endl;
      cross+=obs[i].getDiffCross(kin,2,thick,0,0,0,plevel,phi,maxEval,1);
      pw+=obs[i].getDiffCross(kin,2,thick,0,0,1,plevel,phi,maxEval,1);
    }
    cout << cross << " " << pw << " " << cross/pw << endl;
  }
  cout << endl;

}  