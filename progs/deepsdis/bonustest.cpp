//bonus data test program, test formulas from their MC, my cross section
// generate ratios to compare with data (using scattering parameter values from fit)

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include <constants.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <DeuteronStructure.hpp>
#include <DeuteronMomDistr.hpp>
#include <DeuteronCross.hpp>
#include "bonusdata.h"

int main(int argc, char *argv[])
{
  double Ebeam = data::Ebeam[atoi(argv[1])];
  int Qindex = atoi(argv[2]); 
  int Windex = atoi(argv[3]); 
  int pindex = atoi(argv[4]);
  int offshellset = atoi(argv[5]);
  double sigmain = atof(argv[6]); //sigma in mb
  double betain = atof(argv[7]); //beta in GeV^-2
  
  bool proton=0;
  double phi=0.;

  double Wprime = 0.5*(data::W[Windex]+data::W[Windex+1]);
  double Q2 = 0.5*(data::Q2[Qindex]+data::Q2[Qindex+1]);
  double pr = 0.5*(data::ps[pindex]+data::ps[pindex+1]);
  DeuteronCross test("paris",proton,"CB",sigmain,betain,-0.5,8.,1.2,4,1E03);
//   DeuteronCross test("paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  
  for (int i=0;i<10;i+=1){
    double costhetar=-0.9+i*0.2;
    double res=test.getBonusMCresult(Q2,Wprime,Ebeam,pr,costhetar,proton);
    cout << costhetar << " " << res << endl;
  }
}
