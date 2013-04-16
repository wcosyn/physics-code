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

int main(int argc, char *argv[])
{
  double Wprime=1.25E03;//invariant mass X
  double Q2=1.8E06;
  double Ein=5765.;
  string homedir="/home/wim/Code";
  double pr=300;
  bool proton=0;
  double phi=0.;

  DeuteronCross test("paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  
  for (int i=0;i<40;i+=1){
    double costhetar=-0.975+i*0.05;
  
    double pw=0.,fsi=0.;
    test.getDeepsresult(Q2,Wprime,Ein,pr,costhetar,proton,phi,homedir,pw,fsi);
    cout << costhetar << " " << pw << " " << fsi << endl;
  }
}
