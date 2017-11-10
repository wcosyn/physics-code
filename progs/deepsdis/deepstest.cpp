//program that calculates reduced structure funtions for spectator tagging as they were presented in the deeps data

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <DeuteronStructure.hpp>
#include <DeuteronMomDistr.hpp>
#include <DeuteronCross.hpp>
#include <He3Cross.hpp>
#include <LightConeKin2to3.hpp>
#include <NuclStructure.hpp>

int main(int argc, char *argv[])
{



  int Qindex = atoi(argv[1]); // parse from argv or something
  int Windex = atoi(argv[2]); // parse from argv or something
  int pindex = atoi(argv[3]);
  int offshellset = atoi(argv[4]);
  int looplimit = atoi(argv[5]);
  double sigmain = atof(argv[6]); //sigma in mb
  double betain = atof(argv[7]); //beta in GeV^-2
  
  double Qarray[2]={1.8E06,2.8E06};
  double Warray[5]={1.25E03,1.5E03,1.73E03,2.02E03,2.4E03};
  double prarray[6]={300.,340.,390.,460.,560.,100.};
//   double Wprime=1.25E03;//invariant mass X
//   double Q2=1.8E06;
  double Ein=5765.;
//   string homedir="/home/wim/Code";
//   double pr=300;
  bool proton=0;
  double phi=0.;

  double Wprime = Warray[Windex];
  double Q2 = Qarray[Qindex];
  double pr = prarray[pindex];
  DeuteronCross test("paris",proton,"SLAC",sigmain,betain,-0.5,8.,1.2,offshellset,looplimit);
//   DeuteronCross test("paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  
  
//   for(int i=1;i<50;i++){
//     double x=0.02*i;
//     for(int j=0;j<10;j++){
//       double Q2=1.E06*(1+j);
//       cout << x << " " << Q2/1.E06 << " " << NuclStructure::getG1plusg2_grsv2000(0,x,Q2) << " " << NuclStructure::getG1_grsv2000(0,x,Q2) << endl;
//     }
//   }
//   exit(1);
//   
  
  for (int i=0;i<40;i+=1){
    double costhetar=-0.975+i*0.05;
  
    double pw=0.,fsi=0.,azz=0.,azzfsi=0.;
    test.getDeepsresult(Q2,Wprime,Ein,pr,costhetar,proton,pw,fsi);
//     cout << costhetar << " " << pw << " " << fsi << endl;
//     test.getDeepsAzz(Q2,Wprime,Ein,pr,costhetar,proton,azz,azzfsi,pw,fsi);
//     cout << costhetar << " " << azz << " " << azzfsi << " " << pw << " " << fsi << endl;
  }
//   TVector3 sp1(0.,200.,0.);
//   TVector3 sp2(100.,-100.,100.);
//   TVector3 spcm=sp1+sp2;
//   TVector3 beam(5000.*cos(-1.),0.,5000.*sin(-1.));
//   LightConeKin2to3 kinkin(MASSHE3,2.E06,MASSP,MASSP,0.,4.E03,sp1,sp2,beam,spcm);
//   string a="AV18";
//   He3Cross test2("AV18","/home/wim/Code/share","SLAC",proton);
//   cout << test2.getCross(kinkin) << endl;
}
