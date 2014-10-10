#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <DeuteronCross.hpp>
#include <He3Cross.hpp>
#include <LightConeKin2to3.hpp>
#include <LightConeKin2to2.hpp>

int main(int argc, char *argv[])
{
  string homedir=argv[1]; //sharedir, full path please! For me /home/wim/Code/share
  
  double Ein=12.E03; //[MeV] Incoming electron energy
  double thetain=PI-0.1; // [rad] angle of incoming beam with z-axis 
  TVector3 vecpd(0.,0.,10.E03); //deuteron momentum [MeV]
  TVector3 vecps=TVector3(100.,50.,30.); //spectator momentum [MeV]
  vecps+=0.5*vecpd;
  TVector3 qvec=TVector3(0.,0.,10.E03);
  
  TVector3 veckin(Ein*sin(thetain),0.,Ein*cos(thetain));
  double Q2=0.1E06; //[MeV^2]
  
  
  LightConeKin2to2 kinLC(MASSD,Q2,MASSP,vecpd,qvec,vecps,veckin);   //constants are defined in constants.hpp
  LightConeKin2to2 kinLC2(MASSD,Q2,MASSP,vecpd,qvec,vecps,veckin);   //constants are defined in constants.hpp

  int offshellset = 3; //no offshell, see DeuteronCross constructor
  int looplimit = 1; //fsi loop tries, see DeuteronCross constructor
  double sigmain = 40.; //sigma in mb, rescattering parameter in FSI
  double betain = 8.; //beta in GeV^-2, rescattering parameter in FSI
  double epsin = -0.5; //epsilon [], rescattering parameter in FSI
  const bool proton=0; //nucleon the virtual photon interacts with
  double betaoffshell = 8.; //betaoffshell in GeV^-2, see DeuteronCross constructor
  double lambda = 1.2; // GeV^2 lambda cutoff, see DeuteronCross constructor
  string deuteronwf = "paris"; //see DeuteronCross constructor, deuteron wf name
  string strucfunc = "SLAC"; //see NucleonStructure constructor, structure function parametrization
  
  DeuteronCross test(deuteronwf,proton,strucfunc,sigmain,betain,epsin,betaoffshell,lambda,offshellset,looplimit);
  DeuteronCross test22=test;
  
//   double pwLCavg = test.getavgLCCross(kinLC,1);
// //   double fsiLCavg = test.getavgLCCross(kinLC,0);
//   double pwLCnoavg = test.getLCCross(kinLC,1);
//   double pwVNAavg=test.getavgVNACross(kinLC2,1);
//   double fsiVNAavg=test.getavgVNACross(kinLC,0);
//   double pwVNAnoavg=test22.getVNACross(kinLC2,1);

//   cout << "deuteron " << pwLCavg << " " << pwVNAavg <<" " << pwLCnoavg << " " << pwVNAnoavg << " " << fsiVNAavg << endl;
  
  TVector3 sp1(0.,200.,0.);
  TVector3 sp2(100.,-100.,100.);
  TVector3 spcm=sp1+sp2;
  TVector3 beam(5000.*cos(-1.),0.,5000.*sin(-1.));
  //collinear kinematics, lab frame
  LightConeKin2to3 kinkin(MASSHE3,2.E06,MASSP,MASSP,0.,4.E03,sp1,sp2,beam,spcm);
  LightConeKin2to3 kinkin2=kinkin;
  string he3wf = "AV18"; // he3 wf name, see He3wf constructor
  He3Cross test2(he3wf,homedir,strucfunc,proton);
  He3Cross test33=test2;
  cout << "He3 " << test2.getCross(kinkin) << " " << test33.getCross(kinkin) << endl;
}
