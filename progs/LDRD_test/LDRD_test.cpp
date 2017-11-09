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
  string homedir=argv[1]; //sharedir that contains input files, full path please! For me /home/wim/Code/share
  
  double Ein=12.E03; //[MeV] Incoming electron energy
  double thetain=PI-0.1; // [rad] angle of incoming beam with z-axis 
  TVector3 vecpd(0.,0.,10.E03); //deuteron momentum [MeV]
  TVector3 vecps=TVector3(100.,50.,30.); 
  vecps+=0.5*vecpd; //spectator momentum [MeV] 
  TVector3 qvec=TVector3(0.,0.,10.E03); //qvector virtual photon
  
  TVector3 veckin(Ein*sin(thetain),0.,Ein*cos(thetain)); //incoming beam momentum vector
  double Q2=0.1E06; //[MeV^2]
  
  
  LightConeKin2to2 kinLC(MASSD,Q2,MASSP,vecpd,qvec,vecps,veckin);   //constants are defined in constants.hpp

  int offshellset = 3; //no offshell, see DeuteronCross constructor [offshell effects are pretty small in this model so 3==no offshell is safe choice]
  int looplimit = 1; //fsi loop tries, see DeuteronCross constructor
  double sigmain = 40.; //sigma in mb, rescattering parameter in FSI
  double betain = 8.; //beta in GeV^-2, rescattering parameter in FSI
  double epsin = -0.5; //epsilon [], rescattering parameter in FSI
  const bool proton=0; //nucleon the virtual photon interacts with [0==neutron, 1==proton]
  double betaoffshell = 8.; //betaoffshell in GeV^-2, see DeuteronCross constructor, only used when selecting offshellset=2
  double lambda = 1.2; // GeV^2 lambda cutoff, see DeuteronCross constructor, only used when selecting offshellset=1
  string deuteronwf = "paris"; //see DeuteronCross constructor documentation, deuteron wf name
  string strucfunc = "SLAC"; //see NucleonStructure constructor documentation, structure function parametrization
  
  DeuteronCross test(deuteronwf,proton,strucfunc,sigmain,betain,epsin,betaoffshell,lambda,offshellset,looplimit);
  
  double pwLCavg = test.getavgLCCross(kinLC,1);
//   double fsiLCavg = test.getavgLCCross(kinLC,0); //fsi in LC still under construction
  double pwLCnoavg = test.getLCCross(kinLC,1);
  double pwVNAavg=test.getavgVNACross(kinLC,1);
  double fsiVNAavg=test.getavgVNACross(kinLC,0);
  double pwVNAnoavg=test.getVNACross(kinLC,1);
  double fsiVNAnoavg=test.getVNACross(kinLC,0);

  
  
  cout << "deuteron " << pwLCavg << " " << pwVNAavg <<" " << pwLCnoavg << " " 
      << pwVNAnoavg << " " << fsiVNAavg << " " << fsiVNAnoavg << endl << endl; 

  cout << "changing sigma parameter " << endl;
  //take sigma parameter from deeps fits, this should not be used for masses much bigger than the resonance region! 
  // (as we are doing here, bad example ;-) )
  double sigma_deeps=DeuteronCross::sigmaparam(kinLC.getMassX()*kinLC.getMassX(),kinLC.getQ2())*10.*HBARC*HBARC; //MeV-2 to mb!!
  cout << "new sigma [mb]: " << sigma_deeps << endl;
  test.setScatter(sigma_deeps,betain,epsin);  //set the new scattering parameters
  fsiVNAavg=test.getavgVNACross(kinLC,0);
  fsiVNAnoavg=test.getVNACross(kinLC,0);

  cout << fsiVNAavg << " " << fsiVNAnoavg << endl << endl; 
  //He3 part
  
  TVector3 sp1(0.,200.,0.);
  TVector3 sp2(100.,-100.,100.);
  TVector3 spcm=sp1+sp2;
  TVector3 beam(5000.*cos(-1.),0.,5000.*sin(-1.));
  //collinear kinematics, lab frame
  LightConeKin2to3 kinkin(MASSHE3,2.E06,MASSP,MASSP,0.,4.E03,sp1,sp2,beam,spcm);
  string he3wf = "AV18"; // he3 wf name, see He3wf constructor documentation
  He3Cross test2(he3wf,homedir,strucfunc,proton);
  cout << "He3 " << test2.getCross(kinkin) << endl;
}
