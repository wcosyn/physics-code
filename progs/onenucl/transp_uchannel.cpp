//program used to do some transparency calculations for bill Lee et al. for backward proton production
// this is A(e,e'p pi)A-1 but the pion is unobserved

//cross section is factorized, see notes writings/Notes/QE/electro-photo-transp.xopp
// we integrate over solid angles of pion in (pion,A-1) center of mass system

//elementary process is parametrized using input from Lee & Huber
// see python script  ~/Calculations/A\(ee\'p\)/backwards_T/dsigmaT_fit.py

#include <iostream>
#include <cstdlib>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <constants.hpp>
#include <TKinematics2to3.h>
#include <TElectronKinematics.h>
#include <Cross.hpp>
#include <Utilfunctions.hpp>
#include <vector>
#include <numint/numint.hpp>


void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, double Q2, double omega, int shell);
void adap_intpiCM(numint::vector_d &, double costhetapiCM, double phipiCM, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
		double Q2, double omega, double EN, double thetaN, int current, int thick, double lc_mod, double nkt_mod);


//run ./transp_uchannel  [Q2,GeV^2]  [Eout, GeV]  [theta_e, degr] [pNout, GeV] [theta_p, rad] [precision in integration]
int main(int argc, char *argv[])
{
  
  string homedir=HOMEDIR;
  //kinematics input
  double Q2=atof(argv[2])*1.E06; // input in GeV^2 [8.0,9.4,11.4,14.2]
  double Ein=10.6*1.E03;  //input in GeV 10.6 fixed
  double Eout = atof(argv[3])*1.E03; // scattered electron momentum [MeV]
  double thetae = atof(argv[4])*DEGRTORAD; // electron scatt angle
  double pN = atof(argv[5])*1.E03; // outgoing proton momentum [MeV]
  double EN = sqrt(MASSP*MASSP+pN*pN);
  double thetaN = atof(argv[6])*DEGRTORAD; // final proton scatt angle (hall lab frame)

  //calculate kinematics
  double omega=Ein-Eout;
  double x=Q2/(2.*MASSP*omega);
  double q=sqrt(Q2+omega*omega);
  //cout << Q2 << " " << omega << " " << q << " " << Q2/4./Ein/(Ein-omega) << endl;
  
  //glauber parameters
  int thick=0;
  int current=2;
  double prec=atof(argv[7]);
  bool userset=0;//atoi(argv[8]);
  double screening=0.;//atof(argv[9]);
  double lc_mod = 1.;
  double nkt_mod = 1.;
  
  //create objects for nucleus, kinematics and cross section
  MeanFieldNucleusThick Nucleus(1,homedir);  //1 for Carbon
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
  //Cross obs(*elec,&Nucleus,prec, 2, homedir, userset, screening);

  //arrays to collect results
  //all p levels  + total proton, times 3 (glauber, CT, pw)
  vector<double> totalcross(3*Nucleus.getPLevels()+3,0.); 
//   double cthmax[Nucleus.getTotalLevels()];
//   double cthmin[Nucleus.getTotalLevels()];
//   double max=-1.;

//   //determine integration limits so pm<300 MeV
//   for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
//     cthmax[shell]=1.; cthmin[shell]=-1.;
//     getBound(cthmax[shell],cthmin[shell],Nucleus,Q2,omega,shell);
// //     cout << cthmax[shell] << endl;
//     if(max<cthmax[shell]) max=cthmax[shell];
//   }
//   cout << max << endl;

  //integration structure.
  struct Ftor {

    static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
      Ftor &p = * (Ftor *) param;
      p.f(ret,x[0],x[1],p.pNucleus,p.elec,p.Q2,p.omega,p.EN,p.thetaN,p.current,p.thick,p.lc_mod, p.nkt_mod);
    }
    MeanFieldNucleusThick *pNucleus;
    TElectronKinematics *elec;
    double Q2;
    double omega;
    double EN;
    double thetaN;
    int current;
    int thick;
    double lc_mod;
    double nkt_mod;
    void (*f)(numint::vector_d &, double cospicm, double phiphicm, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
	  double Q2, double omega, double EN, double thetaN, int current, int thick, double lc_mod, double nkt_mod);

  };

  //initialize integration object + parameters needed to calculate integrals
  Ftor F;
  F.pNucleus = &Nucleus;
  F.elec=elec;
  F.Q2=Q2;
  F.omega=omega;
  F.EN=EN;
  F.thetaN=thetaN;
  F.current=current;
  F.thick=thick;
  F.lc_mod = lc_mod;
  F.nkt_mod = nkt_mod;

  numint::mdfunction<numint::vector_d,2> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

  unsigned neval = 0;
  numint::array<double,2> lower = {{-1.,0.}};
  numint::array<double,2> upper = {{1.,2.*PI}};
  
  
  F.f=adap_intpiCM;
  unsigned count=0;
   numint::cube_adaptive(mdf,lower,upper,1.E-20,1.E-02,1.E02,1.E04,totalcross,count,0);
  
  cout << Q2/1.E06 << " " << omega << " " << q  << " ";
  //proton shells transparencies
  //cout << "proton shells T s1/2, p3/2 no CT"
  for(int shell=0;shell<Nucleus.getPLevels();shell++) {
    cout << totalcross[3*shell]/totalcross[3*shell+2] << " ";
  }
  //proton shells T s1/2, p3/2 with CT"
  for(int shell=0;shell<Nucleus.getPLevels();shell++) {
    cout << totalcross[3*shell+1]/totalcross[3*shell+2] << " ";
  }
  //cout << "proton shells total cross no CT, CT, plane-wave"
  cout << totalcross[3*Nucleus.getPLevels()] << " " << totalcross[3*Nucleus.getPLevels()+1] << " " << totalcross[3*Nucleus.getPLevels()+2] << endl;

  //total proton and neutron transparency, proton no CT, with CT; neutron no CT, with CT
  cout << totalcross[3*Nucleus.getPLevels()]/totalcross[3*Nucleus.getPLevels()+2] 
        << " " << totalcross[3*Nucleus.getPLevels()+1]/totalcross[3*Nucleus.getPLevels()+2];
 
  
  delete elec;
  return 0;
}


//integrandum
void adap_intpiCM(numint::vector_d & results, double costhetapiCM, double phipiCM, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
		double Q2, double omega, double EN, double thetaN, int current, int thick, double lc_mod, double nkt_mod){
		  
  results=numint::vector_d(3*pNucleus->getPLevels()+3,0.);
  for(int shell=0;shell<pNucleus->getPLevels();shell++) {
    //7 for user input masses, otherwise fixed by pieter's processes
    // in Pieter's language
    // N is the outgoing nucleon
    // K is the pion
    // Y is the A-1 residual nucleus
    TKinematics2to3 kin("","",7,TKinematics2to3::kN,"qsquared:wlab:pn:thn:kycosthkcm:knphikcm",Q2,omega,EN,thetaN,costhetapiCM,phipiCM,
                pNucleus->getMassA(), MASSP, M_PI0, pNucleus->getMassA_min_proton()+pNucleus->getExcitation()[shell]);

    if(!kin.IsPhysical()){
        double pm=kin.GetPy();
        cout << "bla " << pm << endl;
        }
    numint::vector_d cross=numint::vector_d(thick? 5:3,0.);
    //getAllDiffCross(cross,kin,current,shell,0,0.,20000,0,0, lc_mod, nkt_mod);
    results[3*shell]+=cross[0];
    results[3*shell+1]+=cross[1];
    results[3*shell+2]+=cross[2];
    results[3*pNucleus->getPLevels()]+=cross[0];
    results[3*pNucleus->getPLevels()+1]+=cross[1];
    results[3*pNucleus->getPLevels()+2]+=cross[2];
    // cout << 2*shell << " " << 2*shell+1 << " " << 2*pNucleus->getTotalLevels() << " " << 2*pNucleus->getTotalLevels()+1 << endl;
    // cout << "0 " << shell << " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
    // << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
    // << " " << kin.GetKlab() << " " << kin.GetWlab() <<  " " << results[3*shell] << " " << results[3*shell+2] << endl;
    }

    
}




// // find cm scatt angles so that missing momentum stays below 300 MeV 
// void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, double Q2, double omega, int shell){
//   TKinematics2to2 kin("","",nucleus.getMassA(),
// 		      (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())+nucleus.getExcitation()[shell],
// 		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,(high+low)/2.);
//   double pm=kin.GetPklab();
//   if(pm<300.) low=(high+low)/2.;
//   else high=(high+low)/2.;
//   if((high-low)<1.E-04) { return;}
//   else getBound(high,low,nucleus,Q2,omega,shell);  
// }

