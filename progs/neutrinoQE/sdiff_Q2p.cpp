// Program that computes single differential cross sections for QE neutrino nucleus scattering
// Differential to Q2p, computed from final proton kinematics
// Can be used for different experiments, Miniboone and minerva kinematics currently implemented, T2K will be added as well

//run like
// > sdiff_Q2p [Q2p in GeV^2] [integration algorithm for flux integration (see around line 400 for options)] 
//  [max initial nucl momentun in MeV] [min final nucl momentum in MeV] [plane wave 1 or not 0] [ROMEA 1 or not 0 in fsi]
// [experiment "miniboone", "minerva", "t2k"] [lepton id "electron", "muon", "tau"] [sharedir]

// example for plane wave calculation for miniboone
// > sdiff_Q2 0.1 3 500. 200. 1 0 miniboone ~/physics-code/share/



#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cassert>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <constants.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <TLeptonKinematics.h>
#include <WeakQECross.hpp>
#include <Utilfunctions.hpp>
//headers for integration 
#include <numint/numint.hpp>
#include <numint/numint2Cuba.hpp>

#include "fluxes.h"

//const double massmu = 105.6583715;




//integration struct 
struct Ftor {  //Carbon

  static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
    Ftor &p = * (Ftor *) param;
    p.f(ret,x[0],x[1],x[2],*p.pNucleus,p.lepton_id,p.current,p.Q2p,
        p.prec,p.integrator,p.homedir,p.maxEval,p.charged,p.screening,p.enable_romea,
        p.max_initial_nucl_mom,p.min_final_nucl_mom,p.pw, p.exp
       );
  }
  MeanFieldNucleusThick *pNucleus;
  string lepton_id;
  int current;
  double Q2p;
  double prec;
  int integrator;
  string homedir;
  int maxEval;
  bool charged;
  bool screening;
  bool enable_romea;
  double max_initial_nucl_mom;
  double min_final_nucl_mom;
  int pw;
  int maxEvalweakamp;
  string exp;
  void (*f)(numint::vector_d &, double E_omega, double costheta_mu, double E_in,
            MeanFieldNucleusThick &pNucleus, string lepton_id, int current, double Q2p,
            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening, bool enable_romea,
            double max_initial_nucl_mom, double min_final_nucl_mom, int pw, string exp
           );

};

struct FtorH {  //Hydrogen

  static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
    FtorH &p = * (FtorH *) param;
    p.f(ret,x[0],p.Q2p,p.charged,p.maxbeam,p.exp,p.lepton_id);
  }
  double Q2p;
  bool charged;
  double maxbeam;
  string exp;
  string lepton_id;
  void (*f)(numint::vector_d &, double E_out, double Q2p, bool charged, double maxbeam, string exp, string lepton_id);

}; 



//integration for Hydrogen
void int_hydr(numint::vector_d &, double E_out,double Q2p,bool charged, double maxbeam, string exp, string lepton_id);

//integration for nucleus
void adap_intPm(numint::vector_d &, double omega, double costhetacm, double E_in,
		            MeanFieldNucleusThick &pNucleus, string lepton_id, int current, double Q2p,
		            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening, bool enable_romea,
                            double max_initial_nucl_mom, double min_final_nucl_mom, int pw, string exp
               );



//determine boundaries of kinematics
double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double En, int shell, double max_initial_nucl_mom, double min_final_nucl_mom);
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);



int main(int argc, char *argv[])
{
  double Q2p=atof(argv[1])*1.E06;   // give input in GeV^2
  bool screening=0;//atoi(argv[4]);
  int nucleus=1;//atoi(argv[6]);                     
  double prec=1.E-04;//atof(argv[7]);   //1.E-5
  int integrator=2;//
  int fluxintegrator=atoi(argv[2]);
//int thick=0;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  bool charged=1;
  int current=2;
  double max_initial_nucl_mom=atof(argv[3]);
  double min_final_nucl_mom=atof(argv[4]);
  int pw=atoi(argv[5]); //1 is plane-wave, 0 is with FSI
  bool Pauli = 1;
  bool enable_romea = atoi(argv[6]);
  if(pw==1) enable_romea=0; //no FSI if calculating the pw cross section...
  
  string exp=argv[7]; // possibilities "miniboone", "minerva", "t2k"
  string lepton_id=argv[8]; //possiblities "electron", "muon", "tau"
  string homedir=argv[9];   //"/home/wim/Code/share";

  MeanFieldNucleusThick Nucleus(nucleus,homedir);
  double massAmin1=(Nucleus.getMassA_min_proton()+Nucleus.getMassA_min_neutron())/2.;
  
  vector<double> avgcross(2,0.); //neutrino and antineutrino
  
  double minbeam=0., maxbeam=0.; //different beam energy intervals for the different experiments
  double neutnorm=1.,aneutnorm=1.; // flux normalisations

  if(!exp.compare("miniboone")) { maxbeam=3.E03;}
  else if(!exp.compare("minerva")) { 
    minbeam=1.5E03; 
    maxbeam=10.E03;
    if(!lepton_id.compare("muon")){
      neutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_nu_muon_FHC_flux,500,20);
      aneutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_anu_muon_RHC_flux,500,20);
    }
    else if(!lepton_id.compare("electron")){
      neutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_nu_elec_FHC_flux,500,20);
      aneutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_anu_elec_FHC_flux,500,20);      
    }
  }  
  else if(!exp.compare("t2k")) { maxbeam=2.E02; }  //still to implement!
  else {cerr << "invalid experiment name chosen" << endl << "Choose either miniboone, minerva or t2k" << endl; assert(1==0);}
  
  double leptonmass=0.;
  int flav=999;
  if(!lepton_id.compare("electron")){
    leptonmass=TLeptonKinematics::masse;
    flav=0;
  } 
  else if (!lepton_id.compare("muon")){
    leptonmass=TLeptonKinematics::massmu;
    flav=1;
  } 
  else if (!lepton_id.compare("tau")){
    leptonmass=TLeptonKinematics::masstau;
    flav=2;
  } 
  else {cerr << "invalid lepton name chosen" << endl << "Choose either electron, muon or tau" << endl; assert(1==0);}


  double En_out=(Q2p+pow(MASSn-34.,2.)+MASSn*MASSn)/2./(MASSn-34.); //outgoing nucleon momentum from Q2p definition
  double pn_out=sqrt(En_out*En_out-MASSn*MASSn); //outgoing nucleon momentum

  //setting up variables dat will determine integration limits
  double E_out_max=minbeam;
  double E_out_min=maxbeam;
  double E_low=maxbeam;
  double E_high=minbeam;
  double omega_low=maxbeam;
  double omega_hi=minbeam;
  double max=-1.,min=1.; //overall min and max lepton theta angles


  for(int j=0;j<=100;j++){ //loop over T_mu
    //   cout << j << "/100" << endl;
    double E_out=(maxbeam-leptonmass-1.E-03)/100.*j+leptonmass+1.E-03;
    double p_out=sqrt(E_out*E_out-leptonmass*leptonmass);
    double E_in_lo=E_out+En_out-Nucleus.getMassA()+massAmin1;
    if(E_in_lo<maxbeam){
        for(int i=0;i<=100;i++){
            double E_in=E_in_lo+(maxbeam-E_in_lo)*1.E-02*i;
            double omega=E_in-E_out;
            TLeptonKinematics *lepton=TLeptonKinematics::CreateWithBeamEnergy(static_cast<TLeptonKinematics::Lepton>(flav),E_in);
            for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
                double tempmax=-1., tempmin=1.; //temp limits for costhetamu
                
                //anything above 500 MeV contribution will be negligible
                if(getBound(tempmax,tempmin,Nucleus,*lepton,E_in,E_out,En_out,shell,max_initial_nucl_mom,min_final_nucl_mom)<max_initial_nucl_mom){ 
                    if(E_in<E_low)  E_low=E_in;
                    if(E_in>E_high) E_high=E_in;
                    if(max<tempmax) max=tempmax;
                    if(min>tempmin) min=tempmin;
                    if(E_out_max<E_out) E_out_max=E_out;
                    if(E_out_min>E_out) E_out_min=E_out;
                    if(omega>omega_hi) omega_hi=omega;
                    if(omega<omega_low) omega_low=omega;
                }

                delete lepton;
            }
        }
    }
  }
}

double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double En_out, int shell, double max_initial_nucl_mom, double min_final_nucl_mom){
  
  double omega=E_in-E_out;
//   double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
//       -lepton.GetLeptonMass()*lepton.GetLeptonMass();
//   double tempmin=max_initial_nucl_mom;
//   double costhmin=-1.;
//   for(int i=0;i<=100;i++){ //loop over all com angles to find com scattering angle where initial nucl momentum is minimal
//     TKinematics2to2 kin("","",nucleus.getMassA(),
// 			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
// 			+nucleus.getExcitation()[shell],
// 			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.+i*0.02);
//     double pm=kin.GetPklab(); //initial nucleon momentum
//     if(pm<tempmin&&kin.GetPYlab()>min_final_nucl_mom) {tempmin=pm; costhmin=-1.+i*0.02;}
//   }
//   //   cout << omega << " " << Q2/1.E06 << " " << omega/(2.*MASSP*omega) << " " << high << " " << low << " " << pm << endl;
//   double temptemp=costhmin;//com costheta value where the min initial nucl momentum is
// //   cout << costhmin << " " << tempmin << endl;
//   //if it is below the limit we consider, we find the interval in angles where it is below
//   if(tempmin<max_initial_nucl_mom){
//     TKinematics2to2 kin1("","",nucleus.getMassA(),
// 			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
// 			+nucleus.getExcitation()[shell],
// 			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,1.);
//     high=1;
//     if(kin1.GetPklab()>max_initial_nucl_mom) high=getMax(high,temptemp,nucleus,lepton,Q2,omega,shell,max_initial_nucl_mom);
//     temptemp=costhmin;
//     TKinematics2to2 kinmin1("","",nucleus.getMassA(),
// 			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
// 			+nucleus.getExcitation()[shell],
// 			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
//     low=-1.;
//     if(kinmin1.GetPklab()>max_initial_nucl_mom) getMin(temptemp,low,nucleus,lepton,Q2,omega,shell,max_initial_nucl_mom);
// //     cout << "pm " << kinmin1.GetPklab() << " " << kin1.GetPklab() << " " << costhmin << " " << low << " " << high << endl;
//   }
//   return tempmin;//minimal initial nucleon momentum
}

//recursive function to find min costheta value so that initial nucleon momentum is low enough
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom){
  
  TKinematics2to2 kin("","",nucleus.getMassA(),
		      (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
		      +nucleus.getExcitation()[shell],
		      shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,(high+low)/2.);
  double pm=kin.GetPklab();
//   cout << omega << " " << Q2/1.E06 << " " << omega/(2.*MASSP*omega) << " " << high << " " << low << " " << pm << endl;
  if(pm<max_initial_nucl_mom) high=(high+low)/2.;
  else low=(high+low)/2.;
  if((high-low)<1.E-03) { return high;}
  else return getMin(high,low,nucleus,lepton,Q2,omega,shell,max_initial_nucl_mom);  
}

//recursive function to find max costheta value so that initial nucleon momentum is low enough
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom){
  
  TKinematics2to2 kin("","",nucleus.getMassA(),
		      (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
		      +nucleus.getExcitation()[shell],
		      shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,(high+low)/2.);
  double pm=kin.GetPklab();
//   cout << omega << " " << Q2/1.E06 << " " << omega/(2.*MASSP*omega) << " " << high << " " << low << " " << pm << endl;
  if(pm<max_initial_nucl_mom) low=(high+low)/2.;
  else high=(high+low)/2.;
  if((high-low)<1.E-03) { return low;}
  else return getMax(high,low,nucleus,lepton,Q2,omega,shell,max_initial_nucl_mom);  
}
