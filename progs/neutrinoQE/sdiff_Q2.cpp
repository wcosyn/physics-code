// Program that computes single differential cross sections for QE neutrino nucleus scattering
// Differential to Q2, computed from lepton kinematics
// Can be used for different experiments, Miniboone and minerva kinematics currently implemented, T2K will be added as well

//run like
// > sdiff_Q2 [Q2 in GeV^2] [integration algorithm for flux integration (see around line 400 for options)] 
//  [max initial nucl momentun in MeV] [min final nucl momentum in MeV] [plane wave 1 or not 0] [ROMEA 1 or not 0 in fsi]
// [experiment "miniboone", "minerva", "t2k"] [lepton id "electron", "muon", "tau"] [sharedir]

// example for plane wave calculation for miniboone
// > sdiff_Q2 0.1 3 500. 200. 1 0 miniboone ~/physics-code/share





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
    p.f(ret,x[0],x[1],x[2],*p.pNucleus,p.lepton_id,p.current,p.cthmax,p.Q2,
        p.prec,p.integrator,p.homedir,p.maxEval,p.charged,p.screening,p.enable_romea,
        p.max_initial_nucl_mom,p.min_final_nucl_mom,p.pw, p.exp
       );
  }
  MeanFieldNucleusThick *pNucleus;
  string lepton_id;
  int current;
  double *cthmax;
  double Q2;
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
  void (*f)(numint::vector_d &, double E_omega, double costhetacm, double E_in, 
            MeanFieldNucleusThick &pNucleus, string lepton_id, int current, double *cthmax, double Q2,
            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening, bool enable_romea,
            double max_initial_nucl_mom, double min_final_nucl_mom, int pw, string exp
           );

};

struct FtorH {  //Hydrogen

  static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
    FtorH &p = * (FtorH *) param;
    p.f(ret,x[0],p.Q2,p.charged,p.maxbeam,p.exp,p.lepton_id);
  }
  double Q2;
  bool charged;
  double maxbeam;
  string exp;
  string lepton_id;
  void (*f)(numint::vector_d &, double E_out, double Q2, bool charged, double maxbeam, string exp, string lepton_id);

}; 


struct FtorRFG {  //Relativistic Fermi Gas

  static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
    FtorRFG &p = * (FtorRFG *) param;
    p.f(ret,x[0],x[1],p.lepton_id,p.Q2,p.charged,p.Pauli,p.maxbeam,p.exp);
  }
  string lepton_id;
  double Q2;
  bool charged;
  bool Pauli;
  double maxbeam;
  string exp;
  void (*f)(numint::vector_d &, double omega, double E_out, string lepton_id, double Q2, bool charged, bool Pauli, double maxbeam, string exp);

}; 


//integration for Hydrogen
void int_hydr(numint::vector_d &, double E_out,double Q2,bool charged, double maxbeam, string exp, string lepton_id);

void int_RFG(numint::vector_d &, double omega,double E_out, string lepton_id, double Q2,bool charged, bool Pauli, double maxbeam, string exp);


//determine boundaries of kinematics
double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell, double max_initial_nucl_mom, double min_final_nucl_mom);
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);

// integration over missing momentum + over T_mu(=E_out+C)
void adap_intPm(numint::vector_d &, double omega, double costhetacm, double E_in,
		            MeanFieldNucleusThick &pNucleus, string lepton_id, int current, double *cthmax, double Q2,
		            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening, bool enable_romea,
                            double max_initial_nucl_mom, double min_final_nucl_mom, int pw, string exp
               );


int main(int argc, char *argv[])
{
  double Q2=atof(argv[1])*1.E06;   // give input in GeV^2
//double T_mu=atof(argv[1]); // muon kinetic energy in MeV
//   double costhetamu; //=atof(argv[1]);
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
  
  
  vector<double> avgcross(2,0.); //neutrino and antineutrino
  //estimates for integration bounds
  double cthmax[Nucleus.getTotalLevels()];
  double cthmin[Nucleus.getTotalLevels()];
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cthmax[shell]=-1.;cthmin[shell]=1.;}
  double max=-1.,min=1.; //overall min and max center of mass cosine theta angles
  
  double minbeam=0., maxbeam=0.; //different beam energy intervals for the different experiments
  double neutnorm=1.,aneutnorm=1.; // flux normalisations

  if(!exp.compare("miniboone")) { maxbeam=3.E03;}
  else if(!exp.compare("minerva")) { 
    minbeam=1.5E03; 
    maxbeam=10.E03;
    if(!lepton_id.compare("muon")){
      neutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_nu_muon_FHC_flux,500,round(1500./500.),round(10.E03/500));
      aneutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_anu_muon_RHC_flux,500,round(1500./500.),round(10.E03/500));
//      for(int i=0;i<20;i++) neutrino_flux::Minerva_nu_muon_FHC_flux[i]/=neutnorm;
//      for(int i=0;i<20;i++) neutrino_flux::Minerva_anu_muon_RHC_flux[i]/=aneutnorm;
    }
    else if(!lepton_id.compare("electron")){
      neutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_nu_elec_FHC_flux,500,round(1500./500.),round(10.E03/500));
      aneutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_anu_elec_FHC_flux,500,round(1500./500.),round(10.E03/500));      
    }
  }  
  else if(!exp.compare("t2k")) { maxbeam=2.E02; }  //still to implement!
  else {cerr << "invalid experiment name chosen" << endl << "Choose either miniboone, minerva or t2k" << endl; assert(1==0);}
  
  cout << "norm neut " << neutnorm << endl << "norm aneut " << aneutnorm << endl;

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


  //setting up variables dat will determine integration limits
  double E_out_max=minbeam;
  double E_out_min=maxbeam;
  double E_low=maxbeam;
  double E_high=minbeam;
  double omega_low=maxbeam;
  double omega_hi=minbeam;
  double E_muH_low=maxbeam;
  double E_muH_hi=minbeam;
  double E_RFG_out_max=minbeam;
  double E_RFG_out_min=maxbeam;
  double omega_RFG_low=maxbeam;
  double omega_RFG_hi=minbeam;
  
  
  
  //find reasonable integration limits by scanning kinematics and 
  //making sure nucleons have "reasonable" momenta, lepton kinematics are physical etc.
  for(int j=0;j<=100;j++){ //loop over T_mu
//   cout << j << "/100" << endl;
    double E_out=(maxbeam-leptonmass-1.E-03)/100.*j+leptonmass+1.E-03;
    double p_out=sqrt(E_out*E_out-leptonmass*leptonmass);
    //physical limits so that costhetamu is between -1 and 1
    double E_in_low = (Q2+leptonmass*leptonmass)/2./(E_out+p_out);
    if(E_in_low<E_out) E_in_low=E_out;
    double E_in_hi = (Q2+leptonmass*leptonmass)/2./(E_out-p_out);
    
    //hydrogen piece
    double E_nu_H_n = (Q2+MASSP*MASSP-MASSN*MASSN)/2./MASSN+E_out;
    double E_nu_H_p = (Q2+MASSN*MASSN-MASSP*MASSP)/2./MASSP+E_out;
    if (E_nu_H_p>E_out&&E_nu_H_p<maxbeam&&E_nu_H_p<E_in_hi&&E_nu_H_p>E_in_low){
      if(E_out>E_muH_hi) E_muH_hi = E_out;
      if(E_out<E_muH_low) E_muH_low = E_out;
    }
    if (E_nu_H_n>E_out&&E_nu_H_n<maxbeam&&E_nu_H_n<E_in_hi&&E_nu_H_n>E_in_low){
      if(E_out>E_muH_hi) E_muH_hi = E_out;
      if(E_out<E_muH_low) E_muH_low = E_out;
    }
//     cout << j << " " << E_out << " " << E_nu_H_n << " " << E_nu_H_p << " " << E_in_hi << " " << E_in_low << " " << E_muH_low << " " << E_muH_hi << endl;
    if(E_in_low>maxbeam||E_in_low>E_in_hi){} //outside considered incoming neutrino energies, so don't bother
    else{
      for(int i=0;i<=100;i++){
        if(E_in_hi>maxbeam) E_in_hi=maxbeam;
        double E_in=E_in_low+(E_in_hi-E_in_low)*1.E-02*i;
        double costhetamu=(-Q2-leptonmass*leptonmass+2.*E_in*E_out)/(2.*E_in*p_out);
        double omega=E_in-E_out;
        TLeptonKinematics *lepton=TLeptonKinematics::CreateWithBeamEnergy(static_cast<TLeptonKinematics::Lepton>(flav),E_in);
        //RFG bit
        double qvec = sqrt(Q2+omega*omega);
        double lambda=omega/2./sqrt(MASSP*MASSN);
        double kappa=qvec/2./sqrt(MASSP*MASSN);
        double tau=Q2/4./MASSP/MASSN;
        double k_Fermi=228.;
        
        double xi=sqrt(1.+k_Fermi*k_Fermi/MASSP/MASSN)-1.;
        double psi=pow(lambda-tau,2.)/(xi*((1.+lambda)*tau+kappa*sqrt(tau*(1.+tau))));
      //   cout << "blaaa " << psi << endl;
        if(abs(psi)<1.){
          if(E_out>E_RFG_out_max) E_RFG_out_max = E_out;
          if(E_out<E_RFG_out_min) E_RFG_out_min = E_out;
          if(omega>omega_RFG_hi) omega_RFG_hi = omega;        
          if(omega<omega_RFG_low) omega_RFG_low = omega;
        }

        //Carbon checks
        for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
          double tempmax=-1., tempmin=1.;
          
          //anything above 500 MeV contribution will be negligible
          if(getBound(tempmax,tempmin,Nucleus,*lepton,E_in,E_out,costhetamu,shell,max_initial_nucl_mom,min_final_nucl_mom)<max_initial_nucl_mom){ 
            if(E_in<E_low)  E_low=E_in;
            if(E_in>E_high) E_high=E_in;
            if(max<tempmax) max=tempmax;
            if(min>tempmin) min=tempmin;
            if(cthmax[shell]<tempmax) cthmax[shell]=tempmax;
            if(cthmin[shell]>tempmin) cthmin[shell]=tempmin;
            if(E_out_max<E_out) E_out_max=E_out;
            if(E_out_min>E_out) E_out_min=E_out;
            if(omega>omega_hi) omega_hi=omega;
            if(omega<omega_low) omega_low=omega;
//             cout << E_in << " " << E_out << " " << tempmin << " " << tempmax << " " << pm_min << endl;
          }
        }
        delete lepton;
      }
    }
  }  
  if(E_high>maxbeam) E_high=maxbeam;
  if(E_low<minbeam) E_low=minbeam;
  if(E_muH_low<minbeam) E_muH_low=minbeam;
  if(E_low>E_high) E_low=E_high;
  cout << endl;
  cout << "min=" << min << "   max=" << max << endl;
  cout << "Eoutmin=" << E_out_min << "  Eoutmax=" << E_out_max << endl;
  cout << "omega_low=" << omega_low << " " << "  omega_high=" << omega_hi << endl;
  cout << "E_low=" << E_low << "  E_high=" << E_high << endl << endl;

  cout << "E_muH_low=" << E_muH_low << "  E_muH_hi=" << E_muH_hi << endl << endl;
  
  cout << "E_RFG_out_min=" << E_RFG_out_min << "  E_RFG_out_max=" << E_RFG_out_max << endl;
  cout << "omega_RFG_low=" << omega_RFG_low <<"  omega_RFG_hi=" << omega_RFG_hi << endl << endl;
  
  FtorH FH;
  FH.Q2=Q2;
  FH.charged=charged;
  FH.exp=exp;
  FH.maxbeam=maxbeam;
  FH.lepton_id=lepton_id;

  numint::mdfunction<numint::vector_d,1> mdfH;
  mdfH.func = &FtorH::exec;
  mdfH.param = &FH;

  numint::array<double,1> lowerH = {{E_muH_low}};  
  numint::array<double,1> upperH = {{E_muH_hi}};
  
  FH.f=int_hydr;
  vector<double> avgcrossH(2,0.); 
  unsigned count=0;
  if(!fluxintegrator) numint::cube_romb(mdfH,lowerH,upperH,1.E-30,1.E-03,avgcrossH,count,0); //1.E-20,1.E-03
  else numint::cube_adaptive(mdfH,lowerH,upperH,1.E-30,1.E-03,2E02,4E04,avgcrossH,count,0); 
  
//   cout << "Crosssection H: " << Q2 << " " << avgcrossH[0]*1.E19*2.*PI << " [1E-39 cm2/GeV2]" <<  " "<< count << endl; //2pi because of integration over scattered lepton angle
//   cout << Q2QE << " " << avgcrossH[0]*1.E19*2.*PI << " " << avgcrossH[1]*1.E19*2.*PI  << endl;

 
//RFG calculation
  FtorRFG FRFG;
  
  FRFG.lepton_id=lepton_id;
  FRFG.Q2=Q2;
  FRFG.charged=charged;
  FRFG.Pauli = Pauli;
  FRFG.maxbeam=maxbeam;
  FRFG.exp=exp;
  
  numint::mdfunction<numint::vector_d,2> mdfRFG;
  mdfRFG.func = &FtorRFG::exec;
  mdfRFG.param = &FRFG;

  numint::array<double,2> lowerRFG = {{omega_RFG_low,E_RFG_out_min}};
  numint::array<double,2> upperRFG = {{omega_RFG_hi,E_RFG_out_max}};
  
  FRFG.f=int_RFG;
  vector<double> avgcrossRFG(2,0.); 
  count=0;
  if(!fluxintegrator) numint::cube_romb(mdfRFG,lowerRFG,upperRFG,1.E-30,1.E-03,avgcrossRFG,count,0); //1.E-20,1.E-03
  else numint::cube_adaptive(mdfRFG,lowerRFG,upperRFG,1.E-30,1.E-03,2E02,4E04,avgcrossRFG,count,0); 

//   cout << "Crosssection H: " << Q2 << " " << avgcrossH[0]*1.E19*2.*PI << " [1E-39 cm2/GeV2]" <<  " "<< count << endl; //2pi because of integration over scattered lepton angle
// neutrino (on neutron) first, antineutrino (on proton) second
  cout << Q2*1.E-06 << " " << avgcrossH[0]*1.E19*2.*PI/neutnorm << " " << avgcrossH[1]*1.E19*2.*PI/aneutnorm
            << " " << avgcrossRFG[0]*1.E19*2.*PI/neutnorm << " " << avgcrossRFG[1]*1.E19*2.*PI/aneutnorm  << endl;

  //initialize object -- Carbon
  Ftor F;
  F.pNucleus = &Nucleus;
  F.lepton_id=lepton_id;
  F.current=current;
  F.cthmax=cthmax;
  F.Q2=Q2;
  F.prec=prec;
  F.integrator=integrator;
  F.homedir=homedir;
  F.maxEval=maxEval;
  F.charged=charged;
  F.screening=screening;
  F.enable_romea = enable_romea;
  F.min_final_nucl_mom=min_final_nucl_mom;
  F.max_initial_nucl_mom=max_initial_nucl_mom;
  F.pw=pw;
  F.exp=exp;

  numint::mdfunction<numint::vector_d,3> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

//   numint::array<double,3> lower = {{omega_low,min,Tmin+leptonmass}};
//   numint::array<double,3> upper = {{omega_hi,max,Tmax+leptonmass}};
  numint::array<double,3> lower = {{omega_low,min,E_low}};
  numint::array<double,3> upper = {{omega_hi,max,E_high}};
  
  F.f=adap_intPm;
  count=0;
  string stf=homedir+"/statefile/sdiff_q2";
  ostringstream qstr;
  qstr << Q2*1.E-06 << "." << pw << "." <<enable_romea << "." <<max_initial_nucl_mom << "." <<min_final_nucl_mom;
  string qst=qstr.str();
  stf+=qst;
  int nregions,fail,countt;
  vector<double> err(2,0.);
  vector<double> prob(2,0.);
  int nvec=1;
  double epsrel=1.E-01;
  double epsabs=1.E-25;
  int flags=2;
  int seed=235;
  int minEval=100;
  int maxEvalcuba=20000;
  double xgiven[3]={200.,-1.,500.};
  cout << "start integration C" << endl;
  if(fluxintegrator==0) numint::vegas( mdf, lower,upper,nvec, epsrel,epsabs,flags,seed,minEval, maxEvalcuba,1000, 500, 1000, 0,(stf+"vegas").c_str(),countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==1) numint::cuhre( mdf, lower,upper,nvec, epsrel,epsabs,flags,minEval, maxEvalcuba,11, (stf+"cuhre").c_str(),nregions,countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==2) numint::cube_romb(mdf,lower,upper,1.E-25,1.E-01,avgcross,count,0); //1.E-20,1.E-03
  if(fluxintegrator==3) numint::cube_adaptive(mdf,lower,upper,1.E-25,1.E-03,2E02,2E04,avgcross,count,0);
  if(fluxintegrator==4) numint::divonne(mdf, lower,upper,nvec, epsrel, epsabs,flags,seed,
           minEval, maxEvalcuba,7,7,0,5,0,10,0.25,0,3,NULL, 0, 0,
           (stf+"divonne").c_str(),nregions, countt, fail, avgcross,err,prob );
  
  //factor 2\pi because of integration over muon polar angle
  cout << Q2*1.E-6 << " " << avgcrossH[0]*1.E19*2.*PI/neutnorm << " " << avgcrossH[1]*1.E19*2.*PI/aneutnorm << " " 
      << avgcrossRFG[0]*1.E19*2.*PI/neutnorm << " " << avgcrossRFG[1]*1.E19*2.*PI/aneutnorm << " " <<  
    avgcross[0]*1.E19*2.*PI/Nucleus.getN()/neutnorm << " " << avgcross[1]*1.E19*2.*PI/Nucleus.getZ()/aneutnorm << " " <<
    (avgcross[0]+avgcross[1])*1.E19*2.*PI/Nucleus.getN()/(neutnorm+aneutnorm) << " " <<  count << endl;

}

//integrandum carbon
void adap_intPm(numint::vector_d & results, double omega, double costhetacm, double E_in,
	 MeanFieldNucleusThick &nucleus, string lepton_id, int current, double *cthmax, double Q2,
	 double prec, int integrator, string homedir, int maxEval, bool charged, bool screening,
         bool enable_romea, double max_initial_nucl_mom, 
                double min_final_nucl_mom, int pw, string exp){	  

  results=numint::vector_d(2,0.);
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
  double E_out = E_in-omega;
  double p_out=sqrt(E_out*E_out-leptonmass*leptonmass);
  double costhetamu=(-Q2-leptonmass*leptonmass+2.*E_in*E_out)/(2.*E_in*p_out);
  if(/*E_in>maxbeam||E_in<minbeam*/E_out<0.) { return;}
  if(abs(costhetamu)<=1. /*&& omega>0*/) {
      
    TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(static_cast<TLeptonKinematics::Lepton>(flav),E_in);
    WeakQECross pobs(lepton,&nucleus,prec,integrator,homedir,charged,1.03E03,screening,enable_romea);  
    

    for(int shell=0;shell<nucleus.getTotalLevels();shell++) {
      if(costhetacm<cthmax[shell]){
        TKinematics2to2 kin("","",nucleus.getMassA(),
          (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
          +nucleus.getExcitation()[shell],
          shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,costhetacm);
        if(!kin.IsPhysical()||kin.GetPYlab()<min_final_nucl_mom||kin.GetPklab()>max_initial_nucl_mom){ //final nucleon momentum too low, impose cut!
          for(int i=0;i<2;i++) results[i]+=0.; /*cout << "cut " << kin.GetPYlab() << " " << kin.GetPklab() << endl;*/
        }
        else{
          double result=pobs.getDiffWeakQECross(kin,current,0,0,0,pw,shell,0.,maxEval,WeakQECross::ctheta_cm,1);   // prec..2E04
  //	  			cout << "Result " << result << endl;

          //Jacobian    
          double jcb=2*E_in*p_out;
          result/=jcb;      
    
          results[(shell<nucleus.getPLevels()?1:0)]+= result; //results[0] neutrino, results[1] antineutrino
           cout << shell << " " << E_in <<  " " << costhetacm << " " << p_out << " "  << kin.GetCosthklab()<< " " 
          	<< kin.GetCosthYlab() << " " << kin.GetPklab() << " " << kin.GetPYlab() 
          	<< " " << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetQsquared() << " "
          	<< kin.GetXb()*nucleus.getMassA()/MASSP << " " << result <<  endl;
        }
      }  
    }
    //fold with flux
    if(!exp.compare("miniboone")&&!lepton_id.compare("muon")){
      results[0]*=interpolate(neutrino_flux::MiniBooNE_neut_muon_flux_norm,E_in,25,120,1);
      results[1]*=interpolate(neutrino_flux::MiniBooNE_antineut_muon_flux_norm,E_in,25,120,1);
    }
    else if(!exp.compare("minerva")&&!lepton_id.compare("muon")){
      results[0]*=interpolate(neutrino_flux::Minerva_nu_muon_FHC_flux,E_in,500,40,0);
      results[1]*=interpolate(neutrino_flux::Minerva_anu_muon_RHC_flux,E_in,500,40,0);
    }
    else if(!exp.compare("minerva")&&!lepton_id.compare("electron")){
      results[0]*=interpolate(neutrino_flux::Minerva_nu_elec_FHC_flux,E_in,500,20,0);
      results[1]*=interpolate(neutrino_flux::Minerva_anu_elec_FHC_flux,E_in,500,20,0);
    }
    else if(!exp.compare("t2k")&&!lepton_id.compare("muon")){
      results[0]*=interpolate(neutrino_flux::t2k_neut_muon_flux_norm,E_in,25,80,1);
      results[1]*=interpolate(neutrino_flux::t2k_aneut_muon_flux_norm,E_in,25,80,1);
    }
    else { cerr << "Unsupported combination of experiment and lepton_id:" << exp << " " << lepton_id << endl; assert(1==0);}
  
    delete lepton;
  }
  else {results[0]=0; results[1]=0; /*cout << "unphys " << omega << " " << costhetacm << " " << E_out << " " << E_in << " " << 0. << " " << 0. << endl;*/}

}

//integrandum hydrogen
void int_hydr(numint::vector_d & results, double E_out, 
              double Q2, bool charged, double maxbeam, string exp, string lepton_id){
    
  results=numint::vector_d(2,0.);  
  double leptonmass=0.;
  if(!lepton_id.compare("electron")){
    leptonmass=TLeptonKinematics::masse;
  } 
  else if (!lepton_id.compare("muon")){
    leptonmass=TLeptonKinematics::massmu;
  } 
  else if (!lepton_id.compare("tau")){
    leptonmass=TLeptonKinematics::masstau;
  } 
 //   if(abs(costhetamu)>1.) {/*cout << "costheta" << endl;*/ cout << E_out << " " << result[0] << " " << result[1] << endl; return;}
  double E_in_n=(Q2+MASSP*MASSP-MASSN*MASSN)/2./MASSN+E_out;
  double E_in_p=(Q2+MASSN*MASSN-MASSP*MASSP)/2./MASSP+E_out;
  
  //[0] is neutron [1] is proton
  results[0]=WeakQECross::getElWeakQECross(Q2,E_in_n,0,charged,1.03E03,leptonmass*leptonmass,WeakQECross::realQ2);
  results[1]=WeakQECross::getElWeakQECross(Q2,E_in_p,1,charged,1.03E03,leptonmass*leptonmass,WeakQECross::realQ2);
  
  if(!exp.compare("miniboone")&&!lepton_id.compare("muon")){
    results[0]*=E_in_n>maxbeam? 0. : interpolate(neutrino_flux::MiniBooNE_neut_muon_flux_norm,E_in_n,25,120,1);
    results[1]*=E_in_p>maxbeam? 0. : interpolate(neutrino_flux::MiniBooNE_antineut_muon_flux_norm,E_in_p,25,120,1);
  }
  else if(!exp.compare("minerva")&&!lepton_id.compare("muon")){
    results[0]*=E_in_n>maxbeam? 0. : interpolate(neutrino_flux::Minerva_nu_muon_FHC_flux,E_in_n,500,20,0);
    results[1]*=E_in_p>maxbeam? 0. : interpolate(neutrino_flux::Minerva_anu_muon_RHC_flux,E_in_p,500,20,0);
  }
   else if(!exp.compare("minerva")&&!lepton_id.compare("electron")){
      results[0]*=E_in_n>maxbeam? 0. : interpolate(neutrino_flux::Minerva_nu_elec_FHC_flux,E_in_n,500,20,0);
      results[1]*=E_in_p>maxbeam? 0. : interpolate(neutrino_flux::Minerva_anu_elec_FHC_flux,E_in_p,500,20,0);
    }
  else if(!exp.compare("t2k")&&!lepton_id.compare("muon")){
    results[0]*=E_in_n>maxbeam? 0. : interpolate(neutrino_flux::t2k_neut_muon_flux_norm,E_in_n,25,80,1);
    results[1]*=E_in_p>maxbeam? 0. : interpolate(neutrino_flux::t2k_aneut_muon_flux_norm,E_in_p,25,80,1);
  }
  else {cerr << "Unsupported combination of experiment and lepton_id:" << exp << " " << lepton_id << endl; assert(1==0);}

//   if(isnan(result[1])) {cout << E_out << " " << E_in_n << " " << crossHanu_n << " " << interpolate(MiniBooNE_neut_flux_norm,E_in_n,25,120,1) << " " << Q2 << endl; exit(1);}
//    cout << E_out << " " << result[0] << " " << result[1] << " " << E_in_p << " " << E_in_n << " " <<  interpolate(MiniBooNE_antineut_flux_norm,E_in_p,25,120,1) << " " << interpolate(MiniBooNE_neut_flux_norm,E_in_n,25,120,1) <<endl;
} 


void int_RFG(numint::vector_d & results, double omega,double E_out, string lepton_id, double Q2,bool charged, bool Pauli, double maxbeam, string exp){
    
  results=numint::vector_d(2,0.);  
  double leptonmass=0.;
  if(!lepton_id.compare("electron")){
    leptonmass=TLeptonKinematics::masse;
  } 
  else if (!lepton_id.compare("muon")){
    leptonmass=TLeptonKinematics::massmu;
  } 
  else if (!lepton_id.compare("tau")){
    leptonmass=TLeptonKinematics::masstau;
  } 
  double E_in=omega+E_out;
  if(E_in>maxbeam)  {return;}
  double p_out = sqrt(E_out*E_out-leptonmass*leptonmass);
  double costhetamu = (-Q2-leptonmass*leptonmass+2.*E_in*E_out)/(2.*E_in*p_out);

  if(abs(costhetamu)>1.) {/*cout << "costheta" << endl;*/ return;}

  double kf=228.;
  WeakQECross::getRFGWeakQECross(results[1],results[0],Q2, E_in, omega,kf,leptonmass*leptonmass,1,1.03E03,WeakQECross::realQ2,Pauli);
   //cout << omega << " " << E_out << " " << results[1] << " " << results[0] << endl;
  
  if(!exp.compare("miniboone")&&!lepton_id.compare("muon")){
    results[0]*=interpolate(neutrino_flux::MiniBooNE_neut_muon_flux_norm,E_in,25,120,1);
    results[1]*=interpolate(neutrino_flux::MiniBooNE_antineut_muon_flux_norm,E_in,25,120,1);
  }
  else if(!exp.compare("minerva")&&!lepton_id.compare("muon")){
    results[0]*=interpolate(neutrino_flux::Minerva_nu_muon_FHC_flux,E_in,500,40,0);
    results[1]*=interpolate(neutrino_flux::Minerva_anu_muon_RHC_flux,E_in,500,40,0);
  }
   else if(!exp.compare("minerva")&&!lepton_id.compare("electron")){
      results[0]*=interpolate(neutrino_flux::Minerva_nu_elec_FHC_flux,E_in,500,20,0);
      results[1]*=interpolate(neutrino_flux::Minerva_anu_elec_FHC_flux,E_in,500,20,0);
    }
  else if(!exp.compare("t2k")&&!lepton_id.compare("muon")){
    results[0]*=interpolate(neutrino_flux::t2k_neut_muon_flux_norm,E_in,25,80,1);
    results[1]*=interpolate(neutrino_flux::t2k_aneut_muon_flux_norm,E_in,25,80,1);
  }
  else { cerr << "Unsupported combination of experiment and lepton_id:" << exp << " " << lepton_id << endl; assert(1==0);}
    
   //if(isnan(results[0])) { cout << omega << " " << E_out << " " << results[0] << " " << results[1] << " " 
    //<< interpolate(neutrino_flux::Minerva_nu_elec_FHC_flux,E_in,500,20,0) << " " << interpolate(neutrino_flux::Minerva_anu_elec_FHC_flux,E_in,500,20,0) << endl;}
} 





double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell, double max_initial_nucl_mom, double min_final_nucl_mom){
  
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();
  double tempmin=max_initial_nucl_mom;
  double costhmin=-1.;
  for(int i=0;i<=100;i++){ //loop over all com angles to find com scattering angle where initial nucl momentum is minimal
    TKinematics2to2 kin("","",nucleus.getMassA(),
			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			+nucleus.getExcitation()[shell],
			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.+i*0.02);
    double pm=kin.GetPklab(); //initial nucleon momentum
    if(pm<tempmin&&kin.GetPYlab()>min_final_nucl_mom) {tempmin=pm; costhmin=-1.+i*0.02;}
  }
  //   cout << omega << " " << Q2/1.E06 << " " << omega/(2.*MASSP*omega) << " " << high << " " << low << " " << pm << endl;
  double temptemp=costhmin;//com costheta value where the min initial nucl momentum is
//   cout << costhmin << " " << tempmin << endl;
  //if it is below the limit we consider, we find the interval in angles where it is below
  if(tempmin<max_initial_nucl_mom){
    TKinematics2to2 kin1("","",nucleus.getMassA(),
			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			+nucleus.getExcitation()[shell],
			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,1.);
    high=1;
    if(kin1.GetPklab()>max_initial_nucl_mom) high=getMax(high,temptemp,nucleus,lepton,Q2,omega,shell,max_initial_nucl_mom);
    temptemp=costhmin;
    TKinematics2to2 kinmin1("","",nucleus.getMassA(),
			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			+nucleus.getExcitation()[shell],
			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
    low=-1.;
    if(kinmin1.GetPklab()>max_initial_nucl_mom) getMin(temptemp,low,nucleus,lepton,Q2,omega,shell,max_initial_nucl_mom);
//     cout << "pm " << kinmin1.GetPklab() << " " << kin1.GetPklab() << " " << costhmin << " " << low << " " << high << endl;
  }
  return tempmin;//minimal initial nucleon momentum
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
