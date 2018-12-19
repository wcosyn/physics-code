//run ./neutest [p_\mu_parallel [MeV]] [p_\mu_perp] [integrator for flux] [cut in initial nucl momentum] 
//[cut in final nucl momentum (Pauli blocking motivated)] [plane wave(1) or fsi(0)] [max evaluation in matrix element FT integral=20000 is fine value] [homedir]

//we integrate for miniboone cross section
//integration happens over incoming neutrino energy folded with flux
//and also proton center of mass scattering angle (easier to keep track of physical points and such)
//it takes some effort in finding reasonable integration limits, the rest is basically integrating the QE cross section of course
//pretty straightforward...

//maxEvalweakamp values of 20k make the fluxintegrand converge faster than a value of 2000 (more accurate, so smoother cross sections...).!!!

#include <iostream>
#include <cstdlib>
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
#include "fluxes.h"

  //integration struct
struct Ftor {

  static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
    Ftor &p = * (Ftor *) param;
    p.f(ret,x[0],x[1],*p.pObs,*p.pNucleus,*p.lepton,p.E_out,p.costhetamu,p.current,p.cthmax,
        p.max_initial_nucl_mom,p.min_final_nucl_mom,p.pw,p.maxEvalweakamp,p.exp);
  }
  WeakQECross *pObs;
  MeanFieldNucleusThick *pNucleus;
  TLeptonKinematics *lepton;
  double E_out;
  double costhetamu;
  int current;
  double *cthmax;
  double max_initial_nucl_mom;
  double min_final_nucl_mom;
  int pw;
  int maxEvalweakamp;
  string exp;
  void (*f)(numint::vector_d &, double Ein, double costhetacm,
	    WeakQECross& pObs, MeanFieldNucleusThick &pNucleus, TLeptonKinematics &lepton,
	double E_out, double costhetamu, int current, double *cthmax, double max_initial_nucl_mom, double min_final_nucl_mom,
        int pw, int maxEvalweakamp, string exp);

};


//determine boundaries of kinematics
double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell, double max_initial_nucl_mom, double min_final_nucl_mom);
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);

// integration over missing momentum
void adap_intPm(numint::vector_d &, double E_in, double costhetacm,
		WeakQECross& pObs, MeanFieldNucleusThick &pNucleus, TLeptonKinematics &lepton,
		double E_out, double costhetamu, int current, double *cthmax, double max_initial_nucl_mom, 
                double min_final_nucl_mom, int pw, int maxEvalweakamp, string exp);


int main(int argc, char *argv[])
{
//   double Q2=atof(argv[1])*1.E06;
  double p_mu_z=atof(argv[1]); // muon momentum along beam [MeV]
  double p_mu_perp=atof(argv[2]); // muon momentum perp to beam [MeV]

  bool screening=0;//atoi(argv[4]);
  int nucleus=1;//atoi(argv[6]);                     
  double prec=1.E-04;//atof(argv[7]);
  int integrator=2;//
  int fluxintegrator=atoi(argv[3]);
  int thick=0;//atoi(argv[9]);
  int maxEvalweakamp=atoi(argv[7]);
  bool charged=1;
  int current=2;
  double max_initial_nucl_mom=atof(argv[4]);
  double min_final_nucl_mom=atof(argv[5]);
  int pw=atoi(argv[6]); //1 is plane-wave, 0 is with FSI
  
  bool enable_ROMEA=atoi(argv[8]);
  if(pw==1) enable_ROMEA=0; //no FSI if calculating the pw cross section...
  
//   double omega=Q2/(2.*MASSP*Bjx);
  
  string exp=argv[9]; //"miniboone" "minerva" "t2k"
  string homedir=argv[10];//"/home/wim/Code/share";


  MeanFieldNucleusThick Nucleus(nucleus,homedir);

  double p_mu = sqrt(p_mu_z*p_mu_z+p_mu_perp*p_mu_perp);
  double costhetamu = p_mu_z/p_mu;

  TLeptonKinematics *lepton = TLeptonKinematics::CreateWithCosScatterAngle(TLeptonKinematics::muon,costhetamu); 
  if((p_mu_perp/p_mu_z) > tan(20*DEGRTORAD)){
    cout << p_mu_perp/p_mu_z << " " << tan(20*DEGRTORAD) << endl;
    cout << p_mu_z << " " << p_mu_perp << " " << 0. << " " << 0. << " " << 0 << endl;
    delete lepton;
    return 0;    
  }
  double E_out=sqrt(pow(lepton->GetLeptonMass(),2.)+pow(p_mu,2.));
  double T_mu=E_out-lepton->GetLeptonMass();

  WeakQECross obs(lepton,&Nucleus,prec,integrator,homedir,charged, 1.03E03, screening,enable_ROMEA, 0.);
  
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
    minbeam=0.E03; 
    maxbeam=20.E03;
    neutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_nu_muon_FHC_flux,500,0,40);
    aneutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_anu_muon_RHC_flux,500,0,40);
  }  
  else if(!exp.compare("t2k")) { maxbeam=2.E02; }  //still to implement!
  else {cerr << "invalid experiment name chosen" << endl << "Choose either miniboone, minerva or t2k" << endl; assert(1==0);}

  cout << "norm neut " << neutnorm << endl << "norm aneut " << aneutnorm << endl;


  //find reasonable integration limits
  double E_low=E_out;
  double E_high=0.;
  for(int i=1;i<1000;i++){
    double E_in=E_out+(maxbeam-E_out)*0.001*i; //possible incoming lepton energies, start from lowest value up to max incoming
    double omega=E_in-E_out;
    double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton->GetLeptonMass()*lepton->GetLeptonMass()/(E_out*E_out))*costhetamu)
	-lepton->GetLeptonMass()*lepton->GetLeptonMass();
    double x=Q2/(2.*MASSP*omega);
//     cout << E_in << " " << omega << " " << Q2*1.E-06 << " " << x <<  endl;
    for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
      //lowest p_m is always at theta_cm -1
//       TKinematics2to2 kin("","",Nucleus.getMassA(),
// 			  (shell<Nucleus.getPLevels()? Nucleus.getMassA_min_proton(): Nucleus.getMassA_min_neutron())
// 			  +Nucleus.getExcitation()[shell],
// 			  shell<Nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
	double tempmax=-1., tempmin=1.;
	
      //anything above 500 MeV contribution will be negligible
      if(getBound(tempmax,tempmin,Nucleus,*lepton,E_in,E_out,costhetamu,shell,max_initial_nucl_mom,min_final_nucl_mom)<max_initial_nucl_mom){ 
        if(E_in<E_low) E_low=E_in;
	if(E_in>E_high) E_high=E_in; //update higher neutrino energy integration limit
       //update integration limit for costheta com.
	if(max<tempmax) max=tempmax;
	if(min>tempmin) min=tempmin;
	if(cthmax[shell]<tempmax) cthmax[shell]=tempmax;
	if(cthmin[shell]>tempmin) cthmin[shell]=tempmin;
	
      }
    }
  }
  if(E_low<minbeam) E_low=minbeam;
  if(E_high>maxbeam) E_high=maxbeam;
  cout << E_low << " " << E_high << " " << min << " " << max << endl;
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cout << shell << " " << cthmax[shell] << " " << cthmin[shell] << endl;}
  cout << endl << endl << endl;

  
  //test for comparison with pascal
  //../../bin/neutest 300 -0.9 1 ~/Code/share

//   double E_in=900.;
//   double omega=E_in-E_out;
//   double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton->GetLeptonMass()*lepton->GetLeptonMass()/(E_out*E_out))*costhetamu)
//       -lepton->GetLeptonMass()*lepton->GetLeptonMass();
//    cout << E_in << " " << omega<< " " << Q2*1.E-06 << " " << sqrt(Q2+omega*omega) << endl;
//   obs.getPlepton()->SetBeamEnergy(E_in);
//   double costhetacm=-1.; //for easy comparison, lots of currents are zero here!
//   double costhetacm=-0.9; //for easy comparison, lots of currents are zero here!
//   int shell=3;
//   TKinematics2to2 kin("","",Nucleus.getMassA(),
// 		      Nucleus.getMassA_min_proton()+Nucleus.getExcitation()[shell],
// 		      shell<Nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,costhetacm);
//   double pm=kin.GetPklab();
// //     cout << shell << " " << E_in <<  " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
// //     << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
// //     << " " << kin.GetKlab() << " " << kin.GetWlab() << endl;
//   if(!kin.IsPhysical()){
//     //for(int i=0;i<2;i++) results[i]+=0.;
//     cout << "bla " << E_in << " " << costhetacm << " " << shell << " " << pm << endl;
//   }
//   else{
//     double result=obs.getDiffWeakQECross(kin,current,thick,0,0,1,shell,0.,2E04,1,1);
//     cout << shell << " " << E_in <<  " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
//     << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
//     << " " << kin.GetKlab() << " " << kin.GetWlab() << " " << result << endl;
//   }
//   exit(1);
  
//   for(int j=0;j<20;j++){
//     double E_in=E_low+(E_high-E_low)*j/20.;
//       for(int i=0;i<20;i++){
// 	double costhetacm=-1.+2.*i/20.;
//       adap_intPm(avgcross,E_in,costhetacm,obs,Nucleus,*lepton,E_out,costhetamu,current,cthmax);
//     }
//   }
//   exit(1);
  
  
  //initialize object
  Ftor F;
  F.pObs = &obs;
  F.pNucleus = &Nucleus;
  F.lepton=lepton;
  F.E_out=E_out;
  F.costhetamu=costhetamu;
  F.current=current;
  F.cthmax=cthmax;
  F.min_final_nucl_mom=min_final_nucl_mom;
  F.max_initial_nucl_mom=max_initial_nucl_mom;
  F.pw=pw;
  F.maxEvalweakamp=maxEvalweakamp;
  F.exp=exp;
  
  numint::mdfunction<numint::vector_d,2> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

  numint::array<double,2> lower = {{E_low,min}};
  numint::array<double,2> upper = {{E_high,max}};
  
  
  F.f=adap_intPm;
  unsigned count=0;
  if(!fluxintegrator) numint::cube_romb(mdf,lower,upper,1.E-25,1.E-02,avgcross,count,0); //1.E-20,1.E-03
  else numint::cube_adaptive(mdf,lower,upper,1.E-25,1.E-02,2E02,2E04,avgcross,count,0);
   
  //cross section in 10^-39 cm^2 GeV ^-1 per nucleon!!
  //factor 2\pi because of integration over muon polar angle
//cout << endl << endl << endl;
  cout << p_mu_z << " " << p_mu_perp << " " << avgcross[0]*1.E19*2.*PI/Nucleus.getN()*p_mu_perp/E_out/p_mu/neutnorm
  << " " << avgcross[1]*1.E19*2.*PI/Nucleus.getZ()*p_mu_perp/E_out/p_mu/aneutnorm << " " << count 
  << " " << T_mu << " " << costhetamu << " " <<  avgcross[0]*1.E16*2.*PI/Nucleus.getN()/neutnorm
  << " " << avgcross[1]*1.E16*2.*PI/Nucleus.getZ()/aneutnorm << " "<< p_mu_perp/E_out/p_mu*1000. << endl;
  
  delete lepton;

}

//integrandum
void adap_intPm(numint::vector_d & results, double E_in, double costhetacm,
		WeakQECross& pObs, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
		double E_out, double costhetamu, int current, double *cthmax, double max_initial_nucl_mom, 
                double min_final_nucl_mom, int pw, int maxEvalweakamp, string exp){
		  
  results=numint::vector_d(2,0.);
  //we can fix the vector boson kinematics
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();

  pObs.getPlepton()->SetBeamEnergy(E_in);
//   int shell=0;
  for(int shell=0;shell<nucleus.getTotalLevels();shell++) {
    if(costhetacm<cthmax[shell]){
      TKinematics2to2 kin("","",nucleus.getMassA(),
			  (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			  +nucleus.getExcitation()[shell],
			  shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,costhetacm);
        double pm=kin.GetPklab();
      if(!kin.IsPhysical()||kin.GetPYlab()<min_final_nucl_mom||kin.GetPklab()>max_initial_nucl_mom){ //final nucleon momentum too low, impose cut!
	//for(int i=0;i<2;i++) results[i]+=0.;
// 	cout << "bla " << E_in << " " << costhetacm << " " << shell << " " << pm << " " << kin.GetPYlab() << " " << kin.IsPhysical() << endl;
      }
      else{
	double result=pObs.getDiffWeakQECross(kin,current,0,0,0,pw,shell,0.,maxEvalweakamp,WeakQECross::ctheta_cm,1);
        results[(shell<nucleus.getPLevels()?1:0)]+= result; //results[0] neutrino, results[1] antineutrino

        cout << shell << " " << E_in <<  " " << costhetacm << " "  << kin.GetCosthklab() << " " 
	<< kin.GetCosthYlab() << " " << kin.GetPklab() << " " << kin.GetPYlab() 
	<< " " << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetQsquared() << " "
	<< kin.GetXb()*nucleus.getMassA()/MASSP << " " << result << endl;
      }
    }
    
  }

  if(!exp.compare("miniboone")){
    results[0]*=interpolate(neutrino_flux::MiniBooNE_neut_muon_flux_norm,E_in,25,120,1);
    results[1]*=interpolate(neutrino_flux::MiniBooNE_antineut_muon_flux_norm,E_in,25,120,1);
  }
  else if(!exp.compare("minerva")){
    results[0]*=interpolate(neutrino_flux::Minerva_nu_muon_FHC_flux,E_in,500,20,0);
    results[1]*=interpolate(neutrino_flux::Minerva_anu_muon_RHC_flux,E_in,500,20,0);
  }
  else if(!exp.compare("t2k")){
    results[0]*=interpolate(neutrino_flux::t2k_neut_muon_flux_norm,E_in,25,80,1);
    results[1]*=interpolate(neutrino_flux::t2k_aneut_muon_flux_norm,E_in,25,80,1);
  }
  else assert(1==0);
  
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

