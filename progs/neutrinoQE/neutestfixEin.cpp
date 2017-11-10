//run ./neutestfixEin [E_in [MeV]] [cos(theta_mu)] [fluxintegrator] [max initial nucl mom MeV] [min final nucl mom MeV] [plane wave 1 or FSI 0] [Romea] [homedir]

//we calculate two-fold diff (Anti)neutrino cross sections, not flux folded, fixed neutrino energy
//we integrate over proton center of mass scattering angle (easier to keep track of physical points and such)
//it takes some effort in finding reasonable integration limits, the rest is basically integrating the QE cross section of course
//pretty straightforward...

//version 

#include <iostream>
#include <cstdlib>

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



  //integration struct
struct Ftor {

  static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
    Ftor &p = * (Ftor *) param;
    p.f(ret,x[0],*p.pObs,*p.pNucleus,*p.lepton,p.E_out,p.costhetamu,p.current,p.cthmax,p.E_in,
      p.max_initial_nucl_mom,p.min_final_nucl_mom,p.pw,p.maxEvalweakamp);
  }
  WeakQECross *pObs;
  MeanFieldNucleusThick *pNucleus;
  TLeptonKinematics *lepton;
  double E_out;
  double costhetamu;
  int current;
  double *cthmax;
  double E_in;
  double max_initial_nucl_mom;
  double min_final_nucl_mom;
  int pw;
  int maxEvalweakamp;
  void (*f)(numint::vector_d &, double costhetacm,
	    WeakQECross& pObs, MeanFieldNucleusThick &pNucleus, TLeptonKinematics &lepton,
	double E_out, double costhetamu, int current, double *cthmax, double E_in,
	 double max_initial_nucl_mom, double min_final_nucl_mom, int pw, int maxEvalweakamp);

};


//determine boundaries of kinematics
double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell, double max_initial_nucl_mom, double min_final_nucl_mom);
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);

void adap_intPm(numint::vector_d &, double costhetacm,
		WeakQECross& pObs, MeanFieldNucleusThick &pNucleus, TLeptonKinematics &lepton,
		double E_out, double costhetamu, int current, double *cthmax, double E_in,
		double max_initial_nucl_mom, double min_final_nucl_mom, int pw, int maxEvalweakamp);


int main(int argc, char *argv[])
{
//   double Q2=atof(argv[1])*1.E06;
  double E_in=atof(argv[1]);
  double costhetamu=atof(argv[2]);
  bool screening=0;//atoi(argv[4]);
  int nucleus=1;//atoi(argv[6]);
  double prec=1.E-05;//atof(argv[7]);
  int integrator=2;//
  int fluxintegrator=atoi(argv[3]);
  int thick=0;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  bool charged=1;
  int current=2;
  double max_initial_nucl_mom=atof(argv[4]);
  double min_final_nucl_mom=atof(argv[5]);
  int pw=atoi(argv[6]); //1 is plane-wave, 0 is with FSI
  
  bool enable_ROMEA=atoi(argv[7]);
  if(pw==1) enable_ROMEA=0; //no FSI if calculating the pw cross section...
  
//   double omega=Q2/(2.*MASSP*Bjx);
  
  string homedir=argv[8];//"/home/wim/Code/share";

  MeanFieldNucleusThick Nucleus(nucleus,homedir);
  TLeptonKinematics *lepton = TLeptonKinematics::CreateWithCosScatterAngle(TLeptonKinematics::muon,costhetamu);

  for(int i=1;i<100;i++){
    double omega=E_in/100.*i;

    double E_out=E_in-omega;;
    
    WeakQECross obs(lepton,&Nucleus,prec,integrator,homedir,charged, 1.03E03, screening,enable_ROMEA, 0.);
    
    vector<double> avgcross(2,0.); 
    //estimates for integration bounds
    double cthmax[Nucleus.getTotalLevels()];
    double cthmin[Nucleus.getTotalLevels()];
    for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cthmax[shell]=-1.;cthmin[shell]=1.;}
    double max=-1.,min=1.;

    //find reasonable integration limits
    //double omega=E_in-E_out;
    double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton->GetLeptonMass()*lepton->GetLeptonMass()/(E_out*E_out))*costhetamu)
	-lepton->GetLeptonMass()*lepton->GetLeptonMass();
    double x=Q2/(2.*MASSn*omega);
//     cout << E_in << " " << omega << " " << Q2*1.E-06 << " " << x <<  endl;
    for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
      //lowest p_m is always at theta_cm -1
      double tempmax=-1., tempmin=1.;
	
      //anything above 500 MeV contribution will be negligible
      if(getBound(tempmax,tempmin,Nucleus,*lepton,E_in,E_out,costhetamu,shell,max_initial_nucl_mom,min_final_nucl_mom)<max_initial_nucl_mom){ 
	if(max<tempmax) max=tempmax;
	if(min>tempmin) min=tempmin;
	if(cthmax[shell]<tempmax) cthmax[shell]=tempmax;
	if(cthmin[shell]>tempmin) cthmin[shell]=tempmin;
	
      }
    }
  
//     for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cout << shell << " " << cthmax[shell] << " " << cthmin[shell] << endl;}
//     cout << endl << endl << endl;
//     if(min>=max)   cout << E_in << " " << costhetamu << " " << 0.
//     << " " << 0. << " " << 0 << endl;

  
  
  //initialize object
    Ftor F;
    F.pObs = &obs;
    F.pNucleus = &Nucleus;
    F.lepton=lepton;
    F.E_out=E_out;
    F.costhetamu=costhetamu;
    F.current=current;
    F.cthmax=cthmax;
    F.E_in=E_in;
    F.min_final_nucl_mom=min_final_nucl_mom;
    F.max_initial_nucl_mom=max_initial_nucl_mom;
    F.pw=pw;
    F.maxEvalweakamp=maxEval;

    numint::mdfunction<numint::vector_d,1> mdf;
    mdf.func = &Ftor::exec;
    mdf.param = &F;

    unsigned neval = 0;
    numint::array<double,1> lower = {{min}};
    numint::array<double,1> upper = {{max}};
    
    
    F.f=adap_intPm;
    unsigned count=0;
    if(!fluxintegrator) numint::cube_romb(mdf,lower,upper,1.E-25,1.E-02,avgcross,count,0); //1.E-20,1.E-03
    else numint::cube_adaptive(mdf,lower,upper,1.E-25,1.E-02,2E02,2E04,avgcross,count,0);
    
    double kf=228.;
    double rfgp=0.,rfgn=0.,rfgp_bl=0.,rfgn_bl=0.;
    
    WeakQECross::getRFGWeakQECross(rfgp_bl,rfgn_bl,Q2,E_in,omega,kf,1,1.03E03,0,1);
    WeakQECross::getRFGWeakQECross(rfgp,rfgn,Q2,E_in,omega,kf,1,1.03E03,0,0);
    
    //cross section in 10^-39 cm^2 GeV ^-1 per nucleon!!
    //factor 2\pi because of integration over muon polar angle
    cout << E_in << " " << costhetamu << " " << Q2 << " " << omega << " " << x << " " << avgcross[0]*1.E16*2.*PI/Nucleus.getN()
      << " " << avgcross[1]*1.E16*2.*PI/Nucleus.getZ() << " " << rfgn*1.E16*2.*PI << " " << rfgp*1.E16*2.*PI << " " 
      << rfgn_bl*1.E16*2.*PI << " " << rfgp_bl*1.E16*2.*PI << " " << count << endl;
  }
  delete lepton;

}

//integrandum
void adap_intPm(numint::vector_d & results, double costhetacm,
		WeakQECross& pObs, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
		double E_out, double costhetamu, int current, double *cthmax, double E_in,
		double max_initial_nucl_mom, 
                double min_final_nucl_mom, int pw, int maxEvalweakamp){
		  
  results=numint::vector_d(2,0.);
  //we can fix the vector boson kinematics
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();
//   cout << E_in << " " << omega<< " " << Q2*1.E-06 << " " << sqrt(Q2+omega*omega) << endl;
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
// 	cout << "bla " << E_in << " " << costhetacm << " " << shell << " " << pm << " " << kin.GetPk() << endl;
      }
      else{
	double result=pObs.getDiffWeakQECross(kin,current,0,0,0,pw,shell,0.,maxEvalweakamp,0,1);
	results[(shell<nucleus.getPLevels()?1:0)]+= result; //results[0] neutrino, results[1] antineutrino
// 	cout << shell << " " << E_in <<  " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
// 	<< acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
// 	<< " " << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetQsquared() << " "
// 	<< kin.GetXb()*nucleus.getMassA()/MASSP << " " << result << endl;
      }
    }
    
  }
  //fold with flux
//  results[0]*=interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1);
//  results[1]*=interpolate(MiniBooNE_antineut_flux_norm,E_in,25,120,1);
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

