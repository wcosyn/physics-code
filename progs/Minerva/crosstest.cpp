//run ./neutest [T_\mu [MeV]] [cos(theta_mu)] [integrator] [homedir]

//we integrate for miniboone cross section
//integration happens over incoming neutrino energy folded with flux
//and also proton center of mass scattering angle (easier to keep track of physical points and such)
//it takes some effort in finding reasonable integration limits, the rest is basically integrating the QE cross section of course
//pretty straightforward...

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

const double massmu = 105.6583715;

double Minerva_nu_flux[] = {
2.57268e-06,
6.5321e-06,
1.69721e-05,
2.51453e-05,
3.31235e-05,
4.07319e-05,
4.27611e-05,
3.41954e-05,
2.04086e-05,
1.10596e-05,
6.78507e-06,
4.86896e-06,
3.94903e-06,
3.34018e-06,
2.90956e-06,
2.5465e-06,
2.28787e-06,
2.04961e-06,
1.85345e-06,
1.69827e-06,
};

double Minerva_anu_flux[] = {
2.32864e-06,
6.25619e-06,
1.60002e-05,
2.295e-05,
2.99316e-05,
3.61717e-05,
3.64453e-05,
2.82484e-05,
1.62189e-05,
8.36204e-06,
4.95442e-06,
3.39086e-06,
2.56585e-06,
2.08784e-06,
1.74572e-06,
1.46005e-06,
1.2336e-06,
1.07991e-06,
9.67436e-07,
8.7012e-07,
};


//integration struct
struct Ftor {  //Carbon

  static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
    Ftor &p = * (Ftor *) param;
    p.f(ret,x[0],x[1],*p.pObs,*p.pNucleus,*p.lepton,p.E_out,p.costhetamu,p.current,p.cthmax);
  }
  WeakQECross *pObs;
  MeanFieldNucleusThick *pNucleus;
  TLeptonKinematics *lepton;
  double E_out;
  double costhetamu;
  int current;
  double *cthmax;
  void (*f)(numint::vector_d &, double Ein, double costhetacm,
	    WeakQECross& pObs, MeanFieldNucleusThick &pNucleus, TLeptonKinematics &lepton,
	double E_out, double costhetamu, int current, double *cthmax);

};

void normalize(double flux[], double dx){
  double sum=0;
  for(int i=0;i<20;i++) sum+=flux[i]*dx;
  for(int i=0;i<20;i++) flux[i]/=sum;
  };

//determine boundaries of kinematics
double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell);
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell);
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell);

// integration over missing momentum
void adap_intPm(numint::vector_d &, double E_in, double costhetacm,
		WeakQECross& pObs, MeanFieldNucleusThick &pNucleus, TLeptonKinematics &lepton,
		double E_out, double costhetamu, int current, double *cthmax);


int main(int argc, char *argv[])
{
//   double Q2=atof(argv[1])*1.E06;
  double T_mu=atof(argv[1]); // muon kinetic energy in MeV
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
  
//   double omega=Q2/(2.*MASSP*Bjx);
  
  string homedir=HOMEDIR;//"/home/wim/Code/share";

  MeanFieldNucleusThick Nucleus(nucleus,homedir);
  TLeptonKinematics *lepton = TLeptonKinematics::CreateWithCosScatterAngle(TLeptonKinematics::muon,costhetamu); 

  double E_out=T_mu+lepton->GetLeptonMass();

  WeakQECross obs(lepton,&Nucleus,prec,integrator,homedir,charged, 1.03E03, screening, 0.);
  
  vector<double> avgcross(2,0.); //neutrino and antineutrino
  //estimates for integration bounds
  double cthmax[Nucleus.getTotalLevels()];
  double cthmin[Nucleus.getTotalLevels()];
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cthmax[shell]=-1.;cthmin[shell]=1.;}
  double max=-1.,min=1.; //overall min and max center of mass cosine theta angles

  //find reasonable integration limits
  double E_low=E_out;
  double E_high=E_out;
  bool switchlow=0;
  for(int i=1;i<1000;i++){
    double E_in=E_out+(3.E03-E_out)*0.001*i; //possible incoming lepton energies
    double omega=E_in-E_out;
    double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton->GetLeptonMass()*lepton->GetLeptonMass()/(E_out*E_out))*costhetamu)
	-lepton->GetLeptonMass()*lepton->GetLeptonMass();
    double x=Q2/(2.*MASSP*omega);
//     cout << E_in << " " << omega << " " << Q2*1.E-06 << " " << x <<  endl;
    for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
      //lowest p_m is always at theta_cm -1
      TKinematics2to2 kin("","",Nucleus.getMassA(),
			  (shell<Nucleus.getPLevels()? Nucleus.getMassA_min_proton(): Nucleus.getMassA_min_neutron())
			  +Nucleus.getExcitation()[shell],
			  shell<Nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
	double tempmax=-1., tempmin=1.;
	
      //anything above 300 MeV contribution will be negligible
      if(getBound(tempmax,tempmin,Nucleus,*lepton,E_in,E_out,costhetamu,shell)<300.){ 
	if(!switchlow){
//	  cout<< E_in << " " << kin.GetPklab() << " " << shell << endl;
	  switchlow=1;
	  E_low=E_in;
	}
	else if(E_in>E_high) E_high=E_in;
// 	cout << shell << " " << kin.GetPklab() <<" " << cthmin[shell] << " " << cthmax[shell] << endl;
	if(max<tempmax) max=tempmax;
	if(min>tempmin) min=tempmin;
	if(cthmax[shell]<tempmax) cthmax[shell]=tempmax;
	if(cthmin[shell]>tempmin) cthmin[shell]=tempmin;
	
      }
    }
  }

  //Normalize flux
  normalize(Minerva_nu_flux,500);
  normalize(Minerva_anu_flux,500);
  
  //initialize object
  Ftor F;
  F.pObs = &obs;
  F.pNucleus = &Nucleus;
  F.lepton=lepton;
  F.E_out=E_out;
  F.costhetamu=costhetamu;
  F.current=current;
  F.cthmax=cthmax;

  numint::mdfunction<numint::vector_d,2> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

  unsigned neval = 0;
  numint::array<double,2> lower = {{E_low,min}};
  numint::array<double,2> upper = {{E_high,max}};
  
  
  F.f=adap_intPm;
  unsigned count=0;
  if(!fluxintegrator) numint::cube_romb(mdf,lower,upper,1.E-18,1.E-02,avgcross,count,0); //1.E-20,1.E-03
  else numint::cube_adaptive(mdf,lower,upper,1.E-18,1.E-02,2E02,2E04,avgcross,count,0);
   
  //cross section in 10^-39 cm^2 GeV ^-1 per nucleon!!
  //factor 2\pi because of integration over muon polar angle
//cout << endl << endl << endl;
  cout << T_mu << " " << costhetamu << " " << avgcross[0]*1.E16*2.*PI/Nucleus.getN()
  << " " << avgcross[1]*1.E16*2.*PI/Nucleus.getZ() << " " << count << endl;
  
  delete lepton;

}

//integrandum
void adap_intPm(numint::vector_d & results, double E_in, double costhetacm,
		WeakQECross& pObs, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
		double E_out, double costhetamu, int current, double *cthmax){
		  
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
      if(!kin.IsPhysical()||kin.GetPYlab()<200.){ //final nucleon momentum too low, impose cut!
	//for(int i=0;i<2;i++) results[i]+=0.;
// 	cout << "bla " << E_in << " " << costhetacm << " " << shell << " " << pm << " " << kin.GetPk() << endl;
      }
      else{
	double result=pObs.getDiffWeakQECross(kin,current,1,0,0,1,shell,0.,2E04,0,1);
  results[(shell<nucleus.getPLevels()?1:0)]+= result; //results[0] neutrino, results[1] antineutrino
//	cout << shell << " " << E_in <<  " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
//	<< acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
//	<< " " << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetQsquared() << " "
//	<< kin.GetXb()*nucleus.getMassA()/MASSP << " " << result << " " << result*(shell<nucleus.getPLevels()?
//	interpolate(MiniBooNE_antineut_flux_norm,E_in,25,120,1):interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1)) << endl;
      }
    }
    
  }
  //fold with flux
  results[0]*=interpolate(Minerva_nu_flux,E_in,500,20,0);
  results[1]*=interpolate(Minerva_anu_flux,E_in,500,20,0);
}

double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell){
  
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();
  double tempmin=500.;
  double costhmin=-1.;
  for(int i=0;i<=1000;i++){ //loop over all com angles
    TKinematics2to2 kin("","",nucleus.getMassA(),
			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			+nucleus.getExcitation()[shell],
			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.+i*0.002);
    double pm=kin.GetPklab(); //initial nucleon momentum
    if(pm<tempmin) {tempmin=pm; costhmin=-1.+i*0.002;}
  }
  //   cout << omega << " " << Q2/1.E06 << " " << omega/(2.*MASSP*omega) << " " << high << " " << low << " " << pm << endl;
  double temptemp=costhmin;//com costheta value where the min initial nucl momentum is
//   cout << costhmin << " " << tempmin << endl;
  if(tempmin<300.){
    TKinematics2to2 kin1("","",nucleus.getMassA(),
			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			+nucleus.getExcitation()[shell],
			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,1.);
    high=1;
    if(kin1.GetPklab()>300.) high=getMax(high,temptemp,nucleus,lepton,Q2,omega,shell);
    temptemp=costhmin;
    TKinematics2to2 kinmin1("","",nucleus.getMassA(),
			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			+nucleus.getExcitation()[shell],
			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
    low=-1.;
    if(kinmin1.GetPklab()>300.) getMin(temptemp,low,nucleus,lepton,Q2,omega,shell);
//     cout << "pm " << kinmin1.GetPklab() << " " << kin1.GetPklab() << " " << costhmin << " " << low << " " << high << endl;
  }
  return tempmin;//minimal initial nucleon momentum
}

//recursive function to find min costheta value so that initial nucleon momentum is low enough
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell){
  
  TKinematics2to2 kin("","",nucleus.getMassA(),
		      (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
		      +nucleus.getExcitation()[shell],
		      shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,(high+low)/2.);
  double pm=kin.GetPklab();
//   cout << omega << " " << Q2/1.E06 << " " << omega/(2.*MASSP*omega) << " " << high << " " << low << " " << pm << endl;
  if(pm<300.) high=(high+low)/2.;
  else low=(high+low)/2.;
  if((high-low)<1.E-04) { return high;}
  else return getMin(high,low,nucleus,lepton,Q2,omega,shell);  
}

//recursive function to find max costheta value so that initial nucleon momentum is low enough
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell){
  
  TKinematics2to2 kin("","",nucleus.getMassA(),
		      (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
		      +nucleus.getExcitation()[shell],
		      shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,(high+low)/2.);
  double pm=kin.GetPklab();
//   cout << omega << " " << Q2/1.E06 << " " << omega/(2.*MASSP*omega) << " " << high << " " << low << " " << pm << endl;
  if(pm<300.) low=(high+low)/2.;
  else high=(high+low)/2.;
  if((high-low)<1.E-04) { return low;}
  else return getMax(high,low,nucleus,lepton,Q2,omega,shell);  
}

