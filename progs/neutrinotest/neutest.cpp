//run ./neutest [T_\mu [MeV]] [cos(theta_mu)] 

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

//starts at E=25 MeV with a delta of 25 MeV
const double MiniBooNE_antineut_flux_norm[] = {
0.14760001 , 
0.39853001 , 
0.53648001 , 
0.61304998 , 
0.66586995 , 
0.71846002 , 
0.78075999 , 
0.84440005 , 
0.90257001 , 
0.94977999 , 
0.98400003 , 
1.00589001 , 
1.01067996 , 
1.01136994 , 
1.01204991 , 
1.01820993 , 
1.02638996 , 
1.02365005 , 
1.01616001 , 
1.00451994 , 
0.99015003 , 
0.97646999 , 
0.96209997 , 
0.94156998 , 
0.92031997 , 
0.90051001 , 
0.88066995 , 
0.85808998 , 
0.83551002 , 
0.81224 , 
0.78829002 , 
0.76229 , 
0.73556995 , 
0.70890999 , 
0.68290997 , 
0.65526998 , 
0.62796003 , 
0.60140997 , 
0.57555002 , 
0.54974997 , 
0.52402002 , 
0.49829 , 
0.47262999 , 
0.44738001 , 
0.42289001 , 
0.39935002 , 
0.37665999 , 
0.35499999 , 
0.33372 , 
0.31277999 , 
0.29212001 , 
0.27185997 , 
0.25223002 , 
0.23341 , 
0.21562001 , 
0.19892 , 
0.18325 , 
0.16854 , 
0.15478 , 
0.14192 , 
0.12988001 , 
0.11859 , 
0.10811 , 
0.09853999 , 
0.08971 , 
0.08157 , 
0.07411 , 
0.06730001 , 
0.06102 , 
0.05521 , 
0.04985 , 
0.04503 , 
0.04065 , 
0.03668 , 
0.03308 , 
0.02988 , 
0.02694 , 
0.02424 , 
0.02176 , 
0.0196 , 
0.01763 , 
0.01581 , 
0.01414 , 
0.01268 , 
0.01139 , 
0.01024 , 
0.00921 , 
0.00825 , 
0.0074 , 
0.00667 , 
0.00605 , 
0.00544 , 
0.00488 , 
0.00437 , 
0.00391 , 
0.00351 , 
0.00316 , 
0.00286 , 
0.00259 , 
0.00232 , 
0.00207 , 
0.00185 , 
0.00165 , 
0.00149 , 
0.00135 , 
0.00123 , 
0.00112 , 
0.00101 , 
0.00091 , 
0.0008 , 
0.00071 , 
0.00065 , 
0.0006 , 
0.00054 , 
0.00048 , 
0.00044 , 
0.00042 , 
0.00039 , 
0.00035 , 
0.00031 , 
};

const double MiniBooNE_neut_flux_norm[] = {
0.08814 , 
0.23879001 , 
0.33197999 , 
0.39021999 , 
0.43099001 , 
0.47176 , 
0.51836002 , 
0.57660002 , 
0.64455003 , 
0.66979003 , 
0.70666999 , 
0.72996998 , 
0.75520998 , 
0.77657002 , 
0.79403996 , 
0.81733996 , 
0.83868998 , 
0.85615999 , 
0.86974996 , 
0.87946004 , 
0.88528997 , 
0.88722998 , 
0.88916999 , 
0.88528997 , 
0.88334 , 
0.87946004 , 
0.87558001 , 
0.86781001 , 
0.86004996 , 
0.8484 , 
0.83675003 , 
0.82315999 , 
0.80763 , 
0.79016 , 
0.77267998 , 
0.75326997 , 
0.73580003 , 
0.71638 , 
0.69502997 , 
0.67173004 , 
0.65037 , 
0.62901998 , 
0.60571998 , 
0.58241999 , 
0.55913001 , 
0.53582996 , 
0.51252997 , 
0.48923999 , 
0.46400002 , 
0.43682 , 
0.41545999 , 
0.39217001 , 
0.36887002 , 
0.34557 , 
0.32422 , 
0.30285999 , 
0.28345001 , 
0.26403001 , 
0.24462 , 
0.22714999 , 
0.20967001 , 
0.19395001 , 
0.17861 , 
0.16327 , 
0.15143 , 
0.13901 , 
0.12754999 , 
0.11707 , 
0.10717 , 
0.09804 , 
0.08969 , 
0.08193 , 
0.07494 , 
0.06852999 , 
0.06271 , 
0.05747 , 
0.05261 , 
0.04834 , 
0.04426 , 
0.04058 , 
0.03728 , 
0.03436 , 
0.03165 , 
0.02932 , 
0.02699 , 
0.02504 , 
0.0231 , 
0.02155 , 
0.02 , 
0.01866 , 
0.0174 , 
0.01619 , 
0.01528 , 
0.01439 , 
0.01359 , 
0.01287 , 
0.01223 , 
0.01165 , 
0.01112 , 
0.01062 , 
0.01015 , 
0.00971 , 
0.00936 , 
0.00897 , 
0.00883 , 
0.00878 , 
0.00819 , 
0.00819 , 
0.00775 , 
0.00744 , 
0.00746 , 
0.00726 , 
0.00705 , 
0.00685 , 
0.0067 , 
0.00656 , 
0.00646 , 
0.00635 , 
0.00621 , 
0.00598 , 
};


  //integration struct
struct Ftor {

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

  
void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell);
void adap_intPm(numint::vector_d &, double E_in, double costhetacm,
		WeakQECross& pObs, MeanFieldNucleusThick &pNucleus, TLeptonKinematics &lepton,
		double E_out, double costhetamu, int current, double *cthmax);


int main(int argc, char *argv[])
{
//   double Q2=atof(argv[1])*1.E06;
  double T_mu=atof(argv[1]);
  double costhetamu=atof(argv[2]);
  bool screening=0;//atoi(argv[4]);
  int nucleus=1;//atoi(argv[6]);
  double prec=1.E-05;//atof(argv[7]);
  int integrator=atoi(argv[3]);
  int thick=0;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  bool charged=1;
  int current=2;
  
//   double omega=Q2/(2.*MASSP*Bjx);
  
  string homedir="/home/wim/Code/share";

  MeanFieldNucleusThick Nucleus(nucleus,homedir);
  TLeptonKinematics *lepton = TLeptonKinematics::CreateWithCosScatterAngle(TLeptonKinematics::muon,costhetamu);

  double E_out=T_mu+lepton->GetLeptonMass();

  WeakQECross obs(lepton,&Nucleus,prec,integrator,homedir,charged, 1.03E03, screening, 0.);
  
  vector<double> avgcross(2,0.); 
  //estimates for integration bounds
  double cthmax[Nucleus.getTotalLevels()];
  double cthmin[Nucleus.getTotalLevels()];
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cthmax[shell]=-1.;cthmin[shell]=-1;}
  double max=-1.;

  //find reasonable integration limits
  double E_low=E_out;
  double E_high=E_out;
  bool switchlow=0;
  for(int i=1;i<1000;i++){
    double E_in=E_out+(3.E03-E_out)*0.001*i;
    double omega=E_in-E_out;
    double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton->GetLeptonMass()*lepton->GetLeptonMass()/(E_out*E_out))*costhetamu)
	-lepton->GetLeptonMass()*lepton->GetLeptonMass();
    double x=Q2/(2.*MASSP*omega);
//     cout << E_in << " " << omega << " " << Q2*1.E-06 << " " << x <<  endl;
    for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
      //lowest p_m is always at theta_cm -1
      TKinematics2to2 kin("","",Nucleus.getMassA(),
			  Nucleus.getMassA_min_proton()+Nucleus.getExcitation()[shell],
			  MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
      //anything above 300 MeV contribution will be negligible
      if(kin.GetPklab()<300.){ 
	if(!switchlow){
	  cout<< E_in << " " << kin.GetPklab() << " " << shell << endl;
	  switchlow=1;
	  E_low=E_in;
	}
	else if(E_in>E_high) E_high=E_in;
	double tempmax=1., tempmin=-1.;
	getBound(tempmax,tempmin,Nucleus,*lepton,E_in,E_out,costhetamu,shell);
// 	cout << shell << " " << kin.GetPklab() <<" " << cthmin[shell] << " " << cthmax[shell] << endl;
	if(max<tempmax) max=tempmax;
	if(cthmax[shell]<tempmax) cthmax[shell]=tempmax;
      }
    }
  }
  cout << E_low << " " << E_high << " " << max << endl;
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cout << shell << " " << cthmax[shell] << " " << cthmin[shell] << endl;}
   
   


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
  numint::array<double,2> lower = {{E_low,-1.}};
  numint::array<double,2> upper = {{E_high,max}};
  
  
  F.f=adap_intPm;
  unsigned count=0;
  if(!integrator) numint::cube_romb(mdf,lower,upper,1.E-20,1.E-03,avgcross,count,0);
  else numint::cube_adaptive(mdf,lower,upper,1.E-20,1.E-03,2E02,2E04,avgcross,count,0);
   
  //cross section in 10^-39 cm^2 GeV ^-1
  //factor 2\pi because of integration over muon polar angle
  cout << endl << endl;
  cout << T_mu << " " << costhetamu << " " << avgcross[0]*1.E16*2.*PI << " " << avgcross[1]*1.E16*2.*PI << " " << count << endl;
  
  delete lepton;

}

//integrandum
void adap_intPm(numint::vector_d & results, double E_in, double costhetacm,
		WeakQECross& pObs, MeanFieldNucleusThick &pNucleus, TLeptonKinematics &lepton,
		double E_out, double costhetamu, int current, double *cthmax){
		  
  results=numint::vector_d(2,0.);
  
  //we can fix the vector boson kinematics
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();
  for(int shell=0;shell<pNucleus.getTotalLevels();shell++) {
    if(costhetacm<cthmax[shell]){
      TKinematics2to2 kin("","",pNucleus.getMassA(),
			  pNucleus.getMassA_min_proton()+pNucleus.getExcitation()[shell],
			  MASSP,"qsquared:wlab:costhkcm",Q2,omega,costhetacm);
      double pm=kin.GetPklab();
      if(!kin.IsPhysical()){
	for(int i=0;i<2;i++) results[i]+=0.;
	cout << "bla " << E_in << " " << costhetacm << " " << shell << " " << pm << endl;
      }
      else{
	double result=pObs.getDiffWeakQECross(kin,current,1,0,0,1,shell,0.,2E04,0,1);
	results[(shell<pNucleus.getPLevels()?1:0)]+= result;
	cout << shell << " " << E_in <<  " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
	<< acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
	<< " " << kin.GetKlab() << " " << kin.GetWlab() << " " << result << endl;
      }
    }
    
  }
  //fold with flux
  results[0]*=interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1);
  results[1]*=interpolate(MiniBooNE_antineut_flux_norm,E_in,25,120,1);
}

void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell){
  
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();
  
  TKinematics2to2 kin("","",nucleus.getMassA(),
		      nucleus.getMassA_min_proton()+nucleus.getExcitation()[shell],
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,(high+low)/2.);
  double pm=kin.GetPklab();
//   cout << omega << " " << Q2/1.E06 << " " << omega/(2.*MASSP*omega) << " " << high << " " << low << " " << pm << endl;
  if(pm<300.) low=(high+low)/2.;
  else high=(high+low)/2.;
  if((high-low)<1.E-04) { return;}
  else getBound(high,low,nucleus,lepton,E_in,E_out,costhetamu,shell);  
}

