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
0.14760001E-03 , 
0.39853001E-03 , 
0.53648001E-03 , 
0.61304998E-03 , 
0.66586995E-03 , 
0.71846002E-03 , 
0.78075999E-03 , 
0.84440005E-03 , 
0.90257001E-03 , 
0.94977999E-03 , 
0.98400003E-03 , 
1.00589001E-03 , 
1.01067996E-03 , 
1.01136994E-03 , 
1.01204991E-03 , 
1.01820993E-03 , 
1.02638996E-03 , 
1.02365005E-03 , 
1.01616001E-03 , 
1.00451994E-03 , 
0.99015003E-03 , 
0.97646999E-03 , 
0.96209997E-03 , 
0.94156998E-03 , 
0.92031997E-03 , 
0.90051001E-03 , 
0.88066995E-03 , 
0.85808998E-03 , 
0.83551002E-03 , 
0.81224E-03 , 
0.78829002E-03 , 
0.76229E-03 , 
0.73556995E-03 , 
0.70890999E-03 , 
0.68290997E-03 , 
0.65526998E-03 , 
0.62796003E-03 , 
0.60140997E-03 , 
0.57555002E-03 , 
0.54974997E-03 , 
0.52402002E-03 , 
0.49829E-03 , 
0.47262999E-03 , 
0.44738001E-03 , 
0.42289001E-03 , 
0.39935002E-03 , 
0.37665999E-03 , 
0.35499999E-03 , 
0.33372E-03 , 
0.31277999E-03 , 
0.29212001E-03 , 
0.27185997E-03 , 
0.25223002E-03 , 
0.23341E-03 , 
0.21562001E-03 , 
0.19892E-03 , 
0.18325E-03 , 
0.16854E-03 , 
0.15478E-03 , 
0.14192E-03 , 
0.12988001E-03 , 
0.11859E-03 , 
0.10811E-03 , 
0.09853999E-03 , 
0.08971E-03 , 
0.08157E-03 , 
0.07411E-03 , 
0.06730001E-03 , 
0.06102E-03 , 
0.05521E-03 , 
0.04985E-03 , 
0.04503E-03 , 
0.04065E-03 , 
0.03668E-03 , 
0.03308E-03 , 
0.02988E-03 , 
0.02694E-03 , 
0.02424E-03 , 
0.02176E-03 , 
0.0196E-03 , 
0.01763E-03 , 
0.01581E-03 , 
0.01414E-03 , 
0.01268E-03 , 
0.01139E-03 , 
0.01024E-03 , 
0.00921E-03 , 
0.00825E-03 , 
0.0074E-03 , 
0.00667E-03 , 
0.00605E-03 , 
0.00544E-03 , 
0.00488E-03 , 
0.00437E-03 , 
0.00391E-03 , 
0.00351E-03 , 
0.00316E-03 , 
0.00286E-03 , 
0.00259E-03 , 
0.00232E-03 , 
0.00207E-03 , 
0.00185E-03 , 
0.00165E-03 , 
0.00149E-03 , 
0.00135E-03 , 
0.00123E-03 , 
0.00112E-03 , 
0.00101E-03 , 
0.00091E-03 , 
0.0008E-03 , 
0.00071E-03 , 
0.00065E-03 , 
0.0006E-03 , 
0.00054E-03 , 
0.00048E-03 , 
0.00044E-03 , 
0.00042E-03 , 
0.00039E-03 , 
0.00035E-03 , 
0.00031E-03 , 
};

const double MiniBooNE_neut_flux_norm[] = {
0.08814E-03 , 
0.23879001E-03 , 
0.33197999E-03 , 
0.39021999E-03 , 
0.43099001E-03 , 
0.47176E-03 , 
0.51836002E-03 , 
0.57660002E-03 , 
0.64455003E-03 , 
0.66979003E-03 , 
0.70666999E-03 , 
0.72996998E-03 , 
0.75520998E-03 , 
0.77657002E-03 , 
0.79403996E-03 , 
0.81733996E-03 , 
0.83868998E-03 , 
0.85615999E-03 , 
0.86974996E-03 , 
0.87946004E-03 , 
0.88528997E-03 , 
0.88722998E-03 , 
0.88916999E-03 , 
0.88528997E-03 , 
0.88334E-03 , 
0.87946004E-03 , 
0.87558001E-03 , 
0.86781001E-03 , 
0.86004996E-03 , 
0.8484E-03 , 
0.83675003E-03 , 
0.82315999E-03 , 
0.80763E-03 , 
0.79016E-03 , 
0.77267998E-03 , 
0.75326997E-03 , 
0.73580003E-03 , 
0.71638E-03 , 
0.69502997E-03 , 
0.67173004E-03 , 
0.65037E-03 , 
0.62901998E-03 , 
0.60571998E-03 , 
0.58241999E-03 , 
0.55913001E-03 , 
0.53582996E-03 , 
0.51252997E-03 , 
0.48923999E-03 , 
0.46400002E-03 , 
0.43682E-03 , 
0.41545999E-03 , 
0.39217001E-03 , 
0.36887002E-03 , 
0.34557E-03 , 
0.32422E-03 , 
0.30285999E-03 , 
0.28345001E-03 , 
0.26403001E-03 , 
0.24462E-03 , 
0.22714999E-03 , 
0.20967001E-03 , 
0.19395001E-03 , 
0.17861E-03 , 
0.16327E-03 , 
0.15143E-03 , 
0.13901E-03 , 
0.12754999E-03 , 
0.11707E-03 , 
0.10717E-03 , 
0.09804E-03 , 
0.08969E-03 , 
0.08193E-03 , 
0.07494E-03 , 
0.06852999E-03 , 
0.06271E-03 , 
0.05747E-03 , 
0.05261E-03 , 
0.04834E-03 , 
0.04426E-03 , 
0.04058E-03 , 
0.03728E-03 , 
0.03436E-03 , 
0.03165E-03 , 
0.02932E-03 , 
0.02699E-03 , 
0.02504E-03 , 
0.0231E-03 , 
0.02155E-03 , 
0.02E-03 , 
0.01866E-03 , 
0.0174E-03 , 
0.01619E-03 , 
0.01528E-03 , 
0.01439E-03 , 
0.01359E-03 , 
0.01287E-03 , 
0.01223E-03 , 
0.01165E-03 , 
0.01112E-03 , 
0.01062E-03 , 
0.01015E-03 , 
0.00971E-03 , 
0.00936E-03 , 
0.00897E-03 , 
0.00883E-03 , 
0.00878E-03 , 
0.00819E-03 , 
0.00819E-03 , 
0.00775E-03 , 
0.00744E-03 , 
0.00746E-03 , 
0.00726E-03 , 
0.00705E-03 , 
0.00685E-03 , 
0.0067E-03 , 
0.00656E-03 , 
0.00646E-03 , 
0.00635E-03 , 
0.00621E-03 , 
0.00598E-03 , 
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
  int integrator=2;//
  int fluxintegrator=atoi(argv[3]);
  int thick=0;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  bool charged=1;
  int current=2;
  
//   double omega=Q2/(2.*MASSP*Bjx);
  
  string homedir=argv[4];//"/home/wim/Code/share";

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
  if(!fluxintegrator) numint::cube_romb(mdf,lower,upper,1.E-20,1.E-03,avgcross,count,0);
  else numint::cube_adaptive(mdf,lower,upper,1.E-20,1.E-03,2E02,2E04,avgcross,count,0);
   
  //cross section in 10^-39 cm^2 GeV ^-1 per nucleon!!
  //factor 2\pi because of integration over muon polar angle
  cout << endl << endl;
  cout << T_mu << " " << costhetamu << " " << avgcross[0]*1.E16*2.*PI/Nucleus.getN()
  << " " << avgcross[1]*1.E16*2.*PI/Nucleus.getZ() << " " << count << endl;
  
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
//   cout << E_in << " " << omega<< " " << Q2*1.E-06 << " " << sqrt(Q2+omega*omega) << endl;
  pObs.getPlepton()->SetBeamEnergy(E_in);
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

