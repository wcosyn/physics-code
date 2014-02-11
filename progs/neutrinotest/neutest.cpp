//run ./neutest [T_\mu [MeV]] [cos(theta_mu)] 
#include <iostream>
#include <cstdlib>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <FsiCorrelator.hpp>
#include <FastParticle.hpp>
#include <constants.hpp>
#include <GlauberGrid.hpp>
#include <OneGlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <TLeptonKinematics.h>
#include <WeakQECross.hpp>

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
  int integrator=2;//atoi(argv[8]);
  int thick=1;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  bool charged=1;
  int current=2;
  
//   double omega=Q2/(2.*MASSP*Bjx);
  
  string homedir="/home/wim/Code/share";

  MeanFieldNucleusThick Nucleus(nucleus,homedir);
  //TKinematics2to2 kin("","",Carbon.getMassA(),Carbon.getMassA_min_proton(),MASSP,"qsquared:wlab:pklab",Q2,omega,pm);
//   TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(4627.);  //Monaghan Data
  TLeptonKinematics *lepton = TLeptonKinematics::CreateWithCosScatterAngle(TLeptonKinematics::muon,costhetamu);
  double E_out=T_mu+lepton->GetLeptonMass();

  //cout << lepton->GetBeamEnergy(kin) << " " << lepton->GetBeamEnergy(kin)-omega << " " << lepton->GetCosScatterAngle(kin) << endl;
  WeakQECross obs(lepton,&Nucleus,prec,integrator,homedir,charged, 1.03E03, screening, 0.);
  
  vector<double> avgcross(2,0.); 
  double cthmax[Nucleus.getTotalLevels()];
  double cthmin[Nucleus.getTotalLevels()];
  double max=-1.;

  double E_in=700;
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton->GetLeptonMass()*lepton->GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton->GetLeptonMass()*lepton->GetLeptonMass();
  double x=Q2/(2.*MASSP*omega);
  cout << omega << " " << Q2*1.E-06 << " " << x << endl;
  TKinematics2to2 kin("","",Nucleus.getMassA(),
		      Nucleus.getMassA_min_proton(),
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
  TKinematics2to2 kin2("","",Nucleus.getMassA(),
		      Nucleus.getMassA_min_proton(),
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,1.);
  cout << kin.GetPklab() << " " << kin2.GetPklab() << endl;
  
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
    cthmax[shell]=1.; cthmin[shell]=-1.;
    getBound(cthmax[shell],cthmin[shell],Nucleus,*lepton,E_in,E_out,costhetamu,shell);
    cout << shell << " " << cthmin[shell] << " " << cthmax[shell] << endl;
    if(max<cthmax[shell]) max=cthmax[shell];
  }

  E_in=3.E03;
  omega=E_in-E_out;
  Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton->GetLeptonMass()*lepton->GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton->GetLeptonMass()*lepton->GetLeptonMass();
  x=Q2/(2.*MASSP*omega);
  cout << omega << " " << Q2*1.E-06 << " " << x << endl;
  TKinematics2to2 kin5("","",Nucleus.getMassA(),
		      Nucleus.getMassA_min_proton(),
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
  TKinematics2to2 kin6("","",Nucleus.getMassA(),
		      Nucleus.getMassA_min_proton(),
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,1.);
  cout << kin5.GetPklab() << " " << kin6.GetPklab() << endl;
  
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
    cthmax[shell]=1.; cthmin[shell]=-1.;
    getBound(cthmax[shell],cthmin[shell],Nucleus,*lepton,E_in,E_out,costhetamu,shell);
    cout << shell << " " << cthmin[shell] << " " << cthmax[shell] << endl;
    if(max<cthmax[shell]) max=cthmax[shell];
  }

  E_in=E_out+10.;
  omega=E_in-E_out;
  Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton->GetLeptonMass()*lepton->GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton->GetLeptonMass()*lepton->GetLeptonMass();
  x=Q2/(2.*MASSP*omega);
  cout << omega << " " << Q2*1.E-06 << " " << x << endl;
  TKinematics2to2 kin3("","",Nucleus.getMassA(),
		      Nucleus.getMassA_min_proton(),
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
  TKinematics2to2 kin4("","",Nucleus.getMassA(),
		      Nucleus.getMassA_min_proton(),
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,1.);
  cout << kin3.GetPklab() << " " << kin4.GetPklab() << endl;
  
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
    cthmax[shell]=1.; cthmin[shell]=-1.;
    getBound(cthmax[shell],cthmin[shell],Nucleus,*lepton,E_in,E_out,costhetamu,shell);
    cout << shell << " " << cthmin[shell] << " " << cthmax[shell] << endl;
    if(max<cthmax[shell]) max=cthmax[shell];
  }
  
  
  exit(1);
//   cout << max << endl;
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
  numint::array<double,2> lower = {{E_out,-1.}};
  numint::array<double,2> upper = {{3.E03,max}};
  
  
  F.f=adap_intPm;
  unsigned count=0;
//   numint::cube_romb(mdf,lower,upper,1.E-20,1.E-03,totalcross,count,0);
  numint::cube_adaptive(mdf,lower,upper,1.E-20,1.E-03,2E03,2E06,avgcross,count,0);
   
  //cross section in 10^-39 cm^2 GeV ^-1
  cout << T_mu << " " << costhetamu << " " << avgcross[0]*1.E16 << " " << avgcross[1]*1.E16 << endl;
  
//   double neut,antineut;
//   obs.getDiffWeakQECross(kin,2,thick,0,0,1,0,0.,maxEval,1,1,neut,antineut);
//   cout << neut << " " << antineut << endl;
//   double free=obs.getElCross(kin,2,0.)*HBARC*HBARC;
//   vector<double> cross;
  delete lepton;

}

void adap_intPm(numint::vector_d & results, double E_in, double costhetacm,
		WeakQECross& pObs, MeanFieldNucleusThick &pNucleus, TLeptonKinematics &lepton,
		double E_out, double costhetamu, int current, double *cthmax){
		  
  results=numint::vector_d(2,0.);
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
	cout << "bla " << pm << endl;
      }
  //      numint::vector_d cross=numint::vector_d(5,0.);
  //     pObs->getAllDiffCross(cross,kin,current,shell,1,0.,2000000,0);
      double neut,antineut;
      pObs.getDiffWeakQECross(kin,current,1,0,0,1,shell,0.,2E04,0,1,neut,antineut);
      results[0]+=neut;
      results[1]+=antineut;
      cout << shell << " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
      << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
      << " " << kin.GetKlab() << " " << kin.GetWlab() << " " << neut << endl;
    }
    
  }
  
}

void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell){
  
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();
//   TKinematics2to2 kin1("","",nucleus.getMassA(),
// 		      nucleus.getMassA_min_proton()+nucleus.getExcitation()[shell],
// 		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,1.);
//   double pm1=kin1.GetPklab();
//   TKinematics2to2 kin2("","",nucleus.getMassA(),
// 		      nucleus.getMassA_min_proton()+nucleus.getExcitation()[shell],
// 		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
//   double pmmin1=kin2.GetPklab();
//   cout << pm1 << " " << pmmin1 << endl;
  
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

