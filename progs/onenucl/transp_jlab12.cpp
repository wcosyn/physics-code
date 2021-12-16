//program used to do some transparency calculations for Or et al. for a proposal.
//last used to compute jlab 12Gev A(e,e'p) and n transparencies, no phiaveraging in cross sectionnnnn!!!

#include <iostream>
#include <cstdlib>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <constants.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <Cross.hpp>
#include <Utilfunctions.hpp>
#include <vector>
#include <numint/numint.hpp>


void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, double Q2, double omega, int shell);
void adap_intPm(numint::vector_d &, double costhetacm,Cross* pObs, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
		double Q2, double omega, int current, int thick, double lc_mod, double nkt_mod, double *cthmax);


//run ./transp_jlab12 [Q2,GeV^2] [Ein, MeV] [precision in integration] [userset coherence length factor] [userset initial sigma CT value]
int main(int argc, char *argv[])
{
  
  string homedir=HOMEDIR;
  //kinematics input
  double Q2=atof(argv[3])*1.E06; // input in GeV^2 [8.0,9.4,11.4,14.2]
  double Ein=atof(argv[2])*1.E03;  //input in GeV [6.4,10.6,10.6,10.6]
  
  //calculate kinematics
  double omega=Q2/(2.*MASSn);  //Bjorken x was chosen super close tot 1 in exp kinematics
  double q=sqrt(Q2+omega*omega);
  cout << Q2 << " " << omega << " " << q << " " << Q2/4./Ein/(Ein-omega) << endl;
  
  //glauber parameters
  int thick=0;
  int current=2;
  double prec=atof(argv[4]);
  bool userset=0;//atoi(argv[8]);
  double screening=0.;//atof(argv[9]);
  double lc_mod = atof(argv[5]);
  double nkt_mod = atof(argv[6]);
  
  //create objects for nucleus, kinematics and cross section
  MeanFieldNucleusThick Nucleus(atoi(argv[1]),homedir);
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
  Cross obs(*elec,&Nucleus,prec, 2, homedir, userset, screening);

  //arrays to collect results
  //all levels (glauber, pw) + total proton + total neutron    
  vector<double> totalcross(3*Nucleus.getTotalLevels()+6,0.); 
  double cthmax[Nucleus.getTotalLevels()];
  double cthmin[Nucleus.getTotalLevels()];
  double max=-1.;

  //determine integration limits so pm<300 MeV
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
    cthmax[shell]=1.; cthmin[shell]=-1.;
    getBound(cthmax[shell],cthmin[shell],Nucleus,Q2,omega,shell);
//     cout << cthmax[shell] << endl;
    if(max<cthmax[shell]) max=cthmax[shell];
  }
//   cout << max << endl;

  //integration structure.
  struct Ftor {

    static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
      Ftor &p = * (Ftor *) param;
      p.f(ret,x[0],p.pObs,p.pNucleus,p.elec,p.Q2,p.omega,p.current,p.thick,p.lc_mod, p.nkt_mod, p.cthmax);
    }
    Cross *pObs;
    MeanFieldNucleusThick *pNucleus;
    TElectronKinematics *elec;
    double Q2;
    double omega;
    int current;
    int thick;
    double lc_mod;
    double nkt_mod;
    double *cthmax;
    void (*f)(numint::vector_d &, double pm,Cross* pObs, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
	  double Q2, double omega, int current, int thick, double lc_mod, double nkt_mod, double *cthmax);

  };

  //initialize integration object + parameters needed to calculate integrals
  Ftor F;
  F.pObs = &obs;
  F.pNucleus = &Nucleus;
  F.elec=elec;
  F.Q2=Q2;
  F.omega=omega;
  F.current=current;
  F.thick=thick;
  F.cthmax=cthmax;
  F.lc_mod = lc_mod;
  F.nkt_mod = nkt_mod;

  numint::mdfunction<numint::vector_d,1> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

  unsigned neval = 0;
  numint::array<double,1> lower = {{-1.}};
  numint::array<double,1> upper = {{max}};
  
  
  F.f=adap_intPm;
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
  cout << totalcross[3*Nucleus.getTotalLevels()] << " " << totalcross[3*Nucleus.getTotalLevels()+1] << " " << totalcross[3*Nucleus.getTotalLevels()+2] << endl;

  //neutron shells transparencies
  cout << Q2/1.E06 << " " << omega << " " << q  << " ";
  //neutron s1/2, p3/2 no CT
  for(int shell=Nucleus.getPLevels();shell<Nucleus.getTotalLevels();shell++) {
    cout << totalcross[3*shell]/totalcross[3*shell+2] << " ";
  }
  //neutron s1/2, p3/2 with CT
  for(int shell=Nucleus.getPLevels();shell<Nucleus.getTotalLevels();shell++) {
    cout << totalcross[3*shell+1]/totalcross[3*shell+2] << " ";
  }
  //cout << "neutron shells total cross no CT, CT, plane-wave"
  cout << totalcross[3*Nucleus.getTotalLevels()+3] << " " << totalcross[3*Nucleus.getTotalLevels()+4] 
  << " " << totalcross[3*Nucleus.getTotalLevels()+5] << endl;

  //total proton and neutron transparency, proton no CT, with CT; neutron no CT, with CT
  cout << totalcross[3*Nucleus.getTotalLevels()]/totalcross[3*Nucleus.getTotalLevels()+2] 
        << " " << totalcross[3*Nucleus.getTotalLevels()+1]/totalcross[3*Nucleus.getTotalLevels()+2] 
        << " " << totalcross[3*Nucleus.getTotalLevels()+3]/totalcross[3*Nucleus.getTotalLevels()+5] 
        << " " << totalcross[3*Nucleus.getTotalLevels()+4]/totalcross[3*Nucleus.getTotalLevels()+5] << endl;

  //neutron nominator and denominator cross sections
//   for(int shell=Nucleus.getPLevels();shell<Nucleus.getTotalLevels();shell++) {
//     cout << totalcross[2*shell] << " ";
//   }
//   cout << endl;
//   for(int shell=Nucleus.getPLevels();shell<Nucleus.getTotalLevels();shell++) {
//     cout << totalcross[2*shell+1] << " ";
//   }
//   cout << endl;
  
  delete elec;
  return 0;
}


//integrandum
void adap_intPm(numint::vector_d & results, double costhetacm,Cross* pObs, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
		double Q2, double omega, int current, int thick, double lc_mod, double nkt_mod, double *cthmax){
		  
  results=numint::vector_d(3*pNucleus->getTotalLevels()+6,0.);
  for(int shell=0;shell<pNucleus->getPLevels();shell++) {
    if(costhetacm<cthmax[shell]){
      TKinematics2to2 kin("","",pNucleus->getMassA(),
        pNucleus->getMassA_min_proton()+pNucleus->getExcitation()[shell],
        MASSP,"qsquared:wlab:costhkcm",Q2,omega,costhetacm);
      double pm=kin.GetPklab();
      if(!kin.IsPhysical()){
        // for(int i=0;i<5;i++) results[i]+=0.;
        cout << "bla " << pm << endl;
      }
      numint::vector_d cross=numint::vector_d(thick? 5:3,0.);
      pObs->getAllDiffCross(cross,kin,current,shell,0,0.,20000,0,0, lc_mod, nkt_mod);
      // double pw=pObs->getDiffCross(kin,current,1,0,0,1,shell,0.,20000,0,1);
      results[3*shell]+=cross[0];
      results[3*shell+1]+=cross[thick?2:1];
      results[3*shell+2]+=cross[thick?4:2];
      results[3*pNucleus->getTotalLevels()]+=cross[0];
      results[3*pNucleus->getTotalLevels()+1]+=cross[thick?2:1];
      results[3*pNucleus->getTotalLevels()+2]+=cross[thick?4:2];
      // cout << 2*shell << " " << 2*shell+1 << " " << 2*pNucleus->getTotalLevels() << " " << 2*pNucleus->getTotalLevels()+1 << endl;
      cout << "0 " << shell << " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
      << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
      << " " << kin.GetKlab() << " " << kin.GetWlab() <<  " " << results[3*shell] << " " << results[3*shell+2] << endl;
    }  
  }

  for(int shell=pNucleus->getPLevels();shell<pNucleus->getTotalLevels();shell++) {
    if(costhetacm<cthmax[shell]){
      TKinematics2to2 kin("","",pNucleus->getMassA(),
        pNucleus->getMassA_min_neutron()+pNucleus->getExcitation()[shell],
        MASSN,"qsquared:wlab:costhkcm",Q2,omega,costhetacm);
      double pm=kin.GetPklab();
      if(!kin.IsPhysical()){
        //for(int i=5;i<10;i++) results[i]+=0.;
        cout << "bla " << pm << endl;
      }
      numint::vector_d cross=numint::vector_d(thick?5:3,0.);
      pObs->getAllDiffCross(cross,kin,current,shell,0,0.,20000,0,0, lc_mod, nkt_mod);
      // double pw=pObs->getDiffCross(kin,current,1,0,0,1,shell,0.,20000,0,1);
      //for(int i=5;i<10;++i) results[i]+=cross[i-5];
      results[3*shell]+=cross[0];
      results[3*shell+1]+=cross[thick?2:1];
      results[3*shell+2]+=cross[thick?4:2];
      results[3*pNucleus->getTotalLevels()+3]+=cross[0];
      results[3*pNucleus->getTotalLevels()+4]+=cross[thick?2:1];
      results[3*pNucleus->getTotalLevels()+5]+=cross[thick?4:2];
      // cout << 2*shell << " " << 2*shell+1 << " " << 2*pNucleus->getTotalLevels()+2 << " " << 2*pNucleus->getTotalLevels()+3 << endl;
      cout << "1 " << shell << " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
      << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
      << " " << kin.GetKlab() << " " << kin.GetWlab() << " "  << results[3*shell] << " " << results[3*shell+2] << endl;
    }
    
  }
}




// find cm scatt angles so that missing momentum stays below 300 MeV 
void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, double Q2, double omega, int shell){
  TKinematics2to2 kin("","",nucleus.getMassA(),
		      (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())+nucleus.getExcitation()[shell],
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,(high+low)/2.);
  double pm=kin.GetPklab();
  if(pm<300.) low=(high+low)/2.;
  else high=(high+low)/2.;
  if((high-low)<1.E-04) { return;}
  else getBound(high,low,nucleus,Q2,omega,shell);  
}

