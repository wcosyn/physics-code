//program used to do some transparency calculations for Or et al. for a proposal.
//last used to compute Hall B p and n transparencies, phiaveraging in cross section is on!!!

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


void intPm(const double pm, double *results, va_list ap);
void intcosth(const double costh, double *results, va_list ap);
void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, double Q2, double omega, int shell);
void adap_intPm(numint::vector_d &, double costhetacm,Cross* pObs, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
		double Q2, double omega, int current, double *cthmax);


//run ./transp [nucleus] [Q2,GeV^2] [Ein, MeV] [thetae, degr] [thickness] [current] [precision in integration] [userset FSI sigma] [screening sigma value]
int main(int argc, char *argv[])
{
  
  // string homedir=HOMEDIR;
  // double Q2=atof(argv[2])*1.E06;
  // double Ein=atof(argv[3]);
  // double thetae=atof(argv[4])*DEGRTORAD;
  // double omega=Ein-Q2/(4.*Ein*pow(sin(thetae/2.),2.));
  // int thick=atoi(argv[5]);
  // int current=atoi(argv[6]);
  // double prec=atof(argv[7]);
  // bool userset=atoi(argv[8]);
  // double screening=atof(argv[9]);
  // //cout << Ein-omega << endl;
   string homedir=argv[4];
  double Q2=atof(argv[2])*1.E06;
  double Ein=5010;
  //double q=atof(argv[4]);
  double omega=Q2/(2.*MASSn);
  double q=sqrt(Q2+omega*omega);
  cout << Q2 << " " << omega << " " << q << " " << Q2/4./Ein/(Ein-omega) << endl;
  int thick=1;
  int current=2;
  double prec=atof(argv[3]);
  bool userset=0;//atoi(argv[8]);
  double screening=0.;//atof(argv[9]);
  //cout << Ein-omega << endl;
  
  MeanFieldNucleusThick Nucleus(atoi(argv[1]),homedir);
  //TKinematics2to2 kin("","",.getMassA(),Carbon.getMassA_min_proton(),MASSP,"qsquared:wlab:pklab",atof(argv[1]),atof(argv[2]),atof(argv[3]));
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
  Cross obs(*elec,&Nucleus,prec, 2, homedir, userset, screening);

//   double k=sqrt(Q2+omega*omega);
//   double Eout=Ein-omega;
//   //cout << Q2 << " " << omega << " " << k << " " << Q2/(2.*MASSP*omega) << " " << atan2(-Eout*sin(thetae),Ein-Eout*cos(thetae))*RADTODEGR << endl;
// 
//   for(int i=0;i<=20;i++){
//     double Eoutnew=Eout*(1.+(i-10.)/10.*0.15);
//     omega=Ein-Eoutnew;
//     Q2=4.*Ein*Eoutnew*pow(sin(thetae/2.),2.);
//     cout << Q2*1.E-06 << " " << Q2/(2.*MASSP*omega) << endl;
//     double thetap=0.;
//     for(int shell=0;shell<Nucleus.getPLevels();shell++) {
//       TKinematics2to2 kin("","",Nucleus.getMassA(),
// 			  Nucleus.getMassA_min_proton()+Nucleus.getExcitation()[shell],
// 			  MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
//           cout << kin.IsPhysical() << " " << shell << " " << kin.GetCosthklab() << " " 
//     << kin.GetCosthYlab() << " " << kin.GetPklab() << " " << kin.GetPYlab() 
//     << " " << kin.GetKlab() << " " << kin.GetWlab() << endl;
// 
//     }
//   }
//   exit(1);
  
  

//   double total[5];
//   for(int i=0;i<5;i++) total[i]=0.;
//   
//   for(int shell=0;shell<Nucleus.getPLevels();shell++){
//     double low=-1.,high=1.;
//     getBound(high,low,Nucleus,Q2,omega,shell);
//     double results[5];
//     double pestimate=0.;
//     rombergerN(intPm,-1.,low,5, results,1.E-03,3,10,&pestimate,&obs,elec,&Nucleus,Q2,omega,shell, 
// 	       thick,current);
//     for(int j=0;j<5;j++) total[j]+=results[j];
//     
//     cout << Q2/1.E06 << " " << Ein << " " << shell << " " << results[0]/results[4] << " " << 
// 	    results[1]/results[4] << " " << results[2]/results[4] << " " << results[3]/results[4] << endl;
//   }
//   cout << endl;
//     cout << Q2/1.E06 << " " << Ein << " " << total[0]/total[4] << " " << 
// 	    total[1]/total[4] << " " << total[2]/total[4] << " " << total[3]/total[4] << endl;
//   exit(1);
	
  //all levels (glauber, pw) + total proton + total neutron    
  vector<double> totalcross(2*Nucleus.getTotalLevels()+4,0.); 
  double cthmax[Nucleus.getTotalLevels()];
  double cthmin[Nucleus.getTotalLevels()];
  double max=-1.;
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
    cthmax[shell]=1.; cthmin[shell]=-1.;
    getBound(cthmax[shell],cthmin[shell],Nucleus,Q2,omega,shell);
//     cout << cthmax[shell] << endl;
    if(max<cthmax[shell]) max=cthmax[shell];
  }
//   cout << max << endl;
  struct Ftor {

    static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
      Ftor &p = * (Ftor *) param;
      p.f(ret,x[0],p.pObs,p.pNucleus,p.elec,p.Q2,p.omega,p.current,p.cthmax);
    }
    Cross *pObs;
    MeanFieldNucleusThick *pNucleus;
    TElectronKinematics *elec;
    double Q2;
    double omega;
    int current;
    double *cthmax;
    void (*f)(numint::vector_d &, double pm,Cross* pObs, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
	  double Q2, double omega, int current, double *cthmax);

  };

  Ftor F;
  F.pObs = &obs;
  F.pNucleus = &Nucleus;
  F.elec=elec;
  F.Q2=Q2;
  F.omega=omega;
  F.current=current;
  F.cthmax=cthmax;

  numint::mdfunction<numint::vector_d,1> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

  unsigned neval = 0;
  numint::array<double,1> lower = {{-1.}};
  numint::array<double,1> upper = {{max}};
  
  
  F.f=adap_intPm;
  unsigned count=0;
 //numint::cube_romb(mdf,lower,upper,1.E-20,1.E-03,totalcross,count,0);
   numint::cube_adaptive(mdf,lower,upper,1.E-20,1.E-02,1.E02,1.E04,totalcross,count,0);
  
  cout << Q2/1.E06 << " " << omega << " " << q  << " ";
  //proton shells transparencies
  for(int shell=0;shell<Nucleus.getPLevels();shell++) {
    cout << totalcross[2*shell]/totalcross[2*shell+1] << " ";
  }
  cout << totalcross[2*Nucleus.getTotalLevels()] << " " << totalcross[2*Nucleus.getTotalLevels()+1] << endl;

  //neutron shells transparencies
  cout << Q2/1.E06 << " " << omega << " " << q  << " ";
  for(int shell=Nucleus.getPLevels();shell<Nucleus.getTotalLevels();shell++) {
    cout << totalcross[2*shell]/totalcross[2*shell+1] << " ";
  }
  cout << totalcross[2*Nucleus.getTotalLevels()+2] << " " << totalcross[2*Nucleus.getTotalLevels()+3] << endl;

  //total proton and neutron transparency
  cout << totalcross[2*Nucleus.getTotalLevels()]/totalcross[2*Nucleus.getTotalLevels()+1] << " " << totalcross[2*Nucleus.getTotalLevels()+2]/totalcross[2*Nucleus.getTotalLevels()+3] << endl;

  //neutron nominator and denominator cross sections
  for(int shell=Nucleus.getPLevels();shell<Nucleus.getTotalLevels();shell++) {
    cout << totalcross[2*shell] << " ";
  }
  cout << endl;
  for(int shell=Nucleus.getPLevels();shell<Nucleus.getTotalLevels();shell++) {
    cout << totalcross[2*shell+1] << " ";
  }
  cout << endl;
  
  delete elec;
  return 0;
}

void adap_intPm(numint::vector_d & results, double costhetacm,Cross* pObs, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
		double Q2, double omega, int current, double *cthmax){
		  
  results=numint::vector_d(2*pNucleus->getTotalLevels()+4,0.);
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
      numint::vector_d cross=numint::vector_d(5,0.);
      pObs->getAllDiffCross(cross,kin,current,shell,1,0.,20000,0,1,1.,1.);
      // double pw=pObs->getDiffCross(kin,current,1,0,0,1,shell,0.,20000,0,1);
      results[2*shell]+=cross[0];
      results[2*shell+1]+=cross[4];
      results[2*pNucleus->getTotalLevels()]+=cross[0];
      results[2*pNucleus->getTotalLevels()+1]+=cross[4];
      // cout << 2*shell << " " << 2*shell+1 << " " << 2*pNucleus->getTotalLevels() << " " << 2*pNucleus->getTotalLevels()+1 << endl;
      cout << "0 " << shell << " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
      << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
      << " " << kin.GetKlab() << " " << kin.GetWlab() <<  " " << results[2*shell] << " " << results[2*shell+1] << endl;
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
      numint::vector_d cross=numint::vector_d(5,0.);
      pObs->getAllDiffCross(cross,kin,current,shell,1,0.,20000,0,1,1.,1.);
      // double pw=pObs->getDiffCross(kin,current,1,0,0,1,shell,0.,20000,0,1);
      //for(int i=5;i<10;++i) results[i]+=cross[i-5];
      results[2*shell]+=cross[0];
      results[2*shell+1]+=cross[4];
      results[2*pNucleus->getTotalLevels()+2]+=cross[0];
      results[2*pNucleus->getTotalLevels()+3]+=cross[4];
      // cout << 2*shell << " " << 2*shell+1 << " " << 2*pNucleus->getTotalLevels()+2 << " " << 2*pNucleus->getTotalLevels()+3 << endl;
      cout << "1 " << shell << " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
      << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
      << " " << kin.GetKlab() << " " << kin.GetWlab() << " "  << results[2*shell] << " " << results[2*shell+1] << endl;
    }
    
  }
}


void intPm(const double costhcm, double *results, va_list ap){
  Cross* p_obs = va_arg(ap,Cross*);
  TElectronKinematics* elec = va_arg(ap,TElectronKinematics*);
  MeanFieldNucleusThick* pnucleus = va_arg(ap,MeanFieldNucleusThick*);
  double Q2 = va_arg(ap,double);
  double omega = va_arg(ap,double);
  int shell = va_arg(ap,int);
  int thick = va_arg(ap,int);
  int current = va_arg(ap,int);
//   if(pm<1.e-02){
//     for(int i=0;i<5;i++) results[i]=0.;
//     return;
//   }
  
  TKinematics2to2 kin("","",pnucleus->getMassA(),
		      pnucleus->getMassA_min_proton()+pnucleus->getExcitation()[shell],
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,costhcm);
  double pm=kin.GetPklab();
  if(!kin.IsPhysical()){
    for(int i=0;i<5;i++) results[i]=0.;
    return;
  }
  // kin.Print();
    numint::vector_d cross=numint::vector_d(5,0.);
    p_obs->getAllDiffCross(cross,kin,current,shell,1,0.,2000000,0,1,1.,1.);
    for(int i=0;i<5;++i) results[i]=cross[i];
  cout << pm << " " << kin.IsPhysical() << " " << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetCosthYlab() << " " << acos(kin.GetCosthklab())*RADTODEGR << " " << kin.GetCosthklab() <<" " << results[0] << " " << results[1] << " " << results[4] << endl;
  return;
}

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

