//program that calculates density profiles for Calcium, per request from Eli and Jechiel

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
#include <Cross.hpp>
#include <Utilfunctions.hpp>
#include <NucleonEMOperator.hpp>

//run ./observables [Q2 [MeV^2]] [omega] [missing momentum]
int main(int argc, char *argv[])
{

  string arg_names[argc]={"exec name", 
                            "nucleus", 
                            "nucleus shell",
                            "scattered electron angle [degrees]",
                            "scattered electron energy [MeV]",
                            "phi angle between electron and hadron plane [degrees]",
                            "medium modifications ([0] no, [1] CQM, [2] QSM):"};

  std::cout << "Compiled from file: " << __FILE__ << std::endl;
  Bookkeep(argc,argv,arg_names);  

  double Ein=600.;
  double Eout=atof(argv[4]);
  int shell=atoi(argv[2]);
  double thetae=atof(argv[3])*DEGRTORAD;
  double phi = atof(argv[5])*DEGRTORAD;
  int mmod = atoi(argv[6]);
  //double thetap=atof(argv[3])*DEGRTORAD;
  double omega=Ein-Eout;
//   int medium=atoi(argv[4]);
  
  bool screening=0;//atoi(argv[4]);
  double scr=1.;//atof(argv[5]);
  //string nucleus_name="Ca40";//atoi(argv[6]);
  string nucleus_name=argv[1];//atoi(argv[6]);
  double prec=1.E-05;//atof(argv[7]);
  int integrator=2;//atoi(argv[8]);
  int thick=1;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  
  string homedir=HOMEDIR;
  MeanFieldNucleusThick nucleus(MeanFieldNucleus::TypeNames.at(nucleus_name),homedir);
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
  
  cout << "L of shell: " << nucleus.getL_array()[shell] << endl;

  double qvec=sqrt(Ein*Ein+Eout*Eout-2.*Ein*Eout*cos(thetae));
  double Q2=qvec*qvec-omega*omega;
  double thetaq=atan2(-sin(thetae)*Eout,Ein-Eout*cos(thetae));

  cout << "thetaq: " << thetaq*RADTODEGR << endl; 
  cout << "Q2: " << Q2*1.E-06 << endl;

  Cross obs(*elec,&nucleus,prec,integrator,homedir,screening,scr);

  cout << "q|omega|nucl_p|nucl_thetaq|pm|pm_thetaq|";
  cout << "Px/Pz PWIA" << " " << "Px/Pz RMSGA" << " " << "Px/Pz free" << endl;
  //cout << "xs_fsisrc|xs_pw|Tratio|rhod|rhopw" << endl;
  //cout << "pmiss <rho>_fsi <rho>_pw <r>_fsi <r>_pw rhod_fsi rhod_pw" << endl;
  for(int i=1;i<=30;i+=2){
  double pm=i*10.;
  
  
    TKinematics2to2 kin("","",nucleus.getMassA(),nucleus.getMassA_min_proton()+nucleus.getExcitation()[shell]
		,MASSP,"qsquared:wlab:pklab",Q2,omega,pm);
    // cout << i << endl;
    if(kin.IsPhysical()){
        cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*RADTODEGR 
      << " " << kin.GetPklab() << " " << acos(kin.GetCosthklab())*RADTODEGR << " ";

    numint::vector_d cross=numint::vector_d(5,0.);
    obs.getAllDiffCross(cross,kin,2,shell,1,phi,20000,1,0, 1.,1.);

    numint::vector_d asymm=numint::vector_d(40,0.);
    obs.getAllObs_xyz(asymm, kin, 2, shell, 1, mmod, phi, 20000, 1,1.,1.);
    NucleonEMOperator nuclFF(Q2,1,0);
    double free_ratio = -2.*MASSP/(Ein+Eout)/tan(thetae/2.)*nuclFF.getGE()/nuclFF.getGM();
    
    cout << asymm[37]/asymm[39] << " " << asymm[5]/asymm[7] << " " << free_ratio << endl;
    // obs.getAllObs_tnl(asymm, kin, 2, shell, 1, 0, phi, 20000, 1,1.,1.);
    // cout << asymm[37]/asymm[39] << " " << asymm[5]/asymm[7] << " " << free_ratio << endl;


  //   double rhod=0., rhopw=0.;
  //   for(int ms=-1;ms<=1;ms+=2){
  //     for(int mj = -nucleus.getJ_array()[shell]; mj <=nucleus.getJ_array()[shell]; mj+=2){
  //       // cout << ms << " " << mj << " ";
  //       numint::vector_z phid=numint::vector_z(5,0.);
  //       obs.getPhid(phid,kin,shell,mj,ms,1,maxEval);
  //       rhod+=norm(phid[1]);
  //       rhopw+=norm(phid[4]);
  //     }
  //   }

  //   cout << cross[1] << " " << cross[4] << " " << cross[1]/cross[4] << " " << rhod << " " << rhopw << endl;

  //      obs.printDensity_profile(kin,shell,thick,maxEval);
   }
 }    
  delete elec;
  return 0;




}


