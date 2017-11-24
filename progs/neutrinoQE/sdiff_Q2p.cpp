// Program that computes single differential cross sections for QE neutrino nucleus scattering
// Differential to Q2p, computed from final proton kinematics
// Can be used for different experiments, Miniboone and minerva kinematics currently implemented, T2K will be added as well

//run like
// > sdiff_Q2p [Q2p in GeV^2] [integration algorithm for flux integration (see around line 400 for options)] 
//  [max initial nucl momentun in MeV] [min final nucl momentum in MeV] [plane wave 1 or not 0] [ROMEA 1 or not 0 in fsi]
// [experiment "miniboone", "minerva", "t2k"] [lepton id "electron", "muon", "tau"] [sharedir]

// example for plane wave calculation for miniboone
// > sdiff_Q2 0.1 3 500. 200. 1 0 miniboone ~/physics-code/share/



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
    p.f(ret,x[0],x[1],x[2],*p.pNucleus,p.lepton_id,p.current,p.En_out,
        p.prec,p.integrator,p.homedir,p.maxEval,p.charged,p.screening,p.enable_romea,
        p.max_initial_nucl_mom,p.min_final_nucl_mom,p.pw, p.exp
       );
  }
  MeanFieldNucleusThick *pNucleus;
  string lepton_id;
  int current;
  double En_out;
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
  void (*f)(numint::vector_d &, double E_omega, double costheta_mu, double E_in,
            MeanFieldNucleusThick &pNucleus, string lepton_id, int current, double En_out,
            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening, bool enable_romea,
            double max_initial_nucl_mom, double min_final_nucl_mom, int pw, string exp
           );

};

struct FtorH {  //Hydrogen

  static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
    FtorH &p = * (FtorH *) param;
    p.f(ret,x[0],p.En_out,p.charged,p.maxbeam,p.exp,p.lepton_id);
  }
  double En_out;
  bool charged;
  double maxbeam;
  string exp;
  string lepton_id;
  void (*f)(numint::vector_d &, double E_out, double En_out, bool charged, double maxbeam, string exp, string lepton_id);

}; 



//integration for Hydrogen
void int_hydr(numint::vector_d &, double E_out,double En_out,bool charged, double maxbeam, string exp, string lepton_id);

//integration for nucleus
void adap_intPm(numint::vector_d &, double pm, double costhetamu, double E_in,
		            MeanFieldNucleusThick &pNucleus, string lepton_id, int current, double En_out,
		            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening, bool enable_romea,
                            double max_initial_nucl_mom, double min_final_nucl_mom, int pw, string exp
               );





int main(int argc, char *argv[])
{
  double Q2p=atof(argv[1])*1.E06;   // give input in GeV^2
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
  double massAmin1=(Nucleus.getMassA_min_proton()+Nucleus.getMassA_min_neutron())/2.;
  
  vector<double> avgcross(2,0.); //neutrino and antineutrino
  
  double minbeam=0., maxbeam=0.; //different beam energy intervals for the different experiments
  double neutnorm=1.,aneutnorm=1.; // flux normalisations

  if(!exp.compare("miniboone")) { maxbeam=3.E03;}
  else if(!exp.compare("minerva")) { 
    minbeam=1.5E03; 
    maxbeam=10.E03;
    if(!lepton_id.compare("muon")){
      neutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_nu_muon_FHC_flux,500,round(1500./500.),round(10.E03/500));
      aneutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_anu_muon_RHC_flux,500,round(1500./500.),round(10.E03/500));
    }
    else if(!lepton_id.compare("electron")){
      neutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_nu_elec_FHC_flux,500,round(1500./500.),round(10.E03/500));
      aneutnorm=neutrino_flux::normalize(neutrino_flux::Minerva_anu_elec_FHC_flux,500,round(1500./500.),round(10.E03/500));      
    }
  }  
  else if(!exp.compare("t2k")) { maxbeam=2.E02; }  //still to implement!
  else {cerr << "invalid experiment name chosen" << endl << "Choose either miniboone, minerva or t2k" << endl; assert(1==0);}
  
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


  double En_out=(Q2p+pow(MASSn-34.,2.)+MASSn*MASSn)/2./(MASSn-34.); //outgoing nucleon momentum from Q2p definition
  double pn_out=sqrt(En_out*En_out-MASSn*MASSn); //outgoing nucleon momentum

  //setting up variables dat will determine integration limits
  double E_out_max=minbeam;
  double E_out_min=maxbeam;
  double E_low=maxbeam;
  double E_high=minbeam;
  double E_nuH_low=maxbeam;
  double E_nuH_hi=minbeam;
  double omega_low=maxbeam;
  double omega_hi=0.;
  double max=-1.,min=1.; //overall min and max cm theta angles
  double mumax=-1.,mumin=1;

  // for(int j=0;j<=100;j++){ //loop over T_mu
  //   //   cout << j << "/100" << endl;
  //   double E_out=(maxbeam-leptonmass-1.E-03)/100.*j+leptonmass+1.E-03;
  //   double p_out=sqrt(E_out*E_out-leptonmass*leptonmass);
  //   double E_in_lo=E_out+En_out-Nucleus.getMassA()+massAmin1;
  //   double E_in_hi=E_out+En_out-Nucleus.getMassA()+sqrt(massAmin1*massAmin1+max_initial_nucl_mom*max_initial_nucl_mom);

  //   //hydrogen piece
  //   double E_nu_H_n = E_out+En_out-MASSN;
  //   double E_nu_H_p = E_out+En_out-MASSP;
  //   if (E_nu_H_p>E_out&&E_nu_H_p<maxbeam){
  //     if(E_out>E_muH_hi) E_muH_hi = E_out;
  //     if(E_out<E_muH_low) E_muH_low = E_out;
  //   }
  //   if (E_nu_H_n>E_out&&E_nu_H_n<maxbeam){
  //     if(E_out>E_muH_hi) E_muH_hi = E_out;
  //     if(E_out<E_muH_low) E_muH_low = E_out;
  //   }
    
  //   if(E_in_lo<maxbeam){
  //       if(E_in_hi>maxbeam) E_in_hi=maxbeam;
  //       for(int i=0;i<=100;i++){
  //           double E_in=E_in_lo+(maxbeam-E_in_lo)*1.E-02*i;
  //           double omega=E_in-E_out;
  //           if(omega<0||omega<En_out-Nucleus.getMassA()-massAmin1) cout << "omega err " << omega << endl;
  //           if(E_in_lo<E_low)  E_low=E_in_lo;
  //           if(E_in_hi>E_high) E_high=E_in_hi;
  //           if(E_out_max<E_out) E_out_max=E_out;
  //           if(E_out_min>E_out) E_out_min=E_out;
  //           if(omega>omega_hi) omega_hi=omega;
  //           if(omega<omega_low) omega_low=omega;
            
  //       }
  //   }
  // }

  for(int i=0;i<=100;i++){
    double E_in=minbeam+(maxbeam-minbeam)*i*0.01;
    for(int j=0;j<=100;j++){
      double pm=500.*j*0.01;
      for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
        double resmass=(shell<Nucleus.getPLevels()? Nucleus.getMassA_min_proton(): Nucleus.getMassA_min_neutron())
          +Nucleus.getExcitation()[shell]; //Amin1 mass + excitation energy
      
        double EAmin1=sqrt(resmass*resmass+pm*pm);
        double El_out=E_in+Nucleus.getMassA()-En_out-EAmin1;
        double pl_out=sqrt(El_out*El_out-leptonmass*leptonmass);
        double omega=E_in-El_out;
        if(El_out>leptonmass){
          int phys=0;
          int k=0;
          do{
            double costhetamu=1.-k*0.001;
            double Q2=-leptonmass*leptonmass+2.*E_in*El_out-2.*E_in*pl_out*costhetamu;
            TKinematics2to2 kin("","",Nucleus.getMassA(),Nucleus.getMassA_min_neutron(),
              MASSN,"qsquared:wlab:pklab",Q2,omega,pm);
            // if(kin.IsPhysical()) cout << E_in << " " << pm << " " << costhetamu << " " << kin.GetCosthkcm() << " " << Q2 << " " << Q2p << " " << omega << endl;
            if(kin.IsPhysical()&&phys==0) phys=1;
            if(phys==1&&kin.IsPhysical()){
              if(E_in<E_low) E_low=E_in;
              if(E_in>E_high) E_high=E_in;
              if(El_out<E_out_min) E_out_min=El_out;
              if(El_out>E_out_max) E_out_max=El_out;
              if(omega>omega_hi) omega_hi=omega;
              if(omega<omega_low) omega_low=omega;
              if(kin.GetCosthkcm()<min) min=kin.GetCosthkcm();
              if(kin.GetCosthkcm()>max) max=kin.GetCosthkcm();
              if(costhetamu<mumin) mumin=costhetamu;
              if(costhetamu>mumax) mumax=costhetamu;
            }
            if(phys==1&&(!kin.IsPhysical())) phys=2;
            // cout << i << " " <<j << " " << k << " " << phys << " " << costhetamu << " " << kin.GetCosthkcm() << " " << kin.IsPhysical() << endl;
            k++;       
          }while(phys<2&&k<=2000);
        }
      }
    }    
  }

  if(E_high>maxbeam) E_high=maxbeam;
  if(E_low<minbeam) E_low=minbeam;
  if(E_low>E_high) E_low=E_high;
  if(E_nuH_low<minbeam) E_nuH_low=minbeam;
  if(mumax>0.99) mumax=1.;
  if(mumin<-0.99) mumin=-1.;

  cout << endl;
  cout << "min=" << min << "   max=" << max << endl;
  cout << "mumin=" << mumin << "   mumax=" << mumax << endl;
  cout << "Eoutmin=" << E_out_min << "  Eoutmax=" << E_out_max << endl;
  cout << "omega_low=" << omega_low << " " << "  omega_high=" << omega_hi << endl;
  cout << "E_low=" << E_low << "  E_high=" << E_high << endl << endl;
  cout << "E_nuH_low=" << E_nuH_low << "  E_muH_hi=" << E_nuH_hi << endl << endl;
 

  FtorH FH;
  FH.En_out=En_out;
  FH.charged=charged;
  FH.exp=exp;
  FH.maxbeam=maxbeam;
  FH.lepton_id=lepton_id;

  numint::mdfunction<numint::vector_d,1> mdfH;
  mdfH.func = &FtorH::exec;
  mdfH.param = &FH;

  numint::array<double,1> lowerH = {{minbeam}};  
  numint::array<double,1> upperH = {{maxbeam}};
  
  FH.f=int_hydr;
  vector<double> avgcrossH(2,0.); 
  unsigned count=0;
  if(!fluxintegrator) numint::cube_romb(mdfH,lowerH,upperH,1.E-30,1.E-03,avgcrossH,count,0); //1.E-20,1.E-03
  else numint::cube_adaptive(mdfH,lowerH,upperH,1.E-30,1.E-03,2E02,4E04,avgcrossH,count,0); 
  cout << Q2p*1.E-6 << " " << avgcrossH[0]*1.E19*2.*PI/neutnorm << " " << avgcrossH[1]*1.E19*2.*PI/aneutnorm << endl << endl;


  //initialize object -- Carbon
  Ftor F;
  F.pNucleus = &Nucleus;
  F.lepton_id=lepton_id;
  F.current=current;
  F.En_out=En_out;
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
  
  //we integrate over residual nucleus momentum, lepton costheta, incoming beam energy
  numint::array<double,3> lower = {{0.,E_low,mumin}};
  numint::array<double,3> upper = {{max_initial_nucl_mom,E_high,mumax}};
  
  F.f=adap_intPm;
  count=0;
  string stf=homedir+"/statefile/sdiff_q2p";
  ostringstream qstr;
  qstr << Q2p*1.E-06 << "." << pw << "." <<enable_romea << "." <<max_initial_nucl_mom << "." <<min_final_nucl_mom;
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
  cout << "start integration C" << endl;
  if(fluxintegrator==0) numint::vegas( mdf, lower,upper,nvec, epsrel,epsabs,flags,seed,minEval, maxEvalcuba,1000, 500, 1000, 0,(stf+"vegas").c_str(),countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==1) numint::cuhre( mdf, lower,upper,nvec, epsrel,epsabs,flags,minEval, maxEvalcuba,11, (stf+"cuhre").c_str(),nregions,countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==2) numint::cube_romb(mdf,lower,upper,1.E-25,1.E-01,avgcross,count,0); //1.E-20,1.E-03
  if(fluxintegrator==3) numint::cube_adaptive(mdf,lower,upper,1.E-25,1.E-03,2E02,2E04,avgcross,count,0);
  if(fluxintegrator==4) numint::divonne(mdf, lower,upper,nvec, epsrel, epsabs,flags,seed,
           minEval, maxEvalcuba,7,7,0,5,0,10,0.25,0,3,NULL, 0, 0,
           (stf+"divonne").c_str(),nregions, countt, fail, avgcross,err,prob );
  
  //factor 2\pi because of integration over muon polar angle
  cout << Q2p*1.E-6 << " " << avgcrossH[0]*1.E19*2.*PI/neutnorm << " " << avgcrossH[1]*1.E19*2.*PI/aneutnorm << " " 
     << avgcross[0]*1.E19*2.*PI/Nucleus.getN()/neutnorm << " " << avgcross[1]*1.E19*2.*PI/Nucleus.getZ()/aneutnorm << " " <<  count << endl;



}

//integrandum carbon
void adap_intPm(numint::vector_d & results, double pm, double E_in, double costhetamu, 
	 MeanFieldNucleusThick &nucleus, string lepton_id, int current, double En_out,
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
    

  for(int shell=0;shell<nucleus.getTotalLevels();shell++) {
    double resmass=(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
      +nucleus.getExcitation()[shell]; //Amin1 mass + excitation energy
    double EAmin1=sqrt(resmass*resmass+pm*pm);
    double omega = En_out + EAmin1 -nucleus.getMassA();
    double E_out = E_in-omega;
    double p_out=sqrt(E_out*E_out-leptonmass*leptonmass);
    double Q2=-leptonmass*leptonmass+2.*E_in*E_out-2.*E_in*p_out*costhetamu;
    if(/*E_in>maxbeam||E_in<minbeam*/E_out<0.||Q2<0.) { return;}
      
    TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(static_cast<TLeptonKinematics::Lepton>(flav),E_in);
    WeakQECross pobs(lepton,&nucleus,prec,integrator,homedir,charged,1.03E03,screening,enable_romea);  
      
    TKinematics2to2 kin("","",nucleus.getMassA(),resmass,
      shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:pklab",Q2,omega,pm);
    if(!kin.IsPhysical()||kin.GetPYlab()<min_final_nucl_mom){ //final nucleon momentum too low, impose cut!
      // cout << kin.IsPhysical() << " " << E_out << " " << omega << " " << Q2*1.E-06 << " " << costhetamu << " " << kin.GetCosthkcm() <<  " " << kin.GetEYlab() << " " << En_out << endl;
      for(int i=0;i<2;i++) results[i]+=0.; /*cout << "cut " << kin.GetPYlab() << " " << kin.GetPklab() << endl;*/
      // cout << shell << " " << E_in << " " << costhetamu << " " << p_out << " " << kin.GetCosthklab()<< " " 
      //   << kin.GetCosthYlab() << " " << kin.GetPklab() << " " << kin.GetPYlab() 
      //   << " " << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetQsquared() << " "
      //   << kin.GetXb()*nucleus.getMassA()/MASSP << " 0. 0." << endl; 
    }
    else{
      double result=pobs.getDiffWeakQECross(kin,current,0,0,0,pw,shell,0.,maxEval,WeakQECross::En_lab,1);   // prec..2E04
//	  			cout << "Result " << result << endl;

      //Jacobian to go to dQ2p, pm/E_a-1 is from going from dEleptonout (or domega) to dpm
      double jcb=2*MASSn*EAmin1/pm;
      result/=jcb;      

      results[(shell<nucleus.getPLevels()?1:0)]+= result; //results[0] neutrino, results[1] antineutrino
        cout << shell << " " << E_in <<  " " << costhetamu << " " << p_out << " "  << kin.GetCosthklab()<< " " 
        << kin.GetCosthYlab() << " " << kin.GetPklab() << " " << kin.GetPYlab() 
        << " " << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetQsquared() << " "
        << kin.GetXb()*nucleus.getMassA()/MASSP << " " << result <<  " " 
        << result*interpolate(shell<2? neutrino_flux::Minerva_anu_muon_RHC_flux : neutrino_flux::Minerva_nu_muon_FHC_flux,E_in,500,20,0) << endl; 
    }  
    delete lepton;
  }
  //fold with flux
  if(!exp.compare("miniboone")&&!lepton_id.compare("muon")){
    results[0]*=interpolate(neutrino_flux::MiniBooNE_neut_muon_flux_norm,E_in,25,120,1);
    results[1]*=interpolate(neutrino_flux::MiniBooNE_antineut_muon_flux_norm,E_in,25,120,1);
  }
  else if(!exp.compare("minerva")&&!lepton_id.compare("muon")){
    results[0]*=interpolate(neutrino_flux::Minerva_nu_muon_FHC_flux,E_in,500,20,0);
    results[1]*=interpolate(neutrino_flux::Minerva_anu_muon_RHC_flux,E_in,500,20,0);
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
  return;


}


//integrandum hydrogen
void int_hydr(numint::vector_d & results, double E_in, 
              double En_out, bool charged, double maxbeam, string exp, string lepton_id){
    
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
  double E_out_n=E_in-En_out+MASSN;
  double E_out_p=E_in-En_out+MASSP;
  double p_out_n = sqrt(E_out_n*E_out_n-leptonmass*leptonmass);
  double p_out_p = sqrt(E_out_p*E_out_p-leptonmass*leptonmass);
  double pn_out_n = sqrt(En_out*En_out-MASSP*MASSP);  //initial neutron case
  double pn_out_p = sqrt(En_out*En_out-MASSN*MASSN);  //initial proton case  
  double costhetamu_n = (E_in*E_in+p_out_n*p_out_n-pn_out_n*pn_out_n)/(2.*p_out_n*E_in);
  double costhetamu_p = (E_in*E_in+p_out_p*p_out_p-pn_out_p*pn_out_p)/(2.*p_out_p*E_in);

  //actually this is not needed as in this case the Q2n, Q2p values will be the same for every incoming energy, completely determined by the nucleon kinematics...
  //we leave it here as it's not the bottleneck...
  double Q2_n=-leptonmass*leptonmass+2.*E_in*(E_out_n-p_out_n*costhetamu_n);
  double Q2_p=-leptonmass*leptonmass+2.*E_in*(E_out_p-p_out_p*costhetamu_p);

  //[0] is neutron [1] is proton
  results[0]=WeakQECross::getElWeakQECross(Q2_n,E_in,0,charged,1.03E03,leptonmass*leptonmass,WeakQECross::type_Q2p);
  results[1]=WeakQECross::getElWeakQECross(Q2_p,E_in,1,charged,1.03E03,leptonmass*leptonmass,WeakQECross::type_Q2p);
  
  if(!exp.compare("miniboone")&&!lepton_id.compare("muon")){
    results[0]*=E_in>maxbeam? 0. : interpolate(neutrino_flux::MiniBooNE_neut_muon_flux_norm,E_in,25,120,1);
    results[1]*=E_in>maxbeam? 0. : interpolate(neutrino_flux::MiniBooNE_antineut_muon_flux_norm,E_in,25,120,1);
  }
  else if(!exp.compare("minerva")&&!lepton_id.compare("muon")){
    results[0]*=E_in>maxbeam? 0. : interpolate(neutrino_flux::Minerva_nu_muon_FHC_flux,E_in,500,20,0);
    results[1]*=E_in>maxbeam? 0. : interpolate(neutrino_flux::Minerva_anu_muon_RHC_flux,E_in,500,20,0);
  }
   else if(!exp.compare("minerva")&&!lepton_id.compare("electron")){
      results[0]*=E_in>maxbeam? 0. : interpolate(neutrino_flux::Minerva_nu_elec_FHC_flux,E_in,500,20,0);
      results[1]*=E_in>maxbeam? 0. : interpolate(neutrino_flux::Minerva_anu_elec_FHC_flux,E_in,500,20,0);
    }
  else if(!exp.compare("t2k")&&!lepton_id.compare("muon")){
    results[0]*=E_in>maxbeam? 0. : interpolate(neutrino_flux::t2k_neut_muon_flux_norm,E_in,25,80,1);
    results[1]*=E_in>maxbeam? 0. : interpolate(neutrino_flux::t2k_aneut_muon_flux_norm,E_in,25,80,1);
  }
  else {cerr << "Unsupported combination of experiment and lepton_id:" << exp << " " << lepton_id << endl; assert(1==0);}

//   if(isnan(result[1])) {cout << E_out << " " << E_in << " " << crossHanu_n << " " << interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1) << " " << Q2 << endl; exit(1);}
//    cout << E_out << " " << result[0] << " " << result[1] << " " << E_in << " " << E_in << " " <<  interpolate(MiniBooNE_antineut_flux_norm,E_in,25,120,1) << " " << interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1) <<endl;
} 
