//program to compute d\sigma/dQ^2 for Minerva kinematics, quasi-elastic contribution
//flux averaged
//proton contribution included for antineutrino cross section

//We compute quasi elastic d\sigma^5 and integrate over com scattering angle, energy transfer and muon kinetic energy.
//biggest trouble was finding good integration variables so the integration algorithm didn't have to evaluate a lot of zero volume...
//(phase space has quite some restrictions/cuts)


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
#include <numint/numint2Cuba.hpp>

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

  static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
    Ftor &p = * (Ftor *) param;
    p.f(ret,x[0],x[1],x[2],*p.pNucleus,p.current,p.cthmax,p.Q2,
        p.prec,p.integrator,p.homedir,p.maxEval,p.charged,p.screening,
        p.max_initial_nucl_mom,p.min_final_nucl_mom,p.pw
       );
  }
  MeanFieldNucleusThick *pNucleus;
  int current;
  double *cthmax;
  double Q2;
  double prec;
  int integrator;
  string homedir;
  int maxEval;
  bool charged;
  bool screening;
  double max_initial_nucl_mom;
  double min_final_nucl_mom;
  int pw;
  void (*f)(numint::vector_d &, double E_in, double costhetacm, double E_out, 
            MeanFieldNucleusThick &pNucleus, int current, double *cthmax, double Q2,
            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening,
            double max_initial_nucl_mom, double min_final_nucl_mom, int pw
           );

};

struct FtorH {  //Hydrogen

  static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
    FtorH &p = * (FtorH *) param;
    p.f(ret,x[0],p.current,p.Q2,p.charged);
  }
  int current;
  double Q2;
  bool charged;
  void (*f)(numint::vector_d &, double E_in, int current, double Q2,
             bool charged);

}; 

//determine boundaries of kinematics
double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell, double max_initial_nucl_mom, double min_final_nucl_mom);
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell,double max_initial_nucl_mom);
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell,double max_initial_nucl_mom);

// integration over missing momentum + over T_mu(=E_out+C)
void adap_intPm(numint::vector_d &, double E_in, double costhetacm, double E_out,
		            MeanFieldNucleusThick &pNucleus,int current, double *cthmax, double Q2,
		            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening,
                double max_initial_nucl_mom, double min_final_nucl_mom, int pw);

//integration for Hydrogen
void int_hydr(numint::vector_d &, double E_in,
              int current, double Q2,bool charged);

void normalize(double flux[], double dx){
  double sum=0;
  for(int i=0;i<20;i++) sum+=flux[i]*dx;
  for(int i=0;i<20;i++) flux[i]/=sum;
  };


int main(int argc, char *argv[])
{
  double Q2=atof(argv[1])*1.E06;   // give input in GeV^2
//double T_mu=atof(argv[1]); // muon kinetic energy in MeV
  //double costhetamu; //=atof(argv[1]);
  bool screening=0;//atoi(argv[4]);
  int nucleus=1;//atoi(argv[6]);                     
  double prec=1.E-05;//atof(argv[7]);   //1.E-5
  int integrator=2;//
  int fluxintegrator=atoi(argv[2]);
//int thick=0;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  bool charged=1;
  int current=2;
  int pw=atoi(argv[5]); //1 is plane-wave, 0 is with FSI
  double max_initial_nucl_mom=atof(argv[3]);
  double min_final_nucl_mom=atof(argv[4]);
  
  
  string homedir=/*argv[6];*/   "/home/wim/Code/trunk/share";

  MeanFieldNucleusThick Nucleus(nucleus,homedir);
  
  vector<double> avgcross(2,0.); //neutrino and antineutrino
  //estimates for integration bounds
  double cthmax[Nucleus.getTotalLevels()];
  double cthmin[Nucleus.getTotalLevels()];
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cthmax[shell]=-1.;cthmin[shell]=1.;}
  double max=-1.,min=1.; //overall min and max center of mass cosine theta angles
  double Tmax=0;
  double Tmin=10000;
  double E_low=15000;
  double E_high=0;
  double omega_low=1.E04;
  double omega_hi=0.;
  
  //find reasonable integration limits
//   for(int j=0;j<=100;j++){ //loop over T_mu
// //     cout << j << "/100" << endl;
//     double T_mu=(1.E04-massmu)*0.01*j;    
//     double E_out=T_mu+massmu;
//     double p_out=sqrt(E_out*E_out-massmu*massmu);
//     //physical limits so that costhetamu is between -1 and 1
//     double E_in_low = (Q2+massmu*massmu)/2./(E_out+p_out);
//     double E_in_hi = (Q2+massmu*massmu)/2./(E_out-p_out);
//     if(E_in_low>10000.||E_in_hi<1500.){} //outside considered incoming neutrino energies, so don't bother
//     else{
//       for(int i=0;i<=100;i++){
//         if(E_in_low<1500.) E_in_low=1500.;
//         if(E_in_hi>1.E04) E_in_hi=1.E04;
//         double E_in=E_in_low+(E_in_hi-E_in_low)*1.E-02*i;
//         double costhetamu=(-Q2-massmu*massmu+2.*E_in*E_out)/(2.*E_in*sqrt(E_out*E_out-massmu*massmu));
//         double omega=E_in-E_out;
//         TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(TLeptonKinematics::muon,E_in);
//         
//         for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
//           
//           //anything above 500 MeV contribution will be negligible
//           double tempmax=-1., tempmin=1.;
//           if(getBound(tempmax,tempmin,Nucleus,*lepton,E_in,E_out,costhetamu,shell,max_initial_nucl_mom,min_final_nucl_mom)<max_initial_nucl_mom){ 
//             if(E_in<E_low)  E_low=E_in;
//             if(E_in>E_high) E_high=E_in;
//             if(max<tempmax) max=tempmax;
//             if(min>tempmin) min=tempmin;
//             if(cthmax[shell]<tempmax) cthmax[shell]=tempmax;
//             if(cthmin[shell]>tempmin) cthmin[shell]=tempmin;
//             if(Tmax<T_mu) Tmax=T_mu;
//             if(Tmin>T_mu) Tmin=T_mu;
//             if(omega>omega_hi) omega_hi=omega;
//             if(omega<omega_low) omega_low=omega;
//             
// //               cout << E_in << " " << E_out << " " << tempmin << " " << tempmax << " " << pm_min << endl;
//           }
//         }
//       }
//     }  
//   }
//   if(E_low<1500.) E_low=1500.;  
//   if(E_high<E_low) E_high=E_low;
//   
//   //min=-1;   max=1;
//   //Tmin=100;  Tmax=7100;
//   //E_low=1500;  E_high=7413.36;
// //   for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cout << shell << " " << cthmax[shell] << " " << cthmin[shell] << endl;cthmin[shell]=-1.;cthmax[shell]=1.;}
//   cout << endl;
//   cout << "min=" << min << "   max=" << max << endl;
//   cout << "Tmin=" << Tmin << "  Tmax=" << Tmax << endl;
//   cout << "omega_low=" << omega_low << " " << "  omega_high=" << omega_hi << endl;
//   cout << "E_low=" << E_low << "  E_high=" << E_high << endl << endl;
// //   exit(1);

  //Normalize flux
  normalize(Minerva_nu_flux,500);
  normalize(Minerva_anu_flux,500);
 
  //initialize object -- Hydrogen
  
  for(int i=0;i<200;i++){
    Q2=i*0.01*1.E06;
    
    FtorH FH;
    FH.current=current;
    FH.Q2=Q2;
    FH.charged=charged;

    numint::mdfunction<numint::vector_d,1> mdfH;
    mdfH.func = &FtorH::exec;
    mdfH.param = &FH;

    // find minimum E_in
    double omega=0.5*(MASSN*MASSN-MASSP*MASSP+Q2)/MASSP;
    double qvec=sqrt(Q2+omega*omega);
  //   double EminH=0.25*Q2/MASSP+0.5*sqrt(0.25*Q2*Q2/(MASSP*MASSP)+Q2+massmu*massmu);
    double EminH=omega;
    if (EminH<1500.) EminH=1500.;
//     cout << "EminH: " << EminH << endl;
    numint::array<double,1> lowerH = {{EminH}};  
    numint::array<double,1> upperH = {{10000.}};
    
    FH.f=int_hydr;
    vector<double> avgcrossH(1,0.); 
    unsigned count=0;
    if(!fluxintegrator) numint::cube_romb(mdfH,lowerH,upperH,1.E-30,1.E-03,avgcrossH,count,0); //1.E-20,1.E-03
    else numint::cube_adaptive(mdfH,lowerH,upperH,1.E-30,1.E-03,2E02,4E04,avgcrossH,count,0); 
    
  //   cout << "Crosssection H: " << Q2 << " " << avgcrossH[0]*1.E19*2.*PI << " [1E-39 cm2/GeV2]" <<  " "<< count << endl; //2pi because of integration over scattered lepton angle
    cout << Q2 << " " << avgcrossH[0]*1.E19*2.*PI << endl;
  }
  exit(1);
  
  //initialize object -- Carbon
  Ftor F;
  F.pNucleus = &Nucleus;
  F.current=current;
  F.cthmax=cthmax;
  F.Q2=Q2;
  F.prec=prec;
  F.integrator=integrator;
  F.homedir=homedir;
  F.maxEval=maxEval;
  F.charged=charged;
  F.screening=screening;
  F.min_final_nucl_mom=min_final_nucl_mom;
  F.max_initial_nucl_mom=max_initial_nucl_mom;
  F.pw=pw;

  numint::mdfunction<numint::vector_d,3> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

//unsigned neval = 0;
  numint::array<double,3> lower = {{omega_low,min,Tmin+massmu}};
  numint::array<double,3> upper = {{omega_hi,max,Tmax+massmu}};
  
  F.f=adap_intPm;
  unsigned count=0;
  string stf=homedir+"/statefile/minerva";
  ostringstream qstr;
  qstr << Q2;
  string qst=qstr.str();
  stf+=qst+".FSI";
  ostringstream qstr2;
  qstr2 << (!pw);
  string pwst=qstr2.str();
  stf+=pwst;
  int nregions,fail,countt;
  vector<double> err(2,0.);
  vector<double> prob(2,0.);
  int nvec=1;
  double epsrel=1.E-02;
  double epsabs=1.E-25;
  int flags=2;
  int seed=235;
  int minEval=100;
  int maxEvalcuba=50000;
  double xgiven[3]={200.,-1.,2.5E03};
  cout << "start integration C" << endl;
  if(fluxintegrator==0) numint::vegas( mdf, lower,upper,nvec, epsrel,epsabs,flags,seed,minEval, maxEvalcuba,1000, 500, 1000, 0,(stf+"vegas").c_str(),countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==1) numint::cuhre( mdf, lower,upper,nvec, epsrel,epsabs,flags,minEval, maxEvalcuba,11, (stf+"cuhre").c_str(),nregions,countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==2) numint::cube_romb(mdf,lower,upper,1.E-25,1.E-01,avgcross,count,0); //1.E-20,1.E-03
  if(fluxintegrator==3) numint::cube_adaptive(mdf,lower,upper,1.E-25,1.E-02,2E02,2E03,avgcross,count,0);
  if(fluxintegrator==4) numint::divonne(mdf, lower,upper,nvec, epsrel, epsabs,flags,seed,
           minEval, maxEvalcuba,7,7,0,5,0,10,0.25,0,3,NULL, 0, 0,
           (stf+"divonne").c_str(),nregions, countt, fail, avgcross,err,prob );

   
  //factor 2\pi because of integration over muon polar angle
  cout << "Crosssection C: " << avgcross[1]*1.E19*2.*PI/Nucleus.getZ() << " [1E-39 cm2/GeV2] " << count << " " << countt<< endl << endl;
  
//   cout << Q2 << " " << avgcross[0]*1.E19*2.*PI/Nucleus.getN()
//   << " " << (avgcross[1]+2.*avgcrossH[0])*1.E19*2.*PI/(Nucleus.getZ()+2.) << " " 
//   << avgcrossH[0]*1.E19*2.*PI << " " << avgcross[1]*1.E19*2.*PI/Nucleus.getZ() << " "  << count << " " << countt << endl;

}

//integrandum
void adap_intPm(numint::vector_d & results, double omega, double costhetacm, double E_out,
	 MeanFieldNucleusThick &nucleus, int current, double *cthmax, double Q2,
	 double prec, int integrator, string homedir, int maxEval, bool charged, bool screening, double max_initial_nucl_mom, 
                double min_final_nucl_mom, int pw){	  
  results=numint::vector_d(2,0.);  
  double E_in=E_out+omega;
//   double omega=E_in-E_out;  
  double costhetamu=(-Q2-massmu*massmu+2.*E_in*E_out)/(2.*E_in*sqrt(E_out*E_out-massmu*massmu));
  
  if(E_in<1.5E03 || E_in>10.E03)  {return;}
  if(abs(costhetamu)<=1. /*&& omega>0.*/) {
      
    TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(TLeptonKinematics::muon,E_in);
    WeakQECross pobs(lepton,&nucleus,prec,integrator,homedir,charged,1.03E03,screening,0.);  
    
    //we can fix the vector boson kinematics
//     pobs.getPlepton()->SetBeamEnergy(E_in);

    for(int shell=0;shell<nucleus.getTotalLevels();shell++) {
      if(costhetacm<cthmax[shell]){
        TKinematics2to2 kin("","",nucleus.getMassA(),
          (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
          +nucleus.getExcitation()[shell],
          shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,costhetacm);
        double pm=sqrt(E_out*E_out-massmu*massmu);
        if(!kin.IsPhysical()||kin.GetPYlab()<min_final_nucl_mom||kin.GetPklab()>max_initial_nucl_mom){ //final nucleon momentum too low, impose cut!
          for(int i=0;i<2;i++) results[i]+=0.;/*	cout << "cut " << kin.GetPYlab() << " " << kin.GetPklab() << endl;*/
        }
        else{
          double result=pobs.getDiffWeakQECross(kin,current,0,0,0,pw,shell,0.,maxEval,0,1);   // prec..2E04
// 	  			cout << "Result " << result << endl;

          //Jacobian    
          double Enu=(MASSN*MASSN-(MASSP-30.)*(MASSP-30.)-lepton->GetLeptonMass()*lepton->GetLeptonMass()+2.*(MASSP-30.)*E_out)
                    /(2.*(MASSP-30.-E_out+pm*costhetamu));
          double jcb=2.*pm*Enu*(1.+(E_out-pm*costhetamu)/(MASSP-30.-E_out+pm*costhetamu));
          result/=jcb;      
//           cout << E_in << " " << costhetacm << " " << E_out << " " << omega << " " <<  shell << " " << kin.GetPYlab() << " " << kin.GetPklab() << " " << result << endl;
          
          results[(shell<nucleus.getPLevels()?1:0)]+= result; //results[0] neutrino, results[1] antineutrino
//            cout << shell  << " " << E_in <<  " " << costhetacm << " " << pm << " "  << kin.GetCosthklab() << " " 
//           	<< kin.GetCosthYlab() << " " << kin.GetPklab() << " " << kin.GetPYlab() 
//           	<< " " << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetQsquared() << " "
//           	<< kin.GetXb()*nucleus.getMassA()/MASSP << " " << interpolate(Minerva_nu_flux,E_in,500,20,0) << " " << interpolate(Minerva_anu_flux,E_in,500,20,0) << " " << result << endl;
        }
      }  
    }
//     cout << results[0] << " " << results[1] << endl;
    //fold with flux
    results[0]*=interpolate(Minerva_nu_flux,E_in,500,20,0);
    results[1]*=interpolate(Minerva_anu_flux,E_in,500,20,0);
//     cout << "physical " << omega << " " << costhetacm << " " << E_out << " " << E_in << " " << results[0] << " " << results[1] << endl;
    delete lepton;
  }
  else {results[0]=0; results[1]=0; /*cout << "unphys " << omega << " " << costhetacm << " " << E_out << " " << E_in << " " << 0. << " " << 0. << endl;*/}

}

void int_hydr(numint::vector_d & result, double E_in, 
              int current, double Q2, bool charged){
    
  double omega=0.5*(MASSP*MASSP-MASSN*MASSN+Q2)/MASSP;  
  double E_out=E_in-omega;
  
  TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(TLeptonKinematics::muon,E_in);
  double crossHanu=WeakQECross::getElWeakQECross(Q2,E_in,current,0,charged,1.03E03);
  
  //Jacobian from d\costheta_mu to dQ2    
  double jcb=2.*E_in*sqrt(E_out*E_out-lepton->GetLeptonMass()*lepton->GetLeptonMass());
//	cout << "Jacobian H " << jcb << endl;
  crossHanu/=jcb;
 
  result=numint::vector_d(1,0.); 
  result[0]=crossHanu*interpolate(Minerva_anu_flux,E_in,500,20,0);
  delete lepton;
} 

double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell, double max_initial_nucl_mom, double min_final_nucl_mom){
  
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();
  double tempmin=max_initial_nucl_mom;
  double costhmin=-1.;
  for(int i=0;i<=100;i++){ //loop over all com angles
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
  if((high-low)<1.E-04) { return high;}
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
  if((high-low)<1.E-04) { return low;}
  else return getMax(high,low,nucleus,lepton,Q2,omega,shell,max_initial_nucl_mom);  
}
