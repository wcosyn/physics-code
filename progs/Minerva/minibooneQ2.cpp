#include <iostream>
#include <cstdlib>
#include <sstream>

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
struct Ftor {  //Carbon

  static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
    Ftor &p = * (Ftor *) param;
    p.f(ret,x[0],x[1],x[2],*p.pNucleus,p.current,p.cthmax,p.Q2,
        p.prec,p.integrator,p.homedir,p.maxEval,p.charged,p.screening);
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
  void (*f)(numint::vector_d &, double E_in, double costhetacm, double E_out, 
            MeanFieldNucleusThick &pNucleus, int current, double *cthmax, double Q2,
            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening);

};

//determine boundaries of kinematics
double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell);
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell);
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell);

// integration over missing momentum + over T_mu(=E_out+C)
void adap_intPm(numint::vector_d &, double E_in, double costhetacm, double E_out,
		            MeanFieldNucleusThick &pNucleus,int current, double *cthmax, double Q2,
		            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening);


int main(int argc, char *argv[])
{
  double Q2=atof(argv[1])*1.E06;   // give input in GeV^2
//double T_mu=atof(argv[1]); // muon kinetic energy in MeV
  double costhetamu; //=atof(argv[1]);
  bool screening=0;//atoi(argv[4]);
  int nucleus=1;//atoi(argv[6]);                     
  double prec=1.E-05;//atof(argv[7]);   //1.E-5
  int integrator=2;//
  int fluxintegrator=atoi(argv[2]);
//int thick=0;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  bool charged=1;
  int current=2;
  
  string homedir=argv[3];   //"/home/wim/Code/share";

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
  
  
  //find reasonable integration limits
  double T_mu=0;
  double E_out;  
  for(int j=0;j<=100;j++){ //loop over T_mu
//  cout << j << "/100" << endl;
    T_mu=30*j;    
    E_out=T_mu+massmu;
    for(int i=1;i<1000;i++){
    
//    double E_in=E_out+(3.E03-E_out)*0.001*i; //possible incoming lepton energies
//    costhetamu=(1.-((Q2+massmu*massmu)/(2.*E_in*E_out)))/sqrt(1.-(massmu*massmu)/(E_out*E_out));
//    costhetamu=(-Q2-massmu*massmu+2.*E_in*E_out)/(2.*E_in*sqrt(E_out*E_out-massmu*massmu));
      
			costhetamu=2.E-3*i-1.;
			double E_in=(-Q2-massmu*massmu)/(2.*sqrt(E_out*E_out-massmu*massmu)*costhetamu-2.*E_out);
			double omega=E_in-E_out;
//    cout << "Costhetamu " << costhetamu << endl; 
 
      TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(TLeptonKinematics::muon,E_in);
      WeakQECross obs(lepton,&Nucleus,prec,integrator,homedir,charged, 1.03E03, screening, 0.);
       
//    double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton->GetLeptonMass()*lepton->GetLeptonMass()/(E_out*E_out))*costhetamu)
//    	-lepton->GetLeptonMass()*lepton->GetLeptonMass();
//    double x=Q2/(2.*MASSP*omega);
//    cout << E_in << " " << omega << " " << Q2*1.E-06 << " " << x <<  endl;
      for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
        //lowest p_m is always at theta_cm -1
        TKinematics2to2 kin("","",Nucleus.getMassA(),
        (shell<Nucleus.getPLevels()? Nucleus.getMassA_min_proton(): Nucleus.getMassA_min_neutron()) 
			  +Nucleus.getExcitation()[shell],
			  shell<Nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
	      double tempmax=-1., tempmin=1.;
	
        //anything above 300 MeV contribution will be negligible
        if(getBound(tempmax,tempmin,Nucleus,*lepton,E_in,E_out,costhetamu,shell)<5.E02){ 
          if(E_in<E_low)  E_low=E_in;
	        if(E_in>E_high) E_high=E_in;
	        if(max<tempmax) max=tempmax;
	        if(min>tempmin) min=tempmin;
	        if(cthmax[shell]<tempmax) cthmax[shell]=tempmax;
          if(cthmin[shell]>tempmin) cthmin[shell]=tempmin;
	        if(Tmax<T_mu) Tmax=T_mu;
          if(Tmin>T_mu) Tmin=T_mu;
        }
      }
    }
  }  
  if(E_high>3.E03) E_high=3.E03;
  if(E_low>E_high) E_low=E_high;
 cout << endl;
 cout << "min=" << min << "   max=" << max << endl;
 cout << "Tmin=" << Tmin << "  Tmax=" << Tmax << endl;
 cout << "E_low=" << E_low << "  E_high=" << E_high << endl << endl;
  
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

  numint::mdfunction<numint::vector_d,3> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

  numint::array<double,3> lower = {{E_low,min,Tmin+massmu}};
  numint::array<double,3> upper = {{E_high,max,Tmax+massmu}};
  
  F.f=adap_intPm;
  unsigned count=0;
  string stf=homedir+"/statefile";
  ostringstream qstr;
  qstr << Q2;
  string qst=qstr.str();
  stf+=qst;
  char statefile[1024];
  strncpy(statefile,stf.c_str(),sizeof(statefile));
  statefile[sizeof(statefile)-1]=0;
  int nregions,fail,countt;
  vector<double> err(2,0.);
  vector<double> prob(2,0.);
  if(fluxintegrator==0) numint::vegas( mdf, lower,upper,1, 1.E-01, 1.E-25,0,123,100, 50000,1000, 500, 1000, 0,statefile,countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==1) numint::cuhre( mdf, lower,upper,1, 1.E-01, 1.E-25,0,100, 50000,11, statefile,nregions,countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==2) numint::cube_romb(mdf,lower,upper,1.E-25,1.E-01,avgcross,count,0); //1.E-20,1.E-03
  if(fluxintegrator==3) numint::cube_adaptive(mdf,lower,upper,1.E-25,1.E-03,2E02,2E04,avgcross,count,0);
  count=countt;
  
  cout << Q2*1.E-6 << " " << avgcross[0]*1.E19*2.*PI/Nucleus.getN()
  << " " << avgcross[1]*1.E19*2.*PI/Nucleus.getZ() << " " << count << endl;

}

//integrandum
void adap_intPm(numint::vector_d & results, double E_in, double costhetacm, double E_out,
	 MeanFieldNucleusThick &nucleus, int current, double *cthmax, double Q2,
	 double prec, int integrator, string homedir, int maxEval, bool charged, bool screening){	  

  results=numint::vector_d(2,0.);
  double omega=E_in-E_out;  
  double costhetamu=(-Q2-massmu*massmu+2.*E_in*E_out)/(2.*E_in*sqrt(E_out*E_out-massmu*massmu));
//cout << "Costhetamu " << costhetamu << endl;

  results=numint::vector_d(2,0.);  

  if(costhetamu<=1 && costhetamu>=-1 && omega>0) {
      
    TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(TLeptonKinematics::muon,E_in);
    WeakQECross pobs(lepton,&nucleus,prec,integrator,homedir,charged,1.03E03,screening,0.);  
    
    //we can fix the vector boson kinematics
    pobs.getPlepton()->SetBeamEnergy(E_in);

    for(int shell=0;shell<nucleus.getTotalLevels();shell++) {
      if(costhetacm<cthmax[shell]){
        TKinematics2to2 kin("","",nucleus.getMassA(),
          (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
          +nucleus.getExcitation()[shell],
          shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,costhetacm);
        double pm=sqrt(E_out*E_out-massmu*massmu);
        if(!kin.IsPhysical()||kin.GetPYlab()<200.||kin.GetPklab()>5.E02){ //final nucleon momentum too low, impose cut!
          for(int i=0;i<2;i++) results[i]+=0.;
        }
        else{
          double result=pobs.getDiffWeakQECross(kin,current,1,0,0,1,shell,0.,2E04,0,1);   // prec..2E04
  //	  			cout << "Result " << result << endl;

          //Jacobian    
          double jcb;
          double Enu=(MASSN*MASSN-(MASSP-30.)*(MASSP-30.)-lepton->GetLeptonMass()*lepton->GetLeptonMass()+2.*(MASSP-30.)*E_out)
                    /(2.*(MASSP-30.-E_out+pm*costhetamu));
          jcb=2.*pm*Enu*(1.+(E_out-pm*costhetamu)/(MASSP-30.-E_out+pm*costhetamu));
          result/=jcb;      
  //      cout << "Jacobian C " << jcb << endl;
//           cout << E_in << " " << costhetacm << " " << E_out << " " << omega << " " <<  shell << " " << kin.GetPYlab() << " " << kin.GetPklab() << " " << result << endl;
    
          results[(shell<nucleus.getPLevels()?1:0)]+= result; //results[0] neutrino, results[1] antineutrino
          //  << " " << E_in <<  " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
          //	<< acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
          //	<< " " << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetQsquared() << " "
          //	<< kin.GetXb()*nucleus.getMassA()/MASSP << " " << result << " " << result*(shell<nucleus.getPLevels()?
          //	interpolate(MiniBooNE_antineut_flux_norm,E_in,25,120,1):interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1)) << endl;
        }
      }  
    }
    //fold with flux
    results[0]*=interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1);
    results[1]*=interpolate(MiniBooNE_antineut_flux_norm,E_in,25,120,1);
    delete lepton;
  }
  else {results[0]=0; results[1]=0;}

}



double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell){
  
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();
  double tempmin=500.;
  double costhmin=-1.;
  for(int i=0;i<=100;i++){ //loop over all com angles
    TKinematics2to2 kin("","",nucleus.getMassA(),
			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			+nucleus.getExcitation()[shell],
			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.+i*0.02);
    double pm=kin.GetPklab(); //initial nucleon momentum
    if(pm<tempmin) {tempmin=pm; costhmin=-1.+i*0.02;}
  }
  //   cout << omega << " " << Q2/1.E06 << " " << omega/(2.*MASSP*omega) << " " << high << " " << low << " " << pm << endl;
  double temptemp=costhmin;//com costheta value where the min initial nucl momentum is
  //   cout << costhmin << " " << tempmin << endl;
  if(tempmin<5.E02){
    TKinematics2to2 kin1("","",nucleus.getMassA(),
			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			+nucleus.getExcitation()[shell],
			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,1.);
    high=1;
    if(kin1.GetPklab()>5.E02) high=getMax(high,temptemp,nucleus,lepton,Q2,omega,shell);
    temptemp=costhmin;
    TKinematics2to2 kinmin1("","",nucleus.getMassA(),
			(shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
			+nucleus.getExcitation()[shell],
			shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,-1.);
    low=-1.;
    if(kinmin1.GetPklab()>5.E02) getMin(temptemp,low,nucleus,lepton,Q2,omega,shell);
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
  if(pm<5.E02) high=(high+low)/2.;
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
  if(pm<5.E02) low=(high+low)/2.;
  else high=(high+low)/2.;
  if((high-low)<1.E-04) { return low;}
  else return getMax(high,low,nucleus,lepton,Q2,omega,shell);  
}