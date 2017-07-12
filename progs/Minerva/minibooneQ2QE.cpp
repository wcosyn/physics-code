// Program that computes single differential cross sections for QE neutrino nucleus scattering
//We calculate the single differential dsigma/dQ2_QE cross sections for free nucleon, RFG, carbon.
//Q2_QE is the experimental Q2 that is used to obtain a Q2 from a T_mu cos_mu event.  This gives rise to nasty jacobians that need to be taken into account

// Can be used for different experiments, Miniboone and minerva kinematics currently implemented, T2K will be added as well

//run like
// > minibooneQ2QE [Q2 in GeV^2] [integration algorithm for flux integration (see around line 680 for options)] 
//  [max initial nucl momentun in MeV] [min final nucl momentum in MeV] [plane wave 1 or not 0] [ROMEA 1 or not 0 in fsi]
// [experiment "miniboone", "minerva", "t2k"] [sharedir]

// example for plane wave calculation for miniboone
// > minibooneQ2QE 0.1 3 500. 200. 1 0 miniboone ~/Code/trunk/share  


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

//starts at E=25 MeV with a delta of 25 MeV (up to 2 GeV)
const double t2k_neut_flux_norm[] = {
1.93448341e-05 , 
4.73717519e-05 , 
7.85287411e-05 , 
0.000113062088 , 
0.00015156703 , 
0.000195039029 , 
0.000244842988 , 
0.000302692788 , 
0.000370712689 , 
0.000449333806 , 
0.000550624798 , 
0.000665154366 , 
0.000786867808 , 
0.000913024822 , 
0.0010442721 , 
0.0011845812 , 
0.00134127948 , 
0.00149808032 , 
0.00166478462 , 
0.0018066637 , 
0.00190126372 , 
0.00193840358 , 
0.00192018773 , 
0.00186097308 , 
0.00178742153 , 
0.001663081 , 
0.00157644483 , 
0.00149002427 , 
0.00138141611 , 
0.00124332379 , 
0.00108355703 , 
0.000925001164 , 
0.000805689255 , 
0.000693591835 , 
0.000593686244 , 
0.000506177836 , 
0.000431251246 , 
0.000369091227 , 
0.000319913263 , 
0.00028389183 , 
0.000261221954 , 
0.000233985265 , 
0.000208616329 , 
0.000187916841 , 
0.000172995147 , 
0.000163368924 , 
0.000156893278 , 
0.000149771076 , 
0.000136594026 , 
0.000128640575 , 
0.000120574245 , 
0.000113174974 , 
0.00010683274 , 
0.000101629652 , 
9.72680791e-05 , 
9.31220275e-05 , 
8.81857559e-05 , 
8.37010448e-05 , 
8.03965158e-05 , 
7.77693131e-05 , 
7.54089269e-05 , 
7.29561943e-05 , 
7.01442623e-05 , 
6.67268469e-05 , 
6.25807952e-05 , 
6.06924914e-05 , 
5.85271009e-05 , 
5.65464397e-05 , 
5.49454853e-05 , 
5.36113585e-05 , 
5.21643451e-05 , 
4.98963309e-05 , 
4.58528993e-05 , 
4.5811852e-05 , 
4.46727099e-05 , 
4.31230728e-05 , 
4.16350085e-05 , 
4.04137681e-05 , 
3.94183044e-05 , 
3.83715305e-05
};


const double t2k_aneut_flux_norm[] = {
1.88499998e-05 , 
4.61599993e-05 , 
7.6520002e-05 , 
0.000110170004 , 
0.000147690007 , 
0.000190050006 , 
0.000238580003 , 
0.000294950005 , 
0.000361229992 , 
0.000437840004 , 
0.000536539999 , 
0.000648139976 , 
0.000766740006 , 
0.000889669987 , 
0.00101756002 , 
0.00115428003 , 
0.00130697002 , 
0.00145976001 , 
0.00162220001 , 
0.00176044996 , 
0.00185263006 , 
0.00188881997 , 
0.00187107001 , 
0.00181337004 , 
0.00174169999 , 
0.00162054005 , 
0.00153611996 , 
0.00145191001 , 
0.00134607998 , 
0.00121152005 , 
0.00105584005 , 
0.000901339983 , 
0.000785080018 , 
0.000675850024 , 
0.000578499981 , 
0.000493230007 , 
0.000420219993 , 
0.000359650003 , 
0.000311729993 , 
0.000276629988 , 
0.000254540006 , 
0.000228000004 , 
0.000203279997 , 
0.000183109994 , 
0.000168569997 , 
0.000159190007 , 
0.000152880006 , 
0.000145939994 , 
0.000133099995 , 
0.000125349994 , 
0.000117490003 , 
0.000110280002 , 
0.0001041 , 
9.90300032e-05 , 
9.47799999e-05 , 
9.07400026e-05 , 
8.5929998e-05 , 
8.15599997e-05 , 
7.83400028e-05 , 
7.57800008e-05 , 
7.34799978e-05 , 
7.10899985e-05 , 
6.83499966e-05 , 
6.50200018e-05 , 
6.09800009e-05 , 
5.91399985e-05 , 
5.70299999e-05 , 
5.51000012e-05 , 
5.35399995e-05 , 
5.22400005e-05 , 
5.083e-05 , 
4.86200006e-05 , 
4.46800004e-05 , 
4.46400009e-05 , 
4.35299989e-05 , 
4.20200013e-05 , 
4.05700011e-05 , 
3.93800001e-05 , 
3.84099985e-05 , 
3.73900002e-05
};



void normalize(double flux[], double dx){
  double sum=0;
  for(int i=0;i<20;i++) sum+=flux[i]*dx;
  for(int i=0;i<20;i++) flux[i]/=sum;
  };


//integration struct 
struct Ftor {  //Carbon

  static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
    Ftor &p = * (Ftor *) param;
    p.f(ret,x[0],x[1],x[2],*p.pNucleus,p.current,p.cthmax,p.Q2QE,
        p.prec,p.integrator,p.homedir,p.maxEval,p.charged,p.screening,p.enable_romea,
        p.max_initial_nucl_mom,p.min_final_nucl_mom,p.pw, p.exp
       );
  }
  MeanFieldNucleusThick *pNucleus;
  int current;
  double *cthmax;
  double Q2QE;
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
  void (*f)(numint::vector_d &, double E_omega, double costhetacm, double E_in, 
            MeanFieldNucleusThick &pNucleus, int current, double *cthmax, double Q2QE,
            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening, bool enable_romea,
            double max_initial_nucl_mom, double min_final_nucl_mom, int pw, string exp
           );

};

struct FtorH {  //Hydrogen

  static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
    FtorH &p = * (FtorH *) param;
    p.f(ret,x[0],p.Q2QE,p.charged,p.maxbeam,p.exp);
  }
  int current;
  double Q2QE;
  bool charged;
  double maxbeam;
  string exp;
  void (*f)(numint::vector_d &, double E_out, double Q2QE, bool charged, double maxbeam, string exp);

}; 


struct FtorRFG {  //Hydrogen

  static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
    FtorRFG &p = * (FtorRFG *) param;
    p.f(ret,x[0],x[1],p.Q2QE,p.charged,p.Pauli,p.maxbeam,p.exp);
  }
  int current;
  double Q2QE;
  bool charged;
  bool Pauli;
  double maxbeam;
  string exp;
  void (*f)(numint::vector_d &, double omega, double E_out, double Q2QE, bool charged, bool Pauli, double maxbeam, string exp);

}; 


//integration for Hydrogen
void int_hydr(numint::vector_d &, double E_out,double Q2QE,bool charged, double maxbeam, string exp);

void int_RFG(numint::vector_d &, double omega,double E_out, double Q2QE,bool charged, bool Pauli, double maxbeam, string exp);


//determine boundaries of kinematics
double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell, double max_initial_nucl_mom, double min_final_nucl_mom);
double getMin(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);
double getMax(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double Q2, double omega, int shell, double max_initial_nucl_mom);

// integration over missing momentum + over T_mu(=E_out+C)
void adap_intPm(numint::vector_d &, double omega, double costhetacm, double E_in,
		            MeanFieldNucleusThick &pNucleus,int current, double *cthmax, double Q2QE,
		            double prec, int integrator, string homedir, int maxEval, bool charged, bool screening, bool enable_romea,
                            double max_initial_nucl_mom, double min_final_nucl_mom, int pw, string exp
               );


int main(int argc, char *argv[])
{
  double Q2QE=atof(argv[1])*1.E06;   // give input in GeV^2
//double T_mu=atof(argv[1]); // muon kinetic energy in MeV
//   double costhetamu; //=atof(argv[1]);
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
  string homedir=argv[8];   //"/home/wim/Code/share";

  MeanFieldNucleusThick Nucleus(nucleus,homedir);

   
  vector<double> avgcross(2,0.); //neutrino and antineutrino
  //estimates for integration bounds
  double cthmax[Nucleus.getTotalLevels()];
  double cthmin[Nucleus.getTotalLevels()];
  for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {cthmax[shell]=-1.;cthmin[shell]=1.;}
  double max=-1.,min=1.; //overall min and max center of mass cosine theta angles
  
  double minbeam=0., maxbeam=0.; //different beam energy intervals for the different experiments
  
  if(!exp.compare("miniboone")) { maxbeam=3.E03;}
  else if(!exp.compare("minerva")) { 
    minbeam=1.5E03; 
    maxbeam=10.E03;
    normalize(Minerva_nu_flux,500);
    normalize(Minerva_anu_flux,500);
  }  
  else if(!exp.compare("t2k")) { maxbeam=2.E02; }  //still to implement!
  else {cerr << "invalid experiment name chosen" << endl << "Choose either miniboone, minerva or t2k" << endl; assert(1==0);}
  
  
  
  if(!exp.compare("t2k")) { }  //still to implement!
 //setting up variables dat will determine integration limits
  double E_out_max=minbeam;
  double E_out_min=maxbeam;
  double E_low=maxbeam;
  double E_high=minbeam;
  double omega_low=maxbeam;
  double omega_hi=minbeam;
  double E_muH_low=maxbeam;
  double E_muH_hi=minbeam;
  double E_RFG_out_max=minbeam;
  double E_RFG_out_min=maxbeam;
  double omega_RFG_low=maxbeam;
  double omega_RFG_hi=minbeam;
   
//looking for minimum Emu to have physical costhetamu according to the formula used for Q2QE and EnuQE
  
  double massin=MASSn-34.;
  double massout=MASSn;
  
  double DD=Q2QE-massin*massin+massout*massout;
  
  double A=4.*massin*massin*Q2QE;
  double B=2.*DD*(Q2QE-massmu*massmu)*massin;
  double C=-massmu*massmu*DD*DD-massin*massin*pow(Q2QE+massmu*massmu,2.);
  
  double discr=B*B-4.*A*C;
  
  double minEmu=(-B+sqrt(discr))/2./A;  //other one is always negative, Emu needs to be bigger than this.

//   double Emu=2.*minEmu;
//   double pmu=sqrt(Emu*Emu-massmu*massmu);
//   
//   double costhetamu=(-(Q2QE+massmu*massmu)*(massin-Emu)+(2.*massin*Emu-(massin*massin+massmu*massmu-massout*massout))*Emu)/pmu/(Q2QE+2.*massin*Emu-massin*massin+massout*massout);
//   
//   double EnuQE=(2.*massin*Emu-(massin*massin+massmu*massmu-massout*massout))/2./(massin-Emu+pmu*costhetamu);
//   cout << costhetamu << " " << Q2QE << " " << -massmu*massmu+2.*EnuQE*(Emu-pmu*costhetamu) << endl;
  
  
  //find reasonable integration limits
  for(int j=0;j<=100;j++){ //loop over E_mu
    //reconstruct kinematics according to the QE recipe used in the analysis
    double E_out=(maxbeam-minEmu)/100.*j+minEmu;
    double p_out=sqrt(E_out*E_out-massmu*massmu);
    double costhetamu = (-(Q2QE+massmu*massmu)*(massin-E_out)+(2.*massin*E_out-(massin*massin+massmu*massmu-massout*massout))*E_out)
                          /p_out/(Q2QE+2.*massin*E_out-massin*massin+massout*massout);
    if(abs(costhetamu)<=1.){                              
      //hydrogen limits
      double Enu=(E_out*2.*MASSn-massmu*massmu)/2./(MASSn-E_out+p_out*costhetamu);
      if(Enu>minbeam&&Enu<maxbeam){
        if(E_out>E_muH_hi) E_muH_hi=E_out;
        if(E_out<E_muH_low) E_muH_low=E_out;
      }
                            
                            
      //check for carbon and RFG                     
      for(int i=0;i<=100;i++){
        double E_in=E_out+(maxbeam-E_out)*1.E-02*i;
        double omega=E_in-E_out;
        TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(TLeptonKinematics::muon,E_in);
        
        //RFG checks
        double Q2=-massmu*massmu+2.*E_in*(E_out-p_out*costhetamu);
        double qvec = sqrt(Q2+omega*omega);
        double lambda=omega/2./sqrt(MASSP*MASSN);
        double kappa=qvec/2./sqrt(MASSP*MASSN);
        double tau=Q2/4./MASSP/MASSN;
        double k_Fermi=228.;
        
        double xi=sqrt(1.+k_Fermi*k_Fermi/MASSP/MASSN)-1.;
        double psi=pow(lambda-tau,2.)/(xi*((1.+lambda)*tau+kappa*sqrt(tau*(1.+tau))));
      //   cout << "blaaa " << psi << endl;
        if(abs(psi)<1.){
          if(E_out>E_RFG_out_max) E_RFG_out_max = E_out;
          if(E_out<E_RFG_out_min) E_RFG_out_min = E_out;
          if(omega>omega_RFG_hi) omega_RFG_hi = omega;        
          if(omega<omega_RFG_low) omega_RFG_low = omega;
        }

        //Carbon checks
        for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
          double tempmax=-1., tempmin=1.;
          
          //anything above 500 MeV contribution will be negligible
          if(getBound(tempmax,tempmin,Nucleus,*lepton,E_in,E_out,costhetamu,shell,max_initial_nucl_mom,min_final_nucl_mom)<max_initial_nucl_mom){ 
            if(E_in<E_low)  E_low=E_in;
            if(E_in>E_high) E_high=E_in;
            if(max<tempmax) max=tempmax;
            if(min>tempmin) min=tempmin;
            if(cthmax[shell]<tempmax) cthmax[shell]=tempmax;
            if(cthmin[shell]>tempmin) cthmin[shell]=tempmin;
            if(E_out_max<E_out) E_out_max=E_out;
            if(E_out_min>E_out) E_out_min=E_out;
            if(omega>omega_hi) omega_hi=omega;
            if(omega<omega_low) omega_low=omega;
  //             cout << E_in << " " << E_out << " " << tempmin << " " << tempmax << " " << pm_min << endl;
          }
        }
      }
    }
  }  
  if(E_high>maxbeam) E_high=maxbeam;
  if(E_low<minbeam) E_low=minbeam;
  if(E_muH_low<minbeam) E_muH_low=minbeam;
  if(E_low>E_high) E_low=E_high;
  cout << endl;
  cout << "min=" << min << "   max=" << max << endl;
  cout << "Eoutmin=" << E_out_min << "  Eoutmax=" << E_out_max << endl;
  cout << "omega_low=" << omega_low << " " << "  omega_high=" << omega_hi << endl;
  cout << "E_low=" << E_low << "  E_high=" << E_high << endl << endl;

  cout << "E_muH_low=" << E_muH_low << "  E_muH_hi=" << E_muH_hi << endl << endl;
  
  cout << "E_RFG_out_min=" << E_RFG_out_min << "  E_RFG_out_max=" << E_RFG_out_max << endl;
  cout << "omega_RFG_low=" << omega_RFG_low <<"  omega_RFG_hi=" << omega_RFG_hi << endl << endl;
  
  FtorH FH;
  FH.current=current;
  FH.Q2QE=Q2QE;
  FH.charged=charged;
  FH.exp=exp;
  FH.maxbeam=maxbeam;

  numint::mdfunction<numint::vector_d,1> mdfH;
  mdfH.func = &FtorH::exec;
  mdfH.param = &FH;

  numint::array<double,1> lowerH = {{E_muH_low}};  
  numint::array<double,1> upperH = {{E_muH_hi}};
  
  FH.f=int_hydr;
  vector<double> avgcrossH(2,0.); 
  unsigned count=0;
  if(!fluxintegrator) numint::cube_romb(mdfH,lowerH,upperH,1.E-30,1.E-03,avgcrossH,count,0); //1.E-20,1.E-03
  else numint::cube_adaptive(mdfH,lowerH,upperH,1.E-30,1.E-03,2E02,4E04,avgcrossH,count,0); 
  
//   cout << "Crosssection H: " << Q2 << " " << avgcrossH[0]*1.E19*2.*PI << " [1E-39 cm2/GeV2]" <<  " "<< count << endl; //2pi because of integration over scattered lepton angle
//   cout << Q2QE << " " << avgcrossH[0]*1.E19*2.*PI << " " << avgcrossH[1]*1.E19*2.*PI  << endl;

 
//RFG calculation
  FtorRFG FRFG;
  
  FRFG.current=current;
  FRFG.Q2QE=Q2QE;
  FRFG.charged=charged;
  FRFG.Pauli = Pauli;
  FRFG.maxbeam=maxbeam;
  FRFG.exp=exp;
  
  numint::mdfunction<numint::vector_d,2> mdfRFG;
  mdfRFG.func = &FtorRFG::exec;
  mdfRFG.param = &FRFG;

  numint::array<double,2> lowerRFG = {{omega_RFG_low,E_RFG_out_min}};
  numint::array<double,2> upperRFG = {{omega_RFG_hi,E_RFG_out_max}};
  
  FRFG.f=int_RFG;
  vector<double> avgcrossRFG(2,0.); 
  count=0;
  if(!fluxintegrator) numint::cube_romb(mdfRFG,lowerRFG,upperRFG,1.E-30,1.E-03,avgcrossRFG,count,0); //1.E-20,1.E-03
  else numint::cube_adaptive(mdfRFG,lowerRFG,upperRFG,1.E-30,1.E-03,2E02,4E04,avgcrossRFG,count,0); 

//   cout << "Crosssection H: " << Q2 << " " << avgcrossH[0]*1.E19*2.*PI << " [1E-39 cm2/GeV2]" <<  " "<< count << endl; //2pi because of integration over scattered lepton angle
    cout << Q2QE*1.E-06 << " " << avgcrossH[0]*1.E19*2.*PI << " " << avgcrossH[1]*1.E19*2.*PI << " " << avgcrossRFG[0]*1.E19*2.*PI << " " << avgcrossRFG[1]*1.E19*2.*PI  << endl;

 
  //initialize object -- Carbon
  Ftor F;
  F.pNucleus = &Nucleus;
  F.current=current;
  F.cthmax=cthmax;
  F.Q2QE=Q2QE;
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

//   numint::array<double,3> lower = {{omega_low,min,Tmin+massmu}};
//   numint::array<double,3> upper = {{omega_hi,max,Tmax+massmu}};
  numint::array<double,3> lower = {{omega_low,min,E_low}};
  numint::array<double,3> upper = {{omega_hi,max,E_high}};
  
  F.f=adap_intPm;
  count=0;
  string stf=homedir+"/statefile/minibooneQ2QE";
  ostringstream qstr;
  qstr << Q2QE;
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
  int maxEvalcuba=5000;
  double xgiven[3]={200.,-1.,500.};
  cout << "start integration C" << endl;
  if(fluxintegrator==0) numint::vegas( mdf, lower,upper,nvec, epsrel,epsabs,flags,seed,minEval, maxEvalcuba,1000, 500, 1000, 0,(stf+"vegas").c_str(),countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==1) numint::cuhre( mdf, lower,upper,nvec, epsrel,epsabs,flags,minEval, maxEvalcuba,11, (stf+"cuhre").c_str(),nregions,countt,fail,avgcross,err,prob ); 
  if(fluxintegrator==2) numint::cube_romb(mdf,lower,upper,1.E-25,1.E-01,avgcross,count,0); //1.E-20,1.E-03
  if(fluxintegrator==3) numint::cube_adaptive(mdf,lower,upper,1.E-25,1.E-03,2E02,2E04,avgcross,count,0);
  if(fluxintegrator==4) numint::divonne(mdf, lower,upper,nvec, epsrel, epsabs,flags,seed,
           minEval, maxEvalcuba,7,7,0,5,0,10,0.25,0,3,NULL, 0, 0,
           (stf+"divonne").c_str(),nregions, countt, fail, avgcross,err,prob );
  
  //factor 2\pi because of integration over muon polar angle
  cout << Q2QE*1.E-6 << " " << avgcrossH[0]*1.E19*2.*PI << " " << avgcrossH[1]*1.E19*2.*PI << " " << avgcrossRFG[0]*1.E19*2.*PI << " " << avgcrossRFG[1]*1.E19*2.*PI << " " <<  
    avgcross[0]*1.E19*2.*PI/Nucleus.getN() << " " << avgcross[1]*1.E19*2.*PI/Nucleus.getZ() << " " << count << endl;

}

//integrandum carbon
void adap_intPm(numint::vector_d & results, double omega, double costhetacm, double E_in,
	 MeanFieldNucleusThick &nucleus, int current, double *cthmax, double Q2QE,
	 double prec, int integrator, string homedir, int maxEval, bool charged, bool screening,
         bool enable_romea, double max_initial_nucl_mom, 
                double min_final_nucl_mom, int pw, string exp){	  

  results=numint::vector_d(2,0.);
  double E_out = E_in-omega;
  if(/*E_in>3.E03||E_in<0.*/E_out<0.) { return;}
  double p_out = sqrt(E_out*E_out-massmu*massmu);
  double massin=MASSn-34.;
  double massout=MASSn;
  
  
  double costhetamu=(-(Q2QE+massmu*massmu)*(massin-E_out)+(2.*massin*E_out-(massin*massin+massmu*massmu-massout*massout))*E_out)
                          /p_out/(Q2QE+2.*massin*E_out-massin*massin+massout*massout);
                          
  if(abs(costhetamu)<=1.) {
    
    double Q2=-massmu*massmu+2.*E_in*(E_out-p_out*costhetamu);  //this is the real Q2
//     cout << Q2 << " " << Q2QE << endl;
    TLeptonKinematics *lepton = TLeptonKinematics::CreateWithBeamEnergy(TLeptonKinematics::muon,E_in);
    WeakQECross pobs(lepton,&nucleus,prec,integrator,homedir,charged,1.03E03,screening,enable_romea);  
    
    for(int shell=0;shell<nucleus.getTotalLevels();shell++) {
      if(costhetacm<cthmax[shell]){
        TKinematics2to2 kin("","",nucleus.getMassA(),
          (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())
          +nucleus.getExcitation()[shell],
          shell<nucleus.getPLevels()?MASSN:MASSP,"qsquared:wlab:costhkcm",Q2,omega,costhetacm);
        if(!kin.IsPhysical()||kin.GetPYlab()<min_final_nucl_mom||kin.GetPklab()>max_initial_nucl_mom){ //final nucleon momentum too low, impose cut!
          for(int i=0;i<2;i++) results[i]+=0.; /*cout << "cut " << kin.GetPYlab() << " " << kin.GetPklab() << endl;*/
        }
        else{
          double result=pobs.getDiffWeakQECross(kin,current,0,0,0,pw,shell,0.,maxEval,0,1);   // prec..2E04
  //	  			cout << "Result " << result << endl;

          //Jacobian              
          double EnuQE=(2.*massin*E_out-(massin*massin+massmu*massmu-massout*massout))/2./(massin-E_out+p_out*costhetamu);
          double jcb=-2.*p_out*EnuQE*massin/(massin-E_out+p_out*costhetamu);
          result/=abs(jcb);      
//           cout << E_in << " " << costhetacm << " " << E_out << " " << omega << " " <<  shell << " " << kin.GetPYlab() << " " << kin.GetPklab() << " " << result << endl;
    
          results[(shell<nucleus.getPLevels()?1:0)]+= result; //results[0] neutrino, results[1] antineutrino
           cout << shell << " " << E_in <<  " " << costhetacm << " " << p_out << " "  << kin.GetCosthklab()<< " " 
          	<< kin.GetCosthYlab() << " " << kin.GetPklab() << " " << kin.GetPYlab() 
          	<< " " << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetQsquared() << " "
          	<< kin.GetXb()*nucleus.getMassA()/MASSP << " " << result << endl;
        }
      }  
    }
    //fold with flux
    if(!exp.compare("miniboone")){
      results[0]*=interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1);
      results[1]*=interpolate(MiniBooNE_antineut_flux_norm,E_in,25,120,1);
    }
    else if(!exp.compare("minerva")){
      results[0]*=interpolate(Minerva_nu_flux,E_in,500,20,0);
      results[1]*=interpolate(Minerva_anu_flux,E_in,500,20,0);
    }
    else if(!exp.compare("t2k")){
      results[0]*=interpolate(t2k_neut_flux_norm,E_in,25,80,1);
      results[1]*=interpolate(t2k_aneut_flux_norm,E_in,25,80,1);
    }
    else assert(1==0);
    
    delete lepton;
  }
  else {results[0]=0; results[1]=0; /*cout << "unphys " << omega << " " << costhetacm << " " << E_out << " " << E_in << " " << 0. << " " << 0. << endl;*/}

}

//integrandum hydrogen
void int_hydr(numint::vector_d & results, double E_out, 
              double Q2QE, bool charged, double maxbeam, string exp){
    
  results=numint::vector_d(2,0.);  
  double massin=MASSn-34.;
  double massout=MASSn;
  double p_out=sqrt(E_out*E_out-massmu*massmu);
  double costhetamu = (-(Q2QE+massmu*massmu)*(massin-E_out)+(2.*massin*E_out-(massin*massin+massmu*massmu-massout*massout))*E_out)
                        /p_out/(Q2QE+2.*massin*E_out-massin*massin+massout*massout);
  if(abs(costhetamu)>1.) {/*cout << "costheta" << endl;*/ /*cout << E_out << " " << result[0] << " " << result[1] << endl;*/ return;}
  double EnuQE=(2.*massin*E_out-(massin*massin+massmu*massmu-massout*massout))/2./(massin-E_out+p_out*costhetamu);                          
  double E_in=(E_out*2.*MASSn-massmu*massmu)/2./(MASSn-E_out+p_out*costhetamu);
  double Q2=-massmu*massmu+2.*E_in*(E_out-p_out*costhetamu);
  if(E_in>maxbeam) return;
//  if(E_in>maxbeam)  cout << "int " << Q2QE << " " << EnuQE << " "<< costhetamu << " " << E_out << " " << E_in << " " << E_in-E_out << endl;
 
 
  results[0]=WeakQECross::getElWeakQECross(Q2,E_in,1,charged,1.03E03,2);
  results[1]=WeakQECross::getElWeakQECross(Q2,E_in,0,charged,1.03E03,2);
  
  if(!exp.compare("miniboone")){
    results[0]*=interpolate(MiniBooNE_antineut_flux_norm,E_in,25,120,1);
    results[1]*=interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1);
  }
  else if(!exp.compare("minerva")){
    results[0]*=interpolate(Minerva_nu_flux,E_in,500,20,0);
    results[1]*=interpolate(Minerva_anu_flux,E_in,500,20,0);
  }
  else if(!exp.compare("t2k")){
    results[0]*=interpolate(t2k_neut_flux_norm,E_in,25,80,1);
    results[1]*=interpolate(t2k_aneut_flux_norm,E_in,25,80,1);
  }
  else assert(1==0);

//   cout << E_out << " " << result[0] << " " << result[1] << endl;
} 


void int_RFG(numint::vector_d & results, double omega,double E_out, double Q2QE,bool charged, bool Pauli, double maxbeam, string exp){
    
  results=numint::vector_d(2,0.);  
  double E_in=omega+E_out;
  if(E_in>maxbeam)  {return;}
  double massin=MASSn-34.;
  double massout=MASSn;
  double p_out=sqrt(E_out*E_out-massmu*massmu);
  double costhetamu = (-(Q2QE+massmu*massmu)*(massin-E_out)+(2.*massin*E_out-(massin*massin+massmu*massmu-massout*massout))*E_out)
                        /p_out/(Q2QE+2.*massin*E_out-massin*massin+massout*massout);

  double Q2=-massmu*massmu+2.*E_in*(E_out-p_out*costhetamu);
  if(abs(costhetamu)>1.) {/*cout << "costheta" << endl;*/ return;}
  double EnuQE=(2.*massin*E_out-(massin*massin+massmu*massmu-massout*massout))/2./(massin-E_out+p_out*costhetamu);                          
  double kf=228.;
  WeakQECross::getRFGWeakQECross(results[0],results[1],Q2, E_in, omega,kf,1,1.03E03,2,Pauli);
//   cout << omega << " " << E_out << " " << crossRFG_p << " " << crossRFG_n << endl;
  
  if(!exp.compare("miniboone")){
    results[0]*=interpolate(MiniBooNE_antineut_flux_norm,E_in,25,120,1);
    results[1]*=interpolate(MiniBooNE_neut_flux_norm,E_in,25,120,1);
  }
  else if(!exp.compare("minerva")){
    results[0]*=interpolate(Minerva_nu_flux,E_in,500,20,0);
    results[1]*=interpolate(Minerva_anu_flux,E_in,500,20,0);
  }
  else if(!exp.compare("t2k")){
    results[0]*=interpolate(t2k_neut_flux_norm,E_in,25,80,1);
    results[1]*=interpolate(t2k_aneut_flux_norm,E_in,25,80,1);
  }
  else assert(1==0);
    
} 





double getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, TLeptonKinematics &lepton,
	      double E_in, double E_out, double costhetamu, int shell, double max_initial_nucl_mom, double min_final_nucl_mom){
  
  double omega=E_in-E_out;
  double Q2=2.*E_in*E_out*(1.-sqrt(1.-lepton.GetLeptonMass()*lepton.GetLeptonMass()/(E_out*E_out))*costhetamu)
      -lepton.GetLeptonMass()*lepton.GetLeptonMass();
  double tempmin=max_initial_nucl_mom;
  double costhmin=-1.;
  for(int i=0;i<=100;i++){ //loop over all com angles to find com scattering angle where initial nucl momentum is minimal
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
  //if it is below the limit we consider, we find the interval in angles where it is below
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
  if((high-low)<1.E-03) { return high;}
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
  if((high-low)<1.E-03) { return low;}
  else return getMax(high,low,nucleus,lepton,Q2,omega,shell,max_initial_nucl_mom);  
}
