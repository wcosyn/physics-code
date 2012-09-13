/*! \mainpage Glauber ISI/FSI RMSGA C++ Project
 * \author Wim Cosyn
 * \date 16/08/2011
 * \brief This code implements classes for the Glauber RMSGA formalism, including CT and SRC effects.
 * 
 * \details 
 * - It contains classes for a mean field nucleus, a mean field nucleus with densities (used in thickness calculations). <BR>
 * - A class for correlated FSI calculations containing all the needed functions and a grid for the gamma functions. <BR>
 * - Four abstract classes (one for a general FSI grid, one that adds a CT grid, one for a general thickness FSI grid and one that adds a thickness CT grid). <BR>
 * - Two glauber classes : one without thickness, one with (also adds SRC to the ISI/FSI). <BR>
 * - A special class for the glauber grid of one particle, exploiting the symmetry along the z-axis (no phi dependence).  <BR>
*/
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
void getBound(double &high, double &low, MeanFieldNucleusThick *pnucleus, double Q2, double omega, int shell);
void adap_intPm(numint::vector_d &, double pm,Cross* pObs, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec);


//run ./transp [nucleon] [Q2,GeV^2] [Ein, MeV] [thetae, degr] [thickness] [sharedir]
int main(int argc, char *argv[])
{
  
  string homedir=argv[10];
  double Q2=atof(argv[2])*1.E06;
  double Ein=atof(argv[3]);
  double thetae=atof(argv[4])*DEGRTORAD;
  double omega=Ein-Q2/(4.*Ein*pow(sin(thetae/2.),2.));
  int thick=atoi(argv[5]);
  int current=atoi(argv[6]);
  double prec=atof(argv[7]);
  bool userset=atoi(argv[8]);
  double screening=atof(argv[9]);
  //cout << Ein-omega << endl;
  
  MeanFieldNucleusThick Nucleus(atoi(argv[1]),homedir);
  //TKinematics2to2 kin("","",.getMassA(),Carbon.getMassA_min_proton(),MASSP,"qsquared:wlab:pklab",atof(argv[1]),atof(argv[2]),atof(argv[3]));
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
  Cross obs(*elec,&Nucleus,prec,homedir, userset, screening);

  
//   for(int i=0;i<5;i++) total[i]=0.;
//   
//   for(int shell=0;shell</*Nucleus.getPLevels()*/1;shell++){
//     double low=-1.,high=1.;
//     getBound(high,low,&Nucleus,Q2,omega,shell);
//     double results[5];
//     double pestimate=0.;
//     rombergerN(intPm,-1.,low,5, results,1.E-02,3,10,&pestimate,&obs,elec,&Nucleus,Q2,omega,shell, 
// 	       thick,current);
//     for(int j=0;j<5;j++) total[j]+=results[j];
//     
//     cout << Q2/1.E06 << " " << Ein << " " << shell << " " << results[0]/results[4] << " " << 
// 	    results[1]/results[4] << " " << results[2]/results[4] << " " << results[3]/results[4] << endl;
//   }
//   cout << endl;
//     cout << Q2/1.E06 << " " << Ein << " " << total[0]/total[4] << " " << 
// 	    total[1]/total[4] << " " << total[2]/total[4] << " " << total[3]/total[4] << endl;

  vector<double> totalcross(5,0.) 
  struct Ftor {

    static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
      Ftor &p = * (Ftor *) param;
      p.f(ret,x[0],p.pObs,p.pNucleus,p.elec);
    }
    Cross *pObs;
    MeanFieldNucleusThick *pNucleus;
    TElectronKinematics *elec;
    void (*f)(numint::vector_d &, double pm,Cross* pObs, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec);
  };

  Ftor F;
  F.pObs = &Obs;
  F.pNucleus = &Nucleus;
  F.elec=elec;

  numint::mdfunction<numint::vector_d,1> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

  unsigned neval = 0;
  numint::array<double,1> lower = {{0.}};
  numint::array<double,1> upper = {{1.,1.,2.*PI}};
  
  double SIGNIF = 1.E-05;
  
  
  
  
  delete elec;
  return 0;
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
  results[0]=p_obs->getDiffCross(kin,  current, thick, 0, 0, 0, shell, 0.);
  results[1]=p_obs->getDiffCross(kin,  current, thick, 1, 0, 0, shell, 0.);
  results[2]=p_obs->getDiffCross(kin,  current, thick, 0, 1, 0, shell, 0.);
  results[3]=p_obs->getDiffCross(kin,  current, thick, 1, 1, 0, shell, 0.);
  results[4]=p_obs->getDiffCross(kin,  current, thick, 0, 0, 1, shell, 0.);
  cout << pm << " " << kin.IsPhysical() << " " << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetCosthYlab() << " " << acos(kin.GetCosthklab())*RADTODEGR << " " << kin.GetCosthklab() <<" " << results[0] << " " << results[1] << " " << results[4] << endl;
  return;
}

void getBound(double &high, double &low, MeanFieldNucleusThick *pnucleus, double Q2, double omega, int shell){
  TKinematics2to2 kin("","",pnucleus->getMassA(),
		      pnucleus->getMassA_min_proton()+pnucleus->getExcitation()[shell],
		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,(high+low)/2.);
  double pm=kin.GetPklab();
  if(pm<300.) low=(high+low)/2.;
  else high=(high+low)/2.;
  if((high-low)<1.E-04) { return;}
  else getBound(high,low,pnucleus,Q2,omega,shell);  
}

