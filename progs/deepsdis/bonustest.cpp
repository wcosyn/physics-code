//bonus data test program, test formulas from their MC, my cross section
// generate ratios to compare with data (using scattering parameter values from fit)

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include <constants.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <DeuteronStructure.hpp>
#include <DeuteronMomDistr.hpp>
#include <DeuteronCross.hpp>
#include "bonusdata.h"
#include "bonusfits2.h"
#include "bonusfits.h"

struct Ftor {

  static void exec(const numint::array<double,4> &x, void *param, numint::vector_d &ret) {
    Ftor &p = * (Ftor *) param;
    p.f(ret,x[0],x[1],x[2],x[3],*p.cross,p.Ebeam,p.proton);
  }
  DeuteronCross *cross;
  double Ebeam;
  bool proton;

  void (*f)(numint::vector_d &, double Q2, double W, double ps, double costheta,
	    DeuteronCross &cross,double Ebeam, bool proton);

};

void adap_avg(numint::vector_d &, double Q2, double W, double ps, double costheta,
	    DeuteronCross &cross,double Ebeam, bool proton);

double get_normfit_bonus(int beamindex, int Qindex, int Windex, int psindex, int offshell, bool lc);
double get_normfit_sigma_deeps_bonus(int beamindex, int Qindex, int Windex, int psindex, int offshell, bool lc, bool q2dep);
double get_normrig(int psindex, int beamindex);
double sigmaparam(double W, double Q2, bool Q2dep);

int main(int argc, char *argv[])
{
  double Ebeam = data::Ebeam[atoi(argv[1])]; //[0,1]
  int beamindex=atoi(argv[1]);
  int Qindex = atoi(argv[2]); //[0,2] for 4 GeV, [1,2] for 5 GeV
  int Windex = atoi(argv[3]); //[0,4]
  int pindex = atoi(argv[4]); //[0,3]
  int offshellset = atoi(argv[5]);
//   double sigmain = atof(argv[6]); //sigma in mb
//   double betain = atof(argv[7]); //beta in GeV^-2
  bool lc = 0;//atoi(argv[8]);
  bool q2dep = atoi(argv[6]);
  
  bool proton=0;
  double phi=0.;

  double Wprime = 0.5*(data::W[Windex]+data::W[Windex+1]);
  double Q2 = 0.5*(data::Q2[Qindex]+data::Q2[Qindex+1]);
  double pr = 0.5*(data::ps[pindex]+data::ps[pindex+1]);
  DeuteronCross test("paris",proton,"CB",sigmaparam(Wprime,Q2,q2dep),8.,-0.5,8.,1.2,offshellset,1E03);
//   DeuteronCross test("paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  for (int i=0;i<=36;i+=1){
    double costhetar=-0.9+i*0.05;
    if(1/*abs(costhetar)>1.E-03*/){
        double MCres,modelrespw, modelresfsi;
        test.getBonusMCresult(MCres,modelrespw,modelresfsi,Q2,Wprime,Ebeam,pr,costhetar,proton,0,lc);
        
    //     Ftor F;
    //     F.cross = &test;
    //     F.Ebeam = Ebeam;
    //     F.proton = proton;
    //     
    //     numint::mdfunction<numint::vector_d,4> mdf;
    //     mdf.func = &Ftor::exec;
    //     mdf.param = &F;
    // 
    //     unsigned neval = 0;
    //     numint::array<double,4> lower = {{data::Q2[Qindex],data::W[Windex],data::ps[pindex],costhetar-0.1}};
    //     numint::array<double,4> upper = {{data::Q2[Qindex+1],data::W[Windex+1],data::ps[pindex+1],costhetar+0.1}};
    //     vector<double> avgcross(3,0.); 
    //     F.f=adap_avg;
    //     unsigned count=0;
    // //  numint::cube_romb(mdf,lower,upper,1.E-20,1.E-03,avgcross,count,0);
    //     numint::cube_adaptive(mdf,lower,upper,1.E-20,1.E-03,2E02,1E04,avgcross,count,0);
    //     avgcross[0]/=(data::Q2[Qindex+1]-data::Q2[Qindex])*(data::W[Windex+1]-data::W[Windex])*(data::ps[pindex+1]-data::ps[pindex])*0.2;
    //     avgcross[1]/=(data::Q2[Qindex+1]-data::Q2[Qindex])*(data::W[Windex+1]-data::W[Windex])*(data::ps[pindex+1]-data::ps[pindex])*0.2;
    bool proton =0;
    double massi=proton? MASSP:MASSN;
    double massr=proton? MASSN:MASSP;
    
    double Er=sqrt(massr*massr+pr*pr);
    double Einoff=MASSD-Er;
    double massoff=sqrt(Einoff*Einoff-pr*pr);
    
    double xprime=Q2/(Wprime*Wprime-massoff*massoff+Q2);
    //calc nu
    double prz=pr*costhetar;
    double aaa=Einoff*Einoff-prz*prz;
    double bbb=-Einoff*Q2/xprime;
    double ccc=Q2*Q2/(4.*xprime*xprime)-Q2*prz*prz;
    
    double discr=sqrt(bbb*bbb-4.*aaa*ccc);
    if(std::isnan(discr)) discr=0.;
    double nu1=(-bbb+discr)/(2.*aaa);
    double nu2=(-bbb-discr)/(2.*aaa);
//     cout << nu1 << " " << nu2 << " " << (bbb*bbb-4.*aaa*ccc) << endl;
    double nu=nu2;

    if(costhetar<0.) nu=nu1;
    double qvec=sqrt(Q2+nu*nu);
    /*  double xx=Q2/2./(Einoff*nu+prz*qvec);
    cout << xprime << " " << xx << endl;*/
    double x=Q2/(2.*massi*nu);
        
        cout << Ebeam << " " << Q2 << " " << Wprime << " " << pr << " " << costhetar << " " << x << " " << xprime << " " 
            << MCres/get_normrig(pindex,beamindex) << 
            " " << modelrespw << " " << modelresfsi << " " 
            << get_normrig(pindex,beamindex) 
            <<  /*" " << avgcross[0] << " " << avgcross[1] << " " << count <<*/ endl;
    }
  }
}



void adap_avg(numint::vector_d &result, double Q2, double W, double ps, double costheta,
	    DeuteronCross &cross,double Ebeam, bool proton){
  result=vector<double>(3,0.);  
  cross.getBonusMCresult(result[0],result[1],result[2],Q2,W,Ebeam,ps,costheta,proton,1,0);
//   cout << Q2 << " " << W << " " << ps << " " << costheta << " " << result[0] << endl;
}

double get_normfit_bonus(int beamindex, int Qindex, int Windex, int psindex, int offshell, bool lc){
  if(lc){
    if(offshell==3){
      if(beamindex==0) return bonusfits::normfits_fix3_off3_lc1_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_fix3_off3_lc1_beam5[Qindex-1][Windex][psindex];
    }
    if(offshell==4){
      if(beamindex==0) return bonusfits::normfits_fix3_off4_lc1_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_fix3_off4_lc1_beam5[Qindex-1][Windex][psindex];
    }
  }
  else{
    if(offshell==3){
      if(beamindex==0) return bonusfits::normfits_fix3_off3_lc0_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_fix3_off3_lc0_beam5[Qindex-1][Windex][psindex];
    }
    if(offshell==4){
      if(beamindex==0) return bonusfits::normfits_fix3_off4_lc0_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_fix3_off4_lc0_beam5[Qindex-1][Windex][psindex];
    }
  }
  return 0./0.;
}

double get_normfit_sigma_deeps_bonus(int beamindex, int Qindex, int Windex, int psindex, int offshell, bool lc, bool q2dep){
  if(q2dep){
    if(offshell==3){
      if(beamindex==0){
	if(Qindex<2) return bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep0_beam4[Qindex][Windex][psindex];
	else return bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep1_beam4[0][Windex][psindex];
      }
      else{
	if(Qindex<2) return bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep0_beam5[Qindex-1][Windex][psindex];
	else return bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep1_beam5[0][Windex][psindex];
      }
    }
    if(offshell==4){
      if(beamindex==0) return bonusfits::normfits_sigmadeeps_fix3_off4_lc0_q2dep1_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_sigmadeeps_fix3_off4_lc0_q2dep1_beam5[Qindex-1][Windex][psindex];
    }
  }
  else{
    if(offshell==3){
      if(beamindex==0) return bonusfits::normfits_sigmadeeps_fix3_off3_lc0_q2dep0_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_sigmadeeps_fix3_off3_lc0_q2dep0_beam5[Qindex-1][Windex][psindex];
    }
    if(offshell==4){
      if(beamindex==0) return bonusfits::normfits_sigmadeeps_fix3_off4_lc0_q2dep0_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_sigmadeeps_fix3_off4_lc0_q2dep0_beam5[Qindex-1][Windex][psindex];
    }
  }
  return 0./0.;
}

double get_normrig(int psindex, int beamindex){
  if(beamindex==0) return bonusfits2::normfitsrig_final_q2dep1_beam4_setCB[psindex];
  else return bonusfits2::normfitsrig_final_q2dep1_beam5_setCB[psindex];
}

//sigma in mbarn
double sigmaparam(double W, double Q2, bool Q2dep){
  if(!Q2dep||Q2<1.8E06){
      if(abs(W*W-1.232*1.232E06)<2.5E5) return 65;
  return (25.3*1.E-06*2.3+53*(sqrt(W*W>5.76E06?5.76E06:W*W)-MASSP)*1.E-03)
	/(1.8);
  }
  else{ 
    if(abs(W*W-1.232*1.232E06)<2.5E5) return 65.*1.8E06/Q2;
    return (25.3*1.E-06*2.3+53*(sqrt(W*W>5.76E06?5.76E06:W*W)-MASSP)*1.E-03)
	  /(1.E-06*Q2);
  }
}
