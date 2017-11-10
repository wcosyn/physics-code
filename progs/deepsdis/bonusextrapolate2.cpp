//bonus data test program, test formulas from their MC, my cross section
// generate ratios to compare with data (using scattering parameter values from fit)
//used in extrapolation python script (mainly for plots, not meaningful fits
// this one does not use data, any kinematics

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
#include "bonusfits.h"


double get_normfit_bonus(int beamindex, int Qindex, int Windex, int psindex, int offshell, bool lc);
double get_normfit_sigma_deeps_bonus(int beamindex, int Qindex, int Windex, int psindex, int offshell, bool lc, bool q2dep);
double sigmaparam(double W, double Q2, bool Q2dep);
void get_data(double &data, double& error, int beamindex, int Qindex, int Windex, int psindex, int cosindex);

int main(int argc, char *argv[])
{
  int beamindex = atoi(argv[1]);
  double Ebeam = data::Ebeam[beamindex]; //[0,1]
  double Q2 = atof(argv[2]); //[0,2] for 4 GeV, [1,2] for 5 GeV
  double Wprime = atof(argv[3]); //[0,4]
  double pr = atof(argv[4]); //[0,3]

  int cosindex = atoi(argv[5]);
  double costhetar=-0.9+cosindex*0.2;
  int offshellset = atoi(argv[6]);
  bool q2dep = atoi(argv[7]);
  double xref = atof(argv[8]);
  
  bool lc=0;
  double betain=8.;
  
  bool proton=0;
  double phi=0.;

  DeuteronCross test("paris",proton,"CB",sigmaparam(Wprime,Q2,q2dep),betain,-0.5,8.,1.2,offshellset,1E03);
//     double MCres,modelrespw, modelresfsi;
  
    //not interested in data here
    double result=0., error=0.;
//     get_data(result, error, beamindex,Qindex,Windex,pindex,cosindex);
    test.getBonusextrapolate(Q2,Wprime,Ebeam,pr,costhetar,proton,lc,xref,
      1.,result,error);
    
//     cout << Ebeam << " " << Q2 << " " << Wprime << " " << pr << " " << costhetar << " " << x << " " << xprime << " " 
// 	<< MCres/get_normfit_sigma_deeps_bonus(atoi(argv[1]),Qindex,Windex,pindex,offshellset,lc,q2dep) << " " << modelrespw << " " << modelresfsi << " " << get_normfit_bonus(atoi(argv[1]),Qindex,Windex,pindex,offshellset,lc) <<  /*" " << avgcross[0] << " " << avgcross[1] << " " << count <<*/ endl;
  
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
      if(beamindex==0) return bonusfits::normfits_sigmadeeps_fix3_off3_lc0_q2dep1_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_sigmadeeps_fix3_off3_lc0_q2dep1_beam5[Qindex-1][Windex][psindex];
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

void  get_data(double &data, double& error, int beamindex, int Qindex, int Windex, int psindex, int costhetaindex){
  if(beamindex==0){ 
    error=sqrt(pow(data::bonusdata4[Qindex][Windex][psindex][costhetaindex][2],2.)+
      pow(data::bonusdata4[Qindex][Windex][psindex][costhetaindex][2],2.));
    data=data::bonusdata4[Qindex][Windex][psindex][costhetaindex][1];
  }
  else{
    error=sqrt(pow(data::bonusdata5[Qindex-1][Windex][psindex][costhetaindex][2],2.)+
      pow(data::bonusdata5[Qindex-1][Windex][psindex][costhetaindex][2],2.));
    data=data::bonusdata5[Qindex-1][Windex][psindex][costhetaindex][1];
  }
}
