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



int main(int argc, char *argv[])
{
  double Ebeam = data::Ebeam[atoi(argv[1])]; //[0,1]
  int Qindex = atoi(argv[2]); //[0,2] for 4 GeV, [0,1] for 5 GeV
  int Windex = atoi(argv[3]); //[0,4]
  int pindex = atoi(argv[4]); //[0,3]
  int offshellset = atoi(argv[5]);
  double sigmain = atof(argv[6]); //sigma in mb
  double betain = atof(argv[7]); //beta in GeV^-2
  
  bool proton=0;
  double phi=0.;

  double Wprime = 0.5*(data::W[Windex]+data::W[Windex+1]);
  double Q2 = 0.5*(data::Q2[Qindex]+data::Q2[Qindex+1]);
  double pr = 0.5*(data::ps[pindex]+data::ps[pindex+1]);
  DeuteronCross test("paris",proton,"CB",sigmain,betain,-0.5,8.,1.2,4,1E03);
//   DeuteronCross test("paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  
  for (int i=0;i<10;i+=1){
    double costhetar=-0.9+i*0.2;
    double MCres,modelres;
    test.getBonusMCresult(MCres,modelres,Q2,Wprime,Ebeam,pr,costhetar,proton);
    
    
    Ftor F;
    F.cross = &test;
    F.Ebeam = Ebeam;
    F.proton = proton;
    
    numint::mdfunction<numint::vector_d,4> mdf;
    mdf.func = &Ftor::exec;
    mdf.param = &F;

    unsigned neval = 0;
    numint::array<double,4> lower = {{data::Q2[Qindex],data::W[Windex],data::ps[pindex],costhetar-0.1}};
    numint::array<double,4> upper = {{data::Q2[Qindex+1],data::W[Windex+1],data::ps[pindex+1],costhetar+0.1}};
    vector<double> avgcross(2,0.); 
    F.f=adap_avg;
    unsigned count=0;
//  numint::cube_romb(mdf,lower,upper,1.E-20,1.E-03,avgcross,count,0);
    numint::cube_adaptive(mdf,lower,upper,1.E-20,1.E-03,2E02,1E03,avgcross,count,0);
    avgcross[0]/=(data::Q2[Qindex+1]-data::Q2[Qindex])*(data::W[Windex+1]-data::W[Windex])*(data::ps[pindex+1]-data::ps[pindex])*0.2;
    avgcross[1]/=(data::Q2[Qindex+1]-data::Q2[Qindex])*(data::W[Windex+1]-data::W[Windex])*(data::ps[pindex+1]-data::ps[pindex])*0.2;
    
    cout << costhetar << " " << MCres << " " << modelres << " " << avgcross[0] << " " << avgcross[1] << " " 
	<< modelres/MCres << " " << avgcross[0]/MCres << " " << avgcross[1]/MCres << " " << count << endl;
  }
}



void adap_avg(numint::vector_d &result, double Q2, double W, double ps, double costheta,
	    DeuteronCross &cross,double Ebeam, bool proton){
  result=vector<double>(2,0.);  
  cross.getBonusMCresult(result[0],result[1],Q2,W,Ebeam,ps,costheta,proton);
//   cout << Q2 << " " << W << " " << ps << " " << costheta << " " << result[0] << endl;
}
