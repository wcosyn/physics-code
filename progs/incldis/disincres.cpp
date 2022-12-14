//program calculates F2D inclusive with FSI (three resonance regions), see PhysRevC.89.014612 for details

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <NuclStructure.hpp>
#include "InclusiveCross.hpp"

// make inc
// run ../../bin/inclusive Q2[GeV^2] deuteronwf nucleonstruct
// fi ../../bin/inclusive 5. paris Alekhin

int main(int argc, char *argv[])
{
    double Ein=5000.;  
    double Q2=atof(argv[1])*1.E06;
    string strucname = argv[3];
    string wf = argv[2];
    TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
    std::vector<double> resonances;
    resonances.push_back(1.232E03);
    resonances.push_back(1.5E03);
    resonances.push_back(1.75E03);
    for(int i=1;i<50;i++){      
      double x=0.02*i;
      //TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
      double teller=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffsuppr1=0.,fsioffsuppr2=0.,noemer=1.;
      for(int proton=0;proton<=1;++proton){
	InclusiveCross Dinc(proton,strucname,wf,*elec,resonances,3);
	NuclStructure nucleon(proton, Q2, x, 0, strucname);
	teller+=Dinc.calc_F2Dinc(Q2,x);
	double f1,f2;
  	Dinc.calc_F2DincFSI(f1,f2,Q2,x);
	//cout << "bla " << f1 << " " << f2 << endl;
	fsi1+=f1;
	fsi2+=f2;
	if(!Dinc.getDeutwf()->IsRelativistic()){
// 	  Dinc.calc_F2DincFSI_off(f1,f2,Q2,x);
//   // // 	cout << "bla " << f1 << " " << f2 << endl;
// 	  fsioff1+=f1;
// 	  fsioff2+=f2;
	  //suppressed
	  Dinc.setOffshell(3);
	  Dinc.calc_F2DincFSI_off(f1,f2,Q2,x);
	  fsioffsuppr1+=f1;
	  fsioffsuppr2+=f2;
	}
	//noemer+=nucleon.getF2();
      }
      cout << x << " " << sqrt(Q2*(1./x-1.)+MASSN*MASSN)*1.E-03 << " " << teller/noemer << " " << (teller-fsi1)/noemer << " " << (teller-fsi2)/noemer << " " 
	<<(teller-fsi1-fsioffsuppr1)/noemer << " " << (teller-fsi2-fsioffsuppr2)/noemer << " " << (teller-fsi1-fsioffsuppr1)/noemer << " " 
	<< (teller-fsi2-fsioffsuppr2)/noemer << endl;
      //cout << x << " " << teller << " " << noemer << " " << teller/noemer << endl;
    }
    delete elec;
}
