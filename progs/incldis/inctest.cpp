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

// run ../bin/inclusive 1 5. paris Alekhin

int main(int argc, char *argv[])
{
    double Ein=5000.;  
    double Q2=atof(argv[2])*1.E06;
    string strucname = argv[4];
    string wf = argv[3];
    int symm = atoi(argv[1]);
    TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
    for(int i=1;i<20;i++){      
      double x=0.05*i;
      //TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
      double teller=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffsuppr1=0.,fsioffsuppr2=0.,noemer=0.;
      for(int proton=0;proton<=1;++proton){
	InclusiveCross Dinc(proton,strucname,wf,*elec,symm,4);
	NuclStructure nucleon(proton, Q2, x, 0, strucname);
	teller+=Dinc.calc_F2Dinc(Q2,x);
	double f1,f2;
  	Dinc.calc_F2DincFSI(f1,f2,Q2,x);
	//cout << "bla " << f1 << " " << f2 << endl;
	fsi1+=f1;
	fsi2+=f2;
//   	Dinc.calc_F2DincFSI_off(f1,f2,Q2,x);
	//cout << "bla " << f1 << " " << f2 << endl;
// 	fsioff1+=f1;
// 	fsioff2+=f2;
	//suppressed
	Dinc.setOffshell(3);
//   	Dinc.calc_F2DincFSI_off(f1,f2,Q2,x);
// 	fsioffsuppr1+=f1;
// 	fsioffsuppr2+=f2;
	
	
	noemer+=nucleon.getF2();
      }
      cout << x << " " << teller/noemer << " " << (teller-fsi1)/noemer << " " << (teller-fsi2)/noemer << " " 
	<<(teller-fsi1-fsioff1)/noemer << " " << (teller-fsi2-fsioff2)/noemer << " " << (teller-fsi1-fsioffsuppr1)/noemer << " " 
	<< (teller-fsi2-fsioffsuppr2)/noemer << endl;
      //cout << x << " " << teller << " " << noemer << " " << teller/noemer << endl;
    }
    delete elec;
}
