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

// run ../bin/inclusive 1 3 5. paris Alekhin

int main(int argc, char *argv[])
{
    double Ein=5000.;  
    double Q2=atof(argv[3])*1.E06;
    string strucname = argv[5];
    string wf = argv[4];
    int symm = atoi(argv[1]);
    int offshell = atoi(argv[2]);
    TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
    for(int i=1;i<20;i++){      
      double x=0.05*i;
      //TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
      double teller=0.,fsi1=0.,fsi2=0.,noemer=0.;
      for(int proton=0;proton<=1;++proton){
	InclusiveCross Dinc(proton,strucname,wf,*elec,symm,offshell);
	NuclStructure nucleon(proton, Q2, x, 0, strucname);
	teller+=Dinc.calc_F2Dinc(Q2,x);
	double f1,f2;
  	Dinc.calc_F2DincFSI(f1,f2,Q2,x);
	//cout << "bla " << f1 << " " << f2 << endl;
	fsi1+=f1;
	fsi2+=f2;
	noemer+=nucleon.getF2();
      }
      cout << x << " " << teller/noemer << " " << (teller-fsi1)/noemer << " " << (teller-fsi2)/noemer << endl;
      //cout << x << " " << teller << " " << noemer << " " << teller/noemer << endl;
    }
    delete elec;
}
