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
    double Ein=11000.;  
    double Q2=atof(argv[2])*1.E06;
    string strucname = argv[4];
    string wf = argv[3];
    int symm = atoi(argv[1]);
    TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
    std::vector<double> resonances;
     resonances.push_back(1.232E03);
     resonances.push_back(1.5E03);
    resonances.push_back(1.75E03);
//     for(int i=1;i<100;i++){      
//       double x=0.01*i;
//       cout << x << " " << sqrt(Q2*(1./x-1.)+MASSN*MASSN)*1.E-03 << " " ;
//       //TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
//       double teller=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffsuppr1=0.,fsioffsuppr2=0.,noemer=0.;
//       for(int proton=0;proton<=1;++proton){
// 	InclusiveCross Dinc(proton,strucname,wf,*elec,resonances,symm,4);
// 	NuclStructure nucleon(proton, Q2, x, 0, strucname);
// // 	teller+=Dinc.calc_F2Dinc(Q2,x);
// 	noemer+=nucleon.getF2();
// 	cout << nucleon.getF2() << " ";
//       }
//       cout << endl;
//      }
//     exit(1);
    for(int i=8;i<50;i++){      
      double x=0.02*i;
      //TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
      double teller=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffsuppr1=0.,fsioffsuppr2=0.,noemer=1.;
      double azz=0.,azzfsi1=0.,azzfsi2=0.,azzfsioff1=0.,azzfsioff2=0.,azzfsioffsuppr1=0.,azzfsioffsuppr2=0.;
      for(int proton=0;proton<=1;++proton){
	InclusiveCross Dinc(proton,strucname,wf,*elec,resonances,symm,4);
	NuclStructure nucleon(proton, Q2, x, 0, strucname);
	Dinc.calc_Azzinc(azz,teller,Q2,x,0.);
	double f1,f2,azzf1,azzf2;
//   	Dinc.calc_AzzincFSI(azzf1,azzf2,f1,f2,Q2,x);
	//cout << "bla " << f1 << " " << f2 << endl;
	azzfsi1+=azzf1;
	azzfsi2+=azzf2;
	fsi1+=f1;
	fsi2+=f2;
	if(!Dinc.getDeutwf()->IsRelativistic()){
// 	  Dinc.calc_AzzincFSI_off(azzf1, azzf2,f1,f2,Q2,x,0.);
  // // 	cout << "bla " << f1 << " " << f2 << endl;
// 	  azzfsioff1+=azzf1;
// 	  azzfsioff2+=azzf2;
// 	  fsioff1+=f1;
// 	  fsioff2+=f2;
	  //suppressed
	  Dinc.setOffshell(3);
// 	  Dinc.calc_AzzincFSI_off(azzf1,azzf2,f1,f2,Q2,x,0.);
	  azzfsioffsuppr1+=azzf1;
	  azzfsioffsuppr2+=azzf2;
	  fsioffsuppr1+=f1;
	  fsioffsuppr2+=f2;
	}
	
      }
      cout << x << " " << sqrt(Q2*(1./x-1.)+MASSN*MASSN)*1.E-03 << " " << azz << " " << azzfsi1 << " " << 
      azzfsi2 << " " << azzfsioffsuppr1 << " " << azzfsioffsuppr2 << " " <<  teller << " " << 
      fsi1 << " " << fsi2 << " " << fsioffsuppr1 << " " << fsioffsuppr2 << endl;
      //cout << x << " " << teller << " " << noemer << " " << teller << endl;
    }
    delete elec;
}
