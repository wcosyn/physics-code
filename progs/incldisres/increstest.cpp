#include <iostream>
#include <cstdlib>
#include <cmath>


using namespace std;

#include <constants.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <NuclStructure.hpp>
#include <InclusiveCrossRes.hpp>
#include <Resonance.hpp>

//args are 
//1. symm (-1 asymm delta 0 symm delta 1 no delta)
//2. offshell type
//3. fix propagator masses 
//4. Q2 (GeV)
//5. deuteron wf
//6. Nucleon struct function parametrization
//7. sigma0 (mb)
//8. sigma slope (mb/GeV)
//9. etc resonance masses
//last. t_choice

// run ../../bin/inclusiveres 1 3 0 5. paris Alekhin 40. 0. 1.4 2.4

int main(int argc, char *argv[])
{
    double Ein=5000.;  
    double Q2=atof(argv[4])*1.E06;
    int bla=40;
    TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
    double teller[bla-1], noemer[bla-1], fsi1[bla-1], fsi2[bla-1];
    Resonance res1=Resonance(atof(argv[9])*1.E03,1./sqrt(2),atof(argv[7]), atof(argv[8]));
    Resonance res2=Resonance(atof(argv[10])*1.E03,-1./sqrt(2),atof(argv[7]), atof(argv[8]));
//     double sigma1 = res1.getSigma(10.E06);
//     double sigma2 = res2.getSigma(10.E06);
//     double tttt=sqrt(sigma1*sigma1+sigma2*sigma2);
//     res1.setCoeff(sigma2/tttt);
//     res2.setCoeff(-sigma1/tttt);
    
    for(int i=1;i<bla;i++){
      double x=1./bla*i;
//       cout << x << endl;
      //TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
      teller[i-1]=0.,fsi1[i-1]=0.,fsi2[i-1]=0.,noemer[i-1]=0.;
      for(int proton=0;proton<=1;++proton){
	InclusiveCrossRes Dinc(proton,argv[6],argv[5],*elec,atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),0);
	//Resonance res3=Resonance(atof(argv[11])*1.E03,-1./sqrt(2),atof(argv[7]), atof(argv[8]));
	Dinc.addResonance(res1);
	Dinc.addResonance(res2);
	//Dinc.addResonance(res3);
	NuclStructure nucleon(proton, Q2, x, 0, argv[6]);
	teller[i-1]+=Dinc.calc_F2Dinc(Q2,x);
	double f1,f2;
   	Dinc.calc_F2DincFSI(f1,f2,Q2,x);
	//cout << "bla " << f1 << " " << f2 << endl;
	fsi1[i-1]+=f1;
	fsi2[i-1]+=f2;
	noemer[i-1]+=nucleon.getF2();
// 	cout << endl;
	//cout << noemer << endl;
      }
       cout << x << " " << (Q2*(1./x-1.)+MASSP*MASSP)/1.E06 << " " << teller[i-1]/noemer[i-1] << " " << (teller[i-1]-fsi1[i-1])/noemer[i-1] << " " << (teller[i-1]-fsi2[i-1])/noemer[i-1] << endl;
    }
    delete elec;
//     cout << endl << endl << endl;
//     for(int i=0;i<19;i++){
//       double x = (i+1)*0.05;
//     cout << x << " " << teller[i]/noemer[i] << " " << (teller[i]-fsi1[i])/noemer[i] << " " << (teller[i]-fsi2[i])/noemer[i] << " " << (Q2*(1./x-1.)+MASSP*MASSP)/1.E06 << endl;
//     }
}
