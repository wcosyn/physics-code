// three resonance model for inclusive deuteron DIS, never used in proper calculations (see incldis for that)

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <time.h>

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

// run ../../bin/inclusiveres 1 2 5. paris Alekhin 1.4 2.4

int main(int argc, char *argv[])
{
    double Ein=5000.;  
    double Q2=atof(argv[4])*1.E06;
    srand(atoi(argv[12]));
    int bla=40;
    
    double alpha = float(rand())/RAND_MAX*sqrt(2./3.);
    double beta = 0.;
    if(rand()%2) beta = (-alpha+sqrt(2.-3.*alpha*alpha))/2.;
    else beta = (-alpha-sqrt(2.-3.*alpha*alpha))/2.;
    double gamma = -alpha-beta;
    cout << alpha << " " << beta << " " << gamma << " " << alpha*alpha+beta*beta+gamma*gamma << endl << endl;
    Resonance res1=Resonance(atof(argv[9])*1.E03,alpha,atof(argv[7]), atof(argv[8]));
    Resonance res2=Resonance(atof(argv[10])*1.E03,beta,atof(argv[7]), atof(argv[8]));
    Resonance res3=Resonance(atof(argv[11])*1.E03,gamma,atof(argv[7]), atof(argv[8]));
    for(int i=1;i<bla;i++){
      double x=1./bla*i;
      //TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
      TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
      double teller=0.,fsi1=0.,fsi2=0.,noemer=0.;
      for(int proton=0;proton<=1;++proton){
	InclusiveCrossRes Dinc(proton,argv[6],argv[5],*elec,atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[13]));
	Dinc.addResonance(res1);
	Dinc.addResonance(res2);
	Dinc.addResonance(res3);
	NuclStructure nucleon(proton, Q2, x, 0, argv[6]);
	teller+=Dinc.calc_F2Dinc(Q2,x);
	double f1,f2;
  	Dinc.calc_F2DincFSI(f1,f2,Q2,x);
	//cout << "bla " << f1 << " " << f2 << endl;
	fsi1+=f1;
	fsi2+=f2;
	noemer+=nucleon.getF2();
	//cout << noemer << endl;
      }
      cout << x << " " << (Q2*(1./x-1.)+MASSP*MASSP)/1.E06 << " " <<  teller/noemer << " " << (teller-fsi1)/noemer << " " << (teller-fsi2)/noemer << endl;
    }
    
}
