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


//make disinc
// run ../../bin/disinc x Q2 deuteronwf structfunc
// fi: ../../bin/disinc 0.1 5. paris SLAC

int main(int argc, char *argv[])
{
    double Ein=5000.;  
    double Q2=atof(argv[2])*1.E06;
    string strucname = argv[4];
    string wf = argv[3];
    double x = atof(argv[1]);
    double central = atof(argv[5]);
    double width = atof(argv[6]);
    
    TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
    std::vector<double> resonances;
    resonances.push_back(central);
    double teller=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffsuppr1=0.,fsioffsuppr2=0.,noemer=1.;
    for(int proton=0;proton<=1;++proton){
      InclusiveCross Dinc(proton,strucname,wf,*elec,resonances,4);
      teller+=Dinc.calc_F2Dinc(Q2,x);
      double f1,f2;
      Dinc.calc_F2DincFSI_distr(f1,f2,Q2,x,central,width);
      fsi1+=f1;
      fsi2+=f2;
      if(!Dinc.getDeutwf()->IsRelativistic()){
	Dinc.calc_F2DincFSI_distr_off(f1,f2,Q2,x,central,width);
	fsioff1+=f1;
	fsioff2+=f2;
	//suppressed
	Dinc.setOffshell(3);
	Dinc.calc_F2DincFSI_distr_off(f1,f2,Q2,x,central,width);
	fsioffsuppr1+=f1;
	fsioffsuppr2+=f2;
      }
    }
    cout << x << " " << sqrt(Q2*(1./x-1.)+MASSN*MASSN)*1.E-03 << " " << teller/noemer << " " << (teller-fsi1)/noemer << " " << (teller-fsi2)/noemer << " " 
      <<(teller-fsi1-fsioff1)/noemer << " " << (teller-fsi2-fsioff2)/noemer << " " << (teller-fsi1-fsioffsuppr1)/noemer << " " 
      << (teller-fsi2-fsioffsuppr2)/noemer << endl;
    //cout << x << " " << teller << " " << noemer << " " << teller/noemer << endl;    
    delete elec;
}
