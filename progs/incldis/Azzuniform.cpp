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

//continuum contribution with uniform shape
//make azzuniform
// run ../../bin/disuniform x Q2[GeV^2] deuteronwf structfunc widths(interval = 2*width)[MeV]
// fi: ../../bin/disuniform 0.1 5. paris SLAC 300 500 700

int main(int argc, char *argv[])
{
    double Ein=11000.;  
    double Q2=atof(argv[2])*1.E06;
    string strucname = argv[4];
    string wf = argv[3];
    double x = atof(argv[1]);
    double width1 = atof(argv[5]);
    double width2 = atof(argv[6]);
    double width3 = atof(argv[7]);
    
    TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
    std::vector<double> resonances;
    std::vector<double> widths;
    resonances.push_back(2.025E03);
    widths.push_back(width1);
    resonances.push_back(2.5E03);
    widths.push_back(width2);
    resonances.push_back(3.4E03);
    widths.push_back(width3);
    std::vector<double> centrals = resonances;
    double teller=0.,fsi1=0.,fsi2=0.,fsioffsuppr1=0.,fsioffsuppr2=0.,noemer=1.;
    double azz=0.,azzfsi1=0.,azzfsi2=0.,azzfsioffsuppr1=0.,azzfsioffsuppr2=0.;
    for(int proton=0;proton<=1;++proton){
      InclusiveCross Dinc(proton,strucname,wf,*elec,resonances,3);
      Dinc.calc_Azzinc(azz,teller,Q2,x);
      cout << azz << " " << teller << endl;
      double f1,f2,azz1,azz2;
//       Dinc.calc_AzzincFSI_uniform(azz1,azz2,f1,f2,Q2,x,centrals,widths);
      azzfsi1+=azz1;
      azzfsi2+=azz2;      
      fsi1+=f1;
      fsi2+=f2;
      if(!Dinc.getDeutwf()->IsRelativistic()){
	Dinc.setOffshell(3);
// 	Dinc.calc_AzzincFSI_uniform_off(azz1,azz2,f1,f2,Q2,x,centrals,widths);
	azzfsioffsuppr1+=azz1;
	azzfsioffsuppr2+=azz2;
	fsioffsuppr1+=f1;
	fsioffsuppr2+=f2;
      }
    }
    cout << x << " " << sqrt(Q2*(1./x-1.)+MASSN*MASSN)*1.E-03 << " " << azz << " " << 
    azzfsi1 << " " << azzfsi2 << " " << azzfsioffsuppr1 << " " << azzfsioffsuppr2 << " " <<
    teller << " " << fsi1 << " " << fsi2 << " " << fsioffsuppr1 << " " << fsioffsuppr2 << endl;
    //cout << x << " " << teller << " " << noemer << " " << teller/noemer << endl;    
    delete elec;
}
