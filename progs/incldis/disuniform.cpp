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
//make uniform
// run ../../bin/disuniform x Q2[GeV^2] deuteronwf structfunc widths(interval = 2*width)[MeV]
// fi: ../../bin/disuniform 0.1 5. paris SLAC 300 500 700

int main(int argc, char *argv[])
{
    double Ein=5000.;  
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
    for(int proton=0;proton<=1;++proton){
      InclusiveCross Dinc(proton,strucname,wf,*elec,resonances,3);
      teller+=Dinc.calc_F2Dinc(Q2,x);
      double f1,f2;
      Dinc.calc_F2DincFSI_uniform(f1,f2,Q2,x,centrals,widths);
      fsi1+=f1;
      fsi2+=f2;
      if(!Dinc.getDeutwf()->IsRelativistic()){
	Dinc.setOffshell(3);
	Dinc.calc_F2DincFSI_uniform_off(f1,f2,Q2,x,centrals,widths);
	fsioffsuppr1+=f1;
	fsioffsuppr2+=f2;
      }
    }
    cout << x << " " << sqrt(Q2*(1./x-1.)+MASSN*MASSN)*1.E-03 << " " << teller/noemer << " " << (teller-fsi1)/noemer 
      << " " << (teller-fsi2)/noemer << " " 
     << (teller-fsi1-fsioffsuppr1)/noemer << " " 
      << (teller-fsi2-fsioffsuppr2)/noemer << endl;
    //cout << x << " " << teller << " " << noemer << " " << teller/noemer << endl;    
    delete elec;
}
