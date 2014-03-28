#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <TElectronKinematics.h>
#include <DQEinclusive.hpp>

// run ../bin/DQEinc Q2[GeV] wf ffparam x

int main(int argc, char *argv[])
{
    double Ein=atof(argv[6]);  
    double Q2=atof(argv[1])*1.E06;
    int ffparam=atoi(argv[3]);
    int integrator=atoi(argv[4]);
    double x=atof(argv[5]);
    
    string wf = argv[2];
    TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);

    double qvec=sqrt(Q2+pow(Q2/(2.*MASSP*x),2.));
    double nu=Q2/(2.*MASSP*x);
    //       cout << "q " << qvec << " " << nu << endl;
    //TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
    double teller1=0.,teller2=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffPV1=0.,fsioffPV2=0.,noemer=1.;
    double f2_1,f2_2;
    for(int proton=0;proton<=1;++proton){
      DQEinclusive Dinc(proton,ffparam,wf,*elec,4);
      double f1,f2,f3,f4;
      Dinc.calc_Crossinc(f1,f2,Q2,x,2,integrator);
      cout << "pwresult " << x << " " << proton << " " << f1 << " " << f2 << endl;
      teller1+=f1;
      teller2+=f2;
//       Dinc.calc_CrossincFSI(f1,f2,f3,f4,Q2,x,2,integrator,2E05);
//       cout << "fsiresult " << x << " " << proton << " " << f1 << " " << f2 << " " << f3 << " " << f4 << endl;
//       fsi1+=f1;
//       fsi2+=f2;
//       fsioff1+=f3;
//       fsioff2+=f4;
//       (proton?f2_2:f2_1)=Dinc.getMinpcm();
//       Dinc.calc_CrossincFSI_PVoff(f1,f2,Q2,x,2,integrator);
//       cout << "fsiresultoff " << x << " " << proton << " " << f1 << " " << f2 << " " << endl;
//       fsioffPV1+=f1;
//       fsioffPV2+=f2;
      
    }
    cout << x << " " << teller1 << " " << teller2 << " " << teller1+teller2 << " " << fsi1 << " " << fsi2 << " "
    << fsi1+fsi2 << " " << fsioff1 << " " << fsioff2 << " " << fsioff1+fsioff2 << " " 
    << fsioffPV1 << " " << fsioffPV2 << " " << fsioffPV1+fsioffPV2 << " "<<
    teller1+teller2+fsi1+fsi2+fsioff1+fsioff2 << " " << f2_1 << " " << f2_2 << endl;
//     cout << x << " " << teller << " " << noemer << " " << teller/noemer << endl;

    delete elec;
}
