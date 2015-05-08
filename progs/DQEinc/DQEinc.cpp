#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <TElectronKinematics.h>
#include <DQEinclusive.hpp>

// run ../bin/DQEinc Q2[GeV] wf ffparam x Ebeam[MeV] nopt

int main(int argc, char *argv[])
{
    double Ein=atof(argv[6]);  
    double Q2=atof(argv[1])*1.E06;
    int ffparam=atoi(argv[3]);
    double x=atof(argv[5]);
    int maxEval=atoi(argv[4]);
    bool nopt=atoi(argv[7]);
    int integrator=0;
    bool pv=atoi(argv[8]);
    bool beampol=atoi(argv[9]); //1=polarizatino along bem, 0=polarization along qvec
    string wf = argv[2];
    
    for(int i=0;i<10;i++){
      x=0.8+0.1*i;
    
      
//       double Eout=Ein-nu;
//       double sin2th=Q2/Eout/Ein/4.;
//       cout << nu << " " << Eout << " " << sin2th << " " << qvec << endl;
//       
      TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);

  //     cout << "q " << qvec << " " << nu << " " << Q2/4./Ein/(Ein-nu) << " " << Ein << endl;
      //TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
      double pw_direct=0.,pw_cross=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffPV1=0.,fsioffPV2=0.,noemer=1.;
      //azz have to be divided by the unpolarized cross section afterwards by user!!!!
      double azz_direct=0., azz_cross=0.,azz_fsi1=0.,azz_fsi2=0.,azz_fsioff1=0.,azz_fsioff2=0.;  
      double f2_1,f2_2;
      for(int proton=0;proton<=1;++proton){
	
	double nu=Q2/(2.*(proton?MASSP:MASSN)*x);
	double qvec=sqrt(Q2+nu*nu);
	double Eout=Ein-nu;
	double thetae=2.*asin(sqrt(Q2/(4.*Ein*Eout)));
	double qz=Ein-Eout*cos(thetae);
	double qx=-Eout*sin(thetae);
	double thetaq=atan2(qx,qz); //angle between beam and qvec
	double thetapol=thetaq;
	if(!beampol) thetapol=0.;

	
	DQEinclusive Dinc(proton,ffparam,wf,*elec,4);
	vector<double> pw(4,0.);
	vector<double> fsi(8,0.);
	Dinc.calc_Crossinc(pw,Q2,x,2,integrator,thetapol);
// 	cout << "pwresult " << x << " " << proton << " " << f1 << " " << f2 << endl;
	pw_direct+=pw[0];
	pw_cross+=pw[1];
	azz_direct+=pw[2];
	azz_cross+=pw[3];
// 	Dinc.calc_CrossincFSI(fsi,Q2,x,2,integrator,maxEval,nopt,thetapol);
// 	cout << "fsiresult " << x << " " << proton << " " << f1 << " " << f2 << " " << f3 << " " << f4 << endl;
	fsi1+=fsi[0];
	fsi2+=fsi[1];
	fsioff1+=fsi[2];
	fsioff2+=fsi[3];
	azz_fsi1+=fsi[4];
	azz_fsi2+=fsi[5];
	azz_fsioff1+=fsi[6];
	azz_fsioff2+=fsi[7];
	(proton?f2_2:f2_1)=Dinc.getMinpcm();
//         if(pv) Dinc.calc_CrossincFSI_PVoff2(f1,f2,Q2,x,2,integrator,maxEval);
// 	else Dinc.calc_CrossincFSI_PVoff(f1,f2,Q2,x,2,integrator,maxEval);
//       cout << "fsiresultoff " << x << " " << proton << " " << f1 << " " << f2 << " " << endl;
//         fsioffPV1+=f1;
//         fsioffPV2+=f2;
	
      }
      cout << x << " " << pw_direct << " " << pw_cross << " " << pw_direct+pw_cross << " " << fsi1 << " " << fsi2 << " "
      << fsi1+fsi2 << " " << fsioff1 << " " << fsioff2 << " " << fsioff1+fsioff2 << " " 
      << fsioffPV1 << " " << fsioffPV2 << " " << fsioffPV1+fsioffPV2 << " "<<
      pw_direct+pw_cross+fsi1+fsi2+fsioff1+fsioff2 << " " 
      << azz_direct << " " << azz_cross << " " << azz_direct+azz_cross << " " << azz_fsi1 << " " << azz_fsi2 << " "
      <<  azz_fsi1+ azz_fsi2 << " " <<  azz_fsioff1 << " " <<  azz_fsioff2 << " " <<  azz_fsioff1+ azz_fsioff2 << " " 
      << azz_direct+azz_cross+ azz_fsi1+ azz_fsi2+ azz_fsioff1+ azz_fsioff2 << " " 
      
      << f2_1 << " " << f2_2 << endl;
  //     cout << x << " " << teller << " " << noemer << " " << teller/noemer << endl;
    
      delete elec;
    }
}
