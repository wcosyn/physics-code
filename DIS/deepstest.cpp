#include <iostream>
#include <cstdlib>
#include <cmath>


using namespace std;

#include <constants.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <DeuteronStructure.hpp>
#include "DeuteronMomDistr.hpp"
#include "DeuteronCross.hpp"

int main(int argc, char *argv[])
{
  double Wprime=1.25E03;//invariant mass X
  double Q2=1.8E06;
  double Ein=5765.;
  for (int i=0;i<6;i+=2){
    double costhetar=-0.975+i*0.05;
  
  
    string homedir="/home/wim/Code";
    double pr=300;
    int proton=0;
    double phi=0.;
    
    double massi=proton? MASSP:MASSN;
    double massr=proton? MASSN:MASSP;
    
    double Er=sqrt(massr*massr+pr*pr);
    double Einoff=MASSD-Er;
    double massoff=sqrt(Einoff*Einoff-pr*pr);
    
    double xprime=Q2/(Wprime*Wprime-massoff*massoff+Q2);
    
    //calc nu
    double prz=pr*costhetar;
    double aaa=Einoff*Einoff-prz*prz;
    double bbb=-Einoff*Q2/xprime;
    double ccc=Q2*Q2/(4.*xprime*xprime)-Q2*prz*prz;
    
    double discr=sqrt(bbb*bbb-4.*aaa*ccc);
    double nu1=(-bbb+discr)/(2.*aaa);
    double nu2=(-bbb-discr)/(2.*aaa);
    //cout << nu1 << " " << nu2 << endl;
    double nu=nu2;
    if(costhetar<0.) nu=nu1;
    double qvec=sqrt(Q2+nu*nu);
  /*  double xx=Q2/2./(Einoff*nu+prz*qvec);
    cout << xprime << " " << xx << endl;*/
    double x=Q2/(2.*massi*nu);
      
    TKinematics2to2 kin("","",MASSD,MASSP,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
    TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
    DeuteronCross test(*elec,"paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
    
    double Eout=Ein-nu;
    double thetain=acos((Ein*Ein+qvec*qvec-Eout*Eout)/(2.*Ein*qvec));
    double yprime=(Einoff*nu+prz*qvec)/(Ein*Einoff+Ein*cos(thetain)*pr*costhetar/*+Ein*sin(thetain)*pr*sqrt(1.-costhetar*costhetar)*/);  
    double R=0.18;
    double frontdeeps=(4.*PI*ALPHA*ALPHA)/(xprime*Q2*Q2)*
    (yprime*yprime/(2.*(1.+R))+(1.-yprime)+(pow(massoff*xprime*yprime,2.)*(1.-R))/(Q2*(1+R)));
    double dxprimedx=-2.*xprime*xprime*nu/(x*qvec)*((Einoff+prz)/(nu-qvec)+1./(2.*xprime));
    
    
    cout << costhetar << " " << test.getavgCross(kin,1,Einoff)/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10
     << " " << test.getavgCross(kin,0,Einoff)/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10
     << endl;
    delete elec;
  }
}
