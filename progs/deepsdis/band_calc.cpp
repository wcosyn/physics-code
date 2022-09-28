//program that calculates cross sections to do a FSI estimation for the BAND data

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <DeuteronStructure.hpp>
#include <DeuteronMomDistr.hpp>
#include <DeuteronCross.hpp>
#include <He3Cross.hpp>
#include <LightConeKin2to3.hpp>
#include <NuclStructure.hpp>

int main(int argc, char *argv[])
{

    double xb = atof(argv[1]);
    double Q2 = atof(argv[2])*1.E06;  //input GeV2 -> in code MeV2!!
    double alpha_s = atof(argv[3]);  // ratio of minus components in this convention
    double pT = atof(argv[4])*1.E03; //input GeV -> in code MeV!
    //double xprime = atof(argv[5]);

    int proton = 1;
    int offshellset = 4;
    int looplimit=1;

    double sigmain = 40.;
    double betain = 8.;

    DeuteronCross test("paris",proton,"SLAC",sigmain,betain,-0.5,8.,1.2,offshellset,looplimit);

    double psmin = alpha_s*MASSD/2.;
    double psplus = 2.*(pT*pT+MASSN*MASSN)/alpha_s/MASSD;
    double Es = 0.5*(psmin+psplus);
    double psz = 0.5*(psplus-psmin);

    double theta = atan2(pT,psz);

    double ps = sqrt(psz*psz+pT*pT);

    double nu = Q2/(2.*xb*MASSP);
    cout << "nu " << nu << endl;
    double qvec = sqrt(Q2+nu*nu);

    double Eioff = MASSD-Es;
    double massoff2 = Eioff*Eioff-ps*ps;

    double Wsq = pow(nu+Eioff,2.)-qvec*qvec-ps*ps+2.*qvec*psz;

    double xprime = Q2/(Wsq-massoff2+Q2);

    cout << ps*1.E-03 << " "<< theta*RADTODEGR << " " << sqrt(Wsq)*1.E-03 << " " << xprime << endl;
    sigmain= DeuteronCross::sigmaparam(Wsq,Q2);  // get sigma parameter
    //cout << sigmain*10.*HBARC*HBARC << endl;
    cout << 53.444 << endl;
    test.setScatter(53.444,8.,-0.5);
    for (int i=0;i<40;i+=1){
        double costhetar=-0.975+i*0.05;
        double pw=0.,fsi=0.;
        test.getDeepsresult(Q2,sqrt(Wsq),10.2E03,450,costhetar,proton,pw,fsi);
        cout << costhetar << " " << pw << " " << fsi << " " << fsi/pw << endl;
    }
}
