//test program for the light-front polarized deuteron code
// run like >poldeut [alpha] [D-wave] [Melosh] [wf] [strucfunc]
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <Poldeut.hpp>


int main(int argc, char *argv[])
{


//plotting lf deuteron densities

double S_L=0., S_T=0., deltaS_L=0., deltaS_T=0.,deltaT_S_L=0.,deltaT_S_T=0.;  //MeV^-1
Poldeut deuteron_obj("AV18","SLAC");

for(int i=-60;i<=60;i++){
  for(int j=0;j<=60;j++){
    double alpha_p = 1.+i/100.;
    double pt = double(j)*10.;
    deuteron_obj.getLFdistributions(alpha_p,pt,S_L,S_T,deltaS_L,deltaS_T,deltaT_S_L,deltaT_S_T);

    cout << alpha_p << " " << pt*1.E-03 << " " << S_L*1.E06 << " " << S_T*1.E06 << " " << deltaS_L*1.E06 << " " << deltaS_T*1.E06 << " " << deltaT_S_L*1.E06 << " " << deltaT_S_T*1.E06 << endl;

  }
}
exit(1);




  // unpolarized cross section test
  double s=2.*1000.E06;
  double x=0.05;
  double Q2=20.E06;
  double y=0.4;
 // double alpha_s=atof(argv[1]);
//   double pt=50.;
  bool proton=0;
  double F_LSL,F_U;
  
  
  
  Poldeut test(argv[4],argv[5]);
   test.set_melosh(atoi(argv[3]));
  test.set_D(atoi(argv[2]));
  for(int i=-30;i<30;i++){
    //double pt=0.+i*20;
    double pt=0.;
    double alpha_s=1.+i/60.;
  x=0.05*(2-alpha_s);
  // cout << alpha_s << " " << x << endl;
//     pt *=-1.;
//     test.set_D(1);
    test.calc_Double_Asymm(x,Q2,y,alpha_s,pt,proton, F_LSL, F_U);
//     test.set_D(0);
//     test.calc_Double_Asymm(x,Q2,y,alpha_s,pt,proton, F_LSL, F_U);
// //     test.set_D(0);
// //     test.calc_Double_Asymm(x,Q2,y,alpha_s,pt,proton, F_LSL, F_U);
    
  }
  cout << endl << endl; 

  double alpha_s=1.;
  x=0.05*(2.-alpha_s);
  for(int i=0;i<30;i++){
    double pt=0.+i*20;
//     pt *=-1.;
//     test.set_D(1);
    test.calc_Double_Asymm(x,Q2,y,alpha_s,pt,proton, F_LSL, F_U);
  }

}
