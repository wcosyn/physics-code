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


  //unpolarized cross section test
  // double s=2.*1000.E06;
  // double x=0.05;
  // double Q2=35.E06;
  // double alpha_p=1.;
  // for(int i=1;i<10;i++){
  //   double tprime=-0.005*i*1.E06;
  //   Poldeut cross=Poldeut("AV18","CTEQ");
  //   cout << tprime << " " ;
  //   double result = cross.calc_dsigma_unpol(x,Q2,s,alpha_p,tprime,0);
  //   cout << result << " " << result*PI*MASSD/4./sqrt(((-MASSD*(2.*MASSn-MASSD)-pow(2.*MASSn-MASSD,2.)/2.)-tprime)/2.+MASSn*MASSn) << endl;
  // }


  //ALL test
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
