#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp.in>
#include <Poldeut.hpp>


int main(int argc, char *argv[])
{

  double x=0.05;
  double Q2=20.E06;
  double y=0.4;
  double alpha_s=atof(argv[1]);
  x*=2-alpha_s;
//   double pt=50.;
  bool proton=0;
  double F_LSL,F_U;
  
  
  
  Poldeut test(argv[4],argv[5]);
   test.set_melosh(atoi(argv[3]));
  test.set_D(atoi(argv[2]));
  for(int i=0;i<30;i++){
    double pt=0.+i*20;
//     pt *=-1.;
//     test.set_D(1);
    test.calc_Double_Asymm(x,Q2,y,alpha_s,pt,proton, F_LSL, F_U);
//     test.set_D(0);
//     test.calc_Double_Asymm(x,Q2,y,alpha_s,pt,proton, F_LSL, F_U);
    
  }
}
