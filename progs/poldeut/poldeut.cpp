#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp.in>
#include <Poldeut.hpp>


int main(int argc, char *argv[])
{

  double x=0.3;
  double Q2=10.E06;
  double y=0.5;
  double alpha_s=1.;
//   double pt=50.;
  bool proton=0;
  double F_LSL,F_U;
  
  
  
  Poldeut test("AV18","Alekhin");
  test.set_melosh(0);
  test.set_D(0);
  for(int i=1;i<30;i++){
    double pt=10.+i*20;
    pt *=-1.;
    test.calc_Double_Asymm(x,Q2,y,alpha_s,pt,proton, F_LSL, F_U);
  }
}