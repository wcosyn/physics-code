#include <iostream>
#include <cstdlib>

using namespace std;

#include <RhoTCross.hpp>

int main(int argc, char *argv[])
{
  
  string homedir="/home/wim/Code/share";

  RhoTCross test = RhoTCross(1,400,homedir);
  double results[NROFRES];
//   test.getCrosst(results,5.014,1.4,2.970,-atof(argv[1]));
  double Ebeam = 5.014;
  double Q2=1.4;
  double nu=2.970;
  double z=0.95;
  test.getCrossz(results,Ebeam,Q2,nu,z);
  for(int i=0;i<NROFRES;i++) cout << Q2 << " " << nu << " " << z << " " << results[i] << endl;
  return 0;
}


