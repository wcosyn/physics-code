#include <iostream>
#include <cstdlib>

using namespace std;

#include <RhoTCross.hpp>

int main(int argc, char *argv[])
{
  
  string homedir="/home/wim/Code/share";

  RhoTCross test = RhoTCross(1,300,homedir);
  double results[NROFRES];
//   test.getCrosst(results,5.014,1.4,2.970,-atof(argv[1]));
  test.getCrossz(results,5.014,1.4,2.970,0.95);
  
  return 0;
}


