#include <iostream>
#include <cstdlib>

using namespace std;

#include <RhoTCross.hpp>

int main(int argc, char *argv[])
{

  int torz = atoi(argv[1]);
  int nucleus = atoi(argv[2]);
  double Q2 = atof(argv[3]);
  double nu_min = atof(argv[4]);
  double nu_max = atof(argv[5]);
  double Ebeam = 5.014;
  int j = atoi(argv[6]);
  int k = atoi(argv[7]);
  
  double *results = new double[NROFRES];

  RhoTCross test = RhoTCross(nucleus,400,argv[8]);
  
  double nu = nu_min + j*(nu_max-nu_min)/7.;
  if(torz){
    double t = -0.1 + k*(-0.1);
    test.getCrosst(results,Ebeam,Q2,nu,t);
  }
  else{
    double z = 0.9+k*0.03;
    test.getCrossz(results,Ebeam,Q2,nu,z);
  }
  cout << Q2 << " "  << nu << " " << (torz? -0.1 + k*(-0.1) :  0.9+k*0.03) << " " ;
  for (int l=0;l<NROFRES;l++) cout << results[l] << " ";
  cout << endl;


  delete [] results;
  
  
//   string homedir="/home/wim/Code/share";
// 
//   RhoTCross test = RhoTCross(1,400,homedir);
//   double results[NROFRES];
// //   test.getCrosst(results,5.014,1.4,2.970,-atof(argv[1]));
//   double Ebeam = 5.014;
//   double Q2=1.4;
//   double nu=2.970;
//   double z=0.95;
//   test.getCrossz(results,Ebeam,Q2,nu,z);
//   for(int i=0;i<NROFRES;i++) cout << Q2 << " " << nu << " " << z << " " << results[i] << endl;
//   return 0;
}


