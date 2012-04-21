#include <iostream>
#include <cstdlib>

using namespace std;

#include <RhoTCross.hpp>
#include <TMPI.h>
#include <stdio.h>

int main(int argc, char *argv[])
{

  TMPI mpi(&argc,&argv);

  int n=8; // number of calculations
  double **x = new double*[n];
  for(int i=0;i<n;i++) x[i] = new double[NROFRES];

  for(int i=TMPI::Rank(); i<n; i += TMPI::NumberOfProcesses() ) {
    for(int j=0;j<NROFRES;j++) {
      x[i][j]=i+j;
    }
  }

  TMPI::GatherResults(n,x);

  TMPI::SilenceSlaves();
  for(int i=0; i<n; ++i){
    for (int j=0;j<NROFRES;j++) cout << x[i][j] << " ";
    cout << endl;
  }
  TMPI::SilenceSlaves(false);

  for(int i=0;i<n;i++) delete [] x[i];
  delete [] x;
  
  
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


