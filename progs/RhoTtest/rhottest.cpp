#include <iostream>
#include <cstdlib>

using namespace std;

#include <RhoTCross.hpp>
#include <TMPI.h>
#include <stdio.h>

int main(int argc, char *argv[])
{

  TMPI mpi(&argc,&argv);
  int torz = atoi(argv[1]);
  int nucleus = atoi(argv[2]);
  double Q2 = atof(argv[3]);
  double nu_min = atof(argv[4]);
  double nu_max = atof(argv[5]);
  double Ebeam = 5.014;
  
  int n=32; // number of calculations
  double **results = new double*[n];  
  for(int i=0;i<n;i++) results[i] = new double[NROFRES];

  RhoTCross test = RhoTCross(nucleus,400,argv[6]);
  
  for(int i=TMPI::Rank(); i<n; i += TMPI::NumberOfProcesses() ) {
    int j = i/4;
    int k = i%4;
    double nu = nu_min + j*(nu_max-nu_min)/7.;
    if(torz){
      double t = -0.1 + k*(-0.1);
      test.getCrosst(results[i],Ebeam,Q2,nu,t);
    }
    else{
      double z = 0.9+k*0.03;
      test.getCrossz(results[i],Ebeam,Q2,nu,z);
    }
//     cout << Q2 << " "  << nu << " " << (torz? -0.1 + k*(-0.1) :  0.9+k*0.03) << " " ;
//     for (int l=0;l<NROFRES;l++) cout << results[i][l] << " ";
//     cout << endl;
  }

  TMPI::GatherResults(n,results);

  TMPI::SilenceSlaves();
  for(int i=0; i<n; ++i){
    int j = i/4;
    int k = i%4;
    double nu = nu_min + j*(nu_max-nu_min)/7.;    
    cout << Q2 << " "  << nu << " " << (torz? -0.1 + k*(-0.1) :  0.9+k*0.03) << " " ;
    for (int l=0;l<NROFRES;l++) cout << results[i][l] << " ";
    cout << endl;
  }
  TMPI::SilenceSlaves(false);

  for(int i=0;i<n;i++) delete [] results[i];
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


