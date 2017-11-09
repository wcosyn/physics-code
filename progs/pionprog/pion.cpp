//small program that outputs glauber parameters for pions

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include <FastParticle.hpp>
int main(int argc, char *argv[])
{
  double sigmap[151],sigman[151], beta2p[151], beta2n[151], epsp[151], epsn[151];
  string homedir="/home/wim/Code/share";
  for(int i=0;i<=150;i++){
    double logp = 2.5+i*0.01;
    FastParticle::setPionGlauberData(pow(10.,logp), sigmap[i],beta2p[i],epsp[i],sigman[i],beta2n[i],epsn[i],homedir);
    cout << pow(10.,logp) << " " << logp << " " << sigmap[i] << " " << beta2p[i] << " " << epsp[i] << " " << sigman[i] << " " << beta2n[i] << " " << epsn[i] << endl;
  }
  
//   cout << "{ " ;
//   for(int i=0;i<150;i++) cout << sigmap[i] << (i%10==9? ",\n" : ", ");  
//   cout << sigmap[150] << "} " << endl << endl;
//  
//   cout << "{ " ;
//   for(int i=0;i<150;i++) cout << beta2p[i] << (i%10==9? ",\n" : ", ");  
//   cout << beta2p[150] << "} " << endl << endl;
//  
//   cout << "{ " ;
//   for(int i=0;i<150;i++) cout << epsp[i] << (i%10==9? ",\n" : ", ");  
//   cout << epsp[150] << "} " << endl << endl;
//  
//   cout << "{ " ;
//   for(int i=0;i<150;i++) cout << sigman[i] << (i%10==9? ",\n" : ", ");  
//   cout << sigman[150] << "} " << endl << endl;
//  
//   cout << "{ " ;
//   for(int i=0;i<150;i++) cout << beta2n[i] << (i%10==9? ",\n" : ", ");  
//   cout << beta2n[150] << "} " << endl << endl;
//  
//   cout << "{ " ;
//   for(int i=0;i<150;i++) cout << epsn[i] << (i%10==9? ",\n" : ", ");  
//   cout << epsn[150] << "} " << endl << endl;
 
  
}
      