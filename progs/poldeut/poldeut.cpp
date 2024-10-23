//test program for the light-front polarized deuteron code
// run like >poldeut [alpha] [D-wave] [Melosh] [wf] [strucfunc]
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <Poldeut.hpp>
#include <Utilfunctions.hpp>
#include <NuclStructure.hpp>

int main(int argc, char *argv[])
{


// cout << "output is for neutron SF: Q2[GeV^2]|x|F1|F2|FT|FL|g1|g1+g2 (WW approx)" << endl;
// int neutron=0;
// for(int i=0;i<10;i++){
//   double Q2= 2.+i;
//   for(int j=1;j<=100;j++){
//     double x_nucl=0.002*j;
//     NuclStructure nucl(neutron,Q2*1.E06,x_nucl,0,"SLAC");
//     double F1,F2,g1,g1pg2,FL,FT;
//     g1=NuclStructure::getG1_grsv2000(neutron,x_nucl,Q2*1.E06);
//     g1pg2=NuclStructure::getG1plusg2_grsv2000(neutron,x_nucl,Q2*1.E06);
//     nucl.getF(F1,F2);
//     double gamma_n = 2.*x_nucl*MASSn_G/sqrt(Q2);
//     F1=(1+gamma_n*gamma_n)/(1.+0.18)/(2.*x_nucl)*F2;
//     FT=2.*F1;
//     FL=(1+gamma_n*gamma_n)*F2/x_nucl-2.*F1;
//     cout << Q2 << " " << x_nucl << " " << F1 << " " << F2 << " " << FT << " " << FL << " " << g1 << " " << g1pg2 << endl;
//   }

// }
// exit(1);




std::string wf = argv[1];
std::string arg_names[1]={"wf name"};
std::cout << "Compiled from file: " << __FILE__ << std::endl;
Bookkeep(argc,argv,arg_names);  


// //plotting lf deuteron densities

// double S_L=0., S_T=0., deltaS_L=0., deltaS_T=0.,deltaT_S_L=0.,deltaT_S_T=0.;  //MeV^-1
// Poldeut deuteron_obj("AV18","SLAC");

// for(int i=-60;i<=60;i++){
//   for(int j=0;j<=60;j++){
//     double alpha_p = 1.+i/100.;
//     double pt = double(j)*10.;
//     deuteron_obj.getLFdistributions(alpha_p,pt,S_L,S_T,deltaS_L,deltaS_T,deltaT_S_L,deltaT_S_T);

//     cout << alpha_p << " " << pt*1.E-03 << " " << S_L*1.E06 << " " << S_T*1.E06 << " " << deltaS_L*1.E06 << " " << deltaS_T*1.E06 << " " << deltaT_S_L*1.E06 << " " << deltaT_S_T*1.E06 << endl;

//   }
// }
// exit(1);

//plotting lf deuteron densities

double P_U=0., P_TLL=0., P_TLT=0., P_TTT=0.,P_SLSL=0.,P_STSL=0.,P_SLST=0.,P_STST=0.;  //MeV^-2
Poldeut deuteron_obj(wf,"SLAC");
cout << "// output is (all densities GeV-2) alpha|pt[MeV]|P_U|P_TLL|P_TLT|P_TTT|P_SLSL|P_STSL|P_SLST|P_STST" << endl;

for(int i=-60;i<=60;i++){
  for(int j=0;j<=300;j++){
    double alpha_p = 1.+i/100.;
    double pt = double(j)*2.;
    deuteron_obj.getLFdistributions_new(alpha_p,pt,P_U,P_TLL,P_TLT,P_TTT,P_SLSL,P_STSL,P_SLST,P_STST);

    cout << alpha_p << " " << pt*1.E-03 << " " << P_U*1.E06 << " " << P_TLL*1.E06 << " " << P_TLT*1.E06 << 
            " " << P_TTT*1.E06 << " " << P_SLSL*1.E06 << " " << P_STSL*1.E06 << " " << P_SLST*1.E06 << " " << P_STST*1.E06 << endl;

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
