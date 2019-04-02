#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;


#include <iostream>
#include <cmath>
#include "constants.hpp"

#include "GPD.hpp"
#include "TransGPD_set.hpp"

int main(int argc, char *argv[]){
    std::string wf=argv[2];
    GPD test=GPD("MSTW",wf);
    int model=atoi(argv[1]); //chiral odd GPD model nr
    double xi=atof(argv[4]);
    double t=atof(argv[3])*-1.E06;
    cout << "xi " << xi << " t " << t << " model " << model << endl;
    cout << -4.*MASSD*MASSD*xi*xi/(1-xi*xi)-t << endl;
    for(int i=0;i<5;i++){
        vector< complex<double> > hel(5,0.);
        hel[i]=1.;
        vector< complex<double> > gpd = GPD::helamps_to_gpds_V(xi,t,hel);
        vector< complex<double> > out = GPD::gpds_to_helamps_V(xi,t,gpd);
        for(int j=0;j<5;j++) cout << out[j] << " ";
        cout << endl;
    }
    exit(1);

    for(int i=-99;i<=99;i++){
        double x=i*0.01;
        vector< complex<double> > out = test.gpd_conv(xi,x,t,model);
        vector< complex<double> > gpd = GPD::helamps_to_gpds_T(xi,t,out);
        vector< complex<double> > hel = GPD::gpds_to_helamps_T(xi,t,gpd);
        
        cout << x << " ";
        for(int j=0;j<9;j++) cout << out[j].real() << " ";
        for(int j=0;j<9;j++) cout << gpd[j].real() << " ";
        cout << endl;

        
    }    
}