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
    GPD test=GPD("MSTW","AV18");

    double xi=0.1;
    double t=-0.25E06;
    for(int i=-99;i<=99;i++){
        double x=i*0.01;
        vector< complex<double> > out = test.gpd_conv(xi,x,t,1);
        vector< complex<double> > gpd = GPD::helamps_to_gpds(xi,t,out);
        cout << x << " ";
        for(int j=0;j<9;j++) cout << out[j].real() << " ";
        for(int j=0;j<9;j++) cout << gpd[j].real() << " ";
        cout << endl;
        
    }    
}