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
    for(int i=-10;i<=50;i++){
        double x=i*0.01;
        vector< complex<double> > out = test.gpd_conv(xi,x,t,0);
        cout << x << " ";
        for(int j=0;j<9;j++) cout << out[j] << " ";
        cout << endl;
    }    
}