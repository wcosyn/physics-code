#include <iostream>
#include <cstdlib>

using namespace std;


#include <iostream>
#include <cmath>
#include "constants.hpp"

#include "GPD.hpp"
#include "TransGPD_set.hpp"

int main(int argc, char *argv[]){
    GPD test=GPD("MSTW");

    for(int i=-100;i<=100;i++){
        TransGPD_set out=test.getGK_param(0.01*i,0.201914,-0.158217E06);
        cout << 0.01*i << " " << out.getHTd() << " " << out.getHTu() << " " << out.getEbarTd() << " " << out.getEbarTu() << endl;
    }
    
}