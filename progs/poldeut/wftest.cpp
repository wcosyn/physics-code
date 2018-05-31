#include<OldDeuteron.hpp>



#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;


#include <iostream>
#include <cmath>
#include "constants.hpp"

int main(int argc, char *argv[]){
    OldDeuteron wf1("AV18");
    OldDeuteron wf2("AV18b");
    OldDeuteron wf3("CDBonn");
    OldDeuteron wf4("Paris");
    
    double sum1=0.,sum2=0.,sum3=0.,sum4=0.;

    //normalization test
    for(int i=0; i<=1000; i++){
        sum1+=i*i*(pow(wf1.U(i),2.)+pow(wf1.W(i),2.));
        sum2+=i*i*(pow(wf2.U(i),2.)+pow(wf2.W(i),2.));
        sum3+=i*i*(pow(wf3.U(i),2.)+pow(wf3.W(i),2.));
        sum4+=i*i*(pow(wf4.U(i),2.)+pow(wf4.W(i),2.));
    }

    cout << sum1*4.*PI << " " << sum2*4.*PI << " " << sum3*4.*PI << " " << sum4*4.*PI << endl;
    
}