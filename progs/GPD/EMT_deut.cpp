#include <string>
#include <vector>


#include <constants.hpp>
#include <Deut_Conv_EMTLF.hpp>

using namespace std;


/*
 * Main function.
 */
int main(int argc, char** argv) {

    Deut_Conv_EMTLF test("AV18");
    double t = -1.E-04;
    //test.EMT_conv(-1.E-09);
    //cout << endl;
    //test.EMT_conv(-1.);
//    for(int i=0;i<20;i++){
//	test.EMT_conv_real(t);
//	t*=2.;
//	}
     for(int i=9;i<25;i++){
         test.EMT_conv_real(-1.-i);
    //     t*=2;
     }
    // for(int i=0;i<25;i++){
    //     t=-4-i;
    //     test.EMT_conv_real(t);
    // }
    return 1;
}
