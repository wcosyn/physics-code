#include <vector>
#include <iostream>
#include <vector>
#include <algorithm>

#include <constants.hpp>
#include <Deut_Conv_EMTLF.hpp>
#include <Utilfunctions.hpp>

using namespace std;


/*
 * Main function.
 */
int main(int argc, char** argv) {

    std::cout << "Compiled from file: " << __FILE__ << std::endl;
    string arg_names[1]={""};
    Bookkeep(argc,argv,arg_names);  

    Deut_Conv_EMTLF test("AV18");

    std::vector<double> arr;
    arr.push_back(-1.0E-6);
    // 1) -1e-4, double 18 times
    double x = -1.0e-4;
    for (int i = 0; i < 19; i++) {
        arr.push_back(x);
        x *= 2.0;
    }

    // 2) -0.05 to -0.525 in steps of -0.025
    for (double v = -0.05; v >= -0.525 - 1e-12; v -= 0.025) {
        arr.push_back(v);
    }

    // 3) -6 to -25 in steps of -1
    for (int v = -6; v >= -25; v--) {
        arr.push_back((double)v);
    }

    // 4) -0.75 to -5.25 in steps of -0.25
    for (double v = -0.75; v >= -5.25 - 1e-12; v -= 0.25) {
        arr.push_back(v);
    }

    // Sort from largest to smallest (less negative to more negative)
    std::sort(arr.begin(), arr.end(), std::greater<double>());
    //cout << arr.size() << endl;
    cout << "// t [GeV^2] A_00	A_+0	A_-+	A_0+	A_++	J_00	J_+0	J_-+	J_0+	J_++	tD_00	tD_+0	tD_-+	tD_0+	tD_++" << endl;
    for (size_t i = 0; i < arr.size(); ++i) {
        double t = arr[i];
        //cout << t <<endl;
        test.EMT_conv_real(t);
    }
	return 1;
}
