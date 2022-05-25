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
    test.EMT_conv_real(atof(argv[1]));

    exit(1);


}
