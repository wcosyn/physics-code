//test program for the light-front polarized deuteron code
// run like >azztest
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <Poldeut.hpp>


int main(int argc, char *argv[])
{

  
  
  
  Poldeut test("AV18","SLAC");
  for(int i=0;i<=50;i++){
      for(int j=0;j<=50;j++){
        double alpha_p = 0.8+0.4*i/50;
        double pt=j*6;

        double p_p_plus = alpha_p/2.*MASSD;
        double p_p_minus = (MASSn*MASSn+pt*pt)/p_p_plus;
        double p_z = (p_p_plus - p_p_minus)/2.;
        double pnorm = sqrt(p_z*p_z+pt*pt);

        double lf_ratio=0.,nr_ratio=0.;

        test.Tensor_Compare_nonrel(lf_ratio, nr_ratio, alpha_p,pt);

        cout << alpha_p << " " << pt << " " << p_z << " " << pnorm << " " << lf_ratio << " " << nr_ratio << " " << lf_ratio/nr_ratio << endl;


      }
    
  }
}
