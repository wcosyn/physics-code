#include <iostream>
using std::cout;
using std::endl;

#include <TMFSpinor.hpp>
#include <TSpinor.h>
#include <FourVector.h>
#include <constants.hpp>
#include <MeanFieldNucleus.hpp>
#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
#include "GammaStructure.h"
#include <fstream>
using std::ofstream;


int main()
{
	string homedir = "/home/camille/Wim_Code/share";
	// -- MeanFieldNucleus
	MeanFieldNucleus Carbon = MeanFieldNucleus(4,homedir);
	int shellindex = 7;
	//int m = 1; // m_j times two, invalid values will give nans in calculations
	//double r = 2.; // [fm] i presume
	//double costheta = 0.;
	//double phi = 0.;
	ofstream phi_dat("tmfspinor2D_sq.dat");
	double rstop = 4.;
	for (double r=0.001; r < rstop; r+= 0.05)
	{
		cout << r << " to " << rstop << endl;
		for (double phi=0.; phi<2*PI+0.1; phi+=0.1){
			for (double theta=0.; theta < PI+0.15; theta += 0.15)
			{
				phi_dat << r << "\t" << phi << "\t" << theta;
				for (int m=-Carbon.getJ_array()[shellindex]; m <= Carbon.getJ_array()[shellindex]; m+=2)
				{
					TMFSpinor phi1 = TMFSpinor(Carbon,shellindex,m,r,cos(theta),phi);
					Matrix<1,4> phi1Bar = phi1.H(); // hermitian conjugate
					std::complex<double> phisq = phi1Bar*phi1;
					phi_dat << "\t" << phisq.real();
				}
				phi_dat << endl;
			}
		}
	}
	phi_dat.close();
	return 0;
}

