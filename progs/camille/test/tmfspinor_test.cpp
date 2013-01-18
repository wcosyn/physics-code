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
#include <cassert>
void printArray(MeanFieldNucleus&, const int*);

int main()
{
	string homedir = "/home/camille/Wim_Code/share";
	// -- MeanFieldNucleus
	MeanFieldNucleus nuc = MeanFieldNucleus(7,homedir);
	// print some test info
	cout << nuc.getNucleusName() <<" nucleus, with " << nuc.getZ() << " protons and " << nuc.getN() << " neutrons " << endl;
	cout << "Total number of shells " << nuc.getTotalLevels() << " of which " << nuc.getPLevels() << " are proton shells and " << nuc.getNLevels() << " are neutron shells " << endl;
	cout << "n array " << endl; printArray(nuc,nuc.getN_array());
	cout << "k array " << endl; printArray(nuc,nuc.getKappas());
	cout << "l array " << endl; printArray(nuc,nuc.getL_array());
	cout << "j array " << endl; printArray(nuc,nuc.getJ_array()); // returns j_max times two per level!
	int shellindex = 12;
	//int m = 1; // m_j times two, invalid values will give nans in calculations
	//double r = 1.; // [fm] i presume
	double costheta = 0.;
	double phi = 0.;
	ofstream phi_dat("tmfspinor_sq.dat");
	for (double r=0.; r<12.0; r+=0.01){
		phi_dat << r;
		for (int m=-nuc.getJ_array()[shellindex]; m <= nuc.getJ_array()[shellindex]; m+=2)
		{
			TMFSpinor phi1 = TMFSpinor(nuc,shellindex,m,r,costheta,phi);
			Matrix<1,4> phi1Bar = phi1.H(); // hermitian conjugate
			std::complex<double> phisq = phi1Bar*phi1;
			phi_dat << "\t" << phisq.real();
			if (phisq.real() < 0.)
			{
				std::cerr << " abs squared should not be negative: " << phisq.real() << endl;
				std::cerr << phi1 << endl;		
			}
			assert(phisq.imag() < 1e-13); // imaginary part should be zero
		}
		phi_dat << endl;
	}
	//TMFSpinor phi2 = TMFSpinor(nuc,shellindex,m,r,costheta,phi); // remember pauli exclusion !
	phi_dat.close();
	return 0;
}

void printArray(MeanFieldNucleus &nuc,const int* arr)
{
	cout << " proton shells: " << endl;
	for (int i=0; i<nuc.getPLevels(); i++)
	{
		cout << "   "<< i << ": " << arr[i] << endl;
	}
	cout << " neutron shells: " << endl;
	for (int i=0; i<nuc.getNLevels(); i++)
	{
		cout << "   " << i << ": " << arr[i] << endl;
	}
}
