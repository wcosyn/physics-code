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

void printArray(MeanFieldNucleus&, const int*);

int main()
{
	string homedir = "/home/camille/Wim_Code/share";
	/*
	// -- Vectors
	TVector3 q,pf,pm,Pcm,ps;
	q = TVector3(0.,0.,1); // q along z axis
	pf = 0.8*q;
	pf.RotateY(DEGRTORAD*5); // pf now in xz plane
	cout << pf.X() << "\t" << pf.Y() << "\t" <<  pf.Z() << "\t" << endl;
	pm = pf-q; // also in xz plane, duh
	Pcm = TVector3(0.,0.,400.); // should become grid, one dir along q (z axis), one dir perp on pf_q plane (y axis), and one perp to prev 2 (z-axis)
	ps = Pcm - pm;
	// -- Spinors
	FourVector<double> f_pm = FourVector<double>(sqrt(MASSP*MASSP + pm.Dot(pm)) ,pm.X(),pm.Y(),pm.Z());
	Matrix<1,4> u_pmBar = TSpinor::Bar(TSpinor(f_pm,MASSP,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity)); // with choices above for axis, how to choose polarization of spinors ??
	FourVector<double> f_ps = FourVector<double>(sqrt(MASSP*MASSP + ps.Dot(ps)), ps.X(), ps.Y(), ps.Z());
	Matrix<1,4> u_psBar = TSpinor::Bar(TSpinor(f_ps,MASSP,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity));*/
	// -- MeanFieldNucleus
	MeanFieldNucleus Carbon = MeanFieldNucleus(1,homedir);
	// print some test info
	cout << "Carbon nucleus, with " << Carbon.getZ() << " protons and " << Carbon.getN() << " neutrons " << endl;
	cout << "Total number of shells " << Carbon.getTotalLevels() << " of which " << Carbon.getPLevels() << " are proton shells and " << Carbon.getNLevels() << " are neutron shells " << endl;
	cout << "n array " << endl; printArray(Carbon,Carbon.getN_array());
	cout << "k array " << endl; printArray(Carbon,Carbon.getKappas());
	cout << "l array " << endl; printArray(Carbon,Carbon.getL_array());
	cout << "j array " << endl; printArray(Carbon,Carbon.getJ_array()); // returns j_max times two per level!
	int shellindex = 0;
	//int m = 1; // m_j times two, invalid values will give nans in calculations
	//double r = 1.; // [fm] i presume
	double costheta = 0.;
	double phi = 0.;
	ofstream phi_dat("tmfspinor_sq.dat");
	for (double r=0.; r<12.0; r+=0.01){
		phi_dat << r;
		for (int m=-Carbon.getJ_array()[shellindex]; m <= Carbon.getJ_array()[shellindex]; m+=2)
		{
			TMFSpinor phi1 = TMFSpinor(Carbon,shellindex,m,r,costheta,phi);
			Matrix<1,4> phi1Bar = phi1.H()*G0; // hermitian conjugate * gamma_0
			std::complex<double> phisq = phi1Bar*phi1;
			phi_dat << "\t" << phisq.real();
		}
		phi_dat << endl;
	}
	//TMFSpinor phi2 = TMFSpinor(Carbon,shellindex,m,r,costheta,phi); // remember pauli exclusion !
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
