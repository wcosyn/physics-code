#include <iostream>
#include "Model.hpp"
#include <TKinematics2to2.h>
#include <MeanFieldNucleus.hpp>
#include <string>
using std::string;
using std::cout;
using std::endl;
#include <constants.hpp>

int main()
{
	string homedir = "/home/camille/Wim_Code/share/";

	MeanFieldNucleusThick Carbon(1,homedir);
	cout << "total shell levels of nucleus: " << Carbon.getTotalLevels() << endl;
	double Q2 = 10000; // in MeV
	double omega = 20000; // what to take here?
	double pm = 100; // cutoff @ 300?
	TKinematics2to2 kin("","",Carbon.getMassA(),Carbon.getMassA_min_proton(),MASSP,"qsquared:wlab:pklab",Q2,omega,pm);
	double prec = 1e-05;
	int integrator = 2; // 1 = cube_romb, can converge very slowly
	int max_Eval = 10000; // max evaluations of function to compute integral
	bool user_sigma = false; // true to override default sigma
	Model mod = Model(&Carbon,prec,integrator,homedir,max_Eval,user_sigma);
	Matrix<2,3> res;
	int shellindex=3;
	int m = 1; // m_j times 2!
	int CT = 0; // color transparency
	int pw = 1; // plane wave
	int current = 1; // see T. de Forest Nucl Phys A 392,232 (1983) for CC1, CC2, or CC3
	int SRC = 0; // short range correlations
	int thick = 0; // thickness?
	mod.getMatrixEl(kin,res,shellindex,m,CT,pw,current,SRC,thick,1.,1.);
	cout << res << endl;
	cout << "Matrix element for nucleon helicity -1 " << endl;
	cout << "> photon pol =  0 " << res(0,0) << endl;
	cout << "> photon pol = -1 " << res(0,1) << endl;
	cout << "> photon pol = +1 " << res(0,2) << endl;
	cout << "Matrix element for nucleon helicity +1 " << endl;
	cout << "> photon pol =  0 " << res(1,0) << endl;
	cout << "> photon pol = -1 " << res(1,1) << endl;
	cout << "> photon pol = +1 " << res(1,2) << endl;
	return 0;
}
