/** \file
 * \author Camille Colle
 * This is a simple test file to test the construction of TMFSpinors
 * and the integrator lib. We try to make some TMFSpinors and integrate
 * the wavefunctions to check the normalisation
 */

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
#include <numint/numint.hpp>

void printArray(MeanFieldNucleus&, const int*);
void myIntegrandumFunction(double & res, const numint::array<double,3> &x, MeanFieldNucleus& nuc, int n, int m);
struct myStruct {
	static void exec(const numint::array<double,3> &x, void *param, double &ret) { // static: all instances of myStruct will share this function
		myStruct &p = * (myStruct *) param; // cast the param to a myStruct, a bit weird to pass a struct as a parameter to a struct of the same type?
		p.f(ret,x, *p.nuc,p.n,p.m); // now call the f function with of the passed struct with the parameters of the passes struct
	}
	int n,m;
	MeanFieldNucleus* nuc;
	void (*f)(double & res, const numint::array<double,3> &x, MeanFieldNucleus& nuc, int n, int m); // NOT a void pointer, but a pointer to a void function
};



int main()
{
	string homedir = "/home/camille/Wim_Code/share";
	// -- MeanFieldNucleus
	MeanFieldNucleus nuc = MeanFieldNucleus(2,homedir);
	// print some test info
	cout << nuc.getNucleusName() <<" nucleus, with " << nuc.getZ() << " protons and " << nuc.getN() << " neutrons " << endl;
	cout << "Nucleus range " << nuc.getRange() << endl;
	cout << "Total number of shells " << nuc.getTotalLevels() << " of which " << nuc.getPLevels() << " are proton shells and " << nuc.getNLevels() << " are neutron shells " << endl;
	cout << "n array " << endl; printArray(nuc,nuc.getN_array());
	cout << "k array " << endl; printArray(nuc,nuc.getKappas());
	cout << "l array " << endl; printArray(nuc,nuc.getL_array());
	cout << "j array " << endl; printArray(nuc,nuc.getJ_array()); // returns j_max times two per level!
	int shellindex = 0;
	int m = 1; // m_j times two, invalid values will give nans in calculations, for example, if you want l=0, m=1 because s=1/2 -> m=1/2 ->2*m=m=1
	// now make the weird function struct thingy to test the integration lib with tmf wfns -> e.g. test normalisation
	numint::array<double,3> lowerb = {{0.,-1.,0.}}; // lower r, costheta, phi bounds	
	//numint::array<double,3> upperb = {{1.,1.,2.*PI}}; // SPHERE TEST: upper r, costheta, phi bounds
	numint::array<double,3> upperb = {{nuc.getRange(),1.,2.*PI}}; // upper r, costheta, phi bounds
	numint::mdfunction<double,3> mymdf; // mdf ([m]ulti[d]imensional[f]unction; in the brackets : <double,..> is the return variable type <...,unsigned N> is the dimensionality of the function
	myStruct myStructInstance;
	mymdf.func = &myStruct::exec; // this can be done because exec is static, so exec does not need an actual struct and can be called scopewise
	mymdf.param = &myStructInstance; // because the variables are stored in the struct itself
	myStructInstance.nuc = &nuc; // now init all the parameters of our struct before we try to integrate -> otherwise segfaults will be coming our way
	myStructInstance.n = shellindex; // shellindex
	myStructInstance.m = m; // m_j times two
	myStructInstance.f = myIntegrandumFunction; // function we wish to integrate
	double prec = 1e-10, ret;
	int maxEval = 10000;
	unsigned int count;

	//int succes = numint::cube_adaptive(mymdf,lowerb,upperb,1.E-08,prec,maxEval,ret,count,0);
	//cout << " succes? : " << succes << " result : " << ret << endl;
	
	// now lets try looping over all the levels and checking the normalisation for all the mfspinors
	
	for (int n=0; n<nuc.getTotalLevels(); n++)
	{
		myStructInstance.n = n;
		for (int m = -nuc.getJ_array()[n]; m <= nuc.getJ_array()[n]; m+=2)
		{
			cout << "Integrating TMFSpinor with shellindex : " << n << " 2*m_j : " << m;
			myStructInstance.m = m;
			int succes = numint::cube_adaptive(mymdf,lowerb,upperb,1.E-10,prec,1E02,maxEval,ret,count,0);
			cout << " succes(=0)? : " << succes << " norm. result : " << ret << endl;
		}
	}
	return 0;
}

void myIntegrandumFunction(double &res, const numint::array<double,3> &x, MeanFieldNucleus& nuc, int n, int m)
{
	TMFSpinor phi = TMFSpinor(nuc,n,m,x[0],x[1],x[2]);
	res = (phi.H()*phi).real();	
	//res = x[0]*x[0]; //SPHERE TEST: x[0] = r, because x[1] is cos(theta) the factor sin(theta) in the jacobian is already taken into account
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
