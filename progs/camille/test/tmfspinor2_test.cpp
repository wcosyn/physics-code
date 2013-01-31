/** \file
 * \author Camille Colle
 * This is a simple test file to test the construction of TMFSpinors
 * and the integrator lib. We try to make some TMFSpinors and integrate
 * the wavefunctions to check orthonormalisation
 */

#include <iostream>
using std::cout;
using std::endl;

#include <TMFSpinor.hpp>
#include <TSpinor.h>

#include <constants.hpp>
#include <MeanFieldNucleus.hpp>
#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
#include <numint/numint.hpp>
#include <iomanip>
using std::setw;

void printArray(MeanFieldNucleus&, const int*);
void myIntegrandumFunction(double & res, const numint::array<double,3> &x, MeanFieldNucleus& nuc, int n1, int m1, int n2, int m2);
struct myStruct {
	static void exec(const numint::array<double,3> &x, void *param, double &ret) { // static: all instances of myStruct will share this function
		myStruct &p = * (myStruct *) param; // cast the param to a myStruct, a bit weird to pass a struct as a parameter to a struct of the same type?
		p.f(ret,x, *p.nuc,p.n1,p.m1,p.n2,p.m2); // now call the f function with of the passed struct with the parameters of the passes struct
	}
	int n1,m1,n2,m2;
	MeanFieldNucleus* nuc;
	void (*f)(double & res, const numint::array<double,3> &x, MeanFieldNucleus& nuc, int n1, int m1,int n2, int m2); // NOT a void pointer, but a pointer to a void function
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
	// now make the weird function struct thingy to test the integration lib with tmf wfns -> e.g. test normalisation
	numint::array<double,3> lowerb = {{0.,-1.,0.}}; // lower r, costheta, phi bounds	
	//numint::array<double,3> upperb = {{1.,1.,2.*PI}}; // SPHERE TEST: upper r, costheta, phi bounds
	numint::array<double,3> upperb = {{nuc.getRange(),1.,2.*PI}}; // upper r, costheta, phi bounds
	numint::mdfunction<double,3> mymdf; // mdf ([m]ulti[d]imensional[f]unction; in the brackets : <double,..> is the return variable type <...,unsigned N> is the dimensionality of the function
	myStruct myStructInstance;
	mymdf.func = &myStruct::exec; // this can be done because exec is static, so exec does not need an actual struct and can be called scopewise
	mymdf.param = &myStructInstance; // because the variables are stored in the struct itself
	myStructInstance.nuc = &nuc; // now init all the parameters of our struct before we try to integrate -> otherwise segfaults will be coming our way
	myStructInstance.f = myIntegrandumFunction; // function we wish to integrate
	double prec = 1e-10, ret;
	int maxEval = 10000;
	unsigned int count;

	//int succes = numint::cube_adaptive(mymdf,lowerb,upperb,1.E-08,prec,maxEval,ret,count,0);
	//cout << " succes? : " << succes << " result : " << ret << endl;
	
	// now lets try looping over all the levels and checking the orthonormality for all the mfspinors
	
	for (int n1=0; n1<nuc.getTotalLevels(); n1++)
	{
		myStructInstance.n1 = n1;
		for (int m1 = -nuc.getJ_array()[n1]; m1 <= nuc.getJ_array()[n1]; m1+=2)
		{
			myStructInstance.m1 = m1;
			for (int n2 = n1; n2 < nuc.getTotalLevels(); n2++)
			{
				myStructInstance.n2 = n2;
				for (int m2 = m1; m2 <= nuc.getJ_array()[n2]; m2+=2)
				{
					myStructInstance.m2 = m2;
					cout << "Integrating (n1,l1,m1)*(n2,l2,m2) : (" << setw(2) << nuc.getN_array()[n1] << ", " << setw(2) << nuc.getL_array()[n1] << ", " << setw(2) << m1 << " )*(" <<  nuc.getN_array()[n2] << ", " << setw(2) << nuc.getL_array()[n2] <<  ",  " << setw(2) <<m2 << ")   ";
				int succes = numint::cube_adaptive(mymdf,lowerb,upperb,1.E-10,prec,maxEval,ret,count,0);
				cout << " succes(=0)? : " << setw(2) <<succes << " result : " << ret << endl;
				}
			}
		}
	}
	return 0;
}

void myIntegrandumFunction(double &res, const numint::array<double,3> &x, MeanFieldNucleus& nuc, int n1, int m1,int n2, int m2)
{
	TMFSpinor phi1 = TMFSpinor(nuc,n1,m1,x[0],x[1],x[2]);
	TMFSpinor phi2 = TMFSpinor(nuc,n2,m2,x[0],x[1],x[2]);
	res = (phi2.H()*phi1).real();
	//res = (phi.H()*phi).real();	
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
