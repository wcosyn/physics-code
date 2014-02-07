
#include <stdio.h>
#include "cubatest_integrands.hpp"
#include <vector>
using std::vector;

void test_vegas_comp_arr(); // test multicomponent integration using numint::array to store components
void test_vegas_comp_vec(); // test multicomponent integration using vector<...> to store components

int main(){
	printf("\nRunning using numint::arrays to store multicomponent function\n");
	test_vegas_comp_arr();
	printf("\nRunning using vector<double> to store multicomponent function\n");
	test_vegas_comp_vec();
	printf("4d sphere true integral : %f (first 3 components should be equal to this) \n", 0.5*PI*PI);
	printf("4d cube   true integral : %f (all components above 3 should be equal to this) \n", 16.0); // 16 because each side is 2
}

void test_vegas_comp_arr() {
	const int ndim  = 4;
	const int ncomp = 6;
	numint::mdfunction<numint::array<double,ncomp> ,ndim> mymdf;
	mymdf.func = mixed_sphere_cube_comp<numint::array<double,ncomp>,ndim>;
	//mymdf.func = unit_sphere_comp<numint::array<double,ncomp>,ndim>;
	
	numint::array<double,ndim> start = {{-1,-1,-1,-1}};
	numint::array<double,ndim> stop  = {{ 1, 1, 1, 1}};
	
	// **** INPUT VARIABLES ****** //
	int nvec = 1;// unless you know what SIMD is leave this 1
	double epsrel = 1e-4;
	double epsabs = 1e-4;
	int flags = 0b00000000;
	int seed = 9875661;
	int mineval = 10000;
	int maxeval = 1000000;
	int nstart = 1000;
	int nincrease = 10000;
	int nbatch = 50000;
	int gridno = 0;
	char statefile[] = {""};
	
	// **** OUTPUT VARIABLES ****** //
	int neval=0;
	int fail=0;
	numint::array<double,ncomp> integral;
	numint::array<double,ncomp> error;
	numint::array<double,ncomp> prob;
	
	vegas(mymdf,start,stop,nvec,epsrel,epsabs,flags,seed,
				mineval,maxeval,nstart,nincrease,nbatch,gridno,
				statefile,neval,fail,integral,error,prob);
	for (unsigned i=0; i<integral.size(); i++)
		printf("Vegas  : res = %12f, error = %e, prob = %f, neval = %d ,fail = %d\n",integral[i],error[i],prob[i],neval,fail);
}


void test_vegas_comp_vec() {
	const int ndim  = 4;
	const int ncomp = 6;
	numint::mdfunction<vector<double>, ndim> mymdf;
	mymdf.func = mixed_sphere_cube_comp<vector<double>, ndim>;
	//mymdf.func = unit_sphere_comp<vector<double>, ndim>;
	
	numint::array<double,ndim> start = {{-1,-1,-1,-1}};
	numint::array<double,ndim> stop  = {{ 1, 1, 1, 1}};
	
	// **** INPUT VARIABLES ****** //
	int nvec = 1;// unless you know what SIMD is leave this 1
	double epsrel = 1e-4;
	double epsabs = 1e-4;
	int flags = 0b00000000;
	int seed = 9875661;
	int mineval = 10000;
	int maxeval = 1000000;
	int nstart = 1000;
	int nincrease = 10000;
	int nbatch = 50000;
	int gridno = 0;
	char statefile[] = {""};
	
	// **** OUTPUT VARIABLES ****** //
	int neval=0;
	int fail=0;
	vector<double> integral(ncomp);
	vector<double> error(ncomp);
	vector<double> prob(ncomp);
	
	vegas(mymdf,start,stop,nvec,epsrel,epsabs,flags,seed,
				mineval,maxeval,nstart,nincrease,nbatch,gridno,
				statefile,neval,fail,integral,error,prob);
	for (unsigned i=0; i<integral.size(); i++)
		printf("Vegas  : res = %12f, error = %e, prob = %f, neval = %d ,fail = %d\n",integral[i],error[i],prob[i],neval,fail);
}
