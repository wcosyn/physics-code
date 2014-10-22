#include "cubatest_integrands.hpp"
#include <numint/numint.hpp>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;

void test_cuhre();
void test_divonne();
void test_suave();
void test_vegas();
void test_vegas_comp();

int main(){
	test_cuhre();
	test_divonne();
	test_suave();
	test_vegas();
	printf("4d sphere true integral : %f\n", 0.5*PI*PI);
	return 0;
}

void test_suave() {
	const int ndim = 4;
	numint::mdfunction<double,ndim> mymdf;
	mymdf.func = unit_sphere<ndim>;
	numint::array<double,ndim> start = {{-1,-1,-1,-1}};
	numint::array<double,ndim> stop  = {{ 1, 1, 1, 1}};
	
	// **** INPUT VARIABLES ****** //
	int nvec=1;
	double epsrel=1e-4;
	double epsabs=1e-4;
	int flags = 0x00;
	int seed = 6667;
	int mineval = 10000;
	int maxeval = 1000000;
	int nnew = 5000;
	double flatness = 0.25;
	char* statefile = NULL; //{"suave.state"};
	
	// **** OUTPUT VARIABLES ****** //
	int nregions    = 0;
	int neval       = 0;
	int fail        = 0;
	double integral = 0.;
	double error    = 0.;
	double prob     = 0.;
	
	suave(mymdf,start,stop,nvec,epsrel,epsabs,flags,seed,
				mineval,maxeval,nnew,flatness,
				statefile,nregions,neval,fail,integral,error,prob);
	
	printf("Suave  : res = %f, error = %e, prob = %f, neval = %d ,fail = %d\n",integral,error,prob,neval,fail);
}

void test_vegas() {
	const int ndim = 4;
	numint::mdfunction<double,ndim> mymdf;
	mymdf.func = unit_sphere<ndim>;
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
	char* statefile = NULL; //{"vegas.state"};
	
	// **** OUTPUT VARIABLES ****** //
	int neval=0;
	int fail=0;
	double integral;
	double error;
	double prob;
	
	vegas(mymdf,start,stop,nvec,epsrel,epsabs,flags,seed,
				mineval,maxeval,nstart,nincrease,nbatch,gridno,
				statefile,neval,fail,integral,error,prob);
	
	printf("Vegas  : res = %f, error = %e, prob = %f, neval = %d ,fail = %d\n",integral,error,prob,neval,fail);
}

void test_divonne(){
	const int ndim = 4;
	numint::mdfunction<double,ndim> mymdf;
	mymdf.func = unit_sphere<ndim>;
	numint::array<double,ndim> start = {{-1,-1,-1,-1}};
	numint::array<double,ndim> stop  = {{ 1, 1, 1, 1}};
	
	// **** INPUT VARIABLES ****** //
	int nvec = 1; // unless you know what SIMD is leave this 1
	double epsrel=1e-4;
	double epsabs=1e-4;
	int flags = 0x00;
	int seed  = 1235488;
	int mineval = 10000;
	int maxeval = 1000000;
	int key1 = 1;
	int key2 = 1;
	int key3 = 1;
	int maxpass = 50 ;
	double border = 0.;
	double maxchisq = 5;
	double mindeviation = 0.1;
	int ngiven = 0;
	int ldxgiven = 0;
	double* xgiven = NULL;
	int nextra = 0;
	peakfinder_t peakfinder = 0;
	char* statefile = NULL; //{"divonne.state"};
	
	// **** OUTPUT VARIABLES ****** //
	int nregions=0;
	int neval=0;
	int fail=0;
	double integral = 0.; // init to 0. is not needed but valgrind will warn us that this is pointing to uninitialized bytes if we don't
	double error = 0.;
	double prob = 0.;
	divonne(mymdf,start,stop,nvec,
		epsrel,epsabs,flags,seed,
		mineval,maxeval,key1,key2,key3,
		maxpass,border,maxchisq,mindeviation,
		ngiven,ldxgiven,xgiven,nextra,peakfinder,
		statefile,nregions,neval,fail,
		integral,error,prob);
	
	printf("Divonne: res = %f, error = %e, prob = %f, neval = %d ,fail = %d\n",integral,error,prob,neval,fail);
}

void test_cuhre(){
	const int ndim = 4;
	numint::mdfunction<double,ndim> mymdf;
	mymdf.func = unit_sphere<ndim>;
	numint::array<double,ndim> start = {{-1,-1,-1,-1}};
	numint::array<double,ndim> stop  = {{ 1, 1, 1, 1}};
	
	// **** INPUT VARIABLES ****** //
	int nvec = 1; // unless you know what SIMD is leave this 1
	double epsabs = 1e-4;
	double epsrel = 1e-4;
	int flags = 0x00;
	int minEval = 10000;
	int maxEval = 1000000;
	int key=11;
	char* statefile = NULL; //{"cuhre.state"};
	// **** OUTPUT VARIABLES ****** //
	int nregions =0;
	int neval    =0;
	int fail     =0;
	double res   =0.;
	double err   =0.;
	double prob  =0.;
	cuhre( mymdf,start,stop,nvec,epsrel,epsabs,flags,minEval,maxEval,key,statefile,
			nregions,neval,fail,res,err,prob);
	
	printf("Cuhre  : res = %f, error = %e, prob = %f, neval = %d ,fail = %d\n",res,err,prob,neval,fail);
}

