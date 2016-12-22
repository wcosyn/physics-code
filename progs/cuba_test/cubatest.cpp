#include <Cuba/cuba.h>
#include <cstdio>
#include "cubatest_integrands.hpp"

void test_vegas();
void test_cuhre();
void test_divonne();
void test_suave();

int main()
{
	test_cuhre();
	test_divonne();
	test_suave();
	test_vegas();
	printf("4d sphere true integral : %f\n", 0.5*PI*PI);
	return 0;
}

void test_vegas() {
	int ndim = 4;
	int ncomp =2;
	//integrand_t integrand = &hyper_lin; //hypercube; //&my_integrand;
	integrand_t integrand = &unit_sphere;
	void* userdata = NULL;
	int nvec = 1;
	double epsrel = 1e-4;
	double epsabs = 1e-4;
	int flags = 0b00000000; // see manual what each bit means
	int seed = 9875661;
	int mineval = 10000;
	int maxeval = 100000;
	int nstart = 1000;
	int nincrease = 10000;
	int nbatch = 20000;
	int gridno = 0;
	char statefile[123] = {"vegas.state"};
	int neval = 0;
	int fail = 0;
	double* integral = new double[ncomp];
	double* error = new double[ncomp];
	double* prob = new double[ncomp];
	Vegas(ndim,ncomp,integrand,userdata,nvec,
		epsrel,epsabs,flags,seed,
		mineval,maxeval,nstart,nincrease,
		nbatch,gridno,statefile,&neval,
		&fail,integral,error,prob);
    printf("Vegas  : neval = %d ,fail = %d\n",neval,fail);
    for (int nc=0;nc<ncomp;nc++){
        printf(" >  component[%d] , result = %6.2e, error = %6.2e, fail = %.2f\n",nc,integral[nc],error[nc],prob[nc]);
    }
    printf("\n");
    delete[] integral;
    delete[] error;
    delete[] prob;
}


/** note that cuhre fails miserably for
  * 1D integrals, is OK, was never intended
  * for 1D integrals
  */

void test_cuhre() {
	int ndim = 4;
	int ncomp = 2;
	//integrand_t integrand = &hyper_lin; //hypercube; // &my_integrand;
	//integrand_t integrand = &hyper_gauss;
	integrand_t integrand = &unit_sphere;
	void* userdata = NULL;
	int nvec = 1;
	double epsrel = 1e-4;
	double epsabs = 1e-4;
	int flags = 0x00;
	int mineval = 10000;
	int maxeval = 100000;
	int key = 9;
	char statefile[] = {"cuhre.state"};
	int nregions = 0;
	int neval = 0;
	int fail = 0;
	double* integral = new double[ncomp];
	double* error = new double[ncomp];
	double* prob = new double[ncomp];
	
	Cuhre(ndim,ncomp,integrand,userdata,nvec,
		epsrel,epsabs,flags,
		mineval,maxeval,key,
		statefile,&nregions,&neval,&fail,
		integral,error,prob);
	printf("Cuhre  : neval = %d, fail = %d \n",neval,fail);
    for (int nc=0;nc<ncomp;nc++){
        printf(" >  component[%d] , result = %6.2e, error = %6.2e, fail = %.2f\n",nc,integral[nc],error[nc],prob[nc]);
    }
    printf("\n");
    delete[] integral;
    delete[] error;
    delete[] prob;
}

void test_divonne() {
	int ndim = 4;
	int ncomp = 5;
	integrand_t integrand = &unit_sphere;
	void *userdata=NULL;
	int nvec = 1;
	double epsrel=1e-4;
	double epsabs=1e-4;
	int flags = 0x00;
	int seed  = 1235488;
	int mineval = 10000;
	int maxeval = 100000;
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
	char statefile[] = {"divonne.state"};
	int nregions = 0;
	int neval = 0;
	int fail = 0;
	double* integral = new double[ncomp];
	double* error = new double[ncomp];
	double* prob = new double[ncomp];
	
	Divonne(ndim,ncomp,integrand,userdata,nvec,
		epsrel,epsabs,flags,seed,
		mineval,maxeval,key1,key2,key3,
		maxpass,border,maxchisq,mindeviation,
		ngiven,ldxgiven,xgiven,nextra,peakfinder,
		statefile,&nregions,&neval,&fail,
		integral,error,prob);

	printf("Divonne: neval = %d, fail = %d \n",neval,fail);
    for (int nc=0;nc<ncomp;nc++){
        printf(" >  component[%d] , result = %6.2e, error = %6.2e, fail = %.2f\n",nc,integral[nc],error[nc],prob[nc]);
    }
    printf("\n");
    delete[] integral;
    delete[] error;
    delete[] prob;
}

void test_suave(){
	int ndim = 4;
	int ncomp = 2;
	integrand_t integrand = &unit_sphere; //&hypercube;
	void* userdata = NULL;
	int nvec=1;
	double epsrel=1e-4;
	double epsabs=1e-4;
	int flags = 0x00;
	int seed = 6667;
	int mineval = 10000;
	int maxeval = 100000;
	int nnew = 5000;
	double flatness = 0.25;
	char statefile[100] = {"suave.state"};
	int nregions=0;
	int neval=0;
	int fail=0;
	double* integral = new double[ncomp];
	double* error = new double[ncomp];
	double* prob = new double[ncomp];

	Suave(ndim,ncomp,integrand,userdata,nvec,
		epsrel,epsabs,flags,seed,mineval,
		maxeval,nnew,flatness,statefile,
		&nregions,&neval,&fail,integral,
		error,prob);

	printf("Suave  : neval = %d, fail = %d \n",neval,fail);
    for (int nc=0;nc<ncomp;nc++){
        printf(" >  component[%d] , result = %6.2e, error = %6.2e, fail = %.2f\n",nc,integral[nc],error[nc],prob[nc]);
    }
    printf("\n");
    delete[] integral;
    delete[] error;
    delete[] prob;
}
