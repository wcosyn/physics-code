#include <Cuba/cuba.h>
#include <cstdio>
#include "cubatest_integrands.hpp"

void test_vegas  (int ndim,int ncomp,integrand_t integrand);
void test_cuhre  (int ndim,int ncomp,integrand_t integrand);
void test_divonne(int ndim,int ncomp,integrand_t integrand);
void test_suave  (int ndim,int ncomp,integrand_t integrand);

int main()
{
    integrand_t integrand = unit_sphere; printf("#[Info] Integrand is unit_sphere\n");
	test_cuhre(4,2,integrand);
	test_divonne(4,2,integrand);
	test_suave(4,2,integrand);
	test_vegas(4,2,integrand);
	printf("4d sphere true integral : %f\n", 0.5*PI*PI);
	return 0;
}

void test_vegas(int ndim,int ncomp,integrand_t integrand) {
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
	char* statefile = nullptr; // {"vegas.state"};
	int neval = 0;
	int fail = 0;
	double* integral = new double[ncomp];
	double* error = new double[ncomp];
	double* prob = new double[ncomp];
	Vegas(ndim,ncomp,integrand,userdata,nvec,
		epsrel,epsabs,flags,seed,
		mineval,maxeval,nstart,nincrease,
		nbatch,gridno,statefile,NULL,&neval,
		&fail,integral,error,prob);
    printf("Vegas  : neval = %d ,fail = %d\n",neval,fail);
    for (int nc=0;nc<ncomp;nc++){
        printf(" >  component[%d] , result = %6.2e, error = %6.2e, prob = %.2f\n",nc,integral[nc],error[nc],prob[nc]);
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

void test_cuhre(int ndim,int ncomp,integrand_t integrand) {
	void* userdata = NULL;
	int nvec = 1;
	double epsrel = 1e-4;
	double epsabs = 1e-4;
	int flags = 0x00;
	int mineval = 10000;
	int maxeval = 100000;
	int key = 9;
	char* statefile = nullptr; // {"cuhre.state"};
	int nregions = 0;
	int neval = 0;
	int fail = 0;
	double* integral = new double[ncomp];
	double* error = new double[ncomp];
	double* prob = new double[ncomp];
	
	Cuhre(ndim,ncomp,integrand,userdata,nvec,
		epsrel,epsabs,flags,
		mineval,maxeval,key,
		statefile,NULL,&nregions,&neval,&fail,
		integral,error,prob);
	printf("Cuhre  : neval = %d, fail = %d \n",neval,fail);
    for (int nc=0;nc<ncomp;nc++){
        printf(" >  component[%d] , result = %6.2e, error = %6.2e, prob = %.2f\n",nc,integral[nc],error[nc],prob[nc]);
    }
    printf("\n");
    delete[] integral;
    delete[] error;
    delete[] prob;
}

void test_divonne(int ndim,int ncomp,integrand_t integrand) {
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
	int ngiven = 4;
	int ldxgiven = ndim;
	double* xgiven = new double[ngiven*ldxgiven];
    for (int nx=0;nx<ngiven;nx++){
        for (int nxd=0;nxd<ldxgiven;nxd++){
            xgiven[ldxgiven*nx+nxd] = 0.;
        }
    }
	int nextra = 0;
	peakfinder_t peakfinder = 0;
	char* statefile = nullptr; //{"divonne.state"};
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
		statefile,NULL,&nregions,&neval,&fail,
		integral,error,prob);

	printf("Divonne: neval = %d, fail = %d \n",neval,fail);
    for (int nc=0;nc<ncomp;nc++){
        printf(" >  component[%d] , result = %6.2e, error = %6.2e, prob = %.2f\n",nc,integral[nc],error[nc],prob[nc]);
    }
    printf("\n");
    delete[] integral;
    delete[] error;
    delete[] prob;
    delete[] xgiven;
}

void test_suave(int ndim,int ncomp,integrand_t integrand){
	void* userdata = NULL;
	int nvec=1;
	double epsrel=1e-4;
	double epsabs=1e-4;
	int flags = 0x00;
	int seed = 6667;
	int mineval = 10000;
	int maxeval = 100000;
	int nnew = 5000;
	int nmin = 2;
	double flatness = 0.25;
	char* statefile = nullptr; // {"suave.state"};
	int nregions=0;
	int neval=0;
	int fail=0;
	double* integral = new double[ncomp];
	double* error = new double[ncomp];
	double* prob = new double[ncomp];

	Suave(ndim,ncomp,integrand,userdata,nvec,
		epsrel,epsabs,flags,seed,mineval,
		maxeval,nnew,nmin,flatness,statefile,NULL,
		&nregions,&neval,&fail,integral,
		error,prob);

	printf("Suave  : neval = %d, fail = %d \n",neval,fail);
    for (int nc=0;nc<ncomp;nc++){
        printf(" >  component[%d] , result = %6.2e, error = %6.2e, prob = %.2f\n",nc,integral[nc],error[nc],prob[nc]);
    }
    printf("\n");
    delete[] integral;
    delete[] error;
    delete[] prob;
}
