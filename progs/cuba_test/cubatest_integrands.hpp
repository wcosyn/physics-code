
#ifndef CUBA_TEST_INTEGRANDS
#define CUBA_TEST_INTEGRANDS
#include <cmath>
#include <numint/numint.hpp>

#define SQRT_2PI 2.506628275

int my_integrand(const int* ndim, const double x[], const int *ncomp, double *f , void* userdata) {
	f[0] = x[0] ;//f[0] = 2.*PI*sin(2.*PI*x[0]); // single component function
	return 0;
}

int hypercube(const int* ndim, const double x[], const int*ncomp, double *f, void* userdata) {
	f[0] = 1.;
	return 0;
}

int hyper_lin(const int* ndim, const double x[], const int* ncomp, double *f, void* userdata) {
	f[0] = 0.;
	for (int i=0; i<*ndim; i++)
		f[0] += x[i];
	return 0;
}

int hyper_gauss(const int* ndim, const double x[], const int* ncomp, double *f, void* userdata) {
	double x2 = 0.;
	const double sigma = 0.01; // very peaked if we are in unit hypercube
	for (int i=0; i<*ndim; i++) 
		x2 += (2.*x[i]-1)*(2.*x[i]-1); // shift from 0,1 to -1,1
	f[0] = pow(2,*ndim)*pow( SQRT_2PI*sigma, -(*ndim))*exp(-0.5*x2/(sigma*sigma)); // factor in front is Jacobian of rescaling from 0,1 to -1,1
	return 0;
}

int unit_sphere(const int* ndim, const double x[], const int* ncomp, double *f, void* userdata) {
	double r2 = 0.;
	for (int i=0; i<*ndim; i++)
		r2 += (2.*x[i]-1.)*(2.*x[i]-1);
	f[0] = pow(2,*ndim)*((r2 < 1.)? 1. : 0.);
	return 0;
}

template<unsigned N>
void unit_sphere( const numint::array<double,N> &x, void* param, double &ret){
	double r = 0.;
	for (unsigned i=0; i<N; i++)
		r += x[i]*x[i];
	ret = (r > 1)? 0. : 1. ;
	//cout << "evaluated sphere @ " << x[0] << "," << x[1] << "," << x[2] << " r = " << r << " returned " << ret << endl;
}


#endif // CUBA_TEST_INTEGRANDS
