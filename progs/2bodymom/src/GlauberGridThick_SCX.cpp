#include "GlauberGridThick_SCX.hpp"
#include <gsl/gsl_sf_bessel.h>

/** integration struct 
 * ndim  = 2 (|b|,z)
 * ncomp = 2 real and imaginary piece
 *
 * **/
struct F_SCX {
	static void exec(const int ndim,const double x[],const int* ncomp, double* res,void *param){
		struct F_SCX p = *(struct F_SCX*) param;
		f(res,p.nuc,p.fp);
	}
	MeanFieldNucleusThick* nuc;
	FastParticle* fp;
	double b,z;
	int (*f)(double* res, const double x[], MeanFieldNucleusThick*,Fastparticle*,const double b,const double z);
};

int f_SCX_integrand(double* res, const double x[], MeanFieldNucleusThick* nuc, FastParticle* fp,const double b, const double z){
	const double bi = x[0]; // integration variable b'
	const double zi = x[1]; // integration variable z'
	const double r  = sqrt(bi*bi + zi*zi);
	const double betasq = fp->getBetaSq(...);
	const double dens = (fp->getParticletype()==8) ? nuc->getProtonDensity(r) : nuc->getNeutronDensity(r); // 8 is proton, so SCX with neutron density, vice versa for 9
	// if (r>nuc->getRange()){ ret=0; } else { ???
	*res = bi*exp(-bi*bi/(2.*betasq))*gsl_sf_bessel_I0(bi*b/betasq)*dens;
	return 0; // always "succes"
}
GlauberGridThick_SCX::GlauberGridThick_SCX(MeanFieldNucleusThick* nuc,FastParticle& fp,int bpoints,int zpoints) : 
	_nuc(nuc),
       	_fp(fp), 
	_bpoints(bpoints),
	_zpoints(zpoints),
	_bstep(0.),
	_zstep(0.),
	_grid(NULL)
	{
	_bstep = 2.*_nuc->getRange()/(_bpoints-1.); // 2*nuc.getRange because grid goes from
	_zstep = 2.*_nuc->getRange()/(_zpoints-1.); // -nuc.getRange() to +nuc.getRange()
	constructGlauberGrid();
}

GlauberGridThick_SCX::~GlauberGridThick_SCX(){
	for (int bi=0; bi<_bpoints;bi++)
		delete[] _grid[bi];
	delete[] _grid;
}

void GlauberGridThick_SCX::constructGlauberGrid(){
	assert(_grid==NULL && _bstep > 0. && _zstep > 0.);
	//std::cout << "constructing grid " << _bpoints << " by " << _zpoints << std::endl;
	//std::cout << "steps are " << _bstep << "  " << _zstep << std::endl;
	_grid = new complex<double>*[_bpoints];	
	for (int bi=0; bi<_bpoints; bi++){
		_grid[bi] = new complex<double>[_bpoints];
		for (int zi=0; zi<_zpoints; zi++){
			double b = -_nuc->getRange()+bi*_bstep; 
			double z = -_nuc->getRange()+zi*_zstep;
			_grid[bi][zi] = calcFSI(b,z);
			//std::cout << " [" << bi << ", " << zi << "] corresponding to " << b << ", " << z << " init to " << _grid[bi][zi] << std::endl;
		}
	}
}

void GlauberGridThick_SCX::printGrid(){ // debug function
	for (int bi=0; bi<_bpoints;bi++){
		for (int zi=0; zi<_zpoints;zi++){
			double b = -_nuc->getRange()+bi*_bstep; 
			double z = -_nuc->getRange()+zi*_zstep;
			std::cout << b << "\t" << z << "\t" << _grid[bi][zi] << std::endl;
		}
	}
}

/**calculate the glauber grid
 * will need the proton density for neutron knockout
 * will need the neutron density for proton knockout
 */
complex<double> GlauberGridThick_SCX::calcFSI(double b,double z){
	
	return complex<double>(0.,0.);	
}

