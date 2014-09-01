#include "GlauberGridThick_SCX.hpp"
#include <gsl/gsl_sf_bessel.h>

/** integration struct 
 * ndim  = 2 (|b|,z)
 * ncomp = 2 real and imaginary piece
 *
 * **/
struct F_SCX {
	static int exec(const int* ndim,const double x[],const int* ncomp, double* res,void *param){
		struct F_SCX p = *(struct F_SCX*) param;
		int l = p.f(res,x,p.nuc,p.fp,p.b,p.z);
		//std::cout << " Passing result " << *res << " to integrator " << std::endl;
		return l;
		//return p.f(res,x,p.nuc,p.fp,p.b,p.z);
	}
	MeanFieldNucleusThick* nuc;
	FastParticle* fp;
	double b,z;
	int (*f)(double* res, const double x[], MeanFieldNucleusThick*,FastParticle*,const double b,const double z);
};

int f_SCX_integrand(double* res, const double x[], MeanFieldNucleusThick* nuc, FastParticle* fp,const double b, const double z){
	double bi = x[0]; // integration variable b', here in unit square
	double zi = x[1]; // integration variable z', here in unit square
	bi = -nuc->getRange()+bi*(2.*nuc->getRange()); // map [0,1] -> [-R,R]
	zi = z+zi*(nuc->getRange()-z); // map [0,1] -> [z,R] (z instead of -R because of step function HS(z-z')
	//std::cout << " translated unit coords " << x[0] << " and " << x[1] << " to world coords " << bi << " and " << zi << std::endl;
	const double Jac = 2.*nuc->getRange()*(nuc->getRange()-z); // Jacobian of mapping
	//std::cout << " z is " << z << std::endl;
	//std::cout << " Jacobian of transformation is " << Jac << std::endl;
	const double r  = sqrt(bi*bi + zi*zi);
	const double betasq = (fp->getParticletype()==8) ? fp->getBeta2n() : fp->getBeta2p();
	//std::cout << " betasq is " << betasq << std::endl;
	//if (r>nuc->getRange()){ // outside the range, density is (practically) zero, so return zero
	//	*res=0;
	//} else {
		//double dens   = (fp->getParticletype()==8) ? nuc->getNeutronDensity(r) : nuc->getProtonDensity(r); // 8 is proton, so SCX with neutron density, vice versa for 9
		//dens /= (r>1e-15) ? dens/r/r : dens/(r+1e-12)/(r+1e-12); // second thing to prevent division by zero
		//*res = dens;
		//*res = Jac*bi*exp(-bi*bi/(2.*betasq))*gsl_sf_bessel_I0(bi*b/betasq)*dens; // dens divided by r^2 because dens is rho*r*r
		//double dens = r > 1.0 ? 1.0 : 0.0;
		//*res = Jac*1.; //exp(-(bi-b)*(bi-b)/0.0001);
		//*res = fabs(b) < 4. ?  exp(-pow((bi-b)/4.,2))*Jac*1. : 0.;
		if (fabs(b) < 4.)
			*res = Jac*1.;
		else if (fabs(b) < 8.)
			*res = Jac*exp(-pow(bi-4.,2)/4.);
		else
			*res = 0.;
		//std::cout << " exp factor is " << exp(-bi*bi/(2.*betasq)) << " bessel is " <<  gsl_sf_bessel_I0(bi*b/betasq) << " argument in bessel is " << bi*b/betasq << " dens is " << dens << " result is " << *res << std::endl;
		//std::cout << " additional factor would be approx " << exp(-b*b/(2.*betasq)) << std::endl;
	//}
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
	assert(fp.getParticletype()==8 || fp.getParticletype()==9); // make sure we are getting the glauber params for p/n charge exchange!
	_bstep = 2.*_nuc->getRange()/(_bpoints-1.); // 2*nuc.getRange because grid goes from
	_zstep = 2.*_nuc->getRange()/(_zpoints-1.); // -nuc.getRange() to +nuc.getRange()
	std::cout << "#Particle type is " << fp.getParticletype() << " Beta2p is " << fp.getBeta2p() << " Beta2n is " << fp.getBeta2n() << std::endl;
	std::cout << "#Sigmap is " << fp.getSigmap() << "  Sigman is " << fp.getSigman() << std::endl;
	std::cout << "#Scatterfront is (proton): " << fp.getScatterfront(1) << "  (neutron): " << fp.getScatterfront(0) << std::endl;
	std::cout << "#bsteps, zstep is " << _bstep << ", " << _zstep << " and range is " << nuc->getRange() << std::endl;
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
			std::cout << b << "\t" << z << "\t" << _grid[bi][zi].real() << "\t" << _grid[bi][zi].imag() << std::endl;
		}
	}
}

void GlauberGridThick_SCX::printDensityGrid(){
	for (int bi=0; bi<_bpoints;bi++){
		for (int zi=0; zi<_zpoints;zi++){
			double b = -_nuc->getRange()+bi*_bstep; 
			double z = -_nuc->getRange()+zi*_zstep;
			double r = sqrt(b*b+z*z);
			if (r==0.)
				r=1e-12; // to prevent division by zero!
			std::cout << b << "\t" << z << "\t" ;
			if (r < _nuc->getRange())
				std::cout << _nuc->getProtonDensity(r)/r/r << "\t" << _nuc->getNeutronDensity(r)/r/r << std::endl;
			else
				std::cout << 0. << "\t" << 0. << std::endl;
		}
	}
}

/**calculate the glauber grid
 * will need the proton density for neutron knockout
 * will need the neutron density for proton knockout
 */
complex<double> GlauberGridThick_SCX::calcFSI(double b,double z){
	/** SET UP INTEGRATOR **/
	/** set up integrator struct **/
	struct F_SCX p;
	p.nuc = _nuc;
	p.fp  = &_fp;
	p.b   = b;
	p.z   = z;
	p.f   = f_SCX_integrand;
	/** set up integrator **/
	int ndim = 2; // in b and z
	int ncomp = 1;
	integrand_t integr = F_SCX::exec;
	void* userdata = (void*) &p; // cast the integrator struct to void pointer
	int nvec = 1; // unless you now what SIMD is leave this to 1
	double epsrel = 1e-4;
	double epsabs = 1e-4;
	int flags = 0x00;
	int mineval = 10000;
	int maxeval = 100000;
	int key = 13; // 13 only available in two dimensions
	char* statefile = NULL;
	int nregions,neval,fail;
	double* integral = new double;
	double* error = new double;
	double* prob = new double;
	/** INTEGRATOR SET UP DONE **/

	Cuhre(ndim,ncomp,integr,userdata,nvec,
		epsrel,epsabs,flags,
		mineval,maxeval,key,
		statefile,&nregions,&neval,&fail,
		integral,error,prob);
	printf("# Cuhre  : res = %e, error = %e, prob = %f, neval = %d, fail = %d \n",*integral,*error,*prob,neval,fail);
	complex<double> scatfront = (_fp.getParticletype() == 8) ? _fp.getScatterfront(0) : _fp.getScatterfront(1); // 8 means proton, if so scatterfront of scattering with neutron (0) else _fp is neutron, take scattering with proton (1)
	const double betasq = (_fp.getParticletype()==8) ? _fp.getBeta2n() : _fp.getBeta2p();
	return *integral;
	//return 2.*M_PI*scatfront*exp(-b*b/(2.*betasq))*(*integral);	
}

