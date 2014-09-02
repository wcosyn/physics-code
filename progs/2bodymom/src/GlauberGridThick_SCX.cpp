#include "GlauberGridThick_SCX.hpp"
#include <gsl/gsl_sf_bessel.h>
#include <numint/numint.hpp>
/** integration struct 
 * ndim  = 2 (|b|,z)
 * ncomp = 1 integrand is real imaginary piece
 * comes from prefactor (scatterFront)
 *
 * **/
struct F_SCX {
	static int exec(const int* ndim,const double x[],const int* ncomp, double* res,void *param){
		struct F_SCX p = *(struct F_SCX*) param;
		return p.f(res,x,p.nuc,p.fp,p.b,p.z);
	}
	MeanFieldNucleusThick* nuc;
	FastParticle* fp;
	double b,z;
	int (*f)(double* res, const double x[], MeanFieldNucleusThick*,FastParticle*,const double b,const double z);
};

int f_SCX_integrand(double* res, const double x[], MeanFieldNucleusThick* nuc, FastParticle* fp,const double b, const double z){
	double bi = x[0]; // integration variable b', here in unit square
	double zi = x[1]; // integration variable z', here in unit square
	bi = bi*nuc->getRange(); // map [0,1] -> [0,R]
	zi = z+zi*(nuc->getRange()-z); // map [0,1] -> [z,R] (z instead of -R because of step function HS(z-z')
	assert( bi>=0. && b>=0.); // DEBUG, remove if found to never fail
	//std::cout << " translated unit coords " << x[0] << " and " << x[1] << " to world coords " << bi << " and " << zi << std::endl;
	const double Jac = nuc->getRange()*(nuc->getRange()-z); // Jacobian of mapping
	//std::cout << " z is " << z << std::endl;
	//std::cout << " Jacobian of transformation is " << Jac << std::endl;
	double ri  = sqrt(bi*bi + zi*zi);
	if (ri < 1e-12) // arbitrary value, choose close enough to zero plz (close enough as in distance much smaller than the distance scale of the density
		ri = 1e-12; // prevent division by zero in the density
	const double betasq = (fp->getParticletype()==8) ? fp->getBeta2n() : fp->getBeta2p();
	//std::cout << " betasq is " << betasq << std::endl;
	if (ri>nuc->getRange()){ // outside the range, density is (practically) zero, so return zero
		*res=0;
	} else {
		double dens   = (fp->getParticletype()==8) ? nuc->getNeutronDensity(ri)/ri/ri : nuc->getProtonDensity(ri)/ri/ri; // 8 is proton, so SCX with neutron density, vice versa for 9, dens divided by r^2 because dens is rho*r*r

		*res = Jac*bi*exp(-bi*bi/(2.*betasq))*gsl_sf_bessel_I0(bi*b/betasq)*dens; // !!!WARNING!!! make sure both bi and b are > 0 !!!
	}
	return 0; // always "succes"
}

GlauberGridThick_SCX::GlauberGridThick_SCX(MeanFieldNucleusThick* nuc,FastParticle& fp,int bpoints,int zpoints) : 
	_nuc(nuc),
       	_fp(fp), 
	_bpoints(bpoints),
	_zpoints(zpoints),
	_bstep(0.),
	_zstep(0.),
	_pdens_fctr(1.),
	_ndens_fctr(1.),
	_grid(NULL)
	{
	assert(fp.getParticletype()==8 || fp.getParticletype()==9); // make sure we are getting the glauber params for p/n charge exchange!
	_bstep =    _nuc->getRange()/(_bpoints-1.); // 0 to nuc.getRange
	_zstep = 2.*_nuc->getRange()/(_zpoints-1.); // -nuc.getRange() to +nuc.getRange()
	std::cout << "#Particle type is " << fp.getParticletype() << " Beta2p is " << fp.getBeta2p() << " Beta2n is " << fp.getBeta2n() << std::endl;
	std::cout << "#Sigmap is " << fp.getSigmap() << "  Sigman is " << fp.getSigman() << std::endl;
	std::cout << "#Scatterfront is (proton): " << fp.getScatterfront(1) << "  (neutron): " << fp.getScatterfront(0) << std::endl;
	std::cout << "#bsteps, zstep is " << _bstep << ", " << _zstep << " and range is " << nuc->getRange() << std::endl;
	//constructGlauberGrid();
}

GlauberGridThick_SCX::~GlauberGridThick_SCX(){
	for (int bi=0; bi<_bpoints;bi++){
		delete[] _grid[bi];
		delete[] _errorGrid[bi];
	}
	delete[] _grid;
	delete[] _errorGrid;
}

void GlauberGridThick_SCX::addKnockoutParticle(int level){ // add a knockout particle, on which no FSIs take place, only to correct for density changes
	assert(level < _nuc->getTotalLevels());
	if (level < _nuc->getPLevels()){ // it is a proton
		_pdens_fctr -= 1./_nuc->getZ();
	} else {
		_ndens_fctr -= 1./_nuc->getN();
	}
}

void GlauberGridThick_SCX::clearKnockout(){
	_pdens_fctr = 1.;
	_ndens_fctr = 1.;
}

void GlauberGridThick_SCX::constructGlauberGrid(){
	assert(_grid==NULL && _bstep > 0. && _zstep > 0.);
	std::cout << "#constructing grid " << _bpoints << " by " << _zpoints << std::endl;
	std::cout << "#steps are " << _bstep << "  " << _zstep << std::endl;
	std::cout << "#density correction factors are (p,n) = (" << _pdens_fctr << ", " << _ndens_fctr << ") " << std::endl;
	_grid = new complex<double>*[_bpoints];	
	_errorGrid = new complex<double>*[_bpoints];
	for (int bi=0; bi<_bpoints; bi++){
		_grid[bi] = new complex<double>[_bpoints];
		_errorGrid[bi] = new complex<double>[_bpoints];
		for (int zi=0; zi<_zpoints; zi++){
			double b = bi*_bstep; 
			double z = -_nuc->getRange()+zi*_zstep;
			calcFSI(b,z,_grid[bi][zi],_errorGrid[bi][zi]);
			//std::cout << " [" << bi << ", " << zi << "] corresponding to " << b << ", " << z << " init to " << _grid[bi][zi] << std::endl;
		}
	}
}

void GlauberGridThick_SCX::printGrid(){ // debug function
	std::cout << 0                 << "\t" << _nuc->getRange() << "\t" << _bpoints << std::endl;
	std::cout << -_nuc->getRange() << "\t" << _nuc->getRange() << "\t" << _zpoints << std::endl;
	for (int bi=0; bi<_bpoints;bi++){
		for (int zi=0; zi<_zpoints;zi++){
			//double b = -_nuc->getRange()+bi*_bstep; 
			//double z = -_nuc->getRange()+zi*_zstep;
			std::cout << _grid[bi][zi].real() << "\t" << _grid[bi][zi].imag() << "\t";
			std::cout << _errorGrid[bi][zi].real() << "\t" << _errorGrid[bi][zi].imag() << std::endl;
		}
	}
}

void GlauberGridThick_SCX::printDensityGrid(){
	for (int bi=0; bi<_bpoints;bi++){
		for (int zi=0; zi<_zpoints;zi++){
			double b = bi*_bstep; 
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


/** Interpolation of GlauberGrid.
 *  Transform momentum to frame where mom is // z-axis
 *  to fetch the correct fsi factor
 *  @param x the coordinates in \f$r, \, \cos \theta\ \,\text{and} \phi$
 *  */
complex<double> GlauberGridThick_SCX::getInterp(numint::array<double,3> x){ // x is r,costheta,phi 
	//std::cout << "# GlauberGrid got interpolation request for " << x[0] << " " << x[1] << std::endl;
	
	_fp.setHitcoord(x[0],x[1],sqrt(1.-x[1]*x[1]),cos(x[2]),sin(x[2]) ); // set the hitcoordinates for the fast particle
	//std::cout << "#Interpolation for " << x[0] << ", " << x[1] << ", " << x[2] << "	
	double z = _fp.getHitz();//x[0]*x[1]; // z = r \cos \theta
	double b = _fp.getHitbnorm(); // |b|, always larger than one
	double bi_t,zi_t; 
	double bf = modf(b/_bstep,&bi_t); // (b-bmin(=0))/bstep, split into integer (bi) and fractional (bf) part
	double zf = modf((z+_nuc->getRange())/_zstep,&zi_t); // (z-zmin(=-R))/_zstep, split into integer (zi) and fractional (zf) part
	int bi = (int) bi_t;
	int zi = (int) zi_t;
	if (bi == _bpoints-1) { // take care at high end
		bi-=1; 
		bf =1.;
	}
	if (zi == _zpoints-1){ // take care at high end
		zi-=1;
		zf =1.;
	}
	const double bf_comp = 1.-bf;
	const double zf_comp = 1.-zf;
	//std::cout << "# converted to b,z " << b << ", " << z << " converted to indices " << bi << ", " << zi << " and f parts " << bf << ", " << zf << std::endl;
	// interpolation first twice in z direction
	// then once in b direction
	return bf_comp*( zf_comp*_grid[bi][zi] + zf*_grid[bi][zi+1] ) + bf*( zf_comp*_grid[bi+1][zi] + zf*_grid[bi+1][zi+1] );
}
/** Interpolation of GlauberGrid.
 *  Transform momentum to frame where mom is // z-axis
 *  to fetch the correct fsi factor
 *  @param x the coordinates in \f$r, \, \cos \theta\ \,\text{and} \phi$
 *  */
complex<double> GlauberGridThick_SCX::getInterp(double x[]){ // x is r,costheta,phi 
	//std::cout << "# GlauberGrid got interpolation request for " << x[0] << " " << x[1] << std::endl;
	
	_fp.setHitcoord(x[0],x[1],sqrt(1.-x[1]*x[1]),cos(x[2]),sin(x[2]) ); // set the hitcoordinates for the fast particle
	//std::cout << "#Interpolation for " << x[0] << ", " << x[1] << ", " << x[2] << "	
	double z = _fp.getHitz();//x[0]*x[1]; // z = r \cos \theta
	double b = _fp.getHitbnorm(); // |b|, always larger than one
	double bi_t,zi_t; 
	double bf = modf(b/_bstep,&bi_t); // (b-bmin(=0))/bstep, split into integer (bi) and fractional (bf) part
	double zf = modf((z+_nuc->getRange())/_zstep,&zi_t); // (z-zmin(=-R))/_zstep, split into integer (zi) and fractional (zf) part
	int bi = (int) bi_t;
	int zi = (int) zi_t;
	if (bi == _bpoints-1) { // take care at high end
		bi-=1; 
		bf =1.;
	}
	if (zi == _zpoints-1){ // take care at high end
		zi-=1;
		zf =1.;
	}
	const double bf_comp = 1.-bf;
	const double zf_comp = 1.-zf;
	//std::cout << "# converted to b,z " << b << ", " << z << " converted to indices " << bi << ", " << zi << " and f parts " << bf << ", " << zf << std::endl;
	// interpolation first twice in z direction
	// then once in b direction
	return bf_comp*( zf_comp*_grid[bi][zi] + zf*_grid[bi][zi+1] ) + bf*( zf_comp*_grid[bi+1][zi] + zf*_grid[bi+1][zi+1] );
}

/**calculate the glauber grid
 * will need the proton density for neutron knockout
 * will need the neutron density for proton knockout
 */
void GlauberGridThick_SCX::calcFSI(double b,double z,complex<double>& ret,complex<double>& err){
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
	int maxeval = 1000000;
	char* statefile = NULL;
	int neval,fail;
	double integral;
	double error;
	double prob;
	/** INTEGRATOR SET UP DONE **/
	/** cuhre specific arguments **/
	/**
	int key = 13; // 13 only available in two dimensions
	int nregions;	
	
	Cuhre(ndim,ncomp,integr,userdata,nvec,
		epsrel,epsabs,flags,
		mineval,maxeval,key,
		statefile,&nregions,&neval,&fail,
		&integral,&error,&prob);
	printf("# Cuhre (b,z) = (%6.2f,%6.2f) : res = %e, error = %e, prob = %f, neval = %d, fail = %d \n",b,z,*integral,*error,*prob,neval,fail);
	**/
	
	/** vegas specific arguments **/
	/** vegas found to preform better than cuhre in this case **/
	int seed=1234;
	int nstart=10000;
	int nincrease=20000;
	int nbatch=10000;
	int gridno=0;
	
	Vegas(ndim,ncomp,integr,userdata,nvec,
		epsrel,epsabs,flags,seed,
		mineval,maxeval,nstart,nincrease,
		nbatch,gridno,statefile,&neval,
		&fail,&integral,&error,&prob);
	printf("# Vegas  (b,z) = (%6.2f,%6.2f) : res = %e, error = %e, prob = %f, neval = %d ,fail = %d\n",b,z,integral,error,prob,neval,fail);

	complex<double> scatfront = (_fp.getParticletype() == 8) ? _fp.getScatterfront(0) : _fp.getScatterfront(1); // 8 means proton, if so scatterfront of scattering with neutron (0) else _fp is neutron, take scattering with proton (1)
	const double betasq   = (_fp.getParticletype()==8) ? _fp.getBeta2n() : _fp.getBeta2p();
	const double densfctr = (_fp.getParticletype()==8) ? _ndens_fctr     : _pdens_fctr;
	const complex<double> prefac = 2.*M_PI*scatfront*exp(-b*b/(2.*betasq));
	ret = 1.-densfctr*prefac*integral;
	err = densfctr*prefac*error;
}

