#include "GlauberGridThick_SEL.hpp"
#include <gsl/gsl_sf_bessel.h>
#include <numint/numint.hpp>
#include <fstream>
#include "2bodymom.hpp"
/** integration struct 
 * ndim  = 2 (|b|,z)
 * ncomp = 1 integrand is real imaginary piece
 * comes from prefactor (scatterFront)
 *
 * **/
struct F_SEL {
	static int exec(const int* ndim,const double x[],const int* ncomp, double* res,void *param){
		struct F_SEL p = *(struct F_SEL*) param;
		return p.f(res,x,p.nuc,p.fp,p.b,p.z);
	}
	MeanFieldNucleusThick* nuc;
	FastParticle* fp;
	double b,z;
	int (*f)(double* res, const double x[], MeanFieldNucleusThick*,FastParticle*,const double b,const double z);
};

/** integrand for elastic scattering with protons in the nucleus **/
int f_SEL_integrand_p(double* res, const double x[],MeanFieldNucleusThick* nuc,FastParticle* fp, const double b, const double z){
	double bi    = x[0]; // integration variable |b'|
	double th_bi = x[1]; // integration variable theta(b')
	double zi    = x[2]; // integration variable z'
	bi    = bi*nuc->getRange(); // map [0,1] -> [0,R]
	th_bi = th_bi*M_PI;         // map [0,1] -> [0,\pi], remember to multiply with 2 afterwards!
	zi    = z+zi*(nuc->getRange()-z); // map [0,1] -> [z,R]
	const double Jac    = nuc->getRange()*M_PI*(nuc->getRange()-z); // Jacobian of unit cube transformation
	const double betasq = fp->getBeta2p(); // get betasq of scattering with protons
	double ri = sqrt( bi*bi - 2.*b*bi*cos(th_bi) + b*b + zi*zi);
	if (ri < nuc->getWF_r_step()) // going closer to 0 than this causes divergencies!
		ri = nuc->getWF_r_step();
	if (ri>nuc->getRange()){
		*res = 0;
	} else {
		double dens = nuc->getProtonDensity(ri)/ri/ri; // get proton density
		*res = Jac*bi*exp(-bi*bi/(2.*betasq))*dens;
	}
	return 0;
}

/** integrand for elastic scattering with neutrons in the nucleus **/
int f_SEL_integrand_n(double* res, const double x[],MeanFieldNucleusThick* nuc,FastParticle* fp, const double b, const double z){
	double bi    = x[0]; // integration variable |b'|
	double th_bi = x[1]; // integration variable theta(b')
	double zi    = x[2]; // integration variable z'
	bi    = bi*nuc->getRange(); // map [0,1] -> [0,R]
	th_bi = th_bi*M_PI;         // map [0,1] -> [0,\pi], remember to multiply with two afterwards!
	zi    = z+zi*(nuc->getRange()-z); // map [0,1] -> [z,R]
	const double Jac    = nuc->getRange()*M_PI*(nuc->getRange()-z); // Jacobian of unit cube transformation
	const double betasq = fp->getBeta2n(); // get betasq of scattering with neutrons
	double ri = sqrt( bi*bi - 2.*b*bi*cos(th_bi) + b*b + zi*zi);
	if (ri < nuc->getWF_r_step()) // going closer to 0 than this causes divergencies!
		ri = nuc->getWF_r_step();
	if (ri>nuc->getRange()){
		*res = 0;
	} else {
		double dens = nuc->getNeutronDensity(ri)/ri/ri; // get neutron density
		*res = Jac*bi*exp(-bi*bi/(2.*betasq))*dens;
	}
	return 0;
}

GlauberGridThick_SEL::GlauberGridThick_SEL(MeanFieldNucleusThick* nuc,FastParticle& fp,int bpoints,int zpoints) : 
	_nuc(nuc),
       	_fp(fp), 
	_bpoints(bpoints),
	_zpoints(zpoints),
	_bstep(0.),
	_zstep(0.),
	_pdens_fctr(nuc->getZ()),
	_ndens_fctr(nuc->getN()),
	_grid(NULL)
	{
	assert(fp.getParticletype()==0 || fp.getParticletype()==1); // make sure we are getting the glauber params for p/n elastic!
	_bstep =    _nuc->getRange()/(_bpoints-1.); // 0 to nuc.getRange
	_zstep = 2.*_nuc->getRange()/(_zpoints-1.); // -nuc.getRange() to +nuc.getRange()
	std::cout << "#Particle type is " << fp.getParticletype() << " Beta2p is " << fp.getBeta2p() << " Beta2n is " << fp.getBeta2n() << std::endl;
	std::cout << "#Sigmap is " << fp.getSigmap() << "  Sigman is " << fp.getSigman() << std::endl;
	std::cout << "#Epsilonp is " << fp.getEpsilonp() << " Epsilonn is " << fp.getEpsilonn() << std::endl;
	std::cout << "#Scatterfront is : " << ((fp.getParticletype()==8)? fp.getScatterfront(0) : fp.getScatterfront(1)) << std::endl;
	std::cout << "#bsteps, zstep is " << _bstep << ", " << _zstep << " and range is " << nuc->getRange() << std::endl;
	std::cout << "#Momentum of particle is [mag,theta,phi]: " << fp.getP() << ", " << fp.getTheta() << ", " << fp.getPhi() << std::endl;
	sprintf(_filename,"%s/grids/GlauberGridThick_SEL_nuc%s_type%d_p10_%.0f_beta2p%.2e_beta2n%.2e_bp%d_zp%d.bin",SHAREDIR.c_str(),nuc->getNucleusName().c_str(),fp.getParticletype(),fp.getP()/10.,fp.getBeta2p(),fp.getBeta2n(),bpoints,zpoints); // divide momentum by 10 in filename so that we effectively do a binning in Pmag() space.
	std::cout << "#filename is " << _filename << std::endl;
	//constructGlauberGrid(); // don't call this from constructor yet. Leave up to user so that density correction and arbitrary phase can be set before construction of the grid
}

GlauberGridThick_SEL::~GlauberGridThick_SEL(){
	for (int bi=0; bi<_bpoints;bi++){
		delete[] _grid[bi];
		delete[] _errorGrid[bi];
	}
	delete[] _grid;
	delete[] _errorGrid;
}

void GlauberGridThick_SEL::addKnockoutParticle(int level){ // add a knockout particle, on which no FSIs take place, only to correct for density changes
	assert(level < _nuc->getTotalLevels());
	if (level < _nuc->getPLevels()){ // it is a proton
		_pdens_fctr -= 1.; ///_nuc->getZ();
	} else {
		_ndens_fctr -= 1.; ///_nuc->getN();
	}
	std::cout << "#density correction factors are (p,n) = (" << _pdens_fctr << ", " << _ndens_fctr << ") " << std::endl;
}

void GlauberGridThick_SEL::clearKnockout(){
	_pdens_fctr = _nuc->getZ(); //1.; // because density is normed to one, we need density normed to Z
	_ndens_fctr = _nuc->getN(); //1.; // because density is normed to one, we need density normed to N
}

void GlauberGridThick_SEL::constructGlauberGrid(){
	assert(_grid==NULL && _bstep > 0. && _zstep > 0.);
	std::fstream file;
	file.open(_filename,std::ios::binary | std::ios::in);
	bool input = (bool)file;
	if (!input){ // we have to write the grid!
		file.close();
		std::cout << "#GRID NOT FOUND " << _filename << " PREPARING TO WRITE TO IT " << std::endl;
	} else {
		std::cout << "#GRID EXISTS " << _filename << " PREPARING TO READ FROM IT " << std::endl;
		assert( file.is_open()); // make sure read file has been opened
		// get length of file should be 2*_bpoints*_zpoints*sizeof(complex<double>), if not file is probs. corrupt, recalculate tha shizzle
		file.seekg(0,file.end);
		unsigned int length = file.tellg();
		file.seekg(0,file.beg);
		if (!(length == 2*_bpoints*_zpoints*sizeof(std::complex<double>))){
			std::cout << "#GRID HAS NOT THE RIGHT SIZE, RECALCULATING! " << std::endl;
			input = false;
			file.close();
		}
	}
	// std::cout << "#filename exists? " << file << std::endl;
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
			std::complex<double> res,err;
			if (input){
				file.read((char*)&res,sizeof(std::complex<double>)); // read data
				file.read((char*)&err,sizeof(std::complex<double>)); // read error
			} else {
				calcFSI(b,z,res,err);
			}
			_grid[bi][zi]      =   res;
			_errorGrid[bi][zi] =   err;
			//std::cout << " [" << bi << ", " << zi << "] corresponding to " << b << ", " << z << " init to " << _grid[bi][zi] << std::endl;
		}
	}
	/** 
	 * write all the calculated fsi grid points
	 *  only when they are all calculated, to minimize
	 *  time output file is left open...
	 **/
	if (!input){
		file.open(_filename,std::ios::binary | std::ios::out);
		assert( file.is_open()); // make sure the file has been opened
		for (int bi=0; bi<_bpoints;bi++){
			for (int zi=0; zi<_zpoints;zi++){
				file.write((char*)&_grid[bi][zi],sizeof(std::complex<double>)); // write data
				file.write((char*)&_errorGrid[bi][zi],sizeof(std::complex<double>)); // write error
			}
		}
	}
	assert(file.tellg()==_bpoints*_zpoints*2*sizeof(std::complex<double>)); // this should always be true, just to check, times 2 because data and errors are saved, if found to fail, maybe corrupt file?
	file.close();
}

void GlauberGridThick_SEL::printGrid(){ // debug function
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

void GlauberGridThick_SEL::printDensityGrid(){
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
complex<double> GlauberGridThick_SEL::getInterp(numint::array<double,3> x){ // x is r,costheta,phi 
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
 *  @param x the coordinates in \f$r, \, \cos \theta\ \,\text{and} \phi$\f
 *  */
complex<double> GlauberGridThick_SEL::getInterp(double x[]){ // x is r,costheta,phi 
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
	assert(bi >= 0 && bi < _bpoints && zi >= 0 && zi < _zpoints);
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
 * Warning! Does not include "scatterfront" here
 * to keep the integral prefactor independent (more reusability)
 * so that integral only depends on beta
 * computes the integral
 * \f$ \int \textrm{d}z' \theta(z-z') \int \textrm{d}^{2}b' e^{-\frac{(b-b')^{2}}{2 \beta^{2}} } \rho(b',z') \f$
 * you still have to multiply this with the correct prefactors!
 */
void GlauberGridThick_SEL::calcFSI(double b,double z,std::complex<double>& ret,std::complex<double>& err){
	/** SET UP INTEGRATOR **/
	/** set up integrator struct **/
	struct F_SEL p;
	p.nuc = _nuc;
	p.fp  = &_fp;
	p.b   = b;
	p.z   = z;
	/** set up integrator **/
	int ndim = 3; // in b(2) and z(1)
	int ncomp = 1;
	integrand_t integr = F_SEL::exec;
	void* userdata = (void*) &p; // cast the integrator struct to void pointer
	int nvec = 1; // unless you now what SIMD is leave this to 1
	double epsrel = 1e-6;
	double epsabs = 1e-6;
	int flags = 0x00;
	int mineval = 10000;
	int maxeval = 10000000;
	char* statefile = NULL;
	int neval,fail;
	double integral_p;
	double error_p;
	double prob_p;
	
	double integral_n;
	double error_n;
	double prob_n;
	/** INTEGRATOR SET UP DONE **/
	/** cuhre specific arguments **/
	/*	
	int key = 13; // 13 only available in two dimensions
	int nregions;	
	
	p.f   = f_SEL_integrand_p; // do elastic scattering with protons in the nucleus
	
	Cuhre(ndim,ncomp,integr,userdata,nvec,
		epsrel,epsabs,flags,
		mineval,maxeval,key,
		statefile,&nregions,&neval,&fail,
		&integral_p,&error_p,&prob_p);
	
	p.f = f_SEL_integrand_n; // now do elastic scattering with neutrons
	
	Cuhre(ndim,ncomp,integr,userdata,nvec,
		epsrel,epsabs,flags,
		mineval,maxeval,key,
		statefile,&nregions,&neval,&fail,
		&integral_n,&error_n,&prob_n);
	*/
	/** vegas specific arguments **/
	/** vegas found to preform better than cuhre in this case **/
	
	int seed=1234;
	int nstart=10000;
	int nincrease=20000;
	int nbatch=10000;
	int gridno=0;
	
	p.f = f_SEL_integrand_p; // do elastic scattering with protons in the nucleus

	Vegas(ndim,ncomp,integr,userdata,nvec,
		epsrel,epsabs,flags,seed,
		mineval,maxeval,nstart,nincrease,
		nbatch,gridno,statefile,NULL,&neval,
		&fail,&integral_p,&error_p,&prob_p);

	p.f = f_SEL_integrand_n; // now do elastic scattering with neutrons
	
	Vegas(ndim,ncomp,integr,userdata,nvec,
		epsrel,epsabs,flags,seed,
		mineval,maxeval,nstart,nincrease,
		nbatch,gridno,statefile,NULL,&neval,
		&fail,&integral_n,&error_n,&prob_n);

	//printf("# Vegas  (b,z) = (%6.2f,%6.2f) : res = %e, error = %e, prob = %f, neval = %d ,fail = %d\n",b,z,integral,error,prob,neval,fail);
	// why two 2*dens below you ask? *2 because phi is only integrated over [0,\pi], times two due to symmetry
	// the dens factor comes from the single scattering approximation
	//ret = (1. - 2.*_pdens_fctr*_fp.getScatterfront(true)*integral_p )*( 1. - 2.*_ndens_fctr*_fp.getScatterfront(false)*integral_n); 
	// formula below instead of the above one because of the SINGLE scattering approx (only linear terms in profile functions!)
	ret = 1. - 2.*_pdens_fctr*_fp.getScatterfront(true)*integral_p - 2.*_ndens_fctr*_fp.getScatterfront(false)*integral_n; 
	err =      2.*_fp.getScatterfront(true)*error_p +                2.*_ndens_fctr*_fp.getScatterfront(false)*error_n;

	//printf("# Cuhre (b,z) = (%6.2f,%6.2f) : res = %e, error = %e, prob = %f, neval = %d, fail = %d \n",b,z,*integral,*error,*prob,neval,fail);

}
