#include "ClassGridThick_SCX.hpp"
#include <gsl/gsl_sf_bessel.h>
#include <numint/numint.hpp>
#include <fstream>
#include "2bodymom.hpp"

ClassGridThick_SCX::ClassGridThick_SCX(MeanFieldNucleusThick* nuc,FastParticle* fp,int bpoints,int zpoints) : 
	_nuc(nuc),
	_fp(fp),
	_bpoints(bpoints),
	_zpoints(zpoints),
	_bstep(0.),
	_zstep(0.),
	_pdens_fctr(nuc->getZ()),
	_ndens_fctr(nuc->getN()),
	_pgrid(NULL),
	_ngrid(NULL),
	_perrorGrid(NULL),
	_nerrorGrid(NULL)
	{
	assert( _fp->getParticletype() == FastParticle::Particletype::P_CLASS_SCX || _fp->getParticletype() == FastParticle::Particletype::N_CLASS_SCX);
	_bstep =    _nuc->getRange()/(_bpoints-1.); // 0 to nuc.getRange
	_zstep = 2.*_nuc->getRange()/(_zpoints-1.); // -nuc.getRange() to +nuc.getRange()
	
	sprintf(_filename,"%s/grids/ClassGridThick_SCX_nuc%s_b%dz%d.bin",SHAREDIR.c_str(),nuc->getNucleusName().c_str(),bpoints,zpoints);
	std::cout << "#[ClassGridThick::IO]  IO filename is " << _filename << std::endl;
	std::cout << "#[ClassGridThick::info]  FastParticle (Sigmap,Sigman) is " << "(" << _fp->getSigmap() << ", " << _fp->getSigman() << ")"  << std::endl;
}

ClassGridThick_SCX::~ClassGridThick_SCX(){
	if ( _pgrid != NULL ){ // all the other grids are also NULL if _pgrid is NULL
		for (int bi=0; bi<_bpoints;bi++){
			delete[] _pgrid[bi];
			delete[] _ngrid[bi];
			delete[] _perrorGrid[bi];
			delete[] _nerrorGrid[bi];
		}
		delete[] _pgrid;
		delete[] _ngrid;
		delete[] _perrorGrid;
		delete[] _nerrorGrid;
	}
}

void ClassGridThick_SCX::addKnockoutParticle(int level){ // add a knockout particle, on which no FSIs take place, only to correct for density changes
	assert(level < _nuc->getTotalLevels());
	if (level < _nuc->getPLevels()){ // it is a proton
		_pdens_fctr -= 1.; ///_nuc->getZ();
	} else {
		_ndens_fctr -= 1.; ///_nuc->getN();
	}
	//std::cout << "#[ClassGridThick::info] density correction factors are (p,n) = (" << _pdens_fctr << ", " << _ndens_fctr << ") " << std::endl;
}

void ClassGridThick_SCX::clearKnockout(){
	_pdens_fctr = _nuc->getZ(); //1.; // because density is normed to one, we need density normed to Z
	_ndens_fctr = _nuc->getN(); //1.; // because density is normed to one, we need density normed to N
}

void ClassGridThick_SCX::constructGrid(){
	assert(_pgrid==NULL && _ngrid==NULL && _perrorGrid==NULL && _nerrorGrid==NULL && _bstep > 0. && _zstep > 0.);
	std::fstream file;
	file.open(_filename,std::ios::binary | std::ios::in);
	bool input = (bool)file;
	if (!input){ // we have to write the grid!
		file.close();
		std::cout << "#[ClassGridThick::IO] GRID NOT FOUND " << _filename << " PREPARING TO WRITE TO IT " << std::endl;
	} else {
		std::cout << "#[ClassGridThick::IO] GRID EXISTS " << _filename << " PREPARING TO READ FROM IT " << std::endl;
		assert( file.is_open()); // make sure read file has been opened
		// get length of file should be 4*_bpoints*_zpoints*sizeof(double), if not file is probs. corrupt, recalculate tha shizzle
		file.seekg(0,file.end);
		unsigned int length = file.tellg();
		file.seekg(0,file.beg);
		if (!(length == 4*_bpoints*_zpoints*sizeof(double))){
			std::cout << "#[ClassGridThick::IO ERROR] GRID HAS NOT THE RIGHT SIZE, RECALCULATING! " << std::endl;
			input = false;
			file.close();
		}
	}
	// std::cout << "#[ClassGridThick::IO] filename exists? " << file << std::endl;
	std::cout << "#[ClassGridThick::info] constructing grid " << _bpoints << " by " << _zpoints << std::endl;
	std::cout << "#[ClassGridThick::info] steps are " << _bstep << "  " << _zstep << std::endl;
	std::cout << "#[ClassGridThick::info] density correction factors are (p,n) = (" << _pdens_fctr << ", " << _ndens_fctr << ") " << std::endl;
	_pgrid = new double*[_bpoints];
	_ngrid = new double*[_bpoints];
	_perrorGrid = new double*[_bpoints];
	_nerrorGrid = new double*[_bpoints];
	for (int bi=0; bi<_bpoints; bi++){
		_pgrid[bi] = new double[_zpoints];
		_ngrid[bi] = new double[_zpoints];
		_perrorGrid[bi] = new double[_zpoints];
		_nerrorGrid[bi] = new double[_zpoints];
		for (int zi=0; zi<_zpoints; zi++){
			double b = bi*_bstep; 
			double z = -_nuc->getRange()+zi*_zstep;
			double pres,nres,perr,nerr;
			if (input){
				file.read((char*)&pres,sizeof(double)); // read data
				file.read((char*)&nres,sizeof(double)); // read data
				file.read((char*)&perr,sizeof(double)); // read error
				file.read((char*)&nerr,sizeof(double)); // read error
			} else {
				// dont write here yet, because file will be open verry long
				calcq(b,z,pres,perr,&MeanFieldNucleusThick::getNeutronDensity); // leading proton in pgrid interacts with neutrons
				calcq(b,z,nres,nerr,&MeanFieldNucleusThick::getProtonDensity);  // leading neutron in ngrid interacts with protons
				std::cout << "#[ClassGridThick::info] >> calculated grid point " << bi << ", " << zi << "  ->  " << _bpoints << ", " << _zpoints << std::endl;
			}
			_pgrid[bi][zi] = pres;
			_ngrid[bi][zi] = nres;
			_perrorGrid[bi][zi] = perr;
			_nerrorGrid[bi][zi] = nerr;
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
				file.write((char*)&_pgrid[bi][zi],sizeof(double)); // write data
				file.write((char*)&_ngrid[bi][zi],sizeof(double)); // write data
				file.write((char*)&_perrorGrid[bi][zi],sizeof(double)); // write error
				file.write((char*)&_nerrorGrid[bi][zi],sizeof(double)); // write error
			}
		}
	}	
	assert(file.tellg()==_bpoints*_zpoints*4*sizeof(double)); // this should always be true, just to check, times 2 because data and errors are saved, if found to fail, maybe corrupt file?
	file.close();
	
	
	/**
	 * Now we add the energy dependent factor to the grid. Grid is saved without the energy dependent
	 * sigma (total cross section) or density correction factors
	 * Note that because of SCX either Sigman or Sigmap for (neutron,proton) will be zero!
	 * Bit of unnecessary work, but to keep things flexible....
	 * What is in the grid is the expression ( - \int_{z}^{\+\infty} rho(b,z') \textrm{d}z')
	 * Must be converted to 1. - exp( - sigma * densfctr * (-\int ...) )
	 */
	for (int bi=0;bi<_bpoints; bi++){
		for (int zi=0;zi<_zpoints; zi++){
			//std::cout << "#Transforming energy independent factor pgrid factor " << _pgrid[bi][zi] << " to " <<  1.- pow( _pgrid[bi][zi], _fp->getSigman()*_ndens_fctr) << " power is " << _fp->getSigman()*_ndens_fctr << std::endl;
			_pgrid[bi][zi]      = 1.- exp( _pgrid[bi][zi]*_fp->getSigman()*_ndens_fctr);
			//std::cout << "#Transforming energy independent factor ngrid factor " << _ngrid[bi][zi] << " to " <<  1.- pow( _pgrid[bi][zi], _fp->getSigmap()*_pdens_fctr) << " power is " << _fp->getSigmap()*_pdens_fctr << std::endl;
			_ngrid[bi][zi]      = 1.- exp( _ngrid[bi][zi]*_fp->getSigmap()*_pdens_fctr);
			_perrorGrid[bi][zi]*= _fp->getSigman()*_ndens_fctr;
			_perrorGrid[bi][zi]*= _fp->getSigmap()*_pdens_fctr;
		}
	}
}

void ClassGridThick_SCX::printGrid(){ // debug function
	std::cout << 0                 << "\t" << _nuc->getRange() << "\t" << _bpoints << std::endl;
	std::cout << -_nuc->getRange() << "\t" << _nuc->getRange() << "\t" << _zpoints << std::endl;
	for (int bi=0; bi<_bpoints;bi++){
		for (int zi=0; zi<_zpoints;zi++){
			//double b = -_nuc->getRange()+bi*_bstep; 
			//double z = -_nuc->getRange()+zi*_zstep;
			std::cout << _pgrid[bi][zi] << "\t" << _ngrid[bi][zi]<< "\t";
			std::cout << _perrorGrid[bi][zi] << "\t" << _nerrorGrid[bi][zi] << std::endl;
		}
	}
}

void ClassGridThick_SCX::printDensityGrid(){
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
double ClassGridThick_SCX::getInterp(numint::array<double,3> x){ // x is r,costheta,phi 
	//std::cout << "# GlauberGrid got interpolation request for " << x[0] << " " << x[1] << std::endl;
	assert( _pgrid != NULL); // if one grid is not NULL, all the others will not be NULL as well
	_fp->setHitcoord(x[0],x[1],sqrt(1.-x[1]*x[1]),cos(x[2]),sin(x[2]) ); // set the hitcoordinates for the fast particle
	double** _grid = (_fp->getParticletype() == FastParticle::P_CLASS_SCX) ? _pgrid : _ngrid;
	//std::cout << "#Interpolation for " << x[0] << ", " << x[1] << ", " << x[2] << "
	double z = _fp->getHitz();//x[0]*x[1]; // z = r \cos \theta
	double b = _fp->getHitbnorm(); // |b|, always larger than one
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
double ClassGridThick_SCX::getInterp(double x[]){ // x is r,costheta,phi 
	//std::cout << "# ClassGridThick_SCX got interpolation request for " << x[0] << " " << x[1] << " " << x[2] << std::endl;
	assert( _pgrid != NULL); // if one grid is not NULL, all the others will not be NULL as well
	_fp->setHitcoord(x[0],x[1],sqrt(1.-x[1]*x[1]),cos(x[2]),sin(x[2]) ); // set the hitcoordinates for the fast particle
	double** _grid = (_fp->getParticletype() == FastParticle::P_CLASS_SCX) ? _pgrid : _ngrid;
	//std::cout << "#Interpolation for " << x[0] << ", " << x[1] << ", " << x[2] << "
	double z = _fp->getHitz();//x[0]*x[1]; // z = r \cos \theta
	double b = _fp->getHitbnorm(); // |b|, always larger than one
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

/** integration struct **/
struct F_SCX {
	static double exec(double z,void* param){
		struct F_SCX p = *(struct F_SCX*) param;
		double r = sqrt(p.b*p.b+z*z);
		if (r < p.nuc->getWF_r_step()) // going closer to 0 than this causes divergencies!
			r = p.nuc->getWF_r_step();
		else if ( r > p.nuc->getRange())
			return 0.; // outside RMAX density is for virtually zero
		return -((p.nuc)->*(p.densf))(r)/r/r; // note the minus sign ;), divide by r**2!
		/** A bit explanation of the quite complicated line above
		 * the outer brackets around ( (p.nuc)->*(p.densf) ) are necessary
		 * because (p.nuc) is the meanfieldnucleus pointer
		 * ->* is the member function selection operator
		 * (p.densf) is the member function pointer.
		 * Hence only the full expression ( (p.nuc)->*(p.densf) ) actually
		 * returns the function pointer of the object p.nuc, call this with
		 * parameter r. hence the last (r)
		 */
	}
	MeanFieldNucleusThick* nuc;
	double (MeanFieldNucleusThick::*densf)(const double) const;
	double b;
};

/**calculate the grid factors given (b,z) and density function densf
 * you still have to multiply this with the correct prefactors!
 */
void ClassGridThick_SCX::calcq(double b,double z,double& ret,double& err, double (MeanFieldNucleusThick::*densf)(const double) const){
	/** set up integrator struct **/
	struct F_SCX p;
	p.nuc   = _nuc;
	p.b     = b;
	p.densf = densf;
	/** SET UP GSL 1D INTEGRATOR **/
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (2048);
	gsl_function F;
	F.function = F_SCX::exec;
	F.params = (void*) &p;
	double epsabs=0;
	double epsrel=1e-4;
	int key      =6;
	size_t limit = 2048; // this should never exceed size of the allocated integration workspace
	/** DO THE INTEGRATION **/
	gsl_integration_qag(&F, z, _nuc->getRange(),epsabs,epsrel,limit,key,w,&ret,&err);
}
