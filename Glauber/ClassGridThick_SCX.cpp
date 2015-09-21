#include "ClassGridThick_SCX.hpp"

/**
 * ClassGridThick_SCX constructor
 * @param densf      : Density function pointer. Takes b,z and returns density. 
 *                     Also takes void* pointer for passing extra arguments.
 * @param dens_param : void pointer to pass to the densf function as third argument.
 * @param densName   : The name of the density function you are using.
 *                     This is used to name output files. Use different names
 *                     for different density functions are you will run into trouble!
 * @param rmax       : Maximum radius of the system
 * @param fp         : FastParticle, particle you want to calculate SCX probability for
 * @param bpoints    : number of grid points in impact parameter b
 * @param zpoints    : number of grid points in impact parameter z
 **/
ClassGridThick_SCX::ClassGridThick_SCX(double (*densf)(double,void*),void* dens_param, const char* densName, double rmax, FastParticle* fp,int bpoints,int zpoints,char* SHAREDIR) : 
	_densf(densf),
	_dens_param(dens_param),
	_rmax(rmax),
	_fp(fp),
	_bpoints(bpoints),
	_zpoints(zpoints),
	_bstep(0.),
	_zstep(0.),
	_dens_fctr(1.),
	_grid(NULL),
	_errorGrid(NULL)
	{
	assert( _fp->getParticletype() == FastParticle::Particletype::P_CLASS_SCX || _fp->getParticletype() == FastParticle::Particletype::N_CLASS_SCX);
	_bstep =    _rmax/(_bpoints-1.); // 0 to rmax
	_zstep = 2.*_rmax/(_zpoints-1.); // -rmax to rmax
	
	sprintf(_filename,"%s/grids/ClassGridThick_SCX_%s_b%dz%d.bin",SHAREDIR,densName,bpoints,zpoints);
	std::cout << "#[ClassGridThick::IO]    IO filename is                   " << _filename << std::endl;
    std::cout << "#[ClassGridThick::info]  FastParticle Type is             " << _fp->getParticletype() << std::endl;
	std::cout << "#[ClassGridThick::info]  FastParticle (Sigmap,Sigman) is  " << "(" << _fp->getSigmap() << ", " << _fp->getSigman() << ")"  << std::endl;
}

ClassGridThick_SCX::~ClassGridThick_SCX(){
	if ( _grid != NULL ){ // make sure grid has been constructed
		for (int bi=0; bi<_bpoints;bi++){
			delete[] _grid[bi];
			delete[] _errorGrid[bi];
		}
		delete[] _grid;
		delete[] _errorGrid;
	}
}

/** Constructs the grid with actual SCX probabilities.
 *  Call this once you have set _dens_fctr to your liking...
 *  Remember density should be normed to A.
 * */
void ClassGridThick_SCX::constructGrid(){
	assert(_grid==NULL && _errorGrid==NULL && _bstep > 0. && _zstep > 0.);
	std::fstream file;
	file.open(_filename,std::ios::binary | std::ios::in);
	bool input = (bool)file;
	if (!input){ // we have to write the grid!
		file.close();
		std::cout << "#[ClassGridThick::IO]   GRID NOT FOUND " << _filename << " PREPARING TO WRITE TO IT " << std::endl;
	} else {
		std::cout << "#[ClassGridThick::IO]   GRID EXISTS " << _filename << " PREPARING TO READ FROM IT " << std::endl;
		assert( file.is_open()); // make sure read file has been opened
		// get length of file should be 2*_bpoints*_zpoints*sizeof(double), if not file is probs. corrupt, recalculate tha shizzle
		file.seekg(0,file.end);
		unsigned int length = file.tellg();
		file.seekg(0,file.beg);
		if (!(length == 2*_bpoints*_zpoints*sizeof(double))){
			std::cout << "#[ClassGridThick::IO ERROR] GRID HAS NOT THE RIGHT SIZE, RECALCULATING! " << std::endl;
			input = false;
			file.close();
		}
	}
	std::cout << "#[ClassGridThick::IO]   filename exists? " << &file << std::endl;
	std::cout << "#[ClassGridThick::info] constructing grid " << _bpoints << " by " << _zpoints << std::endl;
	std::cout << "#[ClassGridThick::info] steps are " << _bstep << "  " << _zstep << std::endl;
	std::cout << "#[ClassGridThick::info] density correction factor is = " << _dens_fctr << std::endl;
	_grid = new double*[_bpoints];
	_errorGrid = new double*[_bpoints];
	for (int bi=0; bi<_bpoints; bi++){
		_grid[bi] = new double[_zpoints];
		_errorGrid[bi] = new double[_zpoints];
		for (int zi=0; zi<_zpoints; zi++){
			double b = bi*_bstep;  // b [0,rmax]
			double z = -_rmax + zi*_zstep; // z [-rmax,rmax]
			double res,err;
			if (input){
				file.read((char*)&res,sizeof(double)); // read data
				file.read((char*)&err,sizeof(double)); // read error
			} else {
				// dont write here yet, because file will be open verry long
				calcP(b,z,res,err,_densf);
				//std::cout << "#[ClassGridThick::info] >> calculated grid point " << bi << ", " << zi << "  ->  " << _bpoints << ", " << _zpoints << " : " << res << std::endl;
			}
			_grid[bi][zi]      = res;
			_errorGrid[bi][zi] = err;
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
				file.write((char*)&_grid[bi][zi],sizeof(double)); // write data
				file.write((char*)&_errorGrid[bi][zi],sizeof(double)); // write error
			}
		}
	}	
	assert(file.tellg()==_bpoints*_zpoints*2*sizeof(double)); // this should always be true, just to check, times 2 because data and errors are saved, if found to fail, maybe corrupt file?
	file.close();
	
	
	/**
	 * Now we add the energy dependent factor to the grid. Grid is saved without the energy dependent
	 * sigma (total cross section) or density correction factors
	 * What is in the grid is the expression ( - \int_{z}^{\+\infty} rho(b,z') \textrm{d}z')
	 * Must be converted to 1. - exp( - sigma * densfctr * (-\int ...) )
	 */
	double (FastParticle::*fSigma)() const = (_fp->getParticletype() == FastParticle::Particletype::P_CLASS_SCX) ? &FastParticle::getSigman : &FastParticle::getSigmap ;
	for (int bi=0;bi<_bpoints; bi++){
		for (int zi=0;zi<_zpoints; zi++){
			_grid[bi][zi]      = 1.- exp( _grid[bi][zi]*(_fp->*fSigma)()*_dens_fctr); // This is the SCX probability!
			_errorGrid[bi][zi]*= (_fp->*fSigma)()*_dens_fctr; // error ~ first order taylor of exponential...
		}
	}
}

void ClassGridThick_SCX::printGrid(){ // debug function
	for (int bi=0; bi<_bpoints;bi++){
		for (int zi=0; zi<_zpoints;zi++){
			std::cout << bi << "\t" << zi << "\t" << _grid[bi][zi] << "\t" << _errorGrid[bi][zi] << std::endl;
		}
	}
}

void ClassGridThick_SCX::printDensityGrid(){
	for (int bi=0; bi<_bpoints;bi++){
		for (int zi=0; zi<_zpoints;zi++){
			double b = bi*_bstep; 
			double z = -_rmax + zi*_zstep;
            double r = sqrt(b*b + z*z);
			std::cout << bi << "\t" << zi << "\t" << _densf(r,_dens_param) << std::endl;
		}
	}
}


/** Interpolation of grid.
 *  Transform momentum to frame where mom is // z-axis
 *  to fetch the correct fsi factor
 *  @param x the coordinates in \f$r, \, \cos \theta\ \,\text{and} \phi$\f
 *  */
double ClassGridThick_SCX::getInterp(double x[],double& res,double& err){ // x is r,costheta,phi 
	//std::cout << "# ClassGridThick_SCX got interpolation request for " << x[0] << " " << x[1] << " " << x[2] << std::endl;
	assert( _grid != NULL); // make sure grid has been constructed
	_fp->setHitcoord(x[0],x[1],sqrt(1.-x[1]*x[1]),cos(x[2]),sin(x[2]) ); // set the hitcoordinates for the fast particle
	//std::cout << "#Interpolation for " << x[0] << ", " << x[1] << ", " << x[2] << "
	double z = _fp->getHitz();//x[0]*x[1]; // z = r \cos \theta
	double b = _fp->getHitbnorm(); // |b|
	double bi_t,zi_t; 
	double bf = modf(b/_bstep,&bi_t); // (b-bmin(=0))/bstep, split into integer (bi) and fractional (bf) part
	double zf = modf((z+_rmax)/_zstep,&zi_t); // (z-zmin(=rmax))/_zstep, split into integer (zi) and fractional (zf) part
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
	res =  bf_comp*( zf_comp*_grid[bi][zi] + zf*_grid[bi][zi+1] ) + bf*( zf_comp*_grid[bi+1][zi] + zf*_grid[bi+1][zi+1] );
    err =  bf_comp*( zf_comp*_errorGrid[bi][zi] + zf*_errorGrid[bi][zi+1] ) + bf*( zf_comp*_errorGrid[bi+1][zi] + zf*_errorGrid[bi+1][zi+1] );
}

/** integration struct **/
struct F_SCX {
	static double exec(double z,void* param){
		struct F_SCX p = *(struct F_SCX*) param;
		double r = sqrt( p.b*p.b + z*z);
		return -1.*p.densf(r,p.dens_param); // note the minus sign ;)
	}
	double (*densf)(double r,void* p); /**< density function to use **/
	double b;
    void* dens_param;
};

/**calculate the grid factors given (b,z) and density function densf
 * you still have to multiply this with the correct prefactors!
 */
void ClassGridThick_SCX::calcP(double b,double z,double& ret,double& err, double (*densf)(double,void*) ){
	assert( fabs(b) <= _rmax && fabs(z) <=_rmax ); // range check for b and z
	/** set up integrator struct **/
	struct F_SCX p;
	p.b          = b;
	p.densf      = _densf;
    p.dens_param = _dens_param;
	/** SET UP GSL 1D INTEGRATOR **/
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (2048);
	gsl_function F;
	F.function = F_SCX::exec;
	F.params = (void*) &p;
	double epsabs=0;
	double epsrel=1e-5;
	int key      =6;
	size_t limit = 2048; // this should never exceed size of the allocated integration workspace
	/** DO THE INTEGRATION **/
	gsl_integration_qag(&F, z, _rmax ,epsabs,epsrel,limit,key,w,&ret,&err); // integrate from z' to rmax
}
