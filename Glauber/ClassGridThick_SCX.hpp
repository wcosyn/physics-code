#ifndef CLASS_GRIDTHICK_SCX_HPP
#define CLASS_GRIDTHICK_SCX_HPP

#include <cmath>
#include "FastParticle.hpp"
#include <iostream>
#include <cstdio>
#include <cassert>
#include <gsl/gsl_integration.h>
#include <fstream>

/** A class to calculate charge exchange
 *  using the mean field approximation.
 *  Using a semiclassical formalism (bye bye glauber)
 *  although it resembles the glauber formalism.
 *  The grid contains the probability for interaction.
 *  (The stored grid does NOT, only storing the energy
 *   independent integral).
 *  
 *  NOTE: For the same density function but different normalisations
 *        (e.g. normed to A,A-1,A-2,...) you should not calculate different
 *        grids constructed with different density function (_densf)'s but 
 *        simply set the _dens_fctr accordingly via set_dens_fctr().
 *        
 *        For example MeanFieldNucleusThick densities are normed to 1. in Wims Code.
 *        To get for example a density normed to total protons do
 *           $> get_dens_fctr( MeanFieldNucleusThick::getZ() );
 * 
 *  Class deliberately kept fairly modular. Only uses FastParticle.
 *  Please keep it that way. No hard coded densities etc. plz, thx.
 * 
 *  @AUTHOR Camille Colle
 *  @AFFILIATION Ghent University
 *  @DATE 21/09/2015
 *
 *  */
class ClassGridThick_SCX {
	public:
		/** make sure both bpoints and zpoints are always larger than 2 plz **/
		ClassGridThick_SCX(double (*densf)(double,void*),void* p,const char* densName,double rmax, FastParticle* fp,int bpoints, int zpoints,char* SHAREDIR);
		~ClassGridThick_SCX();
		void constructGrid();
		void printGrid();
		void printDensityGrid();
		void calcP(double b,double z,double& ret,double& err,double (*)(double,void*) );
		double getInterp(double[],double& res,double& err); // remember x is r,cos\theta and \phi
		
		double get_dens_fctr(){ return _dens_fctr; }; 
        void   set_dens_fctr(double d){ assert(_grid==NULL); _dens_fctr = d; }; /**< this should be called before constructGrid */
		
		
	private:
		double (*_densf)(double r,void* p); // density function, spherical symmetry assumed.
		void* _dens_param;    /**< density parameters (void*) pointer */
		double _rmax;         /**< radius of the system */
		FastParticle* _fp;    /**< fastparticle: particle on which to calculate SCX probabilities */
		int _bpoints,_zpoints;/**< number of points in b,z */
		double _bstep,_zstep; /**< step sizes in b and z for interpolation grid */
		double _dens_fctr;    /**< correction for density normalization. Density should be normed to A for A particles...*/
		double** _grid;       /**< SCX probability data goes here */
		double** _errorGrid;  /**< integrator error estimations go here */
		char _filename[256];  /**< input/output file for the grid */
};


#endif // CLASS_GRIDTHICK_SCX_HPP
