#ifndef GLAUBERGRIDTHICK_SCX_HPP
#define GLAUBERGRIDTHICK_SCX_HPP

#include "numint/numint.hpp"
#include <cmath>
#include "FastParticle.hpp"
#include "MeanFieldNucleusThick.hpp"
#include <complex>
using std::complex;
#include <iostream>
#include <cassert>

/** A class to calculate charge exchange
 *  using the mean field approximation.
 *  Doesn't inherit from the other Glauber
 *  classes. Sorry Wim...
 *
 *  */
class GlauberGridThick_SCX {
	public:
		/** make sure both bpoints and zpoints are always larger than 2 plz **/
		GlauberGridThick_SCX(MeanFieldNucleusThick* nuc,FastParticle& fp,int bpoints, int zpoints);
		~GlauberGridThick_SCX();
		void constructGlauberGrid();
		void printGrid();
		void printDensityGrid();
		void addKnockoutParticle(int level); /**< add a knockout particle, on which no FSIs take place, only to correct for density changes */
		void clearKnockout();
		void calcFSI(double b,double z,double& res,double& err);
		complex<double> getFrontFactor();
		complex<double> getInterp(double[]); // remember r,cos\theta and \phi
		complex<double> getInterp(numint::array<double,3>); // remember r,cos\theta and \phi
	private:
		MeanFieldNucleusThick* _nuc; // I use a pointer here because no copy constructor exists yet (needed because dyn. alloc. mem.)
		FastParticle _fp;
		int _bpoints,_zpoints;
		double _bstep,_zstep;
		double _pdens_fctr,_ndens_fctr; // correction for mean density if more than one particle are knocked out
		complex<double>** _grid; /**< data goes here */
		complex<double>** _errorGrid; /**< integrator error estimations go here */
		char _filename[256]; /**< input/output file for the grid */
		double _arbitraryPhase; /**< arbitrary phase **/
};


#endif // GLAUBERGRIDTHICK_SCX_HPP
