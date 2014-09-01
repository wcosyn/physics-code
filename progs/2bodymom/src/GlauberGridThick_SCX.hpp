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
		void calcFSI(double b,double z,complex<double>& res, complex<double>& err);
		complex<double> getInterp(double[]); // remember r,cos\theta and \phi
	private:
		MeanFieldNucleusThick* _nuc; // I use a pointer here because no copy constructor exists yet (needed because dyn. alloc. mem.)
		FastParticle _fp;
		int _bpoints,_zpoints;
		double _bstep,_zstep;
		complex<double>** _grid; /**< data goes here */
		complex<double>** _errorGrid; /**< integrator error estimations go here */
};


#endif // GLAUBERGRIDTHICK_SCX_HPP
