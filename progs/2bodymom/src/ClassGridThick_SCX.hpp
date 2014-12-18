#ifndef CLASS_GRIDTHICK_SCX_HPP
#define CLASS_GRIDTHICK_SCX_HPP

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
 *  Using a semiclassical formalism (bye bye glauber)
 *  although it resembles the glauber formalism.
 *  The grid contains the probability for interaction.
 *  (The stored grid does NOT, only storing the energy
 *   independent integral).
 *
 *  */
class ClassGridThick_SCX {
	public:
		/** make sure both bpoints and zpoints are always larger than 2 plz **/
		ClassGridThick_SCX(MeanFieldNucleusThick* nuc,FastParticle* fp,int bpoints, int zpoints);
		~ClassGridThick_SCX();
		void constructGrid();
		void printGrid();
		void printDensityGrid();
		void addKnockoutParticle(int level); /**< add a knockout particle, on which no FSIs take place, only to correct for density changes */
		void clearKnockout();
		void calcq(double b,double z,double& ret,double& err, double (MeanFieldNucleusThick::*)(const double) const);
		double getInterp(double[]); // remember r,cos\theta and \phi
		double getInterp(numint::array<double,3>); // remember r,cos\theta and \phi
	private:
		MeanFieldNucleusThick* _nuc; // I use a pointer here because no copy constructor exists yet (needed because dyn. alloc. mem.)
		FastParticle* _fp;
		int _bpoints,_zpoints;
		double _bstep,_zstep;
		double _pdens_fctr,_ndens_fctr; // correction for mean density if more than one particle are knocked out
		double** _pgrid; /**< data goes here _pgrid is for a proton projectile, interacting with neutrons in the nucleus*/
		double** _ngrid; /**< data goes here _ngrid is for a neutron projectile, interacting with protons in the nucleus*/
		double** _perrorGrid; /**< integrator error estimations go here */
		double** _nerrorGrid; /**< integrator error estimations go here */
		char _filename[256]; /**< input/output file for the grid */
};


#endif // CLASS_GRIDTHICK_SCX_HPP
