//
// Class of pair ( n1 l1 j1 mj1 t1, n2 l2 j2 mj2 t2)
// 
#ifndef PAIR_H
#define PAIR_H
#include "newcoef.h"
#include <vector>
#include <complex>
#include <string>
using std::vector;
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

using namespace std;

#include "TDeuteron.h"

class Pair
{
public:
	/*
	 * Constructor: creates a pair with given quantum numbers. 
	 * two_t1 and two_t2, -1 or 1, determine if particle is n or p. 
	 */
	Pair(string wfname,  string path, int A, int n1, int l1, int two_ms1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_ms2, int two_j2, int two_mj2, int two_t2 );
	~Pair();

	/*
	 * Gives transformation coefficient <n(lS)jmj,NLML,TMT|n1l1j1mj1t1,n2l2j2mj2t2>
	 */
	double getcoef( int n, int l, int S, int j, int mj, int T, int MT, int N, int L, int ML );

	/*
	 * Gives normalised number of pairs with rel orbital momentum l
	 */
// 	double getRelPair( int l ) const;	


	/*
	 * Gives the com wavefunction \psi_{com,n(lS)jmjT}(\vt{R}) of pair with rel quantum numbers |n(lS)j mj, T >.
	 * Defined as \psi_{com,n(lS)jmjT}(\vt{R}) = \int d \vt{r} \Psi_{n(lS)jmjT}(\vt{r},\vt{R}),
	 * where 
	 * \Psi_{n{lS}jmjT} = \sum_{N L ML \Lambda} C |n(lS)jmj,NLML,TMT >
	 * with C the transformation coefficient given in getcoef().
	 * See getrelwf() for more information about how to get the total wavefunction 
	 */
	complex<double> getcomwf( double r1, double costh1, double phi1, double r2, double costh2, double phi2, int n, int l, int S, int j, int mj, int T) const;
	
	/*
	 * Gives relative wf \psi_{rel, n(lS)jmjT}(\vt{r}).
	 * Total wf with quantum numbers |n(lS)jmjT> and summed over the other quantum numbers (NLML,MT) is the product of the result from getcomwf() and this result.  
	 * REMARK: The transformatiom coefficient is not needed for this calculation
	 *
	 */ 
	complex<double> getrelwf( double r1, double costh1, double phi1, double r2, double costh2, double phi2, int n, int l, int la, int S, int j, int mj, int T) const;
// 	complex<double> getrelwf( double r, double costh, double phi, int n, int l, int S, int j, int mj, int T) const;

	/*
	 * Total wf of the pair. If the rel radial wf is the uncorrelated, this is equal to \phi_{j1}(r1) \phi_{j2}(r2). 
	 * Can easily change this function to for example wf of the pair for certain rel ang momentum l.
	 */
	complex<double> getwf(bool corr, double r1, double costh1, double phi1, double r2, double costh2, double phi2) const;


private:
	vector < Newcoef* > coeflist;
	int n1, l1, two_ms1, two_j1, two_mj1, two_t1;
	int n2, l2, two_ms2, two_j2, two_mj2, two_t2;
	RecMosh mosh;
	int A;
	bool coeflistmade;
	double nu;
	TDeuteron::Wavefunction *deuteronwf;
	void makecoeflist();
	double uncorrelatedradialwf(int n, int l, double r) const;
	double radialwf( int n, int l, int la, int S, int j, int T, double r) const;
	complex<double> angularwf( int l, int ml, double costh, double phi) const;
	complex<double> angularwfsplit( int l, int two_S, int two_ms, int two_j, int two_mj, double costh, double phi) const;
	complex<double> angularwf( int l, int two_S, int two_j, int two_mj, double costh, double phi) const;
	void getcomcoord( double r1, double costh1, double phi1, double r2, double costh2, double phi2, double* R, double* costhR, double* phiR ) const;
	void getrelcoord( double r1, double costh1, double phi1, double r2, double costh2, double phi2, double* r, double* costhr, double* phir ) const;
	double gamma( int two_n) const;
	double log( double x) const;
	double exponent( double x) const;
	double correlation(const double r) const;
};

#endif // PAIR_H
