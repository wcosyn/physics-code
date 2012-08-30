// Calculates the coefficients <n1 l1 j1 mj1 t1, n2 l2 j2 mj2 t2| n (l S) j mj, N L ML, T MT >.
// See notes about two and three body transformation, normalisation and anti-symmetrisation for exact expression.
#ifndef NEWCOEF_H
#define NEWCOEF_H

#include <cmath>
using std::pow;
using std::sqrt;
using std::fabs;

#include "recmosh.h"

#include <gsl/gsl_sf_coupling.h>

// INPUT:
// 	- self explaining but
// 		two_j1 = 2* j1
// 		two_mj1 = 2* mj1
// 		two_t1 = 2* mt1
// 		two_j2 = 2* j2
// 		two_t2 = 2* mt2 
// 		two_mj2 = 2* mj2
// 		RecMosh* mosh: Pointer to corresponding class RecMosh with moshinsky bracket
//
class Newcoef
{
private:
	int N;
	int L;
	int ML;
	int n;
	int l;
	int S;
	int j;
	int mj;
	int T;
	int MT;
	double coeff;
public:
	Newcoef( int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2, RecMosh* mosh, 
		   int N, int Lambda, int mLambda, int n, int l, int S, int j, int mj, int T, int MT, int two_ms1, int two_ms2);
	// Return the calculate coefficients 
	double getCoef() { return coeff;}
	int getN() {return N;}
	int getL() {return L;}
	int getML() {return ML;}
	int getn() { return n;}
	int getl() {return l;}
	int getS() { return S;}
	int getj() { return j;}
	int getmj() {return mj;}
	int getT() {return T;}
	int getMT() {return MT;}

};

#endif // NEWCOEF_H
