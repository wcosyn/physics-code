/*
 * GammaStructure.h
 *
 * This class represents a general complex 4x4 matrix
 * The matrix is decomposed in a complete set of gamma matrices:
 * - unity (U)
 * - gamma_5 (G5)
 * - gamma_mu (Gmu)
 * - gamma_5 gamma_mu (G5Gmu)
 * - gamma_0 gamma_i (Ai)
 * - gamma_5 gamma_0 gamma_i (Si)
 * We choose our gamma matrices in the chiral representation
 *
 * author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 */

#include "GammaStructure.h"
#include "constants.hpp"
#include <iostream>
using std::cout; using std::cerr; using std::endl;
using std::ostream;
#include <complex>
using std::complex;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

/* Multiplication tables multTable[][] and multTableCoeff[][]:
 * hold information about multiplication of two gamma matrices.
 * The product of (gamma matrix i) times (gamma matrix j)
 * is (multTableCoeff[i][j] * gamma matrix multTable[i][j])
 */

const int GammaStructure::multTable[16][16] =
  { {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
    {  1,  0, 10, 11, 12,  6,  5, 13, 14, 15,  2,  3,  4,  7,  8,  9 },
    {  2, 10,  0, 15, 14,  7, 13,  5, 12, 11,  1,  9,  8,  6,  4,  3 },
    {  3, 11, 15,  0, 13,  8, 14, 12,  5, 10,  9,  1,  7,  4,  6,  2 },
    {  4, 12, 14, 13,  0,  9, 15, 11, 10,  5,  8,  7,  1,  3,  2,  6 },
    {  5,  6,  7,  8,  9,  0,  1,  2,  3,  4, 13, 14, 15, 10, 11, 12 },
    {  6,  5, 13, 14, 15,  1,  0, 10, 11, 12,  7,  8,  9,  2,  3,  4 },
    {  7, 13,  5, 12, 11,  2, 10,  0, 15, 14,  6,  4,  3,  1,  9,  8 },
    {  8, 14, 12,  5, 10,  3, 11, 15,  0, 13,  4,  6,  2,  9,  1,  7 },
    {  9, 15, 11, 10,  5,  4, 12, 14, 13,  0,  3,  2,  6,  8,  7,  1 },
    { 10,  2,  1,  9,  8, 13,  7,  6,  4,  3,  0, 15, 14,  5, 12, 11 },
    { 11,  3,  9,  1,  7, 14,  8,  4,  6,  2, 15,  0, 13, 12,  5, 10 },
    { 12,  4,  8,  7,  1, 15,  9,  3,  2,  6, 14, 13,  0, 11, 10,  5 },
    { 13,  7,  6,  4,  3, 10,  2,  1,  9,  8,  5, 12, 11,  0, 15, 14 },
    { 14,  8,  4,  6,  2, 11,  3,  9,  1,  7, 12,  5, 10, 15,  0, 13 },
    { 15,  9,  3,  2,  6, 12,  4,  8,  7,  1, 11, 10,  5, 14, 13,  0 } };

const complex<double> GammaStructure::multTableCoeff[16][16] =
  { { complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0) },

    { complex<double>( 1.0,0.0), complex<double>( 1.0,0.0),
      complex<double>( 1.0,0.0), complex<double>( 1.0,0.0),
      complex<double>( 1.0,0.0), complex<double>(-1.0,0.0),
      complex<double>(-1.0,0.0), complex<double>(-1.0,0.0),
      complex<double>(-1.0,0.0), complex<double>(-1.0,0.0),
      complex<double>( 1.0,0.0), complex<double>( 1.0,0.0),
      complex<double>( 1.0,0.0), complex<double>(-1.0,0.0),
      complex<double>(-1.0,0.0), complex<double>(-1.0,0.0) },
      
    { complex<double>( 1.0,0.0), complex<double>(-1.0, 0.0),
      complex<double>(-1.0,0.0), complex<double>( 0.0,-1.0),
      complex<double>( 0.0,1.0), complex<double>(-1.0, 0.0),
      complex<double>( 1.0,0.0), complex<double>( 1.0, 0.0),
      complex<double>( 0.0,1.0), complex<double>( 0.0,-1.0),
     complex<double>( 1.0,0.0), complex<double>( 0.0,-1.0),
      complex<double>( 0.0,1.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0,1.0), complex<double>( 0.0,-1.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0,-1.0), complex<double>(-1.0, 0.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0,-1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0, 1.0), complex<double>( 1.0, 0.0),
      complex<double>( 0.0,-1.0), complex<double>( 0.0,-1.0),
      complex<double>(-1.0, 0.0), complex<double>( 0.0, 1.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0,0.0),
      complex<double>( 0.0,-1.0), complex<double>( 0.0,1.0),
      complex<double>(-1.0, 0.0), complex<double>(-1.0,0.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0,1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0,0.0),
      complex<double>( 0.0,-1.0), complex<double>( 0.0,1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0,1.0),
      complex<double>( 0.0,-1.0), complex<double>(-1.0,0.0) },

    { complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0),
      complex<double>(1.0,0.0), complex<double>(1.0,0.0) },

    { complex<double>( 1.0,0.0), complex<double>( 1.0,0.0),
      complex<double>( 1.0,0.0), complex<double>( 1.0,0.0),
      complex<double>( 1.0,0.0), complex<double>(-1.0,0.0),
      complex<double>(-1.0,0.0), complex<double>(-1.0,0.0),
      complex<double>(-1.0,0.0), complex<double>(-1.0,0.0),
      complex<double>( 1.0,0.0), complex<double>( 1.0,0.0),
      complex<double>( 1.0,0.0), complex<double>(-1.0,0.0),
      complex<double>(-1.0,0.0), complex<double>(-1.0,0.0) },

    { complex<double>( 1.0,0.0), complex<double>(-1.0, 0.0),
      complex<double>(-1.0,0.0), complex<double>( 0.0,-1.0),
      complex<double>( 0.0,1.0), complex<double>(-1.0, 0.0),
      complex<double>( 1.0,0.0), complex<double>( 1.0, 0.0),// 7
      complex<double>( 0.0,1.0), complex<double>( 0.0,-1.0),
      complex<double>( 1.0,0.0), complex<double>( 0.0,-1.0),
      complex<double>( 0.0,1.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0,1.0), complex<double>( 0.0,-1.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0,-1.0), complex<double>(-1.0, 0.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0,-1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0, 1.0), complex<double>( 1.0, 0.0),
      complex<double>( 0.0,-1.0), complex<double>( 0.0,-1.0),
      complex<double>(-1.0, 0.0), complex<double>( 0.0, 1.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0,0.0),
      complex<double>( 0.0,-1.0), complex<double>( 0.0,1.0),
      complex<double>(-1.0, 0.0), complex<double>(-1.0,0.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0,1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0,0.0),
      complex<double>( 0.0,-1.0), complex<double>( 0.0,1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0,1.0),
      complex<double>( 0.0,-1.0), complex<double>(-1.0,0.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>(-1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0, 0.0),
      complex<double>(-1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0,-1.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 1.0, 0.0),
      complex<double>(-1.0, 0.0), complex<double>( 0.0,-1.0), // 11
      complex<double>(-1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0, 1.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0),
      complex<double>(-1.0, 0.0), complex<double>( 1.0, 0.0),
      complex<double>(-1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0, 0.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>(-1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0, 0.0),
      complex<double>(-1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0,-1.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 1.0, 0.0),
      complex<double>(-1.0, 0.0), complex<double>( 0.0,-1.0),
      complex<double>(-1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0, 1.0) },

    { complex<double>( 1.0, 0.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0),
      complex<double>(-1.0, 0.0), complex<double>( 1.0, 0.0),
      complex<double>(-1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>(-1.0, 0.0),
      complex<double>( 0.0, 1.0), complex<double>( 0.0,-1.0),
      complex<double>( 1.0, 0.0), complex<double>( 0.0, 1.0),
      complex<double>( 0.0,-1.0), complex<double>( 1.0, 0.0) } };



//-------------------------------------------------------------------------

// Constructors
// ------------

/* The most general constructor
 * which specifies all coefficients.
 * Standard values are zero from component 7 onwards
 */
GammaStructure::GammaStructure(const complex<double>& c1,
			       const complex<double>& c2,
			       const complex<double>& c3, 
			       const complex<double>& c4,
			       const complex<double>& c5, 
			       const complex<double>& c6,
			       const complex<double>& c7, 
			       const complex<double>& c8,
			       const complex<double>& c9, 
			       const complex<double>& c10,
			       const complex<double>& c11, 
			       const complex<double>& c12,
			       const complex<double>& c13, 
			       const complex<double>& c14,
			       const complex<double>& c15, 
			       const complex<double>& c16)
  : nrComp(16),
    comp(new int[16]),
    compCoeff(new complex<double>[16])
{
  // fill comp[]
  for(int i=0; i<nrComp; i++)
    comp[i] = i;
  
  // Store coefficients
  compCoeff[0] = c1;
  compCoeff[1] = c3;
  compCoeff[2] = c4;
  compCoeff[3] = c5;
  compCoeff[4] = c6;
  compCoeff[5] = c2;
  compCoeff[6] = c7;
  compCoeff[7] = c8;
  compCoeff[8] = c9;
  compCoeff[9] = c10;
  compCoeff[10] = c11;
  compCoeff[11] = c12;
  compCoeff[12] = c13;
  compCoeff[13] = c14;
  compCoeff[14] = c15;
  compCoeff[15] = c16;
  
  crop();
    
}

//-------------------------------------------------------------------------

// Single component and default constructor

GammaStructure::GammaStructure(const complex<double>& c1)
  : nrComp(1),
    comp(new int[1]),
    compCoeff(new complex<double>[1])
{
  // store coefficients
  comp[0] = 0;
  compCoeff[0] = c1;
}

//-------------------------------------------------------------------------

// 2 component constructor

GammaStructure::GammaStructure(const complex<double>& c1,
			       const complex<double>& c2) 
  : nrComp(2),
    comp(new int[2]),
    compCoeff(new complex<double>[2])
{
  // store coefficients
  comp[0] = 0;
  comp[1] = 5;
  compCoeff[0] = c1;
  compCoeff[1] = c2;
  
  crop();
}

//-------------------------------------------------------------------------

// constructor with maximum 6 components

GammaStructure::GammaStructure(const complex<double>& c1, 
			       const complex<double>& c2,
			       const complex<double>& c3, 
			       const complex<double>& c4,
			       const complex<double>& c5, 
			       const complex<double>& c6) 
  : nrComp(6),
    comp(new int[6]),
    compCoeff(new complex<double>[6])
{
  // store coefficients
  for(int i=0; i<6; ++i)
    comp[i] = i;
    
  compCoeff[0] = c1;
  compCoeff[1] = c3;
  compCoeff[2] = c4;
  compCoeff[3] = c5;
  compCoeff[4] = c6;
  compCoeff[5] = c2;
  
  crop();
}

//-------------------------------------------------------------------------

// Constructor when all member data are given as arguments
// This constructor is private. It is used in operator*()

GammaStructure::GammaStructure(int nrComp_init,
			       int* comp_init,
			       complex<double>* compCoeff_init)
  : nrComp( nrComp_init ),
    comp( comp_init ),
    compCoeff( compCoeff_init ) 
{
    // empty body
}

//-------------------------------------------------------------------------

// Copy constructor

GammaStructure::GammaStructure(const GammaStructure& toCopy) 
  : nrComp( toCopy.nrComp ),
    comp( new int[nrComp] ),
    compCoeff( new complex<double>[nrComp] )
{
  // copy contents of comp[] and compCoeff[]
  for(int i=0; i<nrComp; i++) {
    comp[i] = toCopy.comp[i];
    compCoeff[i] = toCopy.compCoeff[i];
  }
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

/* Crop function
 * -------------
 * it is useless to store components whose coefficients are zero.
 * We check for this and delete them when necessary */

void GammaStructure::crop() 
{
  /* Count the number of coefficients equal to zero (store in countZeros)
   * and store the position of the non-zero components in nonZeroComp_index[]
   */
  int countZeros = 0;
  int nonZeroComp_index[16];
  
  static complex<double> zero(0.0, 0.0);
  
  for(int i=0; i<nrComp; i++) {
    if(compCoeff[i]==zero)
      countZeros++;
    else
      nonZeroComp_index[i-countZeros] = i;
  }
  
  /* if all components are zero, the GammaStructure = zero.
   * This is stored by putting the coeff. of the unity matrix to zero. */
  if(nrComp==countZeros) {
    // Release memory
    delete[] comp;
    delete[] compCoeff;
    
    // reallocate
    comp = new int[1];
    compCoeff = new complex<double>[1];
    
    // store appropriate values
    comp[0] = 0;
    compCoeff[0] = 0.0;
    nrComp = 1;
  }
  
  /* If the GammaStructure is non-zero, but some coefficients
   * are.
   */
  else if(countZeros>0) {
    // allocate memory for a reduced comp[] and compCoeff[]
    int* new_comp = new int[nrComp-countZeros];
    complex<double>* new_compCoeff =
      new complex<double>[nrComp-countZeros];
    
    /* Store non-zero components in new_comp[] and new_compCoeff[]
     * the location of those non-zero components was stored
     * in nonZeroComp_index[] */
    for(int i=0; i<nrComp-countZeros; i++) {
      new_comp[i] = comp[ nonZeroComp_index[i] ];
      new_compCoeff[i] = compCoeff[ nonZeroComp_index[i] ];
    } // end loop over new components
    
    // release memory of old component arrays
    delete[] comp;
    delete[] compCoeff;
    
    // link comp[] and compCoeff[] to the reduced arrays
    comp = new_comp;
    compCoeff = new_compCoeff;
    
    // finaly reset the number of components
    nrComp -= countZeros;
  }
  
  // if countZeros==0, nothing needs to be done
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Destructor
// ----------
GammaStructure::~GammaStructure() {
  delete[] comp;
  delete[] compCoeff;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Assignment
// ----------

GammaStructure& GammaStructure::operator=(const GammaStructure& toCopy) 
{
  if(this != &toCopy) // Avoid self-assignment
    {  
      // if the number of components differs, reallocate memory
      if(nrComp != toCopy.nrComp) {
	delete[] comp;
	delete[] compCoeff;
        
	nrComp = toCopy.nrComp;
	comp = new int[nrComp];
	compCoeff = new complex<double>[nrComp];
      }
      
      // copy contents of comp[] and compCoeff[]
      for(int i=0; i<nrComp; i++) {
	comp[i] = toCopy.comp[i];
	compCoeff[i] = toCopy.compCoeff[i];
      }
    }
  
  return *this;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Addition
// --------

GammaStructure GammaStructure::operator+(const GammaStructure& right) const
{
  /* If both GammaStructures have a single component, addition
   * is straight forward */
  if( nrComp * right.nrComp == 1 ) {
    if(comp[0]==right.comp[0]) {
      GammaStructure sum(1, new int[1], new complex<double>[1]);
      sum.comp[0] = comp[0];
      sum.compCoeff[0] = compCoeff[0] + right.compCoeff[0];
      
      return sum;
    }
    else {
      GammaStructure sum(2, new int[2], new complex<double>[2]);
      
      if(comp[0]>right.comp[0]) {
	sum.comp[0] = right.comp[0];
	sum.compCoeff[0] = right.compCoeff[0];
	sum.comp[1] = comp[0];
	sum.compCoeff[1] = compCoeff[0];
      }
      else {
	sum.comp[0] = comp[0];
	sum.compCoeff[0] = compCoeff[0];
	sum.comp[1] = right.comp[0];
	sum.compCoeff[1] = right.compCoeff[0];
      }
      
      return sum;
    }
  }
  
  else {
    /* initialize the GammaStructure that eventualy will hold the result.
     * It has 16 components that are all zero
     * It's done with the private constructor because the default constructor
     * would call crop() prematurely */
    GammaStructure sum(16, new int[16], new complex<double>[16]);
    
    static complex<double> zero(0.0);
    
    for(int i=0; i<sum.nrComp; i++) {
      sum.comp[i] = i;
      sum.compCoeff[i] = zero;
    }
    
    /* Store all non-zero components of the left side into sum
     * This only works because we know for sure that sum has all
     * 16 components */
    for(int i=0; i<nrComp; i++)
      sum.compCoeff[comp[i]] = compCoeff[i];
    
    /* Finally we add the coefficients of the right side to sum
     * Again this is only possible because sum has all 16 components */
    for(int i=0; i<right.nrComp; i++)
      sum.compCoeff[right.comp[i]] += right.compCoeff[i];
    
    // Remove all redundant components
    sum.crop();
    
    return sum;
  }
  
  // dummy return
  return GammaStructure(1.0);
}

//-------------------------------------------------------------------------

GammaStructure GammaStructure::operator-(const GammaStructure& right) const 
{
  /* If both GammaStructures have a single component, addition
   * is straight forward */
  if( nrComp * right.nrComp == 1 ) {
    if(comp[0]==right.comp[0]) {
      GammaStructure sum(1, new int[1], new complex<double>[1]);
      sum.comp[0] = comp[0];
      sum.compCoeff[0] = compCoeff[0] - right.compCoeff[0];
      
      return sum;
    }
    else {
      GammaStructure sum(2, new int[2], new complex<double>[2]);
      
      if(comp[0]>right.comp[0]) {
	sum.comp[0] = right.comp[0];
	sum.compCoeff[0] = -1.*right.compCoeff[0];
	sum.comp[1] = comp[0];
	sum.compCoeff[1] = compCoeff[0];
      }
      else {
	sum.comp[0] = comp[0];
	sum.compCoeff[0] = compCoeff[0];
	sum.comp[1] = right.comp[0];
	sum.compCoeff[1] = -1.*right.compCoeff[0];
      }
      
      return sum;
    }
  }
  
  else {
    /* initialize the GammaStructure that eventualy will hold the result.
     * It has 16 components that are all zero
     * It's done with the private constructor because the default constructor
     * would call crop() prematurely */
    GammaStructure sum(16, new int[16], new complex<double>[16]);
    
    static complex<double> zero(0.0);
    
    for(int i=0; i<sum.nrComp; i++) {
      sum.comp[i] = i;
      sum.compCoeff[i] = zero;
    }
    
    /* Store all non-zero components of the left side into sum
     * This only works because we know for sure that sum has all
     * 16 components */
    for(int i=0; i<nrComp; i++)
      sum.compCoeff[comp[i]] = compCoeff[i];
    
    /* Finally we add the coefficients of the right side to sum
     * Again this is only possible because sum has all 16 components */
    for(int i=0; i<right.nrComp; i++)
      sum.compCoeff[right.comp[i]] -= right.compCoeff[i];
    
    // Remove all redundant components
    sum.crop();
    
    return sum;
  }
  
  // dummy return
  return GammaStructure(1.0);
}

//-------------------------------------------------------------------------

GammaStructure& GammaStructure::operator+=(const GammaStructure& right)
{  
  /* If both GammaStructures have a single component, addition
   * is straight forward */
  if( nrComp * right.nrComp == 1 ) {
    if(comp[0]==right.comp[0]) {
      // Do the addition
      int *sumcomp = new int[1];
      complex<double> *sumcompCoeff = new complex<double>[1];
      sumcomp[0] = comp[0];
      sumcompCoeff[0] = compCoeff[0] + right.compCoeff[0];
      
      // Store the result in the left side
      nrComp = 1;
      delete[] comp;
      comp = sumcomp;
      delete[] compCoeff;
      compCoeff = sumcompCoeff;
    }
    else {
      int *sumcomp = new int[2];
      complex<double> *sumcompCoeff = new complex<double>[2];
      
      if(comp[0]>right.comp[0]) {
	sumcomp[0] = right.comp[0];
	sumcompCoeff[0] = right.compCoeff[0];
	sumcomp[1] = comp[0];
	sumcompCoeff[1] = compCoeff[0];
      }
      else {
	sumcomp[0] = comp[0];
	sumcompCoeff[0] = compCoeff[0];
	sumcomp[1] = right.comp[0];
	sumcompCoeff[1] = right.compCoeff[0];
      }

      // Store the result in the left side
      nrComp = 2;
      delete[] comp;
      comp = sumcomp;
      delete[] compCoeff;
      compCoeff = sumcompCoeff;
    }
  }
  
  else {
    /* initialize the GammaStructure that eventualy will hold the result.
     * It has 16 components that are all zero
     * It's done with the private constructor because the default constructor
     * would call crop() prematurely */
    int *sumcomp = new int[16];
    complex<double> *sumcompCoeff = new complex<double>[16];
    
    static complex<double> zero(0.0);
    
    for(int i=0; i<16; i++) {
      sumcomp[i] = i;
      sumcompCoeff[i] = zero;
    }
    
    /* Store all non-zero components of the left side into sum
     * This only works because we know for sure that sum has all
     * 16 components */
    for(int i=0; i<nrComp; i++)
      sumcompCoeff[comp[i]] = compCoeff[i];
    
    /* Finally we add the coefficients of the right side to sum
     * Again this is only possible because sum has all 16 components */
    for(int i=0; i<right.nrComp; i++)
      sumcompCoeff[right.comp[i]] += right.compCoeff[i];
    
    // Store the result in the left side
    nrComp = 16;
    delete[] comp;
    comp = sumcomp;
    delete[] compCoeff;
    compCoeff = sumcompCoeff;
  
    // Remove all redundant components
    crop();
  }
  
  return *this;
}

//-------------------------------------------------------------------------

GammaStructure& GammaStructure::operator-=(const GammaStructure& right)
{
  /* If both GammaStructures have a single component, addition
   * is straight forward */
  if( nrComp * right.nrComp == 1 ) {
    if(comp[0]==right.comp[0]) {
      // Do the addition
      int *sumcomp = new int[1];
      complex<double> *sumcompCoeff = new complex<double>[1];
      sumcomp[0] = comp[0];
      sumcompCoeff[0] = compCoeff[0] - right.compCoeff[0];
      
      // Store the result in the left side
      nrComp = 1;
      delete[] comp;
      comp = sumcomp;
      delete[] compCoeff;
      compCoeff = sumcompCoeff;
    }
    else {
      int *sumcomp = new int[2];
      complex<double> *sumcompCoeff = new complex<double>[2];
      
      if(comp[0]>right.comp[0]) {
	sumcomp[0] = right.comp[0];
	sumcompCoeff[0] = -1.*right.compCoeff[0];
	sumcomp[1] = comp[0];
	sumcompCoeff[1] = compCoeff[0];
      }
      else {
	sumcomp[0] = comp[0];
	sumcompCoeff[0] = compCoeff[0];
	sumcomp[1] = right.comp[0];
	sumcompCoeff[1] = -1.*right.compCoeff[0];
      }

      // Store the result in the left side
      nrComp = 2;
      delete[] comp;
      comp = sumcomp;
      delete[] compCoeff;
      compCoeff = sumcompCoeff;
    }
  }
  
  else {
    /* initialize the GammaStructure that eventualy will hold the result.
     * It has 16 components that are all zero
     * It's done with the private constructor because the default constructor
     * would call crop() prematurely */
    int *sumcomp = new int[16];
    complex<double> *sumcompCoeff = new complex<double>[16];
    
    static complex<double> zero(0.0);
    
    for(int i=0; i<16; i++) {
      sumcomp[i] = i;
      sumcompCoeff[i] = zero;
    }
    
    /* Store all non-zero components of the left side into sum
     * This only works because we know for sure that sum has all
     * 16 components */
    for(int i=0; i<nrComp; i++)
      sumcompCoeff[comp[i]] = compCoeff[i];
    
    /* Finally we add the coefficients of the right side to sum
     * Again this is only possible because sum has all 16 components */
    for(int i=0; i<right.nrComp; i++)
      sumcompCoeff[right.comp[i]] -= right.compCoeff[i];
    
    // Store the result in the left side
    nrComp = 16;
    delete[] comp;
    comp = sumcomp;
    delete[] compCoeff;
    compCoeff = sumcompCoeff;
  
    // Remove all redundant components
    crop();
  }
  
  return *this;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Multiplication
// --------------

// (GammaStructure * double)

GammaStructure GammaStructure::operator*(double factor) const
{
  // make a copy of the GammaStructure;
  GammaStructure product = *this;
  
  product *= factor;
  
  return product;
}

//-------------------------------------------------------------------------

// (GammaStructure * complex)

GammaStructure GammaStructure::operator*(const complex<double>& factor) const
{
    // make a copy of the GammaStructure;
    GammaStructure product = *this;
    
    product *= factor;
    
    return product;
}

//-------------------------------------------------------------------------

// (GammaStructure * GammaStructure)

GammaStructure GammaStructure::operator*(const GammaStructure& right) const 
{
  /* If both GammaStructures have a single component, multiplication
   * is straight forward */
  if( nrComp * right.nrComp == 1 ) {
    GammaStructure product; // empty GammaStructure
    
    // product.nrComp is already =1
    product.comp[0] = multTable[comp[0]][right.comp[0]];
    product.compCoeff[0] = multTableCoeff[comp[0]][right.comp[0]]
      * compCoeff[0] * right.compCoeff[0];
    
    return product;
  }
  
  // if the product has multiple components
  else {
    /* initialize the GammaStructure that eventualy will hold the result.
     * It has 16 components that are all zero
     * It's done with the private constructor because the default constructor
     * would call crop() prematurely */
    GammaStructure product(16, new int[16], new complex<double>[16]);
    
    for(int i=0; i<product.nrComp; i++) {
      product.comp[i] = i;
      product.compCoeff[i] = 0.0;
    }
    
    /* We loop over all non-zero components of the left en right side.
     * The resulting component is stored in multTable[][].
     * The coefficient of this component is the product of the coefficients
     * of the left and right side multiplied by the value stored in
     * multTalbeCoeff[][]. */
    for(int i=0; i<nrComp; i++) {
      for(int j=0; j<right.nrComp; j++) {
	product.compCoeff[ multTable[comp[i]][right.comp[j]] ] +=
	  compCoeff[i] * right.compCoeff[j] *
	  multTableCoeff[comp[i]][right.comp[j]];
      }
    }
    
    // Finally, we remove all redundant components
    product.crop();
    
    // and return the result
    return product;
  }
  
  // dummy return
  return GammaStructure(1.0);
}

//-------------------------------------------------------------------------

// (double * GammaStructure)

GammaStructure operator*(double f, const GammaStructure& g) 
{
  return g*f; // this product is commutative
}

//-------------------------------------------------------------------------

// (complex * GammaStructure)

GammaStructure operator*(const complex<double>& f, const GammaStructure& g) 
{
  return g*f; // this product is commutative
}

//-------------------------------------------------------------------------

GammaStructure& GammaStructure::operator*=(double factor)
{
  // Multiply all non-zero components with the complex double
  for(int i=0; i<nrComp; i++)
    compCoeff[i] *= factor;
  
  return *this;
}

//-------------------------------------------------------------------------

GammaStructure& GammaStructure::operator*=(const complex<double>& factor)
{
  // Multiply all non-zero components with the complex double
  for(int i=0; i<nrComp; i++)
    compCoeff[i] *= factor;
  
  return *this;
}

//-------------------------------------------------------------------------

GammaStructure& GammaStructure::operator*=(const GammaStructure& right)
{
  /* If both GammaStructures have a single component, multiplication
   * is straight forward */
  if( nrComp * right.nrComp == 1 ) {

    // product.nrComp is already =1
    compCoeff[0] = multTableCoeff[comp[0]][right.comp[0]]
      * compCoeff[0] * right.compCoeff[0];
    comp[0] = multTable[comp[0]][right.comp[0]];

  }
  
  // if the product has multiple components
  else {
    /* initialize the GammaStructure comp and compCoeff arrays
     * that eventualy will hold the result.
     * It has 16 components that are all zero */
    int *productcomp = new int[16];
    complex<double> *productcompCoeff = new complex<double>[16];

    for(int i=0; i<16; i++) {
      productcomp[i] = i;
      productcompCoeff[i] = 0.0;
    }
    
    /* We loop over all non-zero components of the left en right side.
     * The resulting component is stored in multTable[][].
     * The coefficient of this component is the product of the coefficients
     * of the left and right side multiplied by the value stored in
     * multTalbeCoeff[][]. */
    for(int i=0; i<nrComp; i++) {
      for(int j=0; j<right.nrComp; j++) {
	productcompCoeff[ multTable[comp[i]][right.comp[j]] ] +=
	  compCoeff[i] * right.compCoeff[j] *
	  multTableCoeff[comp[i]][right.comp[j]];
      }
    }
    
    // Store the result in the left side
    nrComp = 16;
    delete[] comp;
    comp = productcomp;
    delete[] compCoeff;
    compCoeff = productcompCoeff;

    // Finally, we remove all redundant components
    crop();
  }
 
  return *this;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Evaluate the structure and return the 4x4 matrix

Matrix<4, 4> GammaStructure::value() const 
{
  /* 4x4 matrix that will hold the result
   * initialized to zero */
  Matrix<4, 4> result;
  
  // Evaluate all non-zero components
  for(int i=0; i<nrComp; i++)
    result.multiAdd(compCoeff[i], gammaBasis[ comp[i] ]);
  
  return result;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Matrix operations
// -----------------

// Hermitian
GammaStructure GammaStructure::Hermitian() const
{
  GammaStructure result = *this;

  // The coefficients of G1,G2,G3 and G5G0 change sign after taking the hermitian.
  for(int i=0; i<nrComp; ++i) {
    if( comp[i]==2 || comp[i]==3 || comp[i]==4 || comp[i]==6 )
      result.compCoeff[i] = complex<double>( -result.compCoeff[i].real(),
					     result.compCoeff[i].imag());
    else
      result.compCoeff[i] = complex<double>( result.compCoeff[i].real(),
					     -result.compCoeff[i].imag());
  }

  return result;
}

//-------------------------------------------------------------------------

// Conjugate
GammaStructure GammaStructure::Conjugate() const
{
  GammaStructure result = *this;

  // The coefficients of G5,G0G1,G0G2 and G0G3 change sign after taking the conjugate
  for(int i=0; i<nrComp; ++i) {
    if( comp[i]==5 || comp[i]==10 || comp[i]==11 || comp[i]==12 )
      result.compCoeff[i] = complex<double>( -result.compCoeff[i].real(),
					     result.compCoeff[i].imag());
    else
      result.compCoeff[i] = complex<double>( result.compCoeff[i].real(),
					     -result.compCoeff[i].imag());
  }
  
  return result;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Output
// ------

// Prints the structure in [..,..,..] format
void GammaStructure::print(void) const 
{
  cout << *this << endl;
}

//-------------------------------------------------------------------------

// Overwrite << operator for standard output
ostream& operator<<(ostream& output, const GammaStructure& gamma) 
{
  output << "[";
  
  int index = 0;
  
  for(int i=0; i<15; i++) {
    if(gamma.comp[index] == i) {
      if (abs(gamma.compCoeff[index]) >= STRANGEUFLOW)
	output << gamma.compCoeff[index];
      else
	output << complex<double>(0.0);
      if(index<(gamma.nrComp-1))
	index++;
    }
    else
      output << complex<double>(0.0);
    output << ",";
  }
  if(gamma.comp[index] == 15 && abs(gamma.compCoeff[index]) >=  STRANGEUFLOW)
    output << gamma.compCoeff[index];
  else
    output << complex<double>(0.0);
  
  output << "]";
  
  return output;
}


