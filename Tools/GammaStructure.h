/*
 * GammaStructure.h
 *
 * This class represents a general complex 4x4 matrix.
 * A complex 4x4 matrix can be decomposed in a complete set 
 * of gamma matrices.
 *
 * Therefore we implement a GammaStructure as the sum of 16 lineary
 * independent gamma matrices. Each of them has a numerical code:
 *  0- unity (=U)
 *  1- gamma^0 (=G0)
 *  2- gamma^1 (=G1)
 *  3- gamma^2 (=G2)
 *  4- gamma^3 (=G3)
 *  5- gamma^5 (=G5)
 *  6- gamma^5 * gamma^0 (=G5G0)
 *  7- gamma^5 * gamma^1 (=G5G1)
 *  8- gamma^5 * gamma^2 (=G5G2)
 *  9- gamma^5 * gamma^3 (=G5G3)
 * 10- gamma^0 * gamma^1 (=A1)
 * 11- gamma^0 * gamma^2 (=A2)
 * 12- gamma^0 * gamma^3 (=A3)
 * 13- gamma^5 * gamma^0 * gamma^1 (=S1)
 * 14- gamma^5 * gamma^0 * gamma^2 (=S2)
 * 15- gamma^5 * gamma^0 * gamma^3 (=S3)
 * We choose the Dirac representation for the gamma matrices.
 * They are are defined below as objects of type Matrix<4,4>.
 *
 * Internally, we only store those components with a non-zero coefficient
 * - nrComp(int) holds the number of non-zero components
 * - comp(int*) is an variable-length array that holds the code of all
 *   gamma matrices with non-zero coefficients
 * - compCoeff(complex<double>*) holds the complex coefficients that
 *   go with the gamma matrices stored in comp[].
 *
 * author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 */

#ifndef GAMMASTRUCTURE_H
#define GAMMASTRUCTURE_H

#include <iostream>
#include <complex>
#include "Matrix.h"

// ********************************************************** //
// First we define the gamma matrices (Dirac representation) //
// ********************************************************** //

const Matrix<4,4> U(1.0, 0.0, 0.0, 0.0,
		    0.0, 1.0, 0.0, 0.0,
		    0.0, 0.0, 1.0, 0.0,
		    0.0, 0.0, 0.0, 1.0);

const Matrix<4,4> G0(1.0, 0.0,  0.0,  0.0,
		     0.0, 1.0,  0.0,  0.0,
		     0.0, 0.0, -1.0,  0.0,
		     0.0, 0.0,  0.0, -1.0);

const Matrix<4,4> G1( 0.0, 0.0, 0.0, 1.0,
		      0.0, 0.0, 1.0, 0.0,
		      0.0,-1.0, 0.0, 0.0,
		     -1.0, 0.0, 0.0, 0.0);

const Matrix<4,4> G2(0.0,0.0,0.0,std::complex<double>(0,-1.0),
		     0.0,0.0,std::complex<double>(0,1.0),0.0,
		     0.0,std::complex<double>(0,1.0),0.0,0.0,
		     std::complex<double>(0,-1.0),0.0,0.0,0.0);

const Matrix<4,4> G3( 0.0, 0.0, 1.0, 0.0,
		      0.0, 0.0, 0.0,-1.0,
		     -1.0, 0.0, 0.0, 0.0,
		      0.0, 1.0, 0.0, 0.0);

const Matrix<4,4> G5(0.0, 0.0, 1.0, 0.0,
		     0.0, 0.0, 0.0, 1.0,
		     1.0, 0.0, 0.0, 0.0,
		     0.0, 1.0, 0.0, 0.0); // i * G0*G1*G2*G3

const Matrix<4,4> G5G0(0.0, 0.0,-1.0, 0.0,
		       0.0, 0.0, 0.0,-1.0,
		       1.0, 0.0, 0.0, 0.0,
		       0.0, 1.0, 0.0, 0.0); // G5*G0

const Matrix<4,4> G5G1( 0.0,-1.0, 0.0, 0.0,
		       -1.0, 0.0, 0.0, 0.0,
		        0.0, 0.0, 0.0, 1.0,
		        0.0, 0.0, 1.0, 0.0); // G5*G1

const Matrix<4,4> G5G2(0.0,std::complex<double>(0,1.0),0.0,0.0,
		       std::complex<double>(0,-1.0),0.0,0.0,0.0,
		       0.0,0.0,0.0,std::complex<double>(0,-1.0),
		       0.0,0.0,std::complex<double>(0,1.0),0.0); // G5*G2

const Matrix<4,4> G5G3(-1.0, 0.0, 0.0, 0.0,
		        0.0, 1.0, 0.0, 0.0,
		        0.0, 0.0, 1.0, 0.0,
		        0.0, 0.0, 0.0,-1.0); // G5*G3

const Matrix<4,4> A1( 0.0, 0.0, 0.0, 1.0,
		      0.0, 0.0, 1.0, 0.0,
		      0.0, 1.0, 0.0, 0.0,
		      1.0, 0.0, 0.0, 0.0); // G0*G1

const Matrix<4,4> A2(0.0,0.0,0.0,std::complex<double>(0,-1.0),
		     0.0,0.0,std::complex<double>(0,1.0),0.0,
		     0.0,std::complex<double>(0,-1.0),0.0,0.0,
		     std::complex<double>(0,1.0),0.0,0.0,0.0); // G0*G2

const Matrix<4,4> A3( 0.0, 0.0, 1.0, 0.0,
		      0.0, 0.0, 0.0,-1.0,
		      1.0, 0.0, 0.0, 0.0,
		      0.0,-1.0, 0.0, 0.0); // G0*G3

const Matrix<4,4> S1( 0.0, 1.0, 0.0, 0.0,
		      1.0, 0.0, 0.0, 0.0,
		      0.0, 0.0, 0.0, 1.0,
		      0.0, 0.0, 1.0, 0.0); // G5*G0*G1

const Matrix<4,4> S2(0.0,std::complex<double>(0,-1.0),0.0,0.0,
		     std::complex<double>(0,1.0),0.0,0.0,0.0,
		     0.0,0.0,0.0,std::complex<double>(0,-1.0),
		     0.0,0.0,std::complex<double>(0,1.0),0.0); // G5*G0*G2

const Matrix<4,4> S3(1.0, 0.0, 0.0, 0.0,
		     0.0,-1.0, 0.0, 0.0,
		     0.0, 0.0, 1.0, 0.0,
		     0.0, 0.0, 0.0,-1.0); // G5*G0*G3


// *********************************************** //
// Our basis of lineary independent gamma matrices //
// *********************************************** //
const Matrix<4,4> gammaBasis[16] =
  { U, G0, G1, G2, G3, G5, G5G0, G5G1, G5G2, G5G3,
    A1, A2, A3, S1, S2, S3 };


// ******************************** //
// This is the GammaStructure class //
// ******************************** //

class GammaStructure
{  
 public:
  // Constructors
  // ------------
  // Syntax of GammaStructure constructor is
  // GammaStructure(unity,gamma^5,gamma^0,gamma^1,
  //                gamma^2,gamma^3,gamma^5*gamma^0,gamma^5*gamma^1,
  //                gamma^5*gamma^2,gamma^5*gamma^3,
  //                A_1,A_2,A_3,S_1,S_2,S_3)
  GammaStructure(const std::complex<double>& =0);
  GammaStructure(const std::complex<double>&, const std::complex<double>&);
  GammaStructure(const std::complex<double>&, const std::complex<double>&,
		 const std::complex<double>&, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0); 
  GammaStructure(const std::complex<double>&, const std::complex<double>&,
		 const std::complex<double>&, const std::complex<double>&,
		 const std::complex<double>&, const std::complex<double>&,
		 const std::complex<double>&, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0); 
  GammaStructure(const GammaStructure&); // copy constructor

  // Destructor
  // ----------
  ~GammaStructure();

  // Assignment
  // ----------
  GammaStructure& operator=(const GammaStructure&);
  
  // Evaluate the structure and return the 4x4 matrix
  // ------------------------------------------------
  Matrix<4,4> value() const;
  operator Matrix<4,4>() const { return value(); }

  // Addition
  // --------
  GammaStructure operator+(const GammaStructure&) const;
  GammaStructure operator-(const GammaStructure&) const;
  
  GammaStructure& operator+=(const GammaStructure&);
  GammaStructure& operator-=(const GammaStructure&);

  /* // Addition with scalar multiplication */
  /* // ----------------------------------- */
  /* GammaStructure& multiAdd(double,const GammaStructure&, */
  /* 			   double,const GammaStructure&, */
  /* 			   double,const GammaStructure&); */
  /* GammaStructure multiAdd(double,const GammaStructure&, */
  /* 			  double,const GammaStructure&, */
  /* 			  double,const GammaStructure&) const; */
  /* static */
  /*   GammaStructure multiAdd(double,const GammaStructure&, */
  /* 			    double,const GammaStructure&, */
  /* 			    double,const GammaStructure&, */
  /* 			    double,const GammaStructure&); */

  // Multiplication
  // --------------
  GammaStructure operator*(double) const;
  GammaStructure operator*(const std::complex<double>&) const; 
  GammaStructure operator*(const GammaStructure&) const;
  template<int K>
  Matrix<4,K>    operator*(const Matrix<4,K>& m) const
    { return (Matrix<4,4>)*this * m; }

  friend
    GammaStructure operator*(double, const GammaStructure&);
  friend
    GammaStructure operator*(const std::complex<double>&, const GammaStructure&);
  template<int K> friend
    Matrix<K,4>    operator*(const Matrix<K,4>& m, const GammaStructure& g)
    { return m * (Matrix<4,4>)g; }

  GammaStructure& operator*=(double);
  GammaStructure& operator*=(const std::complex<double>&);
  GammaStructure& operator*=(const GammaStructure&);

  // Matrix operations
  // -----------------
  GammaStructure Hermitian() const;
  GammaStructure Conjugate() const;
  
  // Output
  // ------
  void print(void) const; // Prints structure in [..,..,..] format
  friend 
    std::ostream& operator<<(std::ostream&, const GammaStructure&); // Overwrite <<

 private:
  // Data members
  // ------------
  int nrComp;
  int* comp;
  std::complex<double>* compCoeff;

  // multiplication tables
  static const int multTable[16][16];
  static const std::complex<double> multTableCoeff[16][16];

  // Private member funtions
  // -----------------------
  GammaStructure(int,int*,std::complex<double>*); // member data constructor
  void crop();                               // Minimize member data
};


/* The 4vector gamma_{mu} (GMU) containing 4 GammaStructures is defined in
 * ../FourVector/FourVector.h for nesting reasons */

#endif
