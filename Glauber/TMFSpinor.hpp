#ifndef TMFSPINOR_H
#define TMFSPINOR_H

#include <Matrix.h>
#include <FourVector.h>
#include <TObject.h>
#include <complex>
#include <GammaStructure.h>
#include <TVector3.h>
#include <iostream>
#include "MeanFieldNucleus.hpp"

class TMFSpinor : public TObject
{
 public:
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double r, double costheta, double phi);
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, double costheta, double phi);
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, const double *Sper_harm, double phi);
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double r, double costheta, double phi, Matrix<4,1> CG);
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, double costheta, double phi, Matrix<4,1> CG);
  TMFSpinor(double Gr, double Fr, int m, const double *Sper_harm, double phi, Matrix<4,1> CG);
  TMFSpinor(TRootIOCtor*); // ROOT I/0 Constructor
  TMFSpinor(const TMFSpinor&);
  virtual ~TMFSpinor();
  TMFSpinor& operator=(const TMFSpinor&);
  virtual TMFSpinor* Clone(const char* ="") const;

  operator const Matrix<4,1>&() const { return *fComponent; }
  template<int K>
    Matrix<4,K> operator*(const Matrix<1,K>& rhs) { return (Matrix<4,1>)(*this) * rhs; }
  
  friend 
    std::ostream& operator<<(std::ostream&, const TMFSpinor&);
  
  Matrix<1,4> H();
   
 protected:
  Matrix<4,1> *fComponent; //! a spinor is a 4-component vector
 
  ClassDef(TMFSpinor,1); // Dirac spinor for on-mass shell spin-1/2 particles
};

Matrix<4,1> operator*(const Matrix<4,4>& lhs, const TMFSpinor& rhs);
std::complex<double> operator*(const Matrix<1,4>& lhs, const TMFSpinor& rhs);
Matrix<4,1> operator*(const GammaStructure& lhs, const TMFSpinor& rhs);


#endif



