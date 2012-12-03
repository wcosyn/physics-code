/*! \file TMFSpinor.hpp 
 * \brief Contains declaration of class TMFSpinor
 * \author Wim Cosyn
 * \date 24/08/2012
 * \addtogroup Glauber
 * @{
 */
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

/*! \brief A class for dirac spinor representing a mean-field wave function of a certain nucleon in a nucleus.
 * 
 * Careful!!!  We include a factor r in the definition here, so dimension is [fm^(-1/2)] <BR>
 * 
 * Four components have the following structure:<BR>
 * 0-1: \f$iG_{n\kappa}(r) Y_{\kappa m}(\Omega, \vec{\sigma})\f$ <BR>
 * 2-3: \f$-F_{n\kappa}(r) Y_{-\kappa m}(\Omega, \vec{\sigma})\f$
 * 
 */


class TMFSpinor : public TObject
{
 public:
  /*! Constructor
   * \param nucleus nucleus we want a wf spinor from
   * \param shellindex shell index of the nucleon 
   * \param m m_j quantum number (times two!!)
   * \param r [fm] radial coordinate
   * \param costheta cos(theta)
   * \param phi azimutal coordinate
   */   
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double r, double costheta, double phi);
  /*! Constructor
   * \param nucleus nucleus we want a wf spinor from
   * \param shellindex shell index of the nucleon 
   * \param m m_j quantum number (times two!!)
   * \param Gr [fm-{1/2}] G(r) already computed
   * \param Fr [fm-{1/2}] F(r) already computed
   * \param costheta cos(theta)
   * \param phi azimutal coordinate
   */   
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, double costheta, double phi);
  /*! Constructor
   * \param nucleus nucleus we want a wf spinor from
   * \param shellindex shell index of the nucleon 
   * \param m m_j quantum number (times two!!)
   * \param Gr [fm-{1/2}] G(r) already computed
   * \param Fr [fm-{1/2}] F(r) already computed
   * \param Spher_harm array with theta part of all 4 components
   * \param phi azimutal coordinate
   */   
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, const double *Spher_harm, double phi);
  /*! Constructor
   * \param nucleus nucleus we want a wf spinor from
   * \param shellindex shell index of the nucleon 
   * \param m m_j quantum number (times two!!)
   * \param r [fm] radial coordinate
   * \param costheta cos(theta)
   * \param phi azimutal coordinate
   * \param CG CG coefficients already computed
   */   
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double r, double costheta, double phi, Matrix<4,1> CG);
  /*! Constructor
   * \param nucleus nucleus we want a wf spinor from
   * \param shellindex shell index of the nucleon 
   * \param m m_j quantum number (times two!!)
   * \param Gr [fm-{1/2}] G(r) already computed
   * \param Fr [fm-{1/2}] F(r) already computed
   * \param costheta cos(theta)
   * \param phi azimutal coordinate
   * \param CG CG coefficients already computed
   */     
  TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, 
	    double costheta, double phi, Matrix<4,1> CG);
  /*! Constructor
   * \param m m_j quantum number (times two!!)
   * \param Gr [fm-{1/2}] G(r) already computed
   * \param Fr [fm-{1/2}] F(r) already computed
   * \param Spher_harm array with theta part of all 4 components
   * \param phi azimutal coordinate
   * \param CG CG coefficients already computed
   */   
  TMFSpinor(double Gr, double Fr, int m, const double *Spher_harm, double phi, Matrix<4,1> CG);
  TMFSpinor(TRootIOCtor*); /*!< Root Constructor */
  TMFSpinor(const TMFSpinor&);  /*!< Copy Constructor */
  virtual ~TMFSpinor();  /*!< Destructor */
  TMFSpinor& operator=(const TMFSpinor&);  /*!< Copy operator */
  virtual TMFSpinor* Clone(const char* ="") const;  /*!< Clone function */

  operator const Matrix<4,1>&() const { return *fComponent; }  /*!< cast to matrix object */
  template<int K>
    Matrix<4,K> operator*(const Matrix<1,K>& rhs) { return (Matrix<4,1>)(*this) * rhs; }  /*!< matrix multiplication */
  
  friend 
    std::ostream& operator<<(std::ostream&, const TMFSpinor&);  /*!< Stream operator */
  
  Matrix<1,4> H();  /*!< gives you the hermitian */
   
 protected:
  Matrix<4,1> *fComponent; /*!< a spinor is a 4-component vector */
 
  ClassDef(TMFSpinor,1); // Dirac spinor for on-mass shell spin-1/2 particles
};

/*! multiplication of a 4*4 matrix with a MF spinor, gives you 4*1 matrix(spinor) */
Matrix<4,1> operator*(const Matrix<4,4>& lhs, const TMFSpinor& rhs);
/*! multiplication of a 1*4 matrix with a MF spinor, gives you scalar */
std::complex<double> operator*(const Matrix<1,4>& lhs, const TMFSpinor& rhs);
/*! multiplication of a 4*4 matrix (of the gamma structure kind) with a MF spinor */
Matrix<4,1> operator*(const GammaStructure& lhs, const TMFSpinor& rhs);

/** @}*/
#endif



