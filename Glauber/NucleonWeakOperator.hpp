/*! \file NucleonEMOperator.hpp 
 * \brief Contains declaration of class NucleonEMoperator
 * \author Wim Cosyn
 * \date 16/10/2012
 * \addtogroup Glauber
 * @{
 */
#ifndef NUCLEONWEAKOPERATOR_HPP
#define NUCLEONWEAKOPERATOR_HPP

#include <iostream>
#include "constants.hpp"
#include "FourVector.h"
#include "GammaStructure.h"
#include "NucleonEMOperator.hpp"


/*! \brief A class that contains the EM coupling to the nucleon.
 * Has Weak force form factors (charged and neutral current), medium modifications (CQM, QSM), 
 * computes dirac matrices with coupling
 *
 */
class NucleonWeakOperator : public NucleonEMOperator{
public:
  NucleonWeakOperator(); /*!< Empty constructor */
  /*! constructor
   * \param [in] Q2 [MeV^2] photon fourmomentum transfer squared
   * \param proton interaction with proton [1] or neutron [0], for charged currents it refers to INITIAL nucleon
   * \param para which parametrization? <BR>
   * - 0=BBA (H. Budd, A. Bodek, and J. Arrington, hep-ex/0308005, J. Arrington, Phys. Rev. C 69, 022201(R) (2004)) <BR>
   * - 1=SLAC dipole <BR>
   * - 2=JJ Kelly (Phys. Rev. C70, 068202, 2004) <BR>
   * -3=BBBA (Bradford, Bodek, Budd, Arrington,  hep-ex/0602017) <BR>
   * \param charged 0: Z-boson, 1: W-boson
   * \param M_A [MeV] axial mass in dipole paramterization
   * \param r_s2 [fm^2] parameter for strange contribution to F1_weak
   * \param mu_s [mu_N] parameter for strange contribution to F2_weak
   * \param gA_s [] parameter for strange contribution to G_A
   */
  NucleonWeakOperator(const double Q2, const bool proton, const int para, const bool charged, const int M_A, const double r_s2,
		    const double mu_s, const double gA_s=-0.19);
  ~NucleonWeakOperator(); /*!< Destructor */
  double getGE_weak() const{ return GE_weak;} /*!< Gives you electric form factor */
  double getGM_weak() const{return GM_weak;} /*!< Gives you magnetic form factor */
  double getGA_weak() const{return GA_weak;} /*!< Gives you axial vector form factor */
  double getGP_weak() const{return GP_weak;} /*!< Gives you axial pseudoscalar form factor */
  double getF1_weak() const{return F1_weak;} /*!< Gives you F1 */
  double getF2_weak() const{return F2_weak;} /*!< Gives you F2 */
  /*! compute GE_weak with medium modification
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param [in] nucleus nucleus class object, contains the densities etc.
   * \return GE_weak with medium modification
   */
  double getGE_weak(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! compute GM_weak with medium modification
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param [in] nucleus nucleus class object, contains the densities etc.
   * \return GM_weak with medium modification
   */
  double getGM_weak(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! compute F1_weak with medium modification
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param [in] nucleus nucleus class object, contains the densities etc.
   * \return F1_weak with medium modification
   */
  double getF1_weak(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! compute F2_weak with medium modification
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param [in] nucleus nucleus class object, contains the densities etc.
   * \return F2_weak with medium modification
   */
  double getF2_weak(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! computes the dirac foton-bound nucleon coupling according to the CC1 description of De Forest et al.
   * \param [in] pi [MeV] initial nucleon momentum fourvector
   * \param [in] pf [MeV] final nucleon momentum fourvector
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param [in] nucleus nucleus class object, contains the densities etc.
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC1_weak(const FourVector<double> &pi, const FourVector<double> &pf, 
				    const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! computes the dirac foton-bound nucleon coupling according to the CC2 description of De Forest et al.
   * \param [in] q [MeV] photon momentum fourvector
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param [in] nucleus nucleus class object, contains the densities etc.
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC2_weak(const FourVector<double> &q, const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! computes the dirac foton-bound nucleon coupling according to the CC3 description of De Forest et al.
   * \param [in] q [MeV] photon momentum fourvector
   * \param [in] pi [MeV] initial nucleon momentum fourvector
   * \param [in] pf [MeV] final nucleon momentum fourvector
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param [in] nucleus nucleus class object, contains the densities etc.
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC3_weak(const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				    const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  
  /*! computes the dirac foton-bound nucleon coupling according to a CC description of De Forest et al.
   * \param current which CC of De Forest [1=CC1, 2=CC2, 3=CC3]
   * \param [in] q [MeV] photon momentum fourvector
   * \param [in] pi [MeV] initial nucleon momentum fourvector
   * \param [in] pf [MeV] final nucleon momentum fourvector
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param [in] nucleus nucleus class object, contains the densities etc.
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC_weak(const int current, const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf, 
				   const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;

  /*! computes the dirac weak axial bound nucleon coupling 
   * \param [in] q [MeV] photon momentum fourvector
   * \return fourvector of dirac matrices with axial coupling to nucleon
   */
  FourVector<GammaStructure> getAxial(const FourVector<double> &q) const;
  
  
  static const FourVector<GammaStructure> gamma_5mu; /*!< gamma_5*gamma_mu matrices */
  static const GammaStructure gamma_5; /*!< gamma_5 matrix */
  
private:
  bool charged; /*!< 0=Z-boson, 1=W-boson */

  
  double GE_weak; /*!< weak electric form factor */
  double GM_weak; /*!< weak magnetic form factor */
  double GA_weak; /*!< weak axial form factor */
  double GP_weak; /*!< weak axial form factor, pseudoscalar one */
  double F1_weak; /*!< weak F1 */
  double F2_weak;  /*!< weak kappa*F2!!!!  with kappa anomalous magneticmoment */
 
  double M_A; /*!< [MeV] axial mass */
  double r_s2; /*!< strangeness contribution parameter to F1_weak */
  double mu_s;/*!< strangeness contribution parameter to F2_weak */
  double gA_s; /*!< strangeness contribution parameter to GA_weak */
 
  NucleonEMOperator isopartner; /*!< we need the isospin partner in EM FF to compute the weak ones! */
 
  void setF1_weak();
  void setF2_weak();
  void setGE_weak();
  void setGM_weak();  
  void setGA_weak();
  /*! gives you a dipole form factor 
   * \param [in] Q2 [MeV^2] fourmomentum sq
   * \param mass [MeV] axial nucleon mass
   * \return axial dipole form factor
   */
  double Get_dipole_mass(const double Q2, const double mass) const; 
  
  
    
};

/** @} */
#endif
