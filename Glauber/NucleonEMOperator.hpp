/*! \file NucleonEMOperator.hpp 
 * \brief Contains declaration of class NucleonEMoperator
 * \author Wim Cosyn
 * \date 16/10/2012
 * \addtogroup Glauber
 * @{
 */
#ifndef NUCLEONEMOPERATOR_HPP
#define NUCLEONEMOPERATOR_HPP

#include <iostream>
#include "constants.hpp"
#include "FourVector.h"
#include "GammaStructure.h"

class MeanFieldNucleusThick; //forward declaration

/*! \brief A class that contains the EM coupling to the nucleon.
 * Has form factors, medium modifications (CQM, QSM), 
 * computes dirac matrices with coupling
 *
 */
class NucleonEMOperator{
public:
  NucleonEMOperator(); /*!< Empty constructor */
  /*! constructor
   * \param Q2 [MeV^2] photon fourmomentum transfer squared
   * \param proton interaction with proton [1] or neutron [0]
   * \param para which parametrization? <BR>
   * - 0=BBA (H. Budd, A. Bodek, and J. Arrington, hep-ex/0308005, J. Arrington, Phys. Rev. C 69, 022201(R) (2004)) <BR>
   * - 1=SLAC dipole <BR>
   * - 2=JJ Kelly (Phys. Rev. C70, 068202, 2004) <BR>
   * - 3=BBBA (Bradford, Bodek, Budd, Arrington,  hep-ex/0602017) <BR>
   */
  NucleonEMOperator(const double Q2, const bool proton, const int para);
  ~NucleonEMOperator(); /*!< Destructor */
  double getGE() const{ return GE;} /*!< Gives you electric form factor */
  double getGM() const{return GM;} /*!< Gives you magnetic form factor */
  double getF1() const{return F1;} /*!< Gives you F1 */
  double getF2() const{return F2;} /*!< Gives you F2 */
  /*! compute GE with medium modification
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return GE with medium modification
   */
  double getGE(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{ return GE*getEmod(r,medium,nucleus);}
  /*! compute GM with medium modification
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return GM with medium modification
   */
  double getGM(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{ return GM*getMmod(r,medium,nucleus);}
  /*! compute F1 with medium modification
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return F1 with medium modification
   */
  double getF1(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! compute F2 with medium modification
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return F2 with medium modification
   */
  double getF2(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! computes the dirac foton-bound nucleon coupling according to the CC1 description of De Forest et al.
   * \param pi [MeV] initial nucleon momentum fourvector
   * \param pf [MeV] final nucleon momentum fourvector
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC1(const FourVector<double> &pi, const FourVector<double> &pf, 
				    const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! computes the dirac foton-bound nucleon coupling according to the CC2 description of De Forest et al.
   * \param q [MeV] photon momentum fourvector
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC2(const FourVector<double> &q, const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! computes the dirac foton-bound nucleon coupling according to the CC3 description of De Forest et al.
   * \param q [MeV] photon momentum fourvector
   * \param pi [MeV] initial nucleon momentum fourvector
   * \param pf [MeV] final nucleon momentum fourvector
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC3(const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				    const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! computes the dirac foton-bound nucleon coupling according to a CC description of De Forest et al.
   * \param current which CC of De Forest [1=CC1, 2=CC2, 3=CC3]
   * \param q [MeV] photon momentum fourvector
   * \param pi [MeV] initial nucleon momentum fourvector
   * \param pf [MeV] final nucleon momentum fourvector
   * \param r [fm] radial coordinate in nucleus
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC(const int current, const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf, 
				   const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! computes the dirac foton-bound nucleon coupling according to the CC1 description of De Forest et al.
   * \param pi [MeV] initial nucleon momentum fourvector
   * \param pf [MeV] final nucleon momentum fourvector
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC1(const FourVector<double> &pi, const FourVector<double> &pf) const;
  /*! computes the dirac foton-bound nucleon coupling according to the CC2 description of De Forest et al.
   * \param q [MeV] photon momentum fourvector
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC2(const FourVector<double> &q) const;
  /*! computes the dirac foton-bound nucleon coupling according to the CC3 description of De Forest et al.
   * \param q [MeV] photon momentum fourvector
   * \param pi [MeV] initial nucleon momentum fourvector
   * \param pf [MeV] final nucleon momentum fourvector
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC3(const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf) const;
  /*! computes the dirac foton-bound nucleon coupling according to a CC description of De Forest et al.
   * \param current which CC of De Forest [1=CC1, 2=CC2, 3=CC3]
   * \param q [MeV] photon momentum fourvector
   * \param pi [MeV] initial nucleon momentum fourvector
   * \param pf [MeV] final nucleon momentum fourvector
   * \return fourvector of dirac matrices with CC1 coupling to nucleon
   */
  FourVector<GammaStructure> getCC(const int current, const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf) const;
  static const FourVector<GammaStructure> gamma_mu; /*!< gamma matrices */
  static const GammaStructure Id; /*!< Unit matrix */
  /*! grid for QMC GE modification, stepsize in Q2 [61] is 0.05 GeV^2, 
   * stepsize in density[15] (norm to A) is 0.016 fm^-3
   */
  static const double QMCGE[61][15]; 
  /*! grid for QMC GM modification, stepsize in Q2 [61] is 0.05 GeV^2, 
   * stepsize in density[15] (norm to A) is 0.016 fm^-3
   */
  static const double QMCGM[61][15];
  /*! grid for QSM GE modification, stepsize in Q2 [61] is 0.05 GeV^2, 
   * stepsize in density[16] (norm to A) is 0.016 fm^-3
   */
  static const double CQSMGE[61][16];
  /*! grid for QSM GM modification, stepsize in Q2 [61] is 0.05 GeV^2, 
   * stepsize in density[16] (norm to A) is 0.016 fm^-3
   */
  static const double CQSMGM[61][16];

  /*! gives you the modification of the electric form factor
   * \param r [fm] radial coordinate
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return modification of GE_null
   */
  double getEmod(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  /*! gives you the modification of the magnetic form factor
   * \param r [fm] radial coordinate
   * \param medium which modification [1=CQM, 2=QSM]
   * \param nucleus nucleus class object, contains the densities etc.
   * \return modification of GM_null
   */
  double getMmod(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;


protected:
  bool proton; /*!< 1=proton, 0=neutron coupling */
  int parametrization; /*!< 0=BBA, 1=dipole */
  double Q2; /*!< [MeV^2] fourmomentum sq */
  int Q2index; /*!< used in medium modification, index in the grid in Q2 */
  double Q2_interp; /*!< interpolation variable for Q2 */
  double tau; /*!< \f$ \tau = \frac{Q^2}{4m_N^2} \f$ */
  double GE_null; /*!< GE value at Q^2=0 */
  double GM_null; /*!< GM value at Q^2=0 */
  double GE; /*!< electric form factor */
  double GM; /*!< magnetic form factor */
  double F1; /*!< F1 */
  double F2;  /*!< kappa*F2!!!!  with kappa anomalous magneticmoment */
  
  
private:
  
  void setGE();
  void setGM();
  /*! gives you a dipole form factor 
   * \param Q2 [MeV^2] fourmomentum sq
   * \return dipole form factor
   */
  double Get_Gdipole(const double Q2); 
  void setF1();
  void setF2();
  
  
  
};

/** @} */
#endif
