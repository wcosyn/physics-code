#ifndef TWOVECTOR_NUCL_HPP
#define TWOVECTOR_NUCL_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <TDeuteron.h>
#include <complex>
#include <numint/numint.hpp>
#include <TInterpolatingWavefunction.h>
#include "GPD_T_Nucl_grid.hpp"


#include <partons/Partons.h>
#include <partons/modules/running_alpha_strong/RunningAlphaStrongStandard.h>
#include <partons/services/GPDService.h>

/**
 * @brief Class implements two vector meson production cross section for the nucleon case, various cross sections: both incoming L and T (virtual) photon; both L and T rho meson
 * 
 */
class TwoVector_Nucl{

public:
/**
 * @brief Construct a new TwoVector_Nucl object,  used to calculate cross sections
 * 
 * @param pGPDService GPD service from PARTONS
 * @param pGPDModel GPD model from PARTONS
 * @param pRunningAlphaStrongModule alpha_S from PARTONS
 * @param scale [GeV] factorization and renorm scale
 * @param xi [] skewness
 * @param Q2 [GeV^2] incoming virtual photon 4mom sq.
 */
TwoVector_Nucl(PARTONS::GPDService* pGPDService, PARTONS::GPDModule* pGPDModel, PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule);

/**
 * @brief Destructor
 * 
 */
~TwoVector_Nucl(){;}

/**
 * @brief calculate the Cross section for gammaL rhoL
 * 
 * @param scale [GeV^2] factorization and renorm scale
 * @param xi [] skewness
 * @param Q2 [GeV^2] incoming virtual photon 4mom sq.
 * @param psq [GeV^2] pomeron hard scale
 * @return double [nb/GeV^4] cross section 
 */
double getCross_gammaL_rhoL(const double scale, const double xi, const double Q2, const double psq);

/**
 * @brief calculate the Cross section for gammaT rhoL
 * 
 * @param scale [GeV^2] factorization and renorm scale
 * @param xi [] skewness
 * @param Q2 [GeV^2] incoming virtual photon 4mom sq.
 * @param psq [GeV^2] pomeron hard scale
 * @return double [nb/GeV^4] cross section 
 */
double getCross_gammaT_rhoL(const double scale, const double xi, const double Q2, const double psq);

/**
 * @brief calculate the Cross section for gammaL rhoT
 * 
 * @param scale [GeV^2] factorization and renorm scale
 * @param xi [] skewness
 * @param Q2 [GeV^2] incoming virtual photon 4mom sq.
 * @param psq [GeV^2] pomeron hard scale
 * @return double [nb/GeV^4] cross section 
 */
double getCross_gammaL_rhoT(const double scale, const double xi, const double Q2, const double psq);

/**
 * @brief calculate the Cross section for gammaT rhoT
 * 
 * @param scale [GeV^2] factorization and renorm scale
 * @param xi [] skewness
 * @param Q2 [GeV^2] incoming virtual photon 4mom sq.
 * @param psq [GeV^2] pomeron hard scale
 * @return double [nb/GeV^4] cross section 
 */
double getCross_gammaT_rhoT(const double scale, const double xi, const double Q2, const double psq);

private:

//needed for integrandum of amplitude
struct Ftor_2vector {

  static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
    Ftor_2vector &p = * (Ftor_2vector *) param;
    p.f(ret,x[0],x[1], p.xi, p.t, p.scale, p.psq, p.Qsq, p.spinout, p.spinin, *p.pobj);
  }
  double xi; ///< [] skewness
  double t; ///< [GeV^2] momentum transfer, taken at t_min for now
  double scale; ///< [GeV] factorization = renormalization scale
  double psq; ///< [GeV^2] pomeron hard scale
  double Qsq; ///< [GeV^2] virtual photon 4mom squared
  int spinin; /// [] spin*2 of incomcing nucleon
  int spinout;/// [] spin*2 of outgoing nucleon
  TwoVector_Nucl* pobj;
   /**
   * @brief integration function, integration over u and z from DA's, see Enberg et al paper EPJC47 87-94
   * 
   */
  void (*f)(numint::vector_d &, double u, double z, double xi, double t, double scale, double psq, double Qsq, int spinout, int spinin, TwoVector_Nucl& twovector);

};

/**
 * @brief  integrandum for longitudinal polarized photon on nucleon, two times rho_L in final state
 * 
 * @param u [] integration variable, momentum fraction in lower rho DA
 * @param z [] integration variable, momentum fraction in upper rho DA
 * @param xi [] skewness
 * @param t  [GeV^2] momentum transfer, taken at t_min for now
 * @param psq [GeV^2] pomeron hard scale
 * @param Qsq [GeV^2] virtual photon 4mom squared
 * @param pGPDService GPD service from PARTONS
 * @param pGPDModel GPD model from PARTONS
 */
static void integrandum_L(numint::vector_d &, double u, double z, double xi, double t, double scale, double psq, double Qsq, int spintout, int spinin, TwoVector_Nucl& twovector);

/**
 * @brief integrandum for transversily polarized photon on nucleon, two times rho_L in final state
 * 
 * @param u [] integration variable, momentum fraction in lower rho DA
 * @param z [] integration variable, momentum fraction in upper rho DA
 * @param xi [] skewness
 * @param t  [GeV^2] momentum transfer, taken at t_min for now
 * @param psq [GeV^2] pomeron hard scale
 * @param Qsq [GeV^2] virtual photon 4mom squared
  * @param pGPDService GPD service from PARTONS
 * @param pGPDModel GPD model from PARTONS
*/
static void integrandum_T(numint::vector_d &, double u, double z, double xi, double t, double scale, double psq, double Qsq, int spintout, int spinin, TwoVector_Nucl& twovector);

/**
 * @brief  integrandum for longitudinal polarized photon on nucleon, rho_L and rho_T in final state
 * 
 * @param u [] integration variable, momentum fraction in lower rho DA
 * @param z [] integration variable, momentum fraction in upper rho DA
 * @param xi [] skewness
 * @param t  [GeV^2] momentum transfer, taken at t_min for now
 * @param psq [GeV^2] pomeron hard scale
 * @param Qsq [GeV^2] virtual photon 4mom squared
 * @param pGPDService GPD service from PARTONS
 * @param pGPDModel GPD model from PARTONS
 */
static void integrandum_T_L(numint::vector_d &, double u, double z, double xi, double t, double scale, double psq, double Qsq, int spintout, int spinin, TwoVector_Nucl& twovector);

/**
 * @brief integrandum for transversily polarized photon on nucleon, rho_L and rho_T in final state
 * 
 * @param u [] integration variable, momentum fraction in lower rho DA
 * @param z [] integration variable, momentum fraction in upper rho DA
 * @param xi [] skewness
 * @param t  [GeV^2] momentum transfer, taken at t_min for now
 * @param psq [GeV^2] pomeron hard scale
 * @param Qsq [GeV^2] virtual photon 4mom squared
  * @param pGPDService GPD service from PARTONS
 * @param pGPDModel GPD model from PARTONS
*/
static void integrandum_T_T(numint::vector_d &, double u, double z, double xi, double t, double scale, double psq, double Qsq, int spintout, int spinin, TwoVector_Nucl& twovector);

PARTONS::GPDService* pGPDService; ///< GPD service from PARTONS
PARTONS::GPDModule* pGPDModel; ///< GPD model from PARTONS
PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule; ///< alpha_S from PARTONS

GPD_T_Nucl_grid gpdTgrid;
};

#endif