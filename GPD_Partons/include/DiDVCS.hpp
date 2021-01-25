#ifndef DiDVCS_HPP
#define DiDVCS_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <complex>
#include <numint/numint.hpp>


#include <partons/Partons.h>
#include <partons/modules/running_alpha_strong/RunningAlphaStrongStandard.h>
#include <partons/services/GPDService.h>

/**
 * @brief Class implements diffractive DVCS amplitude and cross section, see Pire, Szymanowski, Wallon Phys.Rev.D 101 (2020) 7, 074005
 * 
 */
class DiDVCS{

public:



/**
 * @brief Construct a new DiDVCS object,  used to calculate cross sections
 * 
 * @param pGPDService GPD service from PARTONS
 * @param pGPDModel GPD model from PARTONS
 * @param pRunningAlphaStrongModule alpha_S from PARTONS
 */
DiDVCS(PARTONS::GPDService* pGPDService, PARTONS::GPDModule* pGPDModel, PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule);

/**
 * @brief Destructor
 * 
 */
~DiDVCS(){;}


/**
 * @brief calculate the photoproduction Cross section for two vector polarization, general function, pass
 * 
 * @param[out] [nb/GeV^8 nb/GeV^10] results various cross section calculations [0] photoproduction [1] electroproduction dydQ^2dt_n dt_rho ds2 dQ'2
 * @param scale [GeV^2] factorization and renorm scale
 * @param t_rho [GeV^2] (q-q_rho)^2 momentum transfer between virtual photon and diffractive rho
 * @param t_N [GeV^2] (p1-p'1)^2 momentum transfer initial to final nucleon
 * @param Q2in [GeV^2] incoming virtual photon 4mom sq.
 * @param Q2out [GeV^2] incoming virtual photon 4mom sq.
 * @param s2 [GeV^2] invariant mass nucleon outgoing photon
 * @param s [GeV^2] invariant mass incoming electron nucleon system
 * @param y [] inelasticity
 * @param max_integrationsteps number of integration steps in du,dz integral.  1E04 is good value above Q^2=.01, below take 2E05
 * @param no_realpart consider only the (dominating) imaginary part of the CFFs
 * @
 * 
 */
void getCross_DiDVCS(std::vector<double> & results, const double scale, const double t_rho, const double t_N, const double Q2in, const double Q2out, 
                          const double s2, const double s, double y, const int maxintsteps, bool no_realpart);


private:

//needed for integrandum of amplitude, advanced complex version
struct Ftor_DiDVCS_general {

  static void exec(const numint::array<double,1> &x, void *param, numint::vector_z &ret) {
    Ftor_DiDVCS_general &p = * (Ftor_DiDVCS_general *) param;
    p.f(ret,x[0], p.xi, p.mandelstam_t, p.scale, p.Q2out_over_Q2in, p.spinout, p.spinin, *p.pobj);
  }
  double xi; ///< [] skewness
  double mandelstam_t; ///< [GeV^2] momentum transfer, taken at t_min for now
  double scale; ///< [GeV] factorization = renormalization scale
  double Q2out_over_Q2in; ///< [] outgoing over incoming virtual photon virtuality 
  int spinin; ///< [] spin*2 of incomcing nucleon
  int spinout;///< [] spin*2 of outgoing nucleon
  DiDVCS* pobj;
   /**
   * @brief integration function, integration over Bjorken x, see Eq. (23) of PRD paper without the prefactor
   * 
   */
  void (*f)(numint::vector_z &, double x, double xi, double mandelstam_t, double scale, double Q2out_over_Q2in, 
            int spinout, int spinin,DiDVCS& twovector);

};


static void integrandum_DiDVCS(numint::vector_z &, double x, double xi, double mandelstam_t, double scale, 
                                      double Q2out_over_Q2in, int spinout, int spinin,
                                     DiDVCS& twovector);



PARTONS::GPDService* pGPDService; ///< GPD service from PARTONS
PARTONS::GPDModule* pGPDModel; ///< GPD model from PARTONS
PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule; ///< alpha_S from PARTONS

};

#endif