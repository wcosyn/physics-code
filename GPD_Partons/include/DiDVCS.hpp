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
 * @brief enumeration used to pass which photon polarization cross section to calculate
 * 
 */
  enum Photon_pol { kgammaL=1,     ///<   longitudinal photon polarization
		       kgammaT=2 ///< transverse photon polarization
  };


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
 * @param[out] [nb/GeV^4] results various cross section calculations <<ADD indices>>
 * @param scale [GeV^2] factorization and renorm scale
 * @param xi [] skewness
 * @param Q2in [GeV^2] incoming virtual photon 4mom sq.
 * @param Q2out [GeV^2] incoming virtual photon 4mom sq.
 * @param mandelstam_t [GeV^2] momentum transfer sq to the proton
 * @param gammapol [kgammaT] gamma_transverse [kgammaL] gamma_longitudinal
 * @param max_integrationsteps number of integration steps in du,dz integral.  1E04 is good value above Q^2=.01, below take 2E05
 * @
 * 
 */
void getCross_DiDVCS(std::vector<double> & results, const double scale, const double xi, const double Q2in, const double Q2out, 
                          const double mandelstam_t, DiDVCS::Photon_pol gammapol, const int max_integrationsteps);

/**
 * @brief calculate the electroproduction Cross section for two vector polarization, general function, pass
 * 
 * @param[out] results_total [nb/GeV^4] results various electroproduction cross section calculations <<ADD indices>>
 * @param[in] resultsL [nb/GeV^4] results various gammaL cross section calculations <<ADD indices>>
 * @param[in] resultsT [nb/GeV^4] results various gammaT cross section calculations <<ADD indices>>
 * @param y [] fractional energy loss of scattered electron
 * @param Q2 [GeV^2] Qsquared
 * 
 * 
 */
void getElectro_Cross_DiDVCS(std::vector<double> & results_total, std::vector<double> & resultsL, std::vector<double> & resultsT,
                        const double y, const double Q2);


private:

//needed for integrandum of amplitude, advanced complex version
struct Ftor_DiDVCS_general {

  static void exec(const numint::array<double,1> &x, void *param, numint::vector_z &ret) {
    Ftor_DiDVCS_general &p = * (Ftor_DiDVCS_general *) param;
    p.f(ret,x[0], p.xi, p.mandelstam_t, p.scale, p.Qsqin, p.Qsqout, p.spinout, p.spinin, p.gammapol, *p.pobj);
  }
  double xi; ///< [] skewness
  double mandelstam_t; ///< [GeV^2] momentum transfer, taken at t_min for now
  double scale; ///< [GeV] factorization = renormalization scale
  double Qsqin; ///< [GeV^2] incoming virtual photon virtuality 
  double Qsqout; ///< [GeV^2] outgoing virtual photon virtuality 
  int spinin; ///< [] spin*2 of incomcing nucleon
  int spinout;///< [] spin*2 of outgoing nucleon
  DiDVCS::Photon_pol gammapol; ////< [kgammaL/kgammaT] long/transverse photon polarization
  DiDVCS* pobj;
   /**
   * @brief integration function, integration over Bjorken x, see Eq. (23) of PRD paper without the prefactor
   * 
   */
  void (*f)(numint::vector_z &, double x, double xi, double mandelstam_t, double scale, double Qsqin, double Qsqout, 
            int spinout, int spinin, DiDVCS::Photon_pol gammapol, DiDVCS& twovector);

};


static void integrandum_DiDVCS(numint::vector_z &, double x, double xi, double mandelstam_t, double scale, 
                                      double Qsqin, double Qsqout, int spinout, int spinin,
                                      DiDVCS::Photon_pol gamma,  DiDVCS& twovector);



PARTONS::GPDService* pGPDService; ///< GPD service from PARTONS
PARTONS::GPDModule* pGPDModel; ///< GPD model from PARTONS
PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule; ///< alpha_S from PARTONS

};

#endif