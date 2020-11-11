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
 * @brief enumeration to pass rhoT or rhoL polarization (which determins chiral even or odd GPDs)
 * 
 */
  enum Rho_pol { krhoL=1,     ///<   longitudinal rho polarization
		       krhoT=2, ///< transverse rho polarization ///< transverse rho polarization
           kaxial=3 ///< pseudescalar meson
  };


/**
 * @brief enumeration used to pass which photon polarization cross section to calculate
 * 
 */
  enum Photon_pol { kgammaL=1,     ///<   longitudinal photon polarization
		       kgammaT=2 ///< transverse photon polarization
  };


/**
 * @brief Construct a new TwoVector_Nucl object,  used to calculate cross sections
 * 
 * @param pGPDService GPD service from PARTONS
 * @param pGPDModel GPD model from PARTONS
 * @param pRunningAlphaStrongModule alpha_S from PARTONS
 */
TwoVector_Nucl(PARTONS::GPDService* pGPDService, PARTONS::GPDModule* pGPDModel, PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule);

/**
 * @brief Destructor
 * 
 */
~TwoVector_Nucl(){;}


/**
 * @brief calculate the photoproduction Cross section for two vector polarization, general function, pass
 * 
 * @param[out] [nb/GeV^4] results various cross section calculations <<ADD indices>>
 * @param scale [GeV^2] factorization and renorm scale
 * @param xi [] skewness
 * @param Q2 [GeV^2] incoming virtual photon 4mom sq.
 * @param psq [GeV^2] virtuality of hard Pomeron
 * @param gammapol [kgammaT] gamma_transverse [kgammaL] gamma_longitudinal
 * @param rhopol [krhoT] Transverse [krhoL] Longitudinal
 * @param max_integrationsteps number of integration steps in du,dz integral.  1E04 is good value above Q^2=.01, below take 2E05
 * @
 * 
 */
void getCross_twovector(std::vector<double> & results, const double scale, const double xi, const double Q2, const double psq, 
                          TwoVector_Nucl::Photon_pol gammapol, TwoVector_Nucl::Rho_pol rhopol, const int max_integrationsteps);


/**
 * @brief calculate the photoproduction Cross section for two vector polarization, general function, dtrho ds2 dtN
 * 
 * @param[out] [nb/GeV^6] results various cross section calculations <<ADD indices>>
 * @param scale [GeV^2] factorization and renorm scale
 * @param s_eN [GeV^2] invariant mass eN system
 * @param y [] q.pn/pe.pn
 * @param Q2 [GeV^2] incoming virtual photon 4mom sq.
 * @param trho [GeV^2] momentum transfer to rho ~ virtuality of hard Pomeron
 * @param tN [GeV^2] momentum transfer to nucleon
 * @param s2 [GeV^2] invariant mass of meson proton system (GPD part)
 * @param gammapol [kgammaT] gamma_transverse [kgammaL] gamma_longitudinal
 * @param rhopol [krhoT] Transverse [krhoL] Longitudinal
 * @param max_integrationsteps number of integration steps in du,dz integral.  1E04 is good value above Q^2=.01, below take 2E05
 * @
 * 
 */
void getCross_twovector_ds2(std::vector<double> & results, const double scale, const double s_eN, const double y,
                                 const double Q2, const double trho, const double tN, const double s2, 
                          TwoVector_Nucl::Photon_pol gammapol, TwoVector_Nucl::Rho_pol rhopol, const int maxintsteps);


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
void getElectro_Cross_twovector(std::vector<double> & results_total, std::vector<double> & resultsL, std::vector<double> & resultsT,
                        const double y, const double Q2);





private:

//needed for integrandum of amplitude, advanced complex version
struct Ftor_2vector_general {

  static void exec(const numint::array<double,2> &x, void *param, numint::vector_z &ret) {
    Ftor_2vector_general &p = * (Ftor_2vector_general *) param;
    p.f(ret,x[0],x[1], p.xi, p.mandelstam_t, p.scale, p.psq, p.Qsq, p.spinout, p.spinin, p.rhopol,  p.gammapol, *p.pobj);
  }
  double xi; ///< [] skewness
  double mandelstam_t; ///< [GeV^2] momentum transfer, taken at t_min for now
  double scale; ///< [GeV] factorization = renormalization scale
  double psq; ///< [GeV^2] pomeron hard scale
  double Qsq; ///< [GeV^2] virtual photon 4mom squared
  int spinin; ///< [] spin*2 of incomcing nucleon
  int spinout;///< [] spin*2 of outgoing nucleon
  TwoVector_Nucl::Rho_pol rhopol; ////< [krhoL/krhoT] long/transverse rho polarization
  TwoVector_Nucl::Photon_pol gammapol; ////< [kgammaL/kgammaT] long/transverse photon polarization
  TwoVector_Nucl* pobj;
   /**
   * @brief integration function, integration over u and z from DA's, see Enberg et al paper EPJC47 87-94
   * 
   */
  void (*f)(numint::vector_z &, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq, 
            int spinout, int spinin, TwoVector_Nucl::Rho_pol rhopol, TwoVector_Nucl::Photon_pol gammapol, TwoVector_Nucl& twovector);

};


static void integrandum_rho_general(numint::vector_z &, double u, double z, double xi, double mandelstam_t, double scale, 
                                      double psq, double Qsq, int spinout, int spinin, TwoVector_Nucl::Rho_pol rhopol,
                                      TwoVector_Nucl::Photon_pol gamma,  TwoVector_Nucl& twovector);



PARTONS::GPDService* pGPDService; ///< GPD service from PARTONS
PARTONS::GPDModule* pGPDModel; ///< GPD model from PARTONS
PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule; ///< alpha_S from PARTONS

GPD_T_Nucl_grid gpdTgrid;
};

#endif