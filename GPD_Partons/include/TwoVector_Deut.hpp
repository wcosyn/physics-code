#ifndef TWOVECTOR_DEUT_HPP
#define TWOVECTOR_DEUT_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <TDeuteron.h>
#include <complex>
#include <numint/numint.hpp>
#include <TInterpolatingWavefunction.h>
#include "Deut_Conv_GPD_V.hpp"
#include "Deut_Conv_GPD_T.hpp"

#include <partons/Partons.h>
#include <partons/modules/running_alpha_strong/RunningAlphaStrongStandard.h>
#include <partons/services/GPDService.h>

/**
 * @brief Class implements two vector meson production cross section for the deuteron case, various cross sections: both incoming L and T (virtual) photon; both L and T omega meson
 * 
 */
class TwoVector_Deut{

public:


/**
 * @brief enumeration to pass rhoT or rhoL polarization (which determins chiral even or odd GPDs)
 * 
 */
  enum Omega_pol { komegaL=1,     ///<   longitudinal rho polarization
		       komegaT=2, ///< transverse rho polarization ///< transverse rho polarization
  };


/**
 * @brief enumeration used to pass which photon polarization cross section to calculate
 * 
 */
  enum Photon_pol { kgammaL=1,     ///<   longitudinal photon polarization
		       kgammaT=2 ///< transverse photon polarization
  };



/**
 * @brief Construct a new TwoVector_Deut object used to calculate two vector meson prod reactions with deuterons
 * 
 * @param pGPDService GPD service from PARTONS
 * @param pGPDModel GPD model from PARTONS
 * @param pRunningAlphaStrongModule alpha_S from PARTONS
 * @param wfname  deuteron wf used in convolution
 * @param pdfname pdf used in forward limit of transversity nucleon gpd (use "MSTW")
 * @param id classid in Partons of nucleon gpd model
 */
TwoVector_Deut(PARTONS::GPDService* pGPDService, PARTONS::GPDModule* pGPDModel, PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule,
                std::string wfname, std::string pdfname, unsigned int id);

/**
 * @brief Destructor
 * 
 */
~TwoVector_Deut(){;}
/**
 * @brief calculate the Cross section for two vector polarization, general function, pass
 * 
 * @param[out] [nb/GeV^4] results various cross section calculations <<ADD indices>>
 * @param scale [GeV^2] factorization and renorm scale
 * @param xi [] skewness
 * @param Q2 [GeV^2] incoming virtual photon 4mom sq.
 * @param gammapol [kgammaT] gamma_transverse [kgammaL] gamma_longitudinal
 * @param rhopol [krhoT] Transverse [krhoL] Longitudinal
 * @param T_model implementation of gpd_T model
 * @
 * 
 */
void getCross_twovector(std::vector<double> & results, const double scale, const double xi, const double Q2, const double psq, 
                          TwoVector_Deut::Photon_pol gammapol, TwoVector_Deut::Omega_pol omegapol, int T_model);


private:


struct Ftor_2vector_general {

  static void exec(const numint::array<double,2> &x, void *param, numint::vector_z &ret) {
    Ftor_2vector_general &p = * (Ftor_2vector_general *) param;
    p.f(ret,x[0],x[1], p.xi, p.mandelstam_t, p.scale, p.psq, p.Qsq, p.helampindex, p.omegapol,  p.gammapol, p.T_model, *p.pobj);
  }
  double xi; ///< [] skewness
  double mandelstam_t; ///< [GeV^2] momentum transfer, taken at t_min for now
  double scale; ///< [GeV] factorization = renormalization scale
  double psq; ///< [GeV^2] pomeron hard scale
  double Qsq; ///< [GeV^2] virtual photon 4mom squared
  int helampindex; /// [] index of the helicity amplitude
  TwoVector_Deut::Omega_pol omegapol; ///< [krhoL/krhoT] long/transverse rho polarization
  TwoVector_Deut::Photon_pol gammapol; ///< [kgammaL/kgammaT] long/transverse photon polarization
  int T_model; ///< diff implementations of Htilde_T
  TwoVector_Deut* pobj;
   /**
   * @brief integration function, integration over u and z from DA's, see Enberg et al paper EPJC47 87-94
   * 
   */
  void (*f)(numint::vector_z &, const double u, const double z, double xi, double mandelstam_t, double scale, double psq, double Qsq, 
            int helampindex, TwoVector_Deut::Omega_pol rhopol, TwoVector_Deut::Photon_pol gammapol, int T_model, TwoVector_Deut& twovector);

};

/**
 * @brief  integrandum in u and z for amplitude, all cases
 * 
 * @param u [] integration variable, momentum fraction in lower rho DA
 * @param z [] integration variable, momentum fraction in upper rho DA
 * @param xi [] skewness
 * @param mandelstam_t  [GeV^2] momentum transfer, taken at t_min for now
 * @param psq [GeV^2] pomeron hard scale
 * @param Qsq [GeV^2] virtual photon 4mom squared
 * @param helampindex [] index of the helicity amplitude
 * @param omegapol [komegaL/komegaT] long/transverse rho polarization
 * @param gamma [kgammaL/kgammaT] long/transverse photon polarization
 * @param twovector instance that contains all the other functions/variables we need
 */

static void integrandum_omega_general(numint::vector_z &, const double u, const double z, double xi, double mandelstam_t, double scale, 
                                      double psq, double Qsq, int helampindex, TwoVector_Deut::Omega_pol omegapol,
                                      TwoVector_Deut::Photon_pol gamma,  int T_model, TwoVector_Deut& twovector);


PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule; ///< alpha_S from PARTONS

Deut_Conv_GPD_V deut_vector_grid;
Deut_Conv_GPD_T deut_tensor_grid;

};

#endif