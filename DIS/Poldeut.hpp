/*! \file Poldeut.hpp 
 * \brief Contains declaration of class Poldeut, used to calculate observables for tagged spectator DIS with a polarized deuteron
 * Equation numbers refer to the poldeut.pdf draft.
 * \author Wim Cosyn
 * \date 31/03/2016
 * 
 * \addtogroup DIS
 * @{
 */

#ifndef POLDEUT_HPP
#define POLDEUT_HPP

#include <TDeuteron.h>
#include <TInterpolatingWavefunction.h>

#include <string>

/*! \brief A class that computes tagged spectator deuteron cross sections */
class Poldeut{
  
public:
    
    /*! Constructor
   * \param wfname Deuteron wave function name, see TDeuteron constructor for possibilities ["AV18", "Paris", etc.]
   * \param strucname which structure function parametrization ["CB"=Chrisy-Bosted, "SLAC", "Alekhin"], see NucleonStructure 
   */
  Poldeut(std::string wfname, std::string strucname);
  ~Poldeut(); /*!< Destructor */
  Poldeut(const Poldeut&); /*!< Copy Constructor */
  Poldeut& operator=(const Poldeut&); /*!< assignment operator */
  
  void set_S(bool s){ s_wave=s;} /*!< setter for to include s_wave or not */
  void set_D(bool d){ d_wave=d;} /*!< setter for to include d_wave or not */
  void set_melosh(bool m) {melosh=m;} /*!< setter for to include melosh rotation or not */

   /*! calculate the asymmetry F_LSL/(F_T+epsilon F_L), collinear frame, deuteron rest frame
   * \param x [] Bjorken x defined as 2Q^2/(P_d*q)
   * \param Q2 [MeV^2] minus virtual photon fourmomentum squared
   * \param y [] (Pd*q)/(Pd*l_in) kinematical variable
   * \param alpha_s [] lightcone momentum fraction defined as 2p_s^+/p_d^+
   * \param pt [MeV] perp momentum of the tagged nucleon
   * \param proton struck nucleon is proton [1] or neutron [0]
   * \param[out] F_LSL structure function for pol electron and long polarized target that has no phi dep.
   * \param[out] F_U F_U,T+epsilon F_U,L 
   */
  void calc_Double_Asymm(double x, double Q2, double y, double alpha_s, double pt, bool proton, double &F_LSL, double &F_U);

/**
 * @brief compare the nuclear densitiy ratio entering in the Azz asymmetry, between lightcone expression and its nonrelativistic reduction
 * 
 * @param[out] lfratio ''exact'' lightfront expression
 * @param[out] nrratio nonrelativistic reduction
 * @param alpha_p [] detected spectator lf +momentum fraction
 * @param pt  [MeV] detected spectator perp momentum
 */
  void Tensor_Compare_nonrel(double &lfratio, double &nrratio, double alpha_p, double pt);




private:
  std::string strucname; /*!< structure function parametrization */
  std::string wfname; /*!< name of deuteron wf parametrization, see TDeuteron for possibilities */
  TDeuteron::Wavefunction *wfref; /*!< contains instance of deuteron wave function*/
  TInterpolatingWavefunction wf; /*!< array of the wave function that gets interpolated */
  bool s_wave; /*!< include deuteron s wave */
  bool d_wave; /*!< include deuteron d wave */
  bool melosh; /*!< include effect of melosh rotation */
  
  
  /*! get lightcone densities that enter in the structure functions (melosh rotation not included here!!)
   * We do include the factors E_k/alpha_i^2 here for convenience
   * \param knorm [MeV] norm of the relative 3-momentum entering in the lc density
   * \param ktheta [rad] angle theta of the relative 3-momentum with the spin quantization axis
   * \param[out] rhou [MeV-3] unpolarized lightcone deuteron density (Eq C.1)
   * \param[out] rho_l_z [MeV-3] polarized lightcone deuteron density (Eq C.4)
   * \param[out] rho_l_x [MeV-3] polarized lightcone deuteron density (Eq C.5)
   * \param[out] rho_tensor_u [MeV-3] tensor polarized lightcone density (Eq C.16)
  */
  void getDensities(double knorm, double ktheta, double &rhou, double &rho_l_z, double &rho_l_x, double &rho_tensor_u);
  
  /*! obtain deuteron radial s_wave (normalized as int k^2 dk U(k)^2=1) 
   * \param k [MeV] relative momentum
   * \return [MeV^-3/2] radial s-wave U(k)
   */
  double getU(double k){ return (s_wave? wf.GetUp(k):0.);}
  /*! obtain deuteron radial d_wave (normalized as int k^2 dk W(k)^2=1)
   * \param k [MeV] relative momentum
   * \return [MeV^-3/2] radial d-wave W(k)
   */
  double getW(double k){ return (d_wave? wf.GetWp(k):0.);}
  /*! carries out the melosh rotation for the lightcone densities
   * \param k_perp [MeV] x-component of the relative lightcone momentum
   * \param k_plus [MeV} +-component of the relative lightcone momentum
   * \param Ek [MeV] on-shell energy of nucleon with relative lightcone three-momentum 
   * \param[out] rho_l_z [MeV^-3] Melosh rotated polarized lightcone deuteron density (Eq 3.18)
   * \param[out] rho_l_x [MeV^-3] Melosh rotated polarized lightcone deuteron density (Eq 3.19)
   */ 
  void Melosh_rot(double k_perp, double k_plus, double Ek, double &rho_l_z, double &rho_l_x); 
  

};
/*! @} */
#endif 