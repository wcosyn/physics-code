/*! \file DeuteronMomDistr.hpp 
 * \brief Contains declaration of class DeuteronMomDistr, to compute deuteron momentum distributions
 * \author Wim Cosyn
 * \date 28/08/2012
 * 
 * \addtogroup Knockout
 * @{
 */
#ifndef DEUTERONMOMDISTR_HPP
#define DEUTERONMOMDISTR_HPP

#include <TVector3.h>
#include <TDeuteron.h>
#include <TKinematics2to2.h>
#include <TInterpolatingWavefunction.h>


/*! \brief Has methods to compute a ton of deuteron momentum distributions (distorted too) */
class DeuteronMomDistr{
  
public:
    /*! Constructor
   * \param name Deuteron wave function name, see TDeuteron
   * \param massi [MeV] mass of nucleon that interacts with photon
   * \param offshell which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: no off-shell amplitude, fully suppressed
   * - 4: full off-shell amplitude, no suppression
   * \param sigma [mb] eikonal parameter, total cross section of particle in final-state with spectator
   * \param beta [GeV\f$^-2 \f$] eikonal parameter, slope parameter of particle in final-state with spectator
   * \param epsilon [] eikonal parameter, real part of amplitude of particle in final-state with spectator
   * \param betaoff [GeV\f$^-2 \f$] eikonal parameter, off-shell slope parameter of particle in final-state with spectator
   * \param lambda [GeV\f$^2 \f$] cutoff parameter for off-shell part
   */
  DeuteronMomDistr(std::string name, double massi, int offshell, double sigma, double beta, double epsilon, double betaoff, double lambda);
    /*! Constructor
   * \param name Deuteron wave function name, see TDeuteron
   */
  DeuteronMomDistr(std::string name);
  ~DeuteronMomDistr(); /*!<Destructor */
  /*! Computes plane-wave momentum distribution, does not depend on phi
   * \param kin contains the hadron kinematics
   * \return plane-wave momentum distribution [MeV^-3]
   */
  double getMomDistrpw(TKinematics2to2 &kin) const;
  /*! Computes plane-wave momentum distribution, does not depend on phi
   * \param pvec vector of spectator momentum
   * \return plane-wave momentum distribution [MeV^-3]
   */
  double getMomDistrpw(TVector3 &pvec) const;  
  /*! Computes plane-wave momentum distribution, but off-diagnola in deuteron polarization and spectator momentum!
   * \param pvec vector of spectator momentum
   * \param pvec2 vector of other spectator momentum
   * \param M deuteron polarization (-1,0,1)
   * \param M2 other deuteron polarization (-1,0,1)
   * \return plane-wave momentum distribution [MeV^-3]
   */
  std::complex<double> getMomDistrpwcoh(TVector3 &pvec,TVector3 &pvec2,int M, int M2) const;
  /*! Computes distorted momentum distribution for DIS production off deuteron
   * \param kin contains the hadron kinematics
   * \param phi angle between hadron and electron scattering plane
   * \return distorted momentum distribution [MeV^-3]
   */
  double getMomDistrfsi(TKinematics2to2 &kin, double phi);
  /*! Computes distorted momentum distribution for quasi-elastic production off deuteron
   * \param pvec vector of spectator momentum
   * \param nu virtual photon energy [MeV]
   * \param qvec virtual photon momentum [MeV]
   * \param s invariant mass of reaction [MeV^2]
   * \param massother mass of other particle involved in the FSI [MeV]
   * \return distorted momentum distribution [MeV^-3]
   */
  double getMomDistrfsi(TVector3 &pvec, double nu, double qvec, double s, double massother);
  /*! Set scattering parameters for FSI
   * \param sigmain total cross section [mb]
   * \param betain slope parameter [GeV^-2]
   * \param epsin real part of amplitude
   */  
  void setScatter(double sigmain, double betain, double epsin);
private:
  TDeuteron::Wavefunction *wfref; /*!< contains instance of deuteron wave function*/
  TInterpolatingWavefunction wf; /*!< array of the wave function that gets interpolated */
  double sigma;  /*!< [MeV^-2] total cross section, scattering parameter in FSI */
  double beta;  /*!< [MeV^-2] slope parameter, scattering parameter in FSI*/
  double epsilon; /*!< [] real part of scattering amplitude */
  double betaoff; /*!< [MeV^-2] off-shell slope parameter, scattering parameter in FSI*/
  double lambda; /*!< [MeV^2 off-shell cutoff parameter */
  double massi; /*!< mass of nucleon that gets hit by photon */
  double massr; /*!< mass of spectator nucleon */
  
  double przprime; /*!< pole in the fsi amplitude */
  double Wxprime2; /*!< intermediate invariant mass squared of the X in the FSI rescattering */
   /*! \param offshellset which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: no off-shell amplitude, fully suppressed
   * - 4: full off-shell amplitude, no suppression */
  int offshellset; 
  
  /*! function that gets integrated over q_t (perpendicular to photon momentum) 
   * \param qt [MeV] in fsi: momentum transfer component perpendicular to photon momentum 
   * \param result result: contains fsi part of DIS distorted momentum distribution
   * \param ap variable parameter list
   */
  void totdens_qt(const double qt, std::complex<double>* result, va_list ap);
   /*! function that gets integrated over q_phi (perpendicular to photon momentum)
   * \param qphi [] azimuthal angle of momentum transfer in fsi perpendicular to photon momentum 
   * \param result result: contains fsi part of DIS distorted momentum distribution
   * \param ap variable parameter list
   */
  void totdens_qphi(const double qphi, std::complex<double>* result, va_list ap);
  /*! function that gets integrated over q_t (perpendicular to photon momentum)
   * \param qt [MeV] momentum transfer component perpendicular to photon momentum 
   * \param result result: contains fsi part of quasi-elastic momentum distribution
   * \param ap variable parameter list
   */
  void totdens_qt_simple(const double qt, std::complex<double>* result, va_list ap);
   /*! function that gets integrated over q_phi (perpendicular to photon momentum)
   * \param qphi [] azimuthal angle of momentum transfer in fsi perpendicular to photon momentum 
   * \param result result: contains fsi part of QE distorted momentum distribution
   * \param ap variable parameter list
   */
  void totdens_qphi_simple(const double qphi, std::complex<double>* result, va_list ap);
  /*! recursive method to find the pole in the fsi integration, longitudinal part,
   * also determines intermediate mass
   * \param pt [MeV] final transverse spectator momentum
   * \param Er [MeV] final spectator on-shell energy
   * \param pkin contains hadron kinematics
   */
  void get_przprime(double pt, double Er, TKinematics2to2 *pkin);
  /*! gives you the scatter amplitude of the final-state interaction
   * \param t [MeV^2] momentum transfer squared
   * \return \f$ \sigma_{tot} (I+\epsilon) e^{\beta t/2} \f$
   */
  std::complex<double> scatter(double t);
  
};
/** @} */
#endif