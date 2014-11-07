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
#include <numint/numint.hpp>
#include <LightConeKin2to2.hpp>

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
   * \param beta [GeV^-2] eikonal parameter, slope parameter of particle in final-state with spectator
   * \param epsilon [] eikonal parameter, real part of amplitude of particle in final-state with spectator
   * \param betaoff [GeV^-2] eikonal parameter, off-shell slope parameter of particle in final-state with spectator
   * \param lambda [GeV^2] cutoff parameter for off-shell part
   * \param looplimit max number of tries in loop to get prz pole
   */
  DeuteronMomDistr(std::string name, double massi, int offshell, double sigma, double beta, 
		   double epsilon, double betaoff, double lambda, int looplimit);
    /*! Constructor
   * \param name Deuteron wave function name, see TDeuteron
   */
  DeuteronMomDistr(std::string name);
  DeuteronMomDistr(); /*!< Default Constructor */
  DeuteronMomDistr(const DeuteronMomDistr&); /*!< Copy Constructor */
  DeuteronMomDistr& operator=(const DeuteronMomDistr&); /*!< assignment operator */
  ~DeuteronMomDistr(); /*!<Destructor */
  /*! Computes plane-wave momentum distribution in LC formalism according to Frankfurt and Strikman, Phys.Rep.88
   * \param[in] kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \return plane-wave momentum distribution [MeV^-3]
   */
  double getLCMomDistrpw(TKinematics2to2 &kin) const;
  /*! Computes plane-wave momentum distribution in LC formalism according to Frankfurt and Strikman, Phys.Rep.88
  * \param p [MeV] relative momentum in deuteron rest frame
  * \return plane-wave momentum distribution [MeV^-3]
   */
  double getLCMomDistrpw(TVector3 p) const;
  /*! Computes plane-wave momentum distribution, does not depend on phi, includes flux factor for baryon conservation
   * \param[in] kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \return plane-wave momentum distribution [MeV^-3]
   */
  double getMomDistrpw(TKinematics2to2 &kin) const;
  /*! Computes plane-wave momentum distribution in light-cone formalism, does not depend on phi
   * \param[in] kin LC kinematics object containing the gamma+D->X+N kinematics <BR>
   * \return plane-wave momentum distribution [MeV^-3]
   */
  double getMomDistrpwLC(LightConeKin2to2 &kin) const;
  /*! Computes plane-wave Azz distribution, does not depend on phi, includes flux factor for baryon conservation
   * \param[in] kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \return plane-wave Azz deuteron distribution [MeV^-3]
   */
  double getAzzDistrpw(TKinematics2to2 &kin) const;
  /*! Computes distorted Azz deuteron distribution for DIS production off deuteron
   * \param[in] kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param phi angle between hadron and electron scattering plane
   * \return distorted Azz deuteron distribution [MeV^-3]
   */
  double getAzzDistrfsi(TKinematics2to2 &kin, double phi);
  /*! Computes plane-wave momentum distribution, does not depend on phi
   * \param pvec vector of spectator momentum
   * \return plane-wave momentum distribution [MeV^-3]
   */
  double getMomDistrpw(TVector3 pvec) const;  
  /*! Computes plane-wave momentum distribution, but off-diagonal in deuteron polarization and spectator momentum!
   * \param pvec vector of spectator momentum
   * \param pvec2 vector of other spectator momentum
   * \param M deuteron polarization (-1,0,1)
   * \param M2 other deuteron polarization (-1,0,1)
   * \return plane-wave momentum distribution [MeV^-3]
   */
  std::complex<double> getMomDistrpwcoh(TVector3 pvec,TVector3 pvec2,int M, int M2) const;
  /*! Computes distorted momentum distribution for DIS production off deuteron
   * \param[in] kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param phi angle between hadron and electron scattering plane
   * \return distorted momentum distribution [MeV^-3]
   */
  double getMomDistrfsi(TKinematics2to2 &kin, double phi);
  /*! Computes distorted momentum distribution for DIS production off deuteron
   * \param[in] kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * \return distorted momentum distribution [MeV^-3]
   */
  double getMomDistrfsiLC(LightConeKin2to2 &kin);
  /*! Computes distorted momentum distribution for quasi-elastic production off deuteron
   * \param pvec vector of spectator momentum
   * \param nu virtual photon energy [MeV]
   * \param qvec virtual photon momentum [MeV]
   * \param s invariant mass of reaction [MeV^2]
   * \param massother mass of other particle involved in the FSI [MeV]
   * \return distorted momentum distribution [MeV^-3]
   */
  double getMomDistrfsi(TVector3 pvec, double nu, double qvec, double s, double massother);
  /*! Set scattering parameters for FSI
   * \param sigmain total cross section [mb]
   * \param betain slope parameter [GeV^-2]
   * \param epsin real part of amplitude
   */  
  void setScatter(double sigmain, double betain, double epsin);
  TDeuteron::Wavefunction *getDeuteronwf(){return wfref;}
private:
  std::string wfname; /*!< name of deuteron wf parametrization, see TDeuteron for possibilities */
  TDeuteron::Wavefunction *wfref; /*!< contains instance of deuteron wave function*/
  TInterpolatingWavefunction wf; /*!< array of the wave function that gets interpolated */
  double sigma;  /*!< [MeV^-2] total cross section, scattering parameter in FSI */
  double beta;  /*!< [MeV^-2] slope parameter, scattering parameter in FSI*/
  double epsilon; /*!< [] real part of scattering amplitude */
  double betaoff; /*!< [MeV^-2] off-shell slope parameter, scattering parameter in FSI*/
  double lambda; /*!< [MeV^2 off-shell cutoff parameter */
  double massi; /*!< [MeV] mass of nucleon that gets hit by photon */
  double massr; /*!< [MeV] mass of spectator nucleon */
  
  double przprime; /*!< pole in the fsi amplitude */
  double k_z_prime; /*!< pole in the fsi amplitude */
  double Wxprime2; /*!< intermediate invariant mass squared of the X in the FSI rescattering */
   /*! \param offshellset which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: no off-shell amplitude, fully suppressed
   * - 4: full off-shell amplitude, no suppression */
  int offshellset; 
  int looplimit; /*!< max number of tries in loop to get prz pole */
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
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   */
  void get_przprime(double pt, double Er, TKinematics2to2 & kin);
  /*! recursive method to find the pole in the LC fsi integration, longitudinal part,
   * also determines intermediate mass
   * \param k_perp_prime [MeV] threevector of initial spectator LC rescaled perp momentum
   * \param kin LC kinematics object containing the gamma+D->X+N kinematics <BR>
   */
  void get_k_z_prime(TVector3 &k_perp_prime, LightConeKin2to2 & kin);
  /*! gives you the scatter amplitude of the final-state interaction
   * \param t [MeV^2] momentum transfer squared
   * \return \f$ \sigma_{tot} (I+\epsilon) e^{\beta t/2} \f$
   */
  std::complex<double> scatter(double t);
  
  /*! struct that is used for integrators fsi mom distribution*/
  struct Ftor_FSI {

    /*! integrandum function */
    static void exec(const numint::array<double,2> &x, void *param, numint::vector_z &ret) {
      Ftor_FSI &p = * (Ftor_FSI *) param;
      p.f(ret,x[0],x[1],*p.momdistr, *p.kin,p.M,p.spinr,p.Er,p.phi);
    }
    DeuteronMomDistr *momdistr;/*!< pointer to DeuteronMomDistr instance that contains all */
    TKinematics2to2 *kin; /*!< kinematics object */
    int M; /*!< spin projection of the deuteron */
    int spinr; /*!< spin proj of spectator */
    double Er; /*!< [MeV] on-shell energy of spectator */
    double phi; /*!< angle between hadron and electron plane */
    /*! integrandum 
    */
    void (*f)(numint::vector_z & res, double qt, double qphi, DeuteronMomDistr &momdistr,
	      TKinematics2to2 &kin, int M, int spinr, double Er, double phi);
  };
  
  /*! struct that is used for integrators fsi LC mom distribution*/
  struct Ftor_FSILC {

    /*! integrandum function */
    static void exec(const numint::array<double,2> &x, void *param, numint::vector_z &ret) {
      Ftor_FSILC &p = * (Ftor_FSILC *) param;
      p.f(ret,x[0],x[1],*p.momdistr, *p.kin,p.M,p.spinr);
    }
    DeuteronMomDistr *momdistr;/*!< pointer to DeuteronMomDistr instance that contains all */
    LightConeKin2to2 *kin; /*!< kinematics object */
    int M; /*!< spin projection of the deuteron */
    int spinr; /*!< spin proj of spectator */
    /*! integrandum 
    */
    void (*f)(numint::vector_z & res, double qt, double qphi, DeuteronMomDistr &momdistr,
	      LightConeKin2to2  &kin, int M, int spinr);
  };
  
  
  
   /*! integrandum function (clean ones)*/
  static void FSILC_int(numint::vector_z & results, double qt, double qphi, DeuteronMomDistr &momdistr,
		      LightConeKin2to2  &kin, int M, int spinr);
  static void FSI_int(numint::vector_z & results, double qt, double qphi, DeuteronMomDistr &momdistr,
	                      TKinematics2to2 &kin, int M, int spinr, double Er, double phi);
  
  
};
/** @} */
#endif