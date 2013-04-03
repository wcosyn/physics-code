/*! \file DeuteronCross.hpp 
 * \brief Contains declaration of class InclusiveCross, computes deuteron inclusive cross sections
 * \author Wim Cosyn
 * \date 29/08/2012
 * 
 * \addtogroup DIS
 * @{
 */

#ifndef INCLUSIVECROSS_HPP
#define INCLUSIVECROSS_HPP

#include <string>
#include <TElectronKinematics.h>
#include <DeuteronStructure.hpp>
#include <TDeuteron.h>
#include <numint/numint.hpp>


/*! \brief A class that computes inclusive deuteron cross sections */
class InclusiveCross{
  
public:
    /*! Constructor
   * \param elec Electron kinematics, see TElectronKinematics
     * \param proton photon interacting with proton [1] or neutron [0]
   * \param wavename Deuteron wave function name, see TDeuteron
   * \param strucname which structure function parametrization ["CB"=Chrisy-Bosted, "SLAC", "Alekhin"]
   * \param symm parameter that controls symmetric or not rescatterings in the FSI <BR>
   * - 1: no mass terms
   * - 0: mass terms always included
   * - -1: mass terms included like in semi-inclusive DIS (mass must increase)
   * \param offshell which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: no off-shell amplitude, fully suppressed
   * - 4: full off-shell amplitude, no suppression
   * \param sigmain [mb] eikonal parameter, total cross section of particle in final-state with spectator
   * \param betain [GeV^-2] eikonal parameter, slope parameter of particle in final-state with spectator
   * \param epsilonin [] eikonal parameter, real part of amplitude of particle in final-state with spectator
   * \param betaoffin [GeV^-2] eikonal parameter, off-shell slope parameter of particle in final-state with spectator
   * \param lambdain [GeV^2] cutoff parameter for off-shell part
   */
  InclusiveCross(bool proton, std::string strucname, std::string wavename, TElectronKinematics &elec, int symm, int offshell, 
		 double sigmain=40.,double betain=8., double epsilonin=-0.5, double betaoffin=8., double lambdain=1.2);
  ~InclusiveCross(); /*!< Destructor */
  /*! Calculates the plane-wave inclusive cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \return [] inclusive cross section without prefactors (structure functions times momentum distribution integrated)
   */
  double calc_F2Dinc(double Q2,double x);
  /*! Calculates the fsi inclusive cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution
   * \param fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   */
  void calc_F2DincFSI(double &fsi1, double &fsi2, double Q2,double x);
  //void calc_F2DincFSI2(double &fsi1, double &fsi2, double Q2,double x);
  
  
private:
  double sigma;  /*!< [MeV^-2] total cross section, scattering parameter in FSI */
  double beta;  /*!< [MeV^-2] slope parameter, scattering parameter in FSI*/
  double epsilon; /*!< [] real part of scattering amplitude */
  double betaoff; /*!< [MeV^-2] off-shell slope parameter, scattering parameter in FSI*/
  double lambda; /*!< [MeV^2 off-shell cutoff parameter */
  double massi; /*!< [MeV] mass of nucleon that gets hit by photon */
  double massr; /*!< [MeV] mass of spectator nucleon */
    /*! \param symm parameter that controls symmetric or not rescatterings in the FSI <BR>
   * - 1: no mass terms
   * - 0: mass terms always included
   * - -1: mass terms included like in semi-inclusive DIS (mass must increase)
   */
  int symm; 
  
  double przprime; /*!< [MeV] longitudinal spectator momentum after rescattering */
  double prz; /*!< [MeV] lingitudinal spectator momentum before rescattering */
  double Wxprime2; /*!< [MeV^2] invariant mass of the X that is paired with the spectator that gets integrated first */
  double otherWx2; /*!< [MeV^2] invariant mass of the other X that depends on the first ps integration */
  /*! which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: no off-shell amplitude, fully suppressed
   * - 4: full off-shell amplitude, no suppression 
   */
  int offshellset;

  TDeuteron::Wavefunction *wf; /*!< pointer to deuteron wave funtion, see TDeuteron */
  TElectronKinematics electron; /*!< electron kinematis, see TElectronKinematics */
  DeuteronStructure structure;  /*!< deuteron structure functions object, see DeuteronStructure */

  void int_pr(double pr, double *result, va_list ap);
  void int_costheta_incl(double costheta, double *result, va_list ap);
  void int_pr_fsi(double pr, double *result, va_list ap);
  void int_costheta_incl_fsi(double costheta, double *result, va_list ap);
  void int_qt(double qt, double *results, va_list ap);
  void int_qphi(double qphi, double *results, va_list ap); 
  /*! recursive method to find the pole in the fsi integration, longitudinal part,
   * also determines intermediate mass
   * \param pt [MeV] final transverse spectator momentum
   * \param Er [MeV] final spectator on-shell energy
   * \param pkin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param first parameter that controls whether we take the pole before or after the rescattering
   */
  void get_prz(double pt, double Er, TKinematics2to2 *pkin, int first);
  /*! recursive method to find the pole in the fsi integration, longitudinal part,
   * also determines intermediate mass.  Based on new method in bjorken method where pr_z \approx m(x-1)
   * \param pt [MeV] final transverse spectator momentum
   * \param Er [MeV] final spectator on-shell energy
   * \param pkin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param first parameter that controls whether we take the pole before or after the rescattering
   */
  void get_prz2(double pt, double Er, TKinematics2to2 *pkin, int first);
 /*! gives you the scatter amplitude of the final-state interaction
   * \param t [MeV^2] momentum transfer squared
   * \return \f$ \sigma_{tot} (I+\epsilon) e^{\beta t/2} \f$
   */
  std::complex<double> scatter(double t);

//   void int_pperp_fsi(double pperp, double *result, va_list ap);
//   void int_qt_bis(double qt, double *results, va_list ap);
//   void int_qphi_bis(double qphi, double *results, va_list ap); 
//   void get_przs(double pperp2, double pperp2other, double qvec, double nu, int first);
  /*! struct that is used for integrators plane wave ones*/
  struct Ftor_planewave {

    /*! integrandum function */
    static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
      Ftor_planewave &p = * (Ftor_planewave *) param;
      p.f(ret,x[0],x[1],*p.cross,p.Q2,p.x);
    }
    InclusiveCross *cross;/*!< pointer to InclusiveCross instance that contains all */
    double Q2; /*!< [MeV^2] momentum transfer */
    double x; /*!< [] Bjorken x */
    /*! integrandum 
    * \param res results
    * \param pnorm first integration variable
    * \param costheta second integration variable
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    */
    void (*f)(numint::vector_d & res, double pnorm, double costheta, InclusiveCross& cross, double Q2, double x);
  };
  
   /*! integrandum function (clean ones)*/
  static void planewave_int(numint::vector_d & results, double pnorm, double costheta, InclusiveCross& cross, double Q2, double x);
 
  
  /*! struct that is used for integrators fsi on-shell ones*/
  struct Ftor_FSI {

    /*! integrandum function */
    static void exec(const numint::array<double,4> &x, void *param, numint::vector_d &ret) {
      Ftor_FSI &p = * (Ftor_FSI *) param;
      p.f(ret,x[0],x[1],x[2],x[3],*p.cross,p.Q2,p.x);
    }
    InclusiveCross *cross;/*!< pointer to InclusiveCross instance that contains all */
    double Q2; /*!< [MeV^2] momentum transfer */
    double x; /*!< [] Bjorken x */
    /*! integrandum 
    * \param res results
    * \param pnorm first integration variable
    * \param costheta second integration variable
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    */
    void (*f)(numint::vector_d & res, double pnorm, double costheta, double qt, double qphi, InclusiveCross& cross, double Q2, double x);
  };
  
   /*! integrandum function (clean ones)*/
  static void FSI_int(numint::vector_d & results, double pnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x);
 
  
};


/*! @} */
#endif




