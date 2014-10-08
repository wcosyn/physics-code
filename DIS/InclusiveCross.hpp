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
   * \param [in] elec Electron kinematics, see TElectronKinematics
   * \param proton photon interacting with proton [1] or neutron [0]
   * \param wavename Deuteron wave function name, see TDeuteron
   * \param strucname which structure function parametrization ["CB"=Chrisy-Bosted, "SLAC", "Alekhin"]
   * \param [in] res [MeV] vector that contains resonance masses used in FSI
   * \param offshell which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: suppression with mass difference between resonance mass and "true" W (used in paper)
   * - 4: full off-shell amplitude, no suppression
   * \param sigmain [mb] eikonal parameter, total cross section of particle in final-state with spectator
   * \param betain [GeV^-2] eikonal parameter, slope parameter of particle in final-state with spectator
   * \param epsilonin [] eikonal parameter, real part of amplitude of particle in final-state with spectator
   * \param betaoffin [GeV^-2] eikonal parameter, off-shell slope parameter of particle in final-state with spectator
   * \param lambdain [GeV^2] cutoff parameter for off-shell part
   */
  InclusiveCross(bool proton, std::string strucname, std::string wavename, TElectronKinematics &elec,
		 std::vector<double> & res, int offshell, 
		 double sigmain=40.,double betain=8., double epsilonin=-0.5, double betaoffin=8., double lambdain=1.2);
  ~InclusiveCross(); /*!< Destructor */
  /*! Calculates the plane-wave inclusive cross on-shell section without any prefactors, just the deuteron structure functions time
   * momentum distribution
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \return [] inclusive cross section without prefactors (structure functions times momentum distribution integrated)
   */
  double calc_F2Dinc(double Q2,double x);
  /*! Calculates Azz, also gives you the cross section (without prefactors)
   * \param[out] Azz [] Azz observable
   * \param[out] cross [] cross section
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param thetapol [rad] angle between virtual photon momentum and angle of deuteron polarization
   */
  void calc_Azzinc(double &Azz, double &cross, double Q2,double x, double thetapol);
  /*! Calculates the fsi inclusive cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   */
  void calc_F2DincFSI(double &fsi1, double &fsi2, double Q2,double x);
  /*! Calculates the fsi inclusive off-shell cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   */
  void calc_F2DincFSI_off(double &fsi1, double &fsi2, double Q2,double x);
  /*! Calculates the fsi inclusive Azz & cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution
   * \param [out] azz1 [] inclusive fsi Azz
   * calculated in one of two ways (see notes)
   * \param [out] azz2 [] other fsi Azz
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param thetapol [rad] angle between virtual photon momentum and angle of deuteron polarization
   */
  void calc_AzzincFSI(double &azz1, double&azz2, double &fsi1, double &fsi2, double Q2,double x, double thetapol);
  /*! Calculates the fsi inclusive off-shell Azz & cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution
   * \param [out] azz1 [] inclusive fsi Azz
   * calculated in one of two ways (see notes)
   * \param [out] azz2 [] other fsi Azz
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param thetapol [rad] angle between virtual photon momentum and angle of deuteron polarization
   */
  void calc_AzzincFSI_off(double &azz1, double&azz2, double &fsi1, double &fsi2, double Q2,double x, double thetapol);


  /*! Calculates the fsi inclusive cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution <BR>
   * Uses broad resonance with Gaussian distribution
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param [in] centrals [MeV] central value of broad resonance
   * \param [in] width [MeV] widths of broad resonance, Gaussian shape assumed
   */
  void calc_F2DincFSI_Gauss(double &fsi1, double &fsi2, double Q2,double x, 
			    std::vector<double> &centrals, std::vector<double> &width);
  /*! Calculates the fsi inclusive off-shell cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution <BR>
   * Uses broad resonance with Gaussian distribution
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param [in] centrals [MeV] central value of broad resonance
   * \param [in] width [MeV] widths of broad resonance, Gaussian shape assumed
   */
  void calc_F2DincFSI_Gauss_off(double &fsi1, double &fsi2, double Q2,double x,
			    std::vector<double> &centrals, std::vector<double> &width);
  /*! Calculates the fsi inclusive cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution <BR>
   * Uses broad resonance with uniform distribution
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param [in] centrals [MeV] central value of broad resonance
   * \param [in] width [MeV] widths of broad resonance, Gaussian shape assumed
   */
  void calc_F2DincFSI_uniform(double &fsi1, double &fsi2, double Q2,double x,
			    std::vector<double> &centrals, std::vector<double> &width);
  /*! Calculates the fsi inclusive off-shell cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution <BR>
   * Uses broad resonance with uniform distribution
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param [in] centrals [MeV] central value of broad resonance
   * \param [in] width [MeV] widths of broad resonance, Gaussian shape assumed
   */
  void calc_F2DincFSI_uniform_off(double &fsi1, double &fsi2, double Q2,double x, 
			    std::vector<double> &centrals, std::vector<double> &width);
  /*! Calculates the fsi inclusive Azz & cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution <BR>
   * Uses broad resonance with uniform distribution
   * \param [out] azz1 [] inclusive fsi Azz
   * calculated in one of two ways (see notes)
   * \param [out] azz2 [] other fsi Azz
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param [in] centrals [MeV] central value of broad resonance
   * \param [in] width [MeV] widths of broad resonance, Gaussian shape assumed
   */
  void calc_AzzincFSI_uniform(double &azz1, double &azz2, double &fsi1, double &fsi2, double Q2,double x,
			    std::vector<double> &centrals, std::vector<double> &width);
  /*! Calculates the fsi inclusive off-shell Azz & cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution <BR>
   * Uses broad resonance with uniform distribution
   * \param [out] azz1 [] inclusive fsi Azz
   * calculated in one of two ways (see notes)
   * \param [out] azz2 [] other fsi Azz
   * \param [out] fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param [out] fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param [in] centrals [MeV] central value of broad resonance
   * \param [in] width [MeV] widths of broad resonance, Gaussian shape assumed
   */
  void calc_AzzincFSI_uniform_off(double &azz1, double &azz2, double &fsi1, double &fsi2, double Q2,double x, 
			    std::vector<double> &centrals, std::vector<double> &width);
  
  
  
  void setOffshell(const int offshell){offshellset=offshell;} /*!< set offshell parametrization */
  TDeuteron::Wavefunction* getDeutwf() const{return wf;} /*!< get instance to deuteron wave function */
  /*! set resonance value
   * \param it index of resonance
   * \param value [MeV] resonance mass 
   */
  void setResonance(const size_t it, double value) {resonances[it]=value; return;}
  
  
private:
  double sigma;  /*!< [MeV^-2] total cross section, scattering parameter in FSI */
  double beta;  /*!< [MeV^-2] slope parameter, scattering parameter in FSI*/
  double epsilon; /*!< [] real part of scattering amplitude */
  double betaoff; /*!< [MeV^-2] off-shell slope parameter, scattering parameter in FSI*/
  double lambda; /*!< [MeV^2 off-shell cutoff parameter */
  double massi; /*!< [MeV] mass of nucleon that gets hit by photon */
  double massr; /*!< [MeV] mass of spectator nucleon */

  std::vector<double> resonances; /*!<vector with resonance masses included in the FSI */
  
  double Wxprime2; /*!< [MeV^2] invariant mass of the X that is paired with the spectator that gets integrated first */
  double otherWx2; /*!< [MeV^2] invariant mass of the other X that depends on the first ps integration */
  /*! which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: suppressed with mass difference resonance and produced W, diffractive with beta (used in article)
   * - 4: full off-shell amplitude, no suppression 
   */
  int offshellset;

  TDeuteron::Wavefunction *wf; /*!< pointer to deuteron wave funtion, see TDeuteron */
  TElectronKinematics electron; /*!< electron kinematis, see TElectronKinematics */
  DeuteronStructure structure;  /*!< deuteron structure functions object, see DeuteronStructure */

  /*! method to find the pole in the fsi integration, longitudinal part,
   * with resonance mass as input.  Essentially solves Eq. (29) of the paper
  * \param pt2 [MeV^2] final transverse spectator momentum sq
  * \param W_sq [MeV^2] mass of resonance entering
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \return [MeV] on-shell pole p_s,z value
   */
  double get_prz_res(double pt2, double W_sq, TKinematics2to2 & kin);
  /*! method to find the pole in the fsi integration, longitudinal part,
   * with resonance mass as input.  Essentially solves Eq. (29) of the paper
  * \param pt2 [MeV^2] final transverse spectator momentum sq
  * \param W_sq [MeV^2] mass of resonance entering
  * \param Q2 [MeV^2] virtual photon Q^2
   * \param nu [MeV] virtual photon energy
   * \param qvec [MeV] virtual photon momentum
   * \return [MeV] on-shell pole p_s,z value
   */
  double get_prz_res(double pt2, double W_sq, double Q2, double nu, double qvec);
 /*! gives you the scatter amplitude of the final-state interaction
   * \param t [MeV^2] momentum transfer squared
   * \return [MeV^-2] \f$ \sigma_{tot} (I+\epsilon) e^{\beta t/2} \f$
   */
  std::complex<double> scatter(double t);
  /*! gives you sigma in the parametrization we got from deeps, Q^2 dependence included
   * \return [MeV^-2] sigma
   * \param W_sq [MeV^2] invariant mass squared of scatterer
   * \param Q2 [MeV^2] four-momentum transfer squared */
  
  static double sigmaparam(double W_sq, double Q2)  ;
  
  /*! struct that is used for integrators plane wave ones*/
  struct Ftor_planewave {

    /*! integrandum function */
    static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
      Ftor_planewave &p = * (Ftor_planewave *) param;
      p.f(ret,x[0],x[1],*p.cross,p.Q2,p.x,p.thetapol);
    }
    InclusiveCross *cross;/*!< pointer to InclusiveCross instance that contains all */
    double Q2; /*!< [MeV^2] momentum transfer */
    double x; /*!< [] Bjorken x */
    double thetapol; /*!< [rad] angle between q vector and deuteron polarization */
    /*! integrandum 
    * \param res results
    * \param pnorm first integration variable
    * \param costheta second integration variable
    * \param cross instance of InclusiveCross where we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param thetapol [rad] angle between q vector and deuteron polarization 
    */
    void (*f)(numint::vector_d & res, double pnorm, double costheta, InclusiveCross& cross, double Q2, double x, double thetapol);
  };
  
   /*! integrandum function (clean ones), Eq (19) of the paper
    * \param results results
    * \param pnorm first integration variable
    * \param costheta second integration variable
    * \param cross instance of InclusiveCross where we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param thetapol [rad] angle between q vector and deuteron polarization 
    */
  static void planewave_int(numint::vector_d & results, double pnorm, double costheta, 
			    InclusiveCross& cross, double Q2, double x, double thetapol);
   /*! integrandum function (clean ones), Eq (19) of the paper [Azz version]
    * \param results results [0=Azz, 1=cross section]
    * \param pnorm first integration variable
    * \param costheta second integration variable
    * \param cross instance of InclusiveCross where we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param thetapol [rad] angle between q vector and deuteron polarization 
    */
  static void Azzplanewave_int(numint::vector_d & results, double pnorm, double costheta, 
			    InclusiveCross& cross, double Q2, double x, double thetapol);
 
  
  /*! struct that is used for integrators fsi on- and off-shell ones*/
  struct Ftor_FSI {

    /*! integrandum function */
    static void exec(const numint::array<double,4> &x, void *param, numint::vector_d &ret) {
      Ftor_FSI &p = * (Ftor_FSI *) param;
      p.f(ret,x[0],x[1],x[2],x[3],*p.cross,p.Q2,p.x,p.thetapol, p.it,p.it2);
    }
    InclusiveCross *cross;/*!< pointer to InclusiveCross instance that contains all */
    double Q2; /*!< [MeV^2] momentum transfer */
    double x; /*!< [] Bjorken x */
    double thetapol; /*!< [rad] angle between q vector and deuteron polarization */

    size_t it; /*!< iterator for resonance */
    size_t it2; /*!< iterator for resonance */
    /*! integrandum 
    * \param res results
    * \param pnorm first integration variable
    * \param costheta second integration variable
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param thetapol [rad] angle between q vector and deuteron polarization 
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    */
    void (*f)(numint::vector_d & res, double pnorm, double costheta, double qt, double qphi, InclusiveCross& cross, 
	      double Q2, double x, double thetapol, size_t it, size_t it2);
  };
  
   /*! integrandum function for on-shell contribution to the FSI amplitude (Eq. (30) of the paper)
    * \param results results
    * \param pnorm [MeV] norm of spectator momentum
    * \param costheta [] polar cos(theta) of spectator momentum
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param thetapol [rad] angle between q vector and deuteron polarization 
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    */
  static void FSI_int(numint::vector_d & results, double pnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, double thetapol, size_t it, size_t it2);
   /*! integrandum function for off-shell contribution to the FSI amplitude (Eq. (33) of the paper)
    * \param results results
    * \param prt [MeV] norm of transverse spectator momentum (initial)
    * \param W [MeV] invariant mass where the structure function is evaluated
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param thetapol [rad] angle between q vector and deuteron polarization 
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    */
  static void FSI_int_off(numint::vector_d & results, double prt, double W, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, double thetapol, size_t it, size_t it2);
 
   /*! integrandum function for on-shell contribution to the FSI amplitude (Eq. (30) of the paper)
    * \param results results
    * \param pnorm [MeV] norm of spectator momentum
    * \param costheta [] polar cos(theta) of spectator momentum
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param thetapol [rad] angle between q vector and deuteron polarization 
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    */
  static void AzzFSI_int(numint::vector_d & results, double pnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, double thetapol, size_t it, size_t it2);
   /*! integrandum function for off-shell contribution to the FSI amplitude (Eq. (33) of the paper)
    * \param results results
    * \param prt [MeV] norm of transverse spectator momentum (initial)
    * \param W [MeV] invariant mass where the structure function is evaluated
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param thetapol [rad] angle between q vector and deuteron polarization 
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    */
  static void AzzFSI_int_off(numint::vector_d & results, double prt, double W, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, double thetapol, size_t it, size_t it2);
 
  /*! struct that is used for integrators fsi on- and off-shell ones for a DIS broad resonance distribution*/
  struct Ftor_FSI_distr {

    /*! integrandum function */
    static void exec(const numint::array<double,5> &x, void *param, numint::vector_d &ret) {
      Ftor_FSI_distr &p = * (Ftor_FSI_distr *) param;
      p.f(ret,x[0],x[1],x[2],x[3],x[4],*p.cross,p.Q2,p.x,p.it,p.it2,p.centrals,p.widths);
    }
    InclusiveCross *cross;/*!< pointer to InclusiveCross instance that contains all */
    double Q2; /*!< [MeV^2] momentum transfer */
    double x; /*!< [] Bjorken x */
    size_t it; /*!< iterator for resonance */
    size_t it2; /*!< iterator for resonance */
    std::vector<double> centrals; /*!< [MeV] central value of resonance */
    std::vector<double> widths; /*!< [MeV] width of broad resonance*/
    /*! integrandum 
    * \param [out] res results
    * \param pnorm first integration variable
    * \param costheta second integration variable
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    * \param [in] centrals[MeV] central value of broad resonance
    * \param [in] widths [MeV] width of broad resonance, Gaussian shape assumed
    */
    void (*f)(numint::vector_d & res, double mass, double pnorm, double costheta, double qt, double qphi, InclusiveCross& cross, 
	      double Q2, double x, size_t it, size_t it2,  std::vector<double> &centrals, std::vector<double> &widths);
  };

   /*! integrandum function for on-shell contribution to the FSI amplitude with a broad resonance,
    * Gaussian shape for the distribution
    * \param [out] results results
    * \param mass [MeV] pole value for resonance
    * \param pnorm [MeV] norm of spectator momentum
    * \param costheta [] polar cos(theta) of spectator momentum
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    * \param [in] centrals [MeV] central value of broad resonance
    * \param [in] widths [MeV] width of broad resonance, Gaussian shape assumed
    */
  static void FSI_int_Gauss(numint::vector_d & results, double mass, double pnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, size_t it, size_t it2,
		      std::vector<double> &centrals, std::vector<double> &widths);
   /*! integrandum function for on-shell contribution to the FSI amplitude with a broad resonance,
    * Gaussian shape for the distribution
    * \param [out] results results
    * \param mass [MeV] pole value for resonance
    * \param prt [MeV] norm of transverse spectator momentum (initial)
    * \param W [MeV] invariant mass where the structure function is evaluated
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    * \param [in] centrals [MeV] central value of broad resonance
    * \param [in] widths [MeV] width of broad resonance, Gaussian shape assumed
    */
  static void FSI_int_Gauss_off(numint::vector_d & results, double mass, double prt, double W, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, size_t it, size_t it2,
		      std::vector<double> &centrals, std::vector<double> &widths);
 
  
   /*! integrandu0m function for on-shell contribution to the FSI amplitude with a broad resonance,
    * uniform shape for the distribution
    * \param [out] results results
    * \param mass [MeV] pole value for resonance
    * \param pnorm [MeV] norm of spectator momentum
    * \param costheta [] polar cos(theta) of spectator momentum
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    * \param [in] centrals [MeV] central value of broad resonance
    * \param [in] widths [MeV] width of broad resonance, Gaussian shape assumed
    */
  static void FSI_int_uniform(numint::vector_d & results, double mass, double pnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, size_t it, size_t it2,
		      std::vector<double> &centrals, std::vector<double> &widths);
   /*! integrandum function for on-shell contribution to the FSI amplitude with a broad resonance,
    * uniform shape for the distribution
    * \param [out] results results
    * \param mass [MeV] pole value for resonance
    * \param prt [MeV] norm of transverse spectator momentum (initial)
    * \param W [MeV] invariant mass where the structure function is evaluated
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    * \param [in] centrals [MeV] central value of broad resonance
    * \param [in] widths [MeV] width of broad resonance, Gaussian shape assumed
    */
  static void FSI_int_uniform_off(numint::vector_d & results, double mass, double prt, double W, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, size_t it, size_t it2,
		      std::vector<double> &centrals, std::vector<double> &widths);
 
   /*! integrandu0m function for on-shell contribution to the FSI amplitude with a broad resonance,
    * uniform shape for the distribution
    * \param [out] results results
    * \param mass [MeV] pole value for resonance
    * \param pnorm [MeV] norm of spectator momentum
    * \param costheta [] polar cos(theta) of spectator momentum
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    * \param [in] centrals [MeV] central value of broad resonance
    * \param [in] widths [MeV] width of broad resonance, Gaussian shape assumed
    */
  static void AzzFSI_int_uniform(numint::vector_d & results, double mass, double pnorm, double costheta, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, size_t it, size_t it2,
		      std::vector<double> &centrals, std::vector<double> &widths);
   /*! integrandum function for on-shell contribution to the FSI amplitude with a broad resonance,
    * uniform shape for the distribution
    * \param [out] results results
    * \param mass [MeV] pole value for resonance
    * \param prt [MeV] norm of transverse spectator momentum (initial)
    * \param W [MeV] invariant mass where the structure function is evaluated
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  InclusiveCross object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param it iterator for initial resonance
    * \param it2 iterator for final resonance (taken equal to it in our approach, diagonal)
    * \param [in] centrals [MeV] central value of broad resonance
    * \param [in] widths [MeV] width of broad resonance, Gaussian shape assumed
    */
  static void AzzFSI_int_uniform_off(numint::vector_d & results, double mass, double prt, double W, double qt, 
		      double qphi, InclusiveCross& cross, double Q2, double x, size_t it, size_t it2,
		      std::vector<double> &centrals, std::vector<double> &widths);
  
};


/*! @} */
#endif




