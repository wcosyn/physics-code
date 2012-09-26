/*! \file RhoTCross.hpp 
 * \brief Contains declaration of class RhoTCross, used to compute A(e,e'rho) cross sections/transparencies
 * \author Wim Cosyn
 * \date 28/08/2012
 * 
 * \addtogroup Knockout
 * @{
 */
#ifndef RHOTCROSS_HPP
#define RHOTCROSS_HPP


#include <MeanFieldNucleusThick.hpp>
// #include <TKinematics2to2.h>
#include <DistMomDistrGrid.hpp>
#include <GlauberDecayGridThick.hpp>
#include <FastParticle.hpp>

// #define NROFRES 9 /*!< number of different cross sections calculated, 1 pw, 8 glaubers */


/*! \brief A class used to compute A(e,e'rho) cross sections & transparencies*/
class RhoTCross{
public:
  /*! Constructor
   * \param nucleus nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param pmax [MeV] maximum missing momentum in integrations 
   * \param dir share dir where all input is located
   * \param no_cuts apply the experimental cuts in t or z or not 
   * \param userset do we want to set sigma ourselves?
   * \param usersigma [mb] user chosen value of sigma_tot rho-nucleon scattering
   * \param precision precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param maxEval max # of function evaluations in the pm integration
   */
  RhoTCross(const int nucleus, const double pmax, const std::string dir, const bool no_cuts, 
	    const bool userset, const double usersigma, const double precision, const int integrator,
	    const int maxEval
 	  );
  ~RhoTCross(); /*!<Destructor */
  /*! Calculate cross section integrated over z at fixed t, cross section is dsigma/dEe'dOmega_e'dtdphi
   * \param results [GeV^-1] all the different computed cross sections
   * \param Ebeam [GeV] beam energy
   * \param Q2 [GeV^2] Q^2
   * \param nu [GeV] virtual photon energy
   * \param t [GeV^2] momentum transfer squared
   */
  void getCrosst(double *results, const double Ebeam, const double Q2, const double nu, const double t);
  /*! Calculate cross section integrated over t at fixed z, cross section is dsigma/dEe'dOmega_e'dzdphi
   * \param results [GeV] all the different computed cross sections
   * \param Ebeam [GeV] beam energy
   * \param Q2 [GeV^2] Q^2
   * \param nu [GeV] virtual photon energy
   * \param z [] E_rho/nu
   */
  void getCrossz(double *results, const double Ebeam, const double Q2, const double nu, const double z);
  /*! Calculate coherent cross section integrated over z at fixed t, not implemented!!!!
   * \param results [GeV] all the different computed cross sections
   * \param Ebeam [GeV] beam energy
   * \param Q2 [GeV^2] Q^2
   * \param nu [GeV] virtual photon energy
   * \param t [GeV^2] momentum transfer squared
   */
  void getCrosst_coh(double *results, const double Ebeam, const double Q2, const double nu, const double t);
  /*! Calculate coherent cross section integrated over t at fixed z, not implemented!!!
   * \param results [GeV] all the different computed cross sections
   * \param Ebeam [GeV] beam energy
   * \param Q2 [GeV^2] Q^2
   * \param nu [GeV] virtual photon energy
   * \param z [] E_rho/nu
   */
  void getCrossz_coh(double *results, const double Ebeam, const double Q2, const double nu, const double z);
  double getPrec() const{return prec;} /*!< returns the precision */
  int getNrofcross() const{return nrofcross;} /*!<returns number of different cross sections */
  MeanFieldNucleusThick* getNucleusthick() {return &nucleusthick;}
  int getNocuts() const{return nocuts;}
  
private:
  std::string homedir; /*!< share dir where all input is located */
  double pmax; /*!< maximum p_miss in integrations*/
  bool nocuts; /*!< experimental cuts in t or z applied? */
  bool userset; /*!< has sigma parameter been set by user?*/
  double usersigma; /*!< value of user set sigma_rho*/
  double prec; /*!< precision in the integrations */
  int integrator; /*!< choice of integrator */
  int nrofcross; /*!< number of different cross sections (number of FSI grids + 1 for plane-wave */
  double abserror; /*!< absolute error in interations */
  int maxEval; /*!< max # of function evaluations in integration*/
  MeanFieldNucleusThick nucleusthick; /*!< nucleus instance */
  DistMomDistrGrid **pdistgrid; /*!< array of distorted momentum distribution grid, one for each shell level */
  GlauberDecayGridThick **pfsigrid; /*!< array of Glauber FSI grid, one for each shell level */
//   TKinematics2to2 *pkin;
//  FastParticle *prho;
  /*! function that gets integrated over pm, all different fsi outputs / fixed t calculations
   * \param pm [MeV] radial coordinate
   * \param results result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPmt(const double pm, double *results, va_list ap);
  /*! function that gets integrated over cos(theta),  all different fsi outputs / fixed t calculations
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the cross section 
   * \param ap variable parameter list
   */
  void intCosThetat(const double costheta, double *results, va_list ap);
  /*! function that gets integrated over phi,  all different fsi outputs / fixed t calculations
   * \param phi [rad] phi coordinate 
   * \param results result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPhit(const double phi, double *results, va_list ap);
  /*! function that gets integrated over pm, all different fsi outputs / fixed z calculations
   * \param pm [MeV] radial coordinate
   * \param results result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPmz(const double pm, double *results, va_list ap);
  /*! function that gets integrated over cos(theta),  all different fsi outputs  / fixed z calculations
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the cross section 
   * \param ap variable parameter list
   */
  void intCosThetaz(const double costheta, double *results, va_list ap);
  /*! function that gets integrated over phi,  all different fsi outputs  / fixed z calculations
   * \param phi [rad] phi coordinate
   * \param results result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPhiz(const double phi, double *results, va_list ap);
  /*!Calculates an array of momentum distributions for all different FSI situations
   * \param results [fm^3] distorted momentum distributions
   * \param prho [MeV] rho momentum
   * \param thetarho [rad]  theta_rhoq
   * \param t [GeV^2] momentum transfer squared, hard_scale for CT
   * \param shell shell index
   * \param pm [MeV] missing momentum spherical coordinates
   * \param pmcostheta
   * \param pmphi
   */
  void getMomdistr(double *results, double prho, double thetarho, double t, int shell, 
			    double pm, double pmcostheta, double pmphi);
  /*! Calculates a kinematic frontfactor that enters in the integrations
   * \param nu [GeV] virtual photon energy
   * \param qvec [GeV] virtual photon momentum
   * \param Erho [GeV] rho energy
   * \param prho [GeV] rho momentum
   * \param pzrho [GeV] component of rho momentum along q
   * \param pxrho [GeV] transverse rho momentum
   * \param s [GeV^2] invariant mass squared
   * \param Q2 [GeV^2] Q^2
   * \param mN [GeV] mass of nucleon that interacts with photon
   * \param t [GeV^2] momentum transfer squared
   */
  double getfrontfactor(double nu, double qvec, double Erho, double prho, double pzrho, double pxrho,
				 double s, double Q2, double mN, double t, bool torz);
  /*! struct that is used for integrators (clean ones)*/
  struct Ftor_rho {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor_rho &p = * (Ftor_rho *) param;
      p.f(ret,x[0],x[1],x[2],*p.cross,p.Q2,p.nu,p.qvec,p.t,p.Erho,p.prho);
    }
    RhoTCross *cross;/*!< pointer to the grid where the integration is performed */
    double Q2;/*!< [GeV^2] Qsquared */
    double nu; /*!<[GeV] virtual photon energy */
    double qvec; /*!< [GeV] virtual photon momentum */
    double t; /*!< [GeV^2] momentum transfer sq */
    double Erho; /*!< [GeV] rho energy */
    double prho; /*!< [GeV] rho momentum */
    /*! integrandum 
    * \param res results
    * \param pm first integration variable
    * \param costheta second integration variable
    * \param phi third integration variable
    * \param cross the RhoTCross instance
    */
    void (*f)(numint::vector_d & res, double pm, double costheta, double phi, RhoTCross & cross, double Q2, double nu, double qvec,
      double t, double Erho, double prho);
  };
  /*! integrandum function (clean ones)*/
  static void klaas_rho_t(numint::vector_d & res, double pm, double costheta, double phi, RhoTCross & cross, double Q2, double nu, double qvec,
      double t, double dummy, double dummy2);
  /*! integrandum function (clean ones), only CT*/
  static void klaas_rho_z(numint::vector_d & res, double pm, double costheta, double phi, RhoTCross & cross, double Q2, double nu, double qvec,
      double t, double Erho, double prho);
  
  
  
  
  
};

/*! @} */
#endif