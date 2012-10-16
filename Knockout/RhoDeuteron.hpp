/*! \file RhoDeuteron.hpp 
 * \brief Contains declaration of class RhoTDeuteron, used to compute D(e,e'rho) cross sections/transparencies
 * \author Wim Cosyn
 * \date 28/08/2012
 * 
 * \addtogroup Knockout
 * @{
 */
#ifndef RHODEUTERON_HPP
#define RHODEUTERON_HPP


#include <MeanFieldNucleusThick.hpp>
// #include <TKinematics2to2.h>
#include <DistMomDistrGrid.hpp>
#include <GlauberDecayGridThick.hpp>
#include <FastParticle.hpp>
#include <TDeuteron.h>
#include "DeuteronMomDistr.hpp"



/*! \brief A class used to compute D(e,e'rho) cross sections & transparencies*/
class RhoDeuteron{
public:
  /*! Constructor
   * \param wavefunction D wave function name, see TDeuteron
   * \param p_max [MeV] maximum missing momentum in integrations 
   * \param no_cuts apply the experimental cuts in t or z or not 
   * \param integr which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param precision precision you want in the integrations
   * \param max_Eval max # of function evaluations in the pm integration
   * \param fsi fsi or not
   */
  RhoDeuteron(const std::string wavefunction, const double p_max, const bool no_cuts, const int integr,
	      const double precision, const int max_Eval, const bool fsi);
  ~RhoDeuteron(); /*!< Destructor */
  /*! Calculate cross section integrated over z at fixed t
   * \param results [fm^3GeV^4] all the different computed cross sections
   * \param Ebeam [GeV] beam energy
   * \param Q2 [GeV^2] Q^2
   * \param nu [GeV] virtual photon energy
   * \param t [GeV^2] momentum transfer squared
   */
  void getCrosst(double *results, const double Ebeam, const double Q2, const double nu, const double t);
  /*! Calculate cross section integrated over t at fixed z
   * \param results [fm^3GeV^4] all the different computed cross sections
   * \param Ebeam [GeV] beam energy
   * \param Q2 [GeV^2] Q^2
   * \param nu [GeV] virtual photon energy
   * \param z [] E_rho/nu
   */
  void getCrossz(double *results, const double Ebeam, const double Q2, const double nu, const double z);
  /*! Calculate coherent cross section integrated over z at fixed t
   * \param results [fm^3GeV^4] all the different computed cross sections
   * \param Ebeam [GeV] beam energy
   * \param Q2 [GeV^2] Q^2
   * \param nu [GeV] virtual photon energy
   * \param t [GeV^2] momentum transfer squared
   */
  void getCrosst_coh(double *results, const double Ebeam, const double Q2, const double nu, const double t);
  /*! Calculate coherent cross section integrated over t at fixed z
   * \param results [fm^3GeV^4] all the different computed cross sections
   * \param Ebeam [GeV] beam energy
   * \param Q2 [GeV^2] Q^2
   * \param nu [GeV] virtual photon energy
   * \param z [] E_rho/nu
   */
  void getCrossz_coh(double *results, const double Ebeam, const double Q2, const double nu, const double z);
  int getNocuts() const{return nocuts;}
  
private:
  double pmax; /*!< maximum p_miss in integrations*/
  bool nocuts; /*!< experimental cuts in t or z applied? */
  DeuteronMomDistr deuteron; /*!< Momentum distribution instance of deuteron */
  double prec; /*!< precision in the integrations */
  int integrator; /*!< choice of integrator */
  double abserror; /*!< absolute error in interations */
  int maxEval; /*!< max # of function evaluations in integration*/
  bool fsi; /*!< include fsi or not */
  
  /*! Calculates coherent cross section
   * \param pDvec [MeV] final deuteron momentum vec
   * \return cross section []
   */
  double calcCross_coh(TVector3 &pDvec);
  //   TKinematics2to2 *pkin;
  /*! function that gets integrated over pm, all different fsi outputs, incoherent cross section fixed t
   * \param pm [MeV] radial coordinate
   * \param results  result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPmt(const double pm, double *results, va_list ap);
  /*! function that gets integrated over cos(theta), all different fsi outputs, incoherent cross section fixed t
   * \param costheta [] cos of theta coordinate
   * \param results  result: contains the cross section 
   * \param ap variable parameter list
   */
  void intCosThetat(const double costheta, double *results, va_list ap);
  /*! function that gets integrated over phi, all different fsi outputs, incoherent cross section fixed t
   * \param phi [rad] phi coordinate
   * \param results  result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPhit(const double phi, double *results, va_list ap);
  /*! function that gets integrated over pm, all different fsi outputs, incoherent cross section fixed z
   * \param pm [MeV] radial coordinate
   * \param results  result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPmz(const double pm, double *results, va_list ap);
  /*! function that gets integrated over cos(theta), all different fsi outputs, incoherent cross section fixed z
   * \param costheta [] cos of theta coordinate
   * \param results  result: contains the cross section 
   * \param ap variable parameter list
   */
  void intCosThetaz(const double costheta, double *results, va_list ap);
  /*! function that gets integrated over phi, all different fsi outputs, incoherent cross section fixed z
   * \param phi [rad] phi coordinate
   * \param results  result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPhiz(const double phi, double *results, va_list ap);
  
  /*! function that gets integrated over pm, all different fsi outputs, coherent cross section
   * \param pm [MeV] radial coordinate
   * \param results  result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPmcoh(const double pm, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), all different fsi outputs, coherent cross section
   * \param costheta [] cos of theta coordinate
   * \param results  result: contains the cross section 
   * \param ap variable parameter list
   */
  void intCosThetacoh(const double costheta, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, all different fsi outputs, coherent cross section
   * \param phi [rad] phi coordinate
   * \param results  result: contains the cross section 
   * \param ap variable parameter list
   */
  void intPhicoh(const double phi, std::complex<double> *results, va_list ap);
  /*!Calculates an array of momentum distributions for all different FSI situations
   * \param results [MeV^-3] distorted momentum distributions
   * \param Erho [MeV] rho energy
   * \param prho [MeV] rho momentum
   * \param pzrho [MeV] component of rho momentum along q
   * \param pxrho [MeV] transverse rho momentum
   * \param pmvec [MeV] missing momentum
   * \param nu [MeV] virtual photon energy
   * \param qvec [MeV] virtual photon momentum
  */
  void getMomdistr(double *results, double Erho, double prho, double pzrho, 
			      double pxrho, TVector3 &pmvec, double nu, double qvec);
  /*! Calculates a kinematic frontfactor that enters in the integrations
   * \param qvec [GeV] virtual photon momentum
   * \param Erho [GeV] rho energy
   * \param prho [GeV] rho momentum
   * \param pzrho [GeV] component of rho momentum along q
   * \param s [GeV^2] invariant mass squared
   * \param Q2 [GeV^2] Q^2
   * \param En [GeV] energy of recoil nucleon
   * \param t [GeV^2] momentum transfer squared
   * \param Eoff [GeV] initial nucleon off-shell momentum
   * \param torz t [1] or z[0] integrations
   */
   double getfrontfactor(double qvec, double Erho, double prho, double pzrho, double s, double Q2, double t, 
			 double En, double Eoff, bool torz);
  struct Ftor_rhoD {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor_rhoD &p = * (Ftor_rhoD *) param;
      p.f(ret,x[0],x[1],x[2],*p.deuteron,p.Q2,p.nu,p.qvec,p.t,p.Erho,p.prho);
    }
    RhoDeuteron *deuteron;/*!< pointer to rhodeuteron object */
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
    void (*f)(numint::vector_d & res, double pm, double costheta, double phi, RhoDeuteron & deuteron, double Q2, double nu, double qvec,
      double t, double Erho, double prho);
  };
  /*! integrandum function (clean ones)*/
  static void klaas_rhoD_t(numint::vector_d & res, double pm, double costheta, double phi, RhoDeuteron & deuteron, double Q2, double nu, double qvec,
      double t, double dummy, double dummy2);
  /*! integrandum function (clean ones), only CT*/
  static void klaas_rhoD_z(numint::vector_d & res, double pm, double costheta, double phi, RhoDeuteron & deuteron, double Q2, double nu, double qvec,
      double dummy, double Erho, double prho);

};
/*! @} */
#endif