/*! \file PionTCross.hpp 
 * \brief Contains declaration of class PionTCross, used to compute A(e,e'pi) transparencies, based on a factorized model (distorted momentum distributions)
 * \author Wim Cosyn
 * \date 02/04/2025
 * 
 * \addtogroup Knockout
 * @{
 */
#ifndef PIONTCROSS_HPP
#define PIONTCROSS_HPP


#include <MeanFieldNucleusThick.hpp>
// #include <TKinematics2to2.h>
#include <DistMomDistrGrid.hpp>
#include <GlauberGridThick.hpp>
#include <FastParticle.hpp>
#include <map>


/*! \brief A class used to compute A(e,e'rho) cross sections & transparencies*/
class PionTCross{
public:
  /*! Constructor
   * \param nucleus nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param pmax [MeV] maximum missing momentum in integrations 
   * \param dir share dir where all input is located
   * \param userset do we want to set sigma ourselves?
   * \param usersigma [mb] user chosen value of sigma_tot rho-nucleon scattering
   * \param precision precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param maxEval max # of function evaluations in the pm integration
   */
  PionTCross(const int nucleus, const double pmax, const std::string dir, 
	    const bool userset, const double usersigma, const double precision, const int integrator,
	    const int maxEval
 	  );
  ~PionTCross(); /*!<Destructor */
  /*! Calculate cross section integrated over z at fixed t, cross section is dsigma/dEe'dOmega_e'dtdphi
   * \param results [GeV^-1] all the different computed cross sections
   * \param Ebeam [GeV] beam energy
   * \param Q2 [GeV^2] Q^2
   * \param nu [GeV] virtual photon energy
   * \param t [GeV^2] momentum transfer squared
   */
  void getCross(double *results, const double Ebeam, const double Q2, const double nu, const double t);
  double getPrec() const{return prec;} /*!< returns the precision */
  int getNrofcross() const{return nrofcross;} /*!<returns number of different cross sections */
  const MeanFieldNucleusThick& getNucleusthick() const{return nucleusthick;}
  
private:
  std::string homedir; /*!< share dir where all input is located */
  double pmax; /*!< maximum p_miss in integrations*/
  bool userset; /*!< has sigma parameter been set by user?*/
  double usersigma; /*!< value of user set sigma_rho*/
  double prec; /*!< precision in the integrations */
  int integrator; /*!< choice of integrator */
  int nrofcross; /*!< number of different cross sections (number of FSI grids + 1 for plane-wave */
  double abserror; /*!< absolute error in interations */
  int maxEval; /*!< max # of function evaluations in integration*/
  MeanFieldNucleusThick nucleusthick; /*!< nucleus instance */
  DistMomDistrGrid **pdistgrid; /*!< array of distorted momentum distribution grid, one for each shell level */
  std::map<double,DistMomDistrGrid> distgridmap; /*!< map used to store distorted momentum distributions */
  GlauberGridThick **pfsigrid; /*!< array of Glauber FSI grid, one for each shell level */
//   TKinematics2to2 *pkin;
//  FastParticle *prho;
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
   * \param torz t [1] or z[0] integrations
   */
  double getfrontfactor(double nu, double qvec, double Erho, double prho, double pzrho, double pxrho,
				 double s, double Q2, double mN, double t, bool torz);
  /*! struct that is used for integrators (clean ones)*/
  struct Ftor_pion {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor_pion &p = * (Ftor_pion *) param;
      p.f(ret,x[0],x[1],x[2],*p.cross,p.Q2,p.nu,p.qvec,p.t,p.Erho,p.prho);
    }
    PionTCross *cross;/*!< pointer to the class object where the integration is performed */
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
    * \param cross the PionTCross instance
    */
    void (*f)(numint::vector_d & res, double pm, double costheta, double phi, PionTCross & cross, double Q2, double nu, double qvec,
      double t, double Erho, double prho);
  };
  /*! integrandum function (clean ones)*/
  static void klaas_rho_t(numint::vector_d & res, double pm, double costheta, double phi, PionTCross & cross, double Q2, double nu, double qvec,
      double t, double dummy, double dummy2);
  
  
  
  
};

/*! @} */
#endif