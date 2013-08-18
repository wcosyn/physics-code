/*! \defgroup Knockout libKnockout: Library that contains all classes related to knockout cross sections, amplitudes and what not
 * \author Wim Cosyn
 * \date 28/08/2012
 * \brief This code implements classes to compute cross sections and amplitudes for knockout reactions
 * 
 * \details 
 * 
 * - It contains a class to compute plane-wave and distorted momentum distributions for the deuteron <BR>
 * [DeuteronMomDistr] <BR>
 * - It has a class to compute amplitudes and cross sections for A(e,e'p) reactions <BR>
 * [Model, Cross] <BR>
 * - It has a class to computed amplitudes and cross sections for A(e,e'NN) reactions (WORK IN PROGRESS!) <BR>
 * [DoubleNModel, DoubleNCross] <BR>
 * - It has classes to compute rho production on a deuteron and a general nucleus <BR>
 * [RhoDeuteron, RhoTCross] <BR>
 * 
 */



/*! \file Cross.hpp 
 * \brief Contains declaration of class Cross, used to compute A(e,e'N) cross sections
 * \author Wim Cosyn
 * \date 28/08/2012
 * 
 * \addtogroup Knockout
 * @{
 */
#ifndef CROSS_HPP
#define CROSS_HPP

#include <vector>
#include <TElectronKinematics.h>
#include <MeanFieldNucleusThick.hpp>
#include <TKinematics2to2.h>
#include "Model.hpp"

class AbstractFsiGrid;

/*! \brief A class Cross, used to compute A(e,e'N) cross sections */
class Cross{
public:
    /*! Constructor
   * \param elec contains all the electron kinematics
   * \param pnucl pointer to a MF nucleus instance
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param dir std::string that contains dir with all input, should be the ./share subdir of the project!
   * \param user_sigma does the user want to introduce a (anti)screening on sigma
   * \param sigmascreening [%] how much do you want to change the sigma value
   */
  Cross(TElectronKinematics &elec, MeanFieldNucleusThick *pnucl, 
	double prec, int integrator, std::string dir, bool user_sigma, double sigmascreening=0.);
  ~Cross(); /*!< Destructor */
  /*! Computes the differential A(e,e'N) cross section for certain kinematics and a certain shell of the nucleus
   * \param kin contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param SRC do you want to include SRC in the FSI?
   * \param CT do you want to include CT effects?
   * \param pw do you want to compute a plane-wave cross section (nullifies the thick,SRC,CT parameters)
   * \param shellindex selects the shell in the nucleus where the ejected N originates from
   * \param phi angle between electron and hadron plane
   * \param maxEval max # of evaluations in integrations
   * \param lab lab frame of cm frame for hadron part
   * \return differential cross section [fm^2 /MeV/sr^2]
   */
  double getDiffCross(TKinematics2to2 &kin, int current, int thick, int SRC, int CT, int pw, int shellindex, 
		      double phi, int maxEval, bool lab);
  /*! Computes the differential A(e,e'N) cross section for certain kinematics and a certain shell of the nucleus
   * \param cross vector with the different cross sections <BR>
   * differential cross section [fm ^2/MeV/sr^2]
   *  [0]: RMSGA <BR>
   *  [1]: RMSGA+SRC <BR>
   *  [2]: RMSGA+CT <BR>
   *  [3]: RMSGA+SRC+CT <BR>
   *  [4]: plane-wave<BR>
   * if no thickness [1] and [3] are not present (and vector has size 3, change indices accordingly)
   * \param kin contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param shellindex selects the shell in the nucleus where the ejected N originates from
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param phi angle between electron and hadron plane
   * \param maxEval max # of evaluations in integrations
   * \param lab lab frame of cm frame for hadron part
   */
  void getAllDiffCross(std::vector<double> &cross, TKinematics2to2 &kin, int current, 
		       int shellindex, int thick, double phi, int maxEval, bool lab);

  /*! Computes observables for the A(e,e'N) reaction for certain kinematics and a certain shell of the nucleus.
   * polarization axes defined as: z along momentum, y perp to the hadron plane, x perp to z in the hadron plane
   * \param obs (return object) vector with the different differential cross sections [fm^2/MeV/sr^2] <BR>
   *  [0-7]: RMSGA <BR>
   *  [8-15]: RMSGA+SRC <BR>
   *  [16-23]: RMSGA+CT <BR>
   *  [24-31]: RMSGA+SRC+CT <BR>
   *  [32-39]: plane-wave [0=sigma (fm^2/MeV/sr^2), 1=A, 2=Px, 3=Py, 4=Pz, 5=P'x, 6=P'y, 7=P'z]<BR>
   * if no thickness [8-15] and [24-31] are not present (and vector has size 24, change indices accordingly)
   * \param kin contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param shellindex selects the shell in the nucleus where the ejected N originates from
   * \param medium medium modifications in EMcoupling? [0=none, 1=CQM, 2=QSM]
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param phi angle between electron and hadron plane
   * \param maxEval max # of evaluations in integrations
   * \param lab lab frame or cm frame for hadron part
   * \return differential cross section 
   */
  void getAllObs(std::vector<double> &obs, TKinematics2to2 &kin, int current, 
			     int shellindex, int thick, int medium, double phi, int maxEval, bool lab);

  
  /*! Computes the off-shell (e,e'p) cross section for certain kinematics and a certain shell of the nucleus <BR>
   * What is denoted as \f$ K \sigma_{ep} \f$ in the literature
   * \param kin contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param phi angle between electron and hadron plane
   * \return differential cross section [dimensionless]
   */  
  double getElCross(TKinematics2to2 &kin, int current, double phi);
  /*! Computes \f$\delta(r) [fm^2]\f$ like defined in our density papers
   * \param densr different densities [fm^-1] <BR>
   *  [0]: RMSGA <BR>
   *  [1]: RMSGA+SRC <BR>
   *  [2]: RMSGA+CT <BR>
   *  [3]: RMSGA+SRC+CT <BR>
   *  [4]: plane-wave <BR>
   * if no thickness [1] and [3] are not present (and vector has size 3, change indices accordingly)
   * \param kin contains the hadron kinematics
   * \param shellindex selects the shell in the nucleus where the ejected N originates from
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param r [fm] argument of \f$\delta(r)\f$ 
   * \param maxEval max # of evaluations in integrations
   */
  void getDensr(std::vector<double> &densr, const TKinematics2to2 &kin, const int shellindex, 
		const int thick, const double r, const int maxEval);
  /*! Computes \f$\delta(r,\cos{\theta}) [fm^2]\f$ like defined in our density papers
   * \param densr different densities [fm^-1] <BR>
   *  [0]: RMSGA <BR>
   *  [1]: RMSGA+SRC <BR>
   *  [2]: RMSGA+CT <BR>
   *  [3]: RMSGA+SRC+CT <BR>
   *  [4]: plane-wave <BR>
   * if no thickness [1] and [3] are not present (and vector has size 3, change indices accordingly)
   * \param kin contains the hadron kinematics
   * \param shellindex selects the shell in the nucleus where the ejected N originates from
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param r [fm] argument of \f$\delta\f$ 
   * \param costheta argument of \f$\delta\f$ 
   * \param maxEval max # of evaluations in integrations
   */
  void getDensr_ctheta(std::vector<double> &densr, const TKinematics2to2 &kin, const int shellindex, 
		const int thick, const double r, const double costheta, const int maxEval);
  
  /*! Computes \f$\phi_D(\vec{p}_m) [fm^{3/2}]\f$ like defined in our density papers
   * \param phid different distorted fourier transforms [fm^3/2] <BR>
   *  [0]: RMSGA <BR>
   *  [1]: RMSGA+SRC <BR>
   *  [2]: RMSGA+CT <BR>
   *  [3]: RMSGA+SRC+CT <BR>
   *  [4]: plane-wave <BR>
   * if no thickness [1] and [3] are not present (and vector has size 3, change indices accordingly)
   * \param kin contains the hadron kinematics
   * \param shellindex selects the shell in the nucleus where the ejected N originates from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param ms \f$ m_s \f$ times TWO!!! of the final nucleon
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param maxEval max # of evaluations in integrations
   */
  void getPhid(std::vector< std::complex<double> > &phid, const TKinematics2to2 &kin, const int shellindex, const int m,
		const int ms, const int thick, const int maxEval);
  
  double getPrec() const{return prec;} /*!< precision of the integrations */
  double getSigmascreening() const{return sigmascreening;} /*!< [%] screening of sigma */
  bool getUsersigma() const{return usersigma;} /*!< has the user set sigma? */
  MeanFieldNucleusThick * getPnucl() const{return pnucl;}  /*!< pointer to nucleus object */
  
  /*! Prints a density profile for certain kinematics.
   * Grid of 100 points for \f$\delta(r) [fm^{-1}]\f$ like defined in our density papers (getDensr() function) <BR>
   * prints out [missing momentum (MeV)] [r (fm)] [\f$\delta(r) RMSGA+SRC\f$ (fm^-1)] [\f$\delta(r) RPWIA\f$ (fm^-1)]
   * [\f$ r^2\rho_A(r) \f$ [fm ^2] <BR>
   * At the end it prints the avg density \f$<\rho>\f$[fm^-3] for RMSGA+SRC & plane-wave, the avg radius \f$<r>\f$[fm] for RSMGA+SRC
   * and the distorted momentum distribution \f$\rho^D(\vec{p}_m)\f$ [fm^3] for RMSGA+SRC & plane-wave <BR>
   * \param kin contains the hadron kinematics
   * \param shellindex selects the shell in the nucleus where the ejected N originates from
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param maxEval max # of evaluations in integrations
   */
  void printDensity_profile(const TKinematics2to2 &kin, const int shellindex, 
		const int thick, const int maxEval);
  
private:
  std::string homedir; /*!< Contains dir with all input */
  double prec; /*!< precision you want in the integrations */
  double integrator; /*!< choice of integrator */
  TElectronKinematics electron; /*!< electron kinematics */
  MeanFieldNucleusThick *pnucl; /*!< pointer to nucleus */
  double kinfactors[6]; /*!< electron kinematic factors, get multiplied with response functions */
  double response[5][9]; /*!< response functions, from the hadronic tensor */
  double frontfactor; /*!< kinematical front factor */
  double mott; /*!< mott cross section */
  bool usersigma;
  double sigmascreening;
  Model *reacmodel; /*!< class object that computes the amplitudes */

  /*! struct that is used for integration of phi_d */
  struct Ftor_phid {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_phid &p = * (Ftor_phid *) param;
      p.f(ret,x[0],x[1],x[2],*p.cross,p.total,p.grid,p.pm,p.spinor,p.shell,p.m);
    }
    Cross *cross;/*!< pointer to cross instance that contains all */
    int total;/*!< number of results (5 with thicknes,3 without)*/
    AbstractFsiGrid *grid;
    TVector3 pm;
    Matrix<1,4> spinor;
    int shell;
    int m;
    /*! integrandum 
    * \param res results
    * \param r first integration variable
    * \param costheta second integration variable
    * \param phi third integration variable
    * \param model Cross instance
    * \param grid glauber grid
    * \param pm missing momentum vector
    * \param spinor outgoing nucleon spinor
    * \param shell shell index
    * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
    */
    void (*f)(numint::vector_z & res, double r, double costheta, double phi, Cross & cross, int total,
	      AbstractFsiGrid *grid, TVector3 &pm, Matrix<1,4> &spinor, int shell, int m);
  };
  static void klaas_phid(numint::vector_z & res, double r, double costheta, double phi, Cross & cross, int total,
	      AbstractFsiGrid *grid, TVector3 &pm, Matrix<1,4> &spinor, int shell, int m);
  
  
  /*! struct that is used for integration of phi_d */
  struct Ftor_densr {

    /*! integrandum function */
    static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
      Ftor_densr &p = * (Ftor_densr *) param;
      p.f(ret,x[0],x[1],*p.cross,p.total,p.grid,p.pm,p.spinordown,p.spinorup,p.phid,p.r,p.shell);
    }
    Cross *cross;/*!< pointer to cross instance that contains all */
    int total;/*!< number of results (5 with thicknes,3 without)*/
    AbstractFsiGrid *grid;
    TVector3 pm;
    Matrix<1,4> spinordown;
    Matrix<1,4> spinorup;
    std::complex<double>  **phid;
    int shell;
    double r;
    /*! integrandum 
    * \param res results
    * \param r first integration variable
    * \param costheta second integration variable
    * \param phi third integration variable
    * \param model Cross instance
    * \param grid glauber grid
    * \param pm missing momentum vector
    * \param spinorup outgoing nucleon spinor
    * \param spinordown outgoing nucleon spinor
    * \param phid distorted fourier transform
    * \param r [fm] r coordinate
    * \param shell shellindex
    */
    void (*f)(numint::vector_d & res, double costheta, double phi, Cross & cross, int total,
	      AbstractFsiGrid *grid, TVector3 &pm, Matrix<1,4> &spinordown, Matrix<1,4> &spinorup,
	      std::complex<double>  **phid, double r, int shell);
  };
  static void klaas_densr(numint::vector_d & res, double costheta, double phi, Cross & cross, int total,
	      AbstractFsiGrid *grid, TVector3 &pm, Matrix<1,4> &spinordown, Matrix<1,4> &spinorup,
	      std::complex<double>  **phid, double r, int shell);
 /*! struct that is used for integration of phi_d */
  struct Ftor_densr_ctheta {

    /*! integrandum function */
    static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
      Ftor_densr_ctheta &p = * (Ftor_densr_ctheta *) param;
      p.f(ret,x[0],*p.cross,p.total,p.grid,p.pm,p.spinordown,p.spinorup,p.phid,p.r,p.ctheta, p.shell);
    }
    Cross *cross;/*!< pointer to cross instance that contains all */
    int total;/*!< number of results (5 with thicknes,3 without)*/
    AbstractFsiGrid *grid;
    TVector3 pm;
    Matrix<1,4> spinordown;
    Matrix<1,4> spinorup;
    std::complex<double>  **phid;
    int shell;
    double r;
    double ctheta;
    /*! integrandum 
    * \param res results
    * \param r first integration variable
    * \param costheta second integration variable
    * \param phi third integration variable
    * \param model Cross instance
    * \param grid glauber grid
    * \param pm missing momentum vector
    * \param spinorup outgoing nucleon spinor
    * \param spinordown outgoing nucleon spinor
    * \param phid distorted fourier transform
    * \param r [fm] r coordinate
    * \param shell shellindex
    */
    void (*f)(numint::vector_d & res, double phi, Cross & cross, int total,
	      AbstractFsiGrid *grid, TVector3 &pm, Matrix<1,4> &spinordown, Matrix<1,4> &spinorup,
	      std::complex<double>  **phid, double r, double costheta, int shell);
  };
  static void klaas_densr_ctheta(numint::vector_d & res, double phi, Cross & cross, int total,
	      AbstractFsiGrid *grid, TVector3 &pm, Matrix<1,4> &spinordown, Matrix<1,4> &spinorup,
	      std::complex<double>  **phid, double r, double costheta, int shell);
};
/** @} */  
#endif