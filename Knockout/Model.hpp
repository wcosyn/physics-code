/*! \file Model.hpp 
 * \brief Contains declaration of class Model, used to compute A(e,e'N) amplitudes
 * \author Wim Cosyn
 * \date 28/08/2012
 * 
 * \addtogroup Knockout
 * @{
 */

#ifndef MODEL_HPP
#define MODEL_HPP

#include <FourVector.h>
#include <MeanFieldNucleusThick.hpp>
#include <AbstractFsiCTGrid.hpp>
#include <TKinematics2to2.h>
#include <NucleonEMOperator.hpp>
#include <TVector3.h>
#include <GlauberGridThick.hpp>
#include <OneGlauberGrid.hpp>
#include <numint/numint.hpp>

#include <string>
#include <cstdarg>



/*! \brief A class Model, used to compute A(e,e'N) amplitudes */
class Model{
public:
    /*! Constructor
   * \param pnucleus pointer to a MF nucleus instance
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   * \param user_sigma does the user want to change sigma?
   * \param sigma_screening [%] screening change of sigma
   */
  Model(MeanFieldNucleusThick *pnucleus, double prec, int integrator, std::string dir, 
	bool user_sigma, double sigma_screening=0.);
  ~Model(); /*!< Destructor */
//   void setSRC(int setSRC) {SRC=setSRC;} /*!< sets SRC in the Glauber FSI */
//   /*! Computes an amplitude \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$  <BR>
//    * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
//    * \param tk contains the hadron kinematics
//    * \param spinout spin \f$ m_f \f$ of the outgoing nucleon times TWO!!! (-1 or +1)
//    * \param photonpol polarization of the photon (-1,0 or +1)
//    * \param shellindex which \f$ \alpha \f$ shell do we eject from
//    * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
//    * \param CT do you want CT effects in the Glauber FSI?
//    * \param pw no FSI?
//    * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
//    * \param SRC SRC or not in the FSI
//    * \param thick thickness in the glauber
//    * \return matrix element in [ fm \f$^{3/2}\f$]
//    */
//   std::complex<double> getMatrixEl(TKinematics2to2 &tk,int spinout, int photonpol, int shellindex, int m, int CT, 
// 				   int pw, int current, int SRC, int thick);
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all three photon polarizations (0,-1,+1) and both final nucleon helicities (-1,+1) <BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param results [ fm \f$^{3/2}\f$] contains amplitudes <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param CT do you want CT effects in the Glauber FSI?
   * \param pw no FSI?
   * \param SRC SRC or not in the FSI
   * \param thick thickness in the glauber
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   */
  void getMatrixEl(TKinematics2to2 &tk, Matrix<2,3> &results, int shellindex, int m, int CT, int pw, int current,
    int SRC, int thick);
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all three photon polarizations (0,-1,+1) and both final nucleon helicities (-1,+1) and all glauber varieties<BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param results [ fm \f$^{3/2}\f$] contains amplitudes (0=RMSGA,1=+SRC,2=+CT,3+CT+SRC,4=plane-wave) <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick thickness in the glauber
   */
  void getAllMatrixElMult(TKinematics2to2 &tk, Matrix<2,3> *results, int shellindex, int m, int current, int thick);
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all three photon polarizations (0,-1,+1) and both final nucleon helicities (-1,+1) and all glauber varieties<BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param results [ fm \f$^{3/2}\f$] contains amplitudes<BR>
   * (0=RMSGA,1=+SRC,2=+CT,3+CT+SRC,4=plane-wave) if thickness <BR>
   * (0=RMSGA,1=+SRC,2=plane-wave) if no thickness <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick thickness in the glauber
   */
  void getAllMatrixEl(TKinematics2to2 &tk, Matrix<2,3> *results, int shellindex, int m, int current, int thick);
   /*! Computes the off-shell amplitude \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu u(\vec{p}_m,m_i) \f$
   * \param tk contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param spinin spin \f$ m_i \f$ of the initial nucleon times TWO!!! (-1 or +1)
   * \param spinout spin \f$ m_f \f$ of the outgoing nucleon times TWO!!! (-1 or +1)
   * \param photonpol polarization of the photon (-1,0 or +1)
   * \return [dimensionless] off-shell amplitude
   */
  std::complex<double> getFreeMatrixEl(TKinematics2to2 &tk, int current, int spinin, int spinout, int photonpol);
  double getPrec() const{return prec;} /*!< returns precision in the integrations */
  bool getUsersigma() const{return usersigma;} /*!< returns 1 if user has changed sigma with some screening */
  double getSigmascreening() const{return sigmascreening;} /*!< [%] returns screening change to sigma */
  AbstractFsiCTGrid *getGrid() const{return grid;}
  GlauberGridThick* getThickGrid() {return &gridthick;}
  int getShell() const{return shell;}
  int getM() const{return mm;}
  int getIntegrator() const{return integrator;}
  MeanFieldNucleusThick* getPnucleus() {return pnucl;}
  Matrix<1,4> & getBarcontract0up() {return barcontract0up;}
  Matrix<1,4> & getBarcontractminup() {return barcontractminup;}
  Matrix<1,4> & getBarcontractplusup() {return barcontractplusup;}
  Matrix<1,4> & getBarcontract0down() {return barcontract0down;}
  Matrix<1,4> & getBarcontractmindown() {return barcontractmindown;}
  Matrix<1,4> & getBarcontractplusdown() {return barcontractplusdown;}
  Matrix<1,4> & getBarcontract() {return barcontract;}
  TVector3 & getPm() {return pm;}
  
private:
  double prec; /*!< precision in the integrations */
  int integrator; /*!< choice of integrator */
  AbstractFsiCTGrid *grid; /*!< pointer to Glauber Grid */
  MeanFieldNucleusThick *pnucl; /*!<  pointer to nucleus */
  NucleonEMOperator *J; /*!<  pointer to nucleon formfactor instance */
  Matrix<1,4> barcontract; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontract0up; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractminup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractplusup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontract0down; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractmindown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractplusdown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  TVector3 pm; /*!< missing momentum */
  std::string homedir; /*!< contains the share dir with all input */
  bool usersigma; /*!< has the user set some screening to sigma */
  double sigmascreening; /*!< [%] screening effect to sigma */
  GlauberGridThick gridthick;
  OneGlauberGrid onegrid;
  
  int shell; /*!< what shell has the ejected nucleon */
  int mm; /*!< \f$ m_j \f$ quantum number of ejected nucleon */
    /*! function that gets integrated over r, integral over r of matrix element
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJR(const double r, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), integral of matrix element
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJCosTheta(const double costheta, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, integral of matrix element
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJPhi(const double phi, std::complex<double> *results, va_list ap);
    /*! function that gets integrated over r, computes 12 different amplitudes:
     * 3 photon polarizations, 2 final nucleon helicities, including and excluding CT
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJR12(const double r, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), computes 12 different amplitudes:
     * 3 photon polarizations, 2 final nucleon helicities, including and excluding CT
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJCosTheta12(const double costheta, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over phi,computes 12 different amplitudes:
     * 3 photon polarizations, 2 final nucleon helicities, including and excluding CT
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJPhi12(const double phi, std::complex<double> *results, va_list ap);
  
  /*! struct that is used for integrators (clean ones)*/
  struct Ftor_one {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_one &p = * (Ftor_one *) param;
      p.f(ret,x[0],x[1],x[2],*p.model,p.SRC,p.pw);
    }
    Model *model;/*!< pointer to model instance that contains all */
    int SRC;/*!< SRC or not */
    int pw; /*!< plane-wave? */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param model Model instance
    * \param SRC SRC in FSI?
    * \param pw plane-wave?
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, Model & model, int SRC, int pw);
  };
  /*! struct that is used for integrators (radial ones)*/
  struct Ftor_r {

    /*! integrandum function */
    static void exec(const numint::array<double,1> &x, void *param, numint::vector_z &ret) {
      Ftor_r &p = * (Ftor_r *) param;
      p.f(ret,x[0],*p.model,p.thick,p.total);
    }
    Model *model;/*!< pointer to model instance that contains all */
    int thick;/*!< SRC or not */
    int total; /*!< plane-wave? */
    /*! integrandum 
    * \param res results
    * \param r first integration variable
    * \param model Model instance
    * \param thick thick in FSI?
    * \param total number of FSI situations?
    */
    void (*f)(numint::vector_z & res, double r, Model & model, int thick, int total);
  };
  /*! struct that is used for integrators (angular ones)*/
  struct Ftor_angles {

    /*! integrandum function */
    static void exec(const numint::array<double,2> &x, void *param, numint::vector_z &ret) {
      Ftor_angles &p = * (Ftor_angles *) param;
      p.f(ret,x[0],x[1],*p.model,p.thick,p.total,p.r);
    }
    Model *model;/*!< pointer to model instance that contains all */
    int thick;/*!< SRC or not */
    int total; /*!< plane-wave? */
    double r; /*!<radial coordinate*/
    /*! integrandum 
    * \param res results
    * \param r first integration variable
    * \param model Model instance
    * \param thick thick in FSI?
    * \param total number of FSI situations?
    */
    void (*f)(numint::vector_z & res, double costheta, double phi, Model & model, int thick, int total, double r);
  };
  /*! integrandum function (clean ones)*/
  static void klaas_one_amp(numint::vector_z & results, double r, double costheta, double phi, Model & model, int SRC, int pw);
  /*! integrandum function (clean ones)*/
  static void klaas_all_amp(numint::vector_z & results, double r, double costheta, double phi, Model & model, int SRC, int pw);
  /*! integrandum function (clean ones)*/
  static void klaas_mult_amp(numint::vector_z & results, double r, double costheta, double phi, Model & model, int SRC, int pw);
  /*! integrandum function (radial ones)*/
  static void klaas_all_radial(numint::vector_z & results, double r, Model & model, int SRC, int pw);
  /*! integrandum function (angular ones)*/
  static void klaas_all_angular(numint::vector_z & results, double costheta, double phi, Model & model, int SRC, int pw, double r);
  

};
/** @} */
#endif