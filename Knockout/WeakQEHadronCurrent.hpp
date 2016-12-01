/*! \file WeakQEHadronCurrent.hpp 
 * \brief Contains declaration of class WeakQEHadronCurrent, used to compute A(nu,lN) amplitudes
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
#include <NucleonWeakOperator.hpp>
#include <TVector3.h>
#include <GlauberGridThick.hpp>
#include <OneGlauberGrid.hpp>
#include <numint/numint.hpp>
#include <TSpinor.h>

#include <string>
#include <cstdarg>



/*! \brief A class WeakQEHadronCurrent, used to compute A(nu,lN) amplitudes */
class WeakQEHadronCurrent{
public:
    /*! Constructor
   * \param pnucleus pointer to a MF nucleus instance
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   * \param max_Eval max # of function evaluations in the coordinate integration
   * \param charged neutral weak current [0] or charged weak current [1]
   * \param M_A_in [MeV] value of axial mass
   * \param user_sigma does the user want to change sigma?
   * \param gA_s [] strange contrib to G_A
   * \param r_s2 [fm^2] strange contrib to F1weak
   * \param mu_s [] strange contrib to F2weak
   * \param sigma_screening [%] screening change of sigma
   */
  WeakQEHadronCurrent(MeanFieldNucleusThick *pnucleus, double prec, int integrator, std::string dir, int max_Eval,
	bool charged, double M_A_in, bool user_sigma, double gA_s=-0.19, double r_s2=0., 
	double mu_s=0., double sigma_screening=0.);
  ~WeakQEHadronCurrent(); /*!< Destructor */
  
  
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all four vector boson polarizations (0,-1,+1,z) and both final nucleon helicities (-1,+1) <BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param [out] results [ fm^{3/2}] contains amplitudes <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1, 3 is z) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param CT do you want CT effects in the Glauber FSI?
   * \param pw no FSI?
   * \param SRC SRC or not in the FSI
   * \param thick thickness in the glauber
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param vector_or_axial [1] calc vector part of the current [0] calc axial vector part of the current
   */
  void getMatrixEl(TKinematics2to2 &tk, Matrix<2,4> &results, int shellindex, int m, int CT, int pw, int current,
    int SRC, int thick, bool vector_or_axial);
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all four vector boson polarizations (0,-1,+1,z) and both final nucleon helicities (-1,+1) and all glauber varieties<BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param [out] results [ fm^{3/2}] contains amplitudes (0=RMSGA,1=+SRC,2=+CT,3+CT+SRC,4=plane-wave) <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1, 3 is z) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick thickness in the glauber
   */
  void getAllMatrixElMult(TKinematics2to2 &tk, Matrix<2,4> *results, int shellindex, int m, int current, int thick);
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all four vector boson polarizations (0,-1,+1,z) and both final nucleon helicities (-1,+1) and all glauber varieties<BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param [out] results [ fm^{3/2} contains amplitudes<BR>
   * (0=RMSGA,1=+SRC,2=+CT,3+CT+SRC,4=plane-wave) if thickness <BR>
   * (0=RMSGA,1=+SRC,2=plane-wave) if no thickness <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1,3 is z) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick thickness in the glauber
   * \param medium medium modifications (0=none, 1=QCM, 2=CQSM)
   */
  void getAllMatrixEl(TKinematics2to2 &tk, Matrix<2,4> *results, int shellindex, int m, int current, int thick, int medium);
   /*! Computes the on-shell amplitude \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu u(0,m_i) \f$
   * \param Q2 [MeV^2] virtual photon Q-squared
   * \param proton [1] proton or [0] neutron
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param spinin spin \f$ m_i \f$ of the initial nucleon times TWO!!! (-1 or +1)
   * \param spinout spin \f$ m_f \f$ of the outgoing nucleon times TWO!!! (-1 or +1)
   * \param photonpol polarization of the photon (-1,0, +1, or z(3))
   * \return [dimensionless] on-shell amplitude
   */
  std::complex<double> getFreeMatrixEl(double Q2, bool proton, int current, int spinin, int spinout, int photonpol);
  double getPrec() const{return prec;} /*!< returns precision in the integrations */
  bool getUsersigma() const{return usersigma;} /*!< returns 1 if user has changed sigma with some screening */
  double getSigmascreening() const{return sigmascreening;} /*!< [%] returns screening change to sigma */
  AbstractFsiCTGrid *getGrid() const{return grid;}
  GlauberGridThick *getThickGrid() {return &gridthick;}
  int getShell() const{return shell;}
  int getM() const{return mm;}
  int getIntegrator() const{return integrator;}
  MeanFieldNucleusThick* getPnucleus() {return pnucl;}
  const Matrix<1,4> & getBarcontract0up() const{return barcontract0up;}
  const Matrix<1,4> & getBarcontractminup() const{return barcontractminup;}
  const Matrix<1,4> & getBarcontractplusup() const{return barcontractplusup;}
  const Matrix<1,4> & getBarcontractzup() const{return barcontractzup;}
  const Matrix<1,4> & getBarcontract0down() const{return barcontract0down;}
  const Matrix<1,4> & getBarcontractmindown() const{return barcontractmindown;}
  const Matrix<1,4> & getBarcontractplusdown() const{return barcontractplusdown;}
  const Matrix<1,4> & getBarcontractzdown() const{return barcontractzdown;}
  const Matrix<1,4> & getBarcontract() const{return barcontract;}
  const NucleonWeakOperator & getJ() const{return *J;}
  TVector3 & getPm() {return pm;}
  
private:
  double prec; /*!< precision in the integrations */
  int integrator; /*!< choice of integrator */
  AbstractFsiCTGrid *grid; /*!< pointer to Glauber Grid */
  MeanFieldNucleusThick *pnucl; /*!<  pointer to nucleus */
  NucleonWeakOperator *J; /*!<  pointer to nucleon formfactor instance */
  Matrix<1,4> barcontract; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontract0up; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractminup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractplusup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractzup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontract0down; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractmindown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractplusdown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractzdown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  TVector3 pm; /*!< missing momentum */
  std::string homedir; /*!< contains the share dir with all input */
  bool charged; /*!< neutral weak current [0] or charged [1] */
  bool usersigma; /*!< has the user set some screening to sigma */
  double sigmascreening; /*!< [%] screening effect to sigma */
  GlauberGridThick gridthick; /*!< grid for calc with thickness */
  OneGlauberGrid onegrid; /*!< grid for calc without thickness */
  int maxEval; /*!< max Evals in integrations */
  double gA_s; /*!< [] strange contrib to G_A */
  double mu_s; /*!< [] strange contrib to F2_weak */
  double r_s2; /*!< [fm^2] strange contrib to F1_weak */
  double M_A; /*!< [MeV] axial mass value */
  
  int shell; /*!< what shell has the ejected nucleon */
  int mm; /*!< \f$ m_j \f$ quantum number of ejected nucleon */

  
  /*! struct that is used for integrators (clean ones)*/
  struct Ftor_one {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_one &p = * (Ftor_one *) param;
      p.f(ret,x[0],x[1],x[2],*p.model,p.SRC,p.pw);
    }
    WeakQEHadronCurrent *model;/*!< pointer to model instance that contains all */
    int SRC;/*!< SRC or not */
    int pw; /*!< plane-wave? */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param model WeakQEHadronCurrent instance
    * \param SRC SRC in FSI?
    * \param pw plane-wave?
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, WeakQEHadronCurrent & model, int SRC, int pw);
  };

  /*! struct that is used for integrators (medium modif ones)*/
  struct Ftor_medium {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_medium &p = * (Ftor_medium *) param;
      p.f(ret,x[0],x[1],x[2],*p.model,p.thick,p.total,p.spinup, p.spindown, p.medium, p.q, p.pi,p.pf, 
	  p.polmin, p.polplus, p.pol0, p.polz, p.current);
    }
    WeakQEHadronCurrent *model;/*!< pointer to model instance that contains all */
    int thick;/*!< thickness or not */
    int total; /*!< number of different FSI? */
    Matrix<1,4> spinup;
    Matrix<1,4> spindown;
    int medium; /*!< which medium modif */
    FourVector<double> q;
    FourVector<double> pi;
    FourVector<double> pf;
    FourVector<std::complex<double> > polmin;
    FourVector<std::complex<double> > polplus;
    FourVector<std::complex<double> > pol0;
    FourVector<std::complex<double> > polz;
    int current; /*!<which current operator */
    
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param model WeakQEHadronCurrent instance
    * \param SRC SRC in FSI?
    * \param pw plane-wave?
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, WeakQEHadronCurrent & model, int thick, int total,
	      Matrix<1,4> &spinup, Matrix<1,4> &spindown, const int medium,
	      const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				   const FourVector<std::complex<double> > &polmin, const FourVector<std::complex<double> > &polplus, 
				   const FourVector<std::complex<double> > &pol0, const FourVector<std::complex<double> > &polz,
				   int current);
    
  };
  /*! integrandum function (clean ones)*/
  static void klaas_one_amp(numint::vector_z & results, double r, double costheta, double phi, WeakQEHadronCurrent & model, int SRC, int pw);
  /*! integrandum function (clean ones)*/
  static void klaas_all_amp(numint::vector_z & results, double r, double costheta, double phi, WeakQEHadronCurrent & model, int thick, int total);
  /*! integrandum function (medium modification)*/
  static void klaas_all_amp_medium(numint::vector_z & results, double r, double costheta, double phi, WeakQEHadronCurrent & model,
				   int thick, int total,Matrix<1,4> &spinup, Matrix<1,4> &spindown, const int medium,
				   const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				   const FourVector<std::complex<double> > &polmin, const FourVector<std::complex<double> > &polplus, 
				   const FourVector<std::complex<double> > &pol0,const FourVector<std::complex<double> > &polz,
				   int current);
  /*! integrandum function (clean ones)*/
  static void klaas_mult_amp(numint::vector_z & results, double r, double costheta, double phi, WeakQEHadronCurrent & model, int SRC, int pw);
   

};
/** @} */
#endif