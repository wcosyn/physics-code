/*! \file Cross.hpp 
 * \brief Contains declaration of class , used to compute D(e,e') inclusive cross sections (QE part)
 * \author Wim Cosyn
 * \date 12/01/2014
 * 
 * \addtogroup Knockout
 * @{
 */
#ifndef DQEINCLUSIVE_HPP
#define DQEINCLUSIVE_HPP

#include <string>
#include <TElectronKinematics.h>
#include <NucleonEMOperator.hpp>
#include <TDeuteron.h>
#include <numint/numint.hpp>
#include <GammaStructure.h> 
#include <vector>


/*! \brief A class that computes inclusive deuteron cross sections */
class DQEinclusive{
  
public:
    /*! Constructor
   * \param elec Electron kinematics, see TElectronKinematics
   * \param proton photon interacting with proton [1] or neutron [0]
   * \param wavename Deuteron wave function name, see TDeuteron
   * \param FFparam which  form factor parametrization [see NucleonEMOperator]
   * \param offshell which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: suppression with mass difference between resonance mass and "true" W (used in paper)
   * - 4: full off-shell amplitude, no suppression
   * \param betaoffin [GeV^-2] eikonal parameter, off-shell slope parameter of particle in final-state with spectator
   * \param lambdain [GeV^2] cutoff parameter for off-shell part
   */
  DQEinclusive(bool proton, int FFparam, std::string wavename, TElectronKinematics &elec,
		 int offshell, double betaoffin=8., double lambdain=1.2);
  ~DQEinclusive(); /*!< Destructor */
  /*! Calculates the plane-wave inclusive cross on-shell section without any prefactors, just the deuteron incl F_L structure function
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param contrib1 [] direct born part
   * \param contrib2 [] crossed born part
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   */
  void calc_F2Dinc(double &contrib1, double &contrib2, double Q2,double x, int current);
  /*! Calculates the fsi inclusive cross section without any prefactors, just the deuteron incl F_L structure function
   * \param fsi1 [] inclusive fsi cross section without prefactors, direct term 
   * \param fsi2 [] crossed term
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   */
  void calc_F2DincFSI(double &fsi1, double &fsi2, double Q2,double x, int current);
  /*! Calculates the fsi inclusive off-shell cross section without any prefactors, just the deuteron incl F_L structure function
   * \param fsi1 [] inclusive fsi cross section without prefactors, direct term
   * \param fsi2 [] other fsi formula, crossed term
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   */
  void calc_F2DincFSI_off(double &fsi1, double &fsi2, double Q2,double x, int current);

  
  
  void setOffshell(const int offshell){offshellset=offshell;} /*!< set offshell parametrization */
  TDeuteron::Wavefunction* getDeutwf() const{return wf;} /*!< get instance to deuteron wave function */
  NucleonEMOperator* getFFactorseq() const{return ffactorseq;} /*!< get instance to form factors */
  NucleonEMOperator* getFFactorsdiff() const{return ffactorsdiff;} /*!< get instance to form factors */
  double getMassr() const{return massr;}
  double getMassi() const{return massi;}
  
  /*! method to find the pole in the fsi integration, longitudinal part,
   * \param sol [MeV] vector with solution (max 2) of on-shell ps,z values
  * \param pt [MeV] transverse spectator momentum 
  * \param W_sq [MeV^2] mass of resonance entering
   * \param nu [MeV] virtual photon energy
   * \param qvec [MeV] virtual photon momentum
   * \return [0] no solutions [1] solutions
   */
  bool get_prz(std::vector<double> &sol, double pt, double Q2, double nu, double qvec);
  
  
private:
  double betaoff; /*!< [MeV^-2] off-shell slope parameter, scattering parameter in FSI*/
  double lambda; /*!< [MeV^2 off-shell cutoff parameter */
  double massi; /*!< [MeV] mass of nucleon that gets hit by photon */
  double massr; /*!< [MeV] mass of spectator nucleon */
  bool proton; /*!< photon interaction with [0] neutron [1] proton */
  
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
  NucleonEMOperator *ffactorseq;  /*!< nucleon form factors, see NucleonEMOperator, nucleon is equal to "proton" from constructor */
  NucleonEMOperator *ffactorsdiff;  /*!< nucleon form factors, see NucleonEMOperator, nucleon is diff to "proton", needed in crossed diagrams */
  int ffparam; /*!< parametrization of formfactors, see NucleonEMOperator */
  
  /*! struct that is used for integrators plane wave ones*/
  struct Ftor_planewave {

    /*! integrandum function */
    static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
      Ftor_planewave &p = * (Ftor_planewave *) param;
      p.f(ret,x[0],*p.cross,p.Q2,p.x,p.current);
    }
    DQEinclusive *cross;/*!< pointer to  instance that contains all */
    double Q2; /*!< [MeV^2] momentum transfer */
    double x; /*!< [] Bjorken x */
    int current; /*!< CC1,CC2,CC3*/
    /*! integrandum 
    * \param res results
    * \param pperp [MeV] integration variable, transverse spectator momentum
    * \param cross instance of  where we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).  
    */
    void (*f)(numint::vector_d & res, double pperp, DQEinclusive& cross, double Q2, double x, int current);
  };
  
   /*! integrandum function
    * \param results results
    * \param pperp [MeV] integration variable, transverse spectator momentum
    * \param cross instance of DQEinclusive where we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
    */
  static void planewave_int(numint::vector_d & results, double pperp, 
			    DQEinclusive& cross, double Q2, double x, int current);
 
  
  /*! struct that is used for integrators fsi on- and off-shell ones*/
  struct Ftor_FSI {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor_FSI &p = * (Ftor_FSI *) param;
      p.f(ret,x[0],x[1],x[2],*p.cross,p.Q2,p.x,p.current);
    }
    DQEinclusive *cross;/*!< pointer to  instance that contains all */
    double Q2; /*!< [MeV^2] momentum transfer */
    double x; /*!< [] Bjorken x */
    int current; /*!< CC1,CC2,CC3*/
    /*! integrandum 
    * \param res results
    * \param pperp [MeV] integration variable, transverse spectator momentum
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  where we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
    */
    void (*f)(numint::vector_d & res, double pperp, double qt, double qphi, DQEinclusive& cross, 
	      double Q2, double x, int current);
  };
  
   /*! integrandum function for on-shell contribution to the FSI amplitude
    * \param results results
    * \param pperp [MeV] integration variable, transverse spectator momentum
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  DQEinclusive object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
    */
  static void FSI_int(numint::vector_d & results, double pperp, double qt, 
		      double qphi, DQEinclusive& cross, double Q2, double x, int current);
   /*! integrandum function for off-shell contribution to the FSI amplitude
    * \param results results
    * \param pperp [MeV] integration variable, transverse spectator momentum
    * \param qt [MeV] norm of transverse momentum transfer in FSI
    * \param qphi [] radial angle of transverse momentum transfer in FSI
    * \param cross instance of  DQEinclusive object we perform the integration on
    * \param Q2 [MeV^2] momentum transfer 
    * \param x [] Bjorken x
    * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
    */
  static void FSI_int_off(numint::vector_d & results, double pperp, double qt, 
		      double qphi, DQEinclusive& cross, double Q2, double x, int current);
 
  
};




/*! @} */
#endif

