

/*! \file WeakQECross.hpp 
 * \brief Contains declaration of class WeakQECross, used to compute A(nu,lN) cross sections, both charged
 * and neutral weak currents
 * \author Wim Cosyn
 * \date 13/12/2013
 * 
 * \addtogroup Knockout
 * @{
 */
#ifndef WEAKQECROSS_HPP
#define WEAKQECROSS_HPP

#include <vector>
#include <TLeptonKinematics.h>
#include <TElectronKinematics.h>
#include <MeanFieldNucleusThick.hpp>
#include <TKinematics2to2.h>
#include "WeakQEHadronCurrent.hpp"

class AbstractFsiGrid;

/*! \brief A class WeakQECross, used to compute A(nu,lN) cross sections (CC and NC weak quasi-elastic knockout) */
class WeakQECross{
public:
    /*! Constructor for weak charged current processes
   * \param lepton contains all the lepton kinematics
   * \param pnucl pointer to a MF nucleus instance
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param dir std::string that contains dir with all input, should be the ./share subdir of the project!
  * \param charged neutral weak current [0] or charged weak current [1]
   * \param M_A_in [MeV] value of axial mass
   * \param user_sigma does the user want to change sigma?
   * \param gA_s [] strange contrib to G_A
   * \param r_s2 [fm^2] strange contrib to F1weak
   * \param mu_s [] strange contrib to F2weak
   * \param sigma_screening [%] how much do you want to change the sigma value
   */
  WeakQECross(TLeptonKinematics *lepton, MeanFieldNucleusThick *pnucl, 
	double prec, int integrator, std::string dir, bool charged, double M_A_in, bool user_sigma,
	double gA_s=-0.19, double r_s2=0., 
	double mu_s=0., double sigma_screening=0.);
    /*! Constructor for weak neutral current processes
   * \param elec contains all the electron kinematics
   * \param pnucl pointer to a MF nucleus instance
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param dir std::string that contains dir with all input, should be the ./share subdir of the project!
  * \param charged neutral weak current [0] or charged weak current [1]
   * \param M_A_in [MeV] value of axial mass
   * \param user_sigma does the user want to change sigma?
   * \param gA_s [] strange contrib to G_A
   * \param r_s2 [fm^2] strange contrib to F1weak
   * \param mu_s [] strange contrib to F2weak
   * \param sigma_screening [%] how much do you want to change the sigma value
   */
  WeakQECross(TElectronKinematics *elec, MeanFieldNucleusThick *pnucl, 
	double prec, int integrator, std::string dir, bool charged, double M_A_in, bool user_sigma,
	double gA_s=-0.19, double r_s2=0., 
	double mu_s=0., double sigma_screening=0.);
  ~WeakQECross(); /*!< Destructor */
  /*! Computes the differential A(nu,lN) cross section for certain kinematics and a certain shell of the nucleus,
   * automatic care is taken in the final state of isospin change for CC
   * \param kin contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param SRC do you want to include SRC in the FSI?
   * \param CT do you want to include CT effects?
   * \param pw do you want to compute a plane-wave cross section (nullifies the thick,SRC,CT parameters)
   * \param shellindex selects the shell in the nucleus where the ejected N originates from. <BR>
   * This also determines reaction type!! initial proton means antineutrino beam, initial neutron means neutrino beam.
   * \param phi angle between lepton and hadron plane
   * \param maxEval max # of evaluations in integrations
   * \param lab lab frame or cm frame for hadron part
   * \param phi_int integrated over phi (angle between lepton and hadron plane)?
   * \return [fm^2/MeV/sr^2] differential cross section 
   */
  double getDiffWeakQECross(TKinematics2to2 &kin, int current, int thick, int SRC, int CT, int pw, int shellindex, 
		      double phi, int maxEval, bool lab, bool phi_int);
  /*! Computes the differential A(nu,lN) cross section for certain kinematics and a certain shell of the nucleus
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
   * \param phi angle between lepton and hadron plane
   * \param maxEval max # of evaluations in integrations
   * \param lab lab frame or cm frame for hadron part
   * \param phi_int integrated over phi (angle between lepton and hadron plane)?
   * \param neutrino [1] neutrino or [0] antineutrino
   */
  void getAllDiffWeakQECross(std::vector<double> &cross, TKinematics2to2 &kin, int current, 
		       int shellindex, int thick, double phi, int maxEval, bool lab, bool phi_int, bool neutrino);

  /*! Computes the off-shell (nu,lN) cross section for certain kinematics <BR>
   * What is denoted as \f$ K \sigma_{ep} \f$ in the literature
   * \param kin contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param phi angle between lepton and hadron plane
   * \param maxEval max # of evaluations in integrations
   * \param phi_int integrated over phi (angle between lepton and hadron plane)?
   * \param neutrino [1] neutrino or [0] antineutrino
   * \param proton reaction on a proton [1] or neutron [0]
   * \return differential cross section [dimensionless]
   */  
  double getElWeakQECross(double Q2, double E_in, int current, double phi, int maxEval, bool phi_int,
			   bool neutrino, bool proton);
  
  double getPrec() const{return prec;} /*!< precision of the integrations */
  double getSigmascreening() const{return sigmascreening;} /*!< [%] screening of sigma */
  bool getUsersigma() const{return usersigma;} /*!< has the user set sigma? */
  MeanFieldNucleusThick * getPnucl() const{return pnucl;}  /*!< pointer to nucleus object */
  TLeptonKinematics * getPlepton() const{return lepton;} /*!< get pointer to lepton object */
  TElectronKinematics * getPelectron() const{return electron;} /*!< get pointer to lepton object */
  
  
private:
  std::string homedir; /*!< Contains dir with all input */
  double prec; /*!< precision you want in the integrations */
  double integrator; /*!< choice of integrator */
  TLeptonKinematics *lepton; /*!<  kinematics */
  TElectronKinematics *electron; /*!< NC lepton kinematics */
  MeanFieldNucleusThick *pnucl; /*!< pointer to nucleus */
//   double kinfactors[6]; /*!< electron kinematic factors, get multiplied with response functions */
//   double response[5][9]; /*!< response functions, from the hadronic tensor */
//   double frontfactor; /*!< kinematical front factor */
//   double mott; /*!< mott cross section */
  bool charged; /*!< neutral weak current [0] or charged weak current [1] */
  bool usersigma; /*!< userset value for sigma in the glauber parameters */
  double sigmascreening; /*!< screening percentage for sigma */
  double gA_s; /*!< [] strange contrib to G_A */
  double mu_s; /*!< [] strange contrib to F2_weak */
  double r_s2; /*!< [fm^2] strange contrib to F1_weak */
  double M_A; /*!< [MeV] axial mass value */

  WeakQEHadronCurrent *reacmodel; /*!< class object that computes the amplitudes */

};
/** @} */  
#endif