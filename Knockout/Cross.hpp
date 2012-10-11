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
  /*! Computes the differential \f$ A(e,e'N)\f$ cross section for certain kinematics and a certain shell of the nucleus
   * \param kin contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param SRC do you want to include SRC in the FSI?
   * \param CT do you want to include CT effects?
   * \param pw do you want to compute a plane-wave cross section (nullifies the thick,SRC,CT parameters)
   * \param shellindex selects the shell in the nucleus where the ejected N originates from
   * \param phi angle between electron and hadron plane
   * \return differential cross section [fm \f$ ^2 \f$/MeV/sr \f$ ^2 \f$]
   */
  double getDiffCross(TKinematics2to2 &kin, int current, int thick, int SRC, int CT, int pw, int shellindex, double phi);
  /*! Computes the differential \f$ A(e,e'N)\f$ cross section for certain kinematics and a certain shell of the nucleus
   * \param cross vector with the different cross sections <BR>
   *  [0]: plane-wave<BR>
   *  [1]: RMSGA <BR>
   *  [2]: RMSGA+SRC <BR>
   *  [3]: RMSGA+CT <BR>
   *  [4]: RMSGA+SRC+CT <BR>
   * \param kin contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param shellindex selects the shell in the nucleus where the ejected N originates from
   * \param thick do you want thickness in the Glauber FSI or not?
   * \param phi angle between electron and hadron plane
   * \param maxEval max # of evaluations in integrations
   * \param lab lab frame of cm frame for hadron part
   * \return differential cross section [fm \f$ ^2 \f$/MeV/sr \f$ ^2 \f$]
   */
  void getAllDiffCross(std::vector<double> &cross, TKinematics2to2 &kin, int current, 
		       int shellindex, int thick, double phi, int maxEval, bool lab);

  void getAllObs(std::vector<double> &obs, TKinematics2to2 &kin, int current, 
			     int shellindex, int thick, double phi, int maxEval, bool lab);

  
  /*! Computes the off-shell \f$ (e,e'p)\f$ cross section for certain kinematics and a certain shell of the nucleus <BR>
   * What is denoted as \f$ K \sigma_{ep} \f$ in the literature
   * \param kin contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param phi angle between electron and hadron plane
   * \return differential cross section [dimensionless]
   */  
  double getElCross(TKinematics2to2 &kin, int current, double phi);
  double getPrec() const{return prec;} /*!< precision of the integrations */
  double getSigmascreening() const{return sigmascreening;} /*!< [%] screening of sigma */
  bool getUsersigma() const{return usersigma;} /*!< has the user set sigma? */
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

};
/** @} */  
#endif