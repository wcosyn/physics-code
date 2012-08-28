/*! \file DoubleNCross.hpp 
 * \brief Contains declaration of class DoubleNCross, used to compute A(e,e'NN) cross sections
 * \author Wim Cosyn
 * \date 28/08/2012
 * 
 * \addtogroup Knockout
 * @{
 */
#ifndef DOUBLENCROSS_HPP
#define DOUBLENCROSS_HPP

#include <TElectronKinematics.h>
#include <MeanFieldNucleusThick.hpp>
#include <TKinematics2to3WithLabAngles.h>
#include "DoubleNModel.hpp"

/*! \brief A class used to compute A(e,e'NN) cross sections */
class DoubleNCross{
public:
    /*! Constructor
   * \param elec contains all the electron kinematics
   * \param pnucl pointer to a MF nucleus instance
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   */
  DoubleNCross(TElectronKinematics &elec, MeanFieldNucleusThick *pnucl, string dir);
  ~DoubleNCross();/*!< Destructor */
  /*! Computes the differential \f$ A(e,e'NN)\f$ cross section for certain kinematics and a certain shell of the nucleus
   * \param kin contains the hadron kinematics
   * \param SRC do you want to include SRC in the FSI?
   * \param CT do you want to include CT effects?
   * \param pw do you want to compute a plane-wave cross section (nullifies the thick,SRC,CT parameters)
   * \param corr use correlated wave functions (from Maarten)
   * \param shellindex1 selects the shell in the nucleus where the first ejected N originates from
   * \param shellindex2 selects the shell in the nucleus where the first ejected N originates from
   * \param phi angle between electron and hadron plane
   * \return differential cross section [fm^4 for now, am not including frontfactor]
   */
  double getDiffCross(const TKinematics2to3 &kin, bool SRC, bool CT, bool pw, bool corr, int shellindex1, int shellindex2,
		      double phi);
private:
  string homedir; /*!< points to the share directory with all input */
  TElectronKinematics electron; /*!< electron kinematics */
  MeanFieldNucleusThick *pnucl; /*!< pointer to instance of nucleus */
  double kinfactors[6]; /*!< electron kinematic factors [dimensionless] */
  double response[9]; /*!< response functions [] */
  double frontfactor; /*!< kinematical front factor */
  double mott; /*!< mott cross section [MeV^-2] */
  DoubleNModel *reacmodel; /*!< pointer to reaction model that computes amplitueds */
  
};

/*! @} */

#endif