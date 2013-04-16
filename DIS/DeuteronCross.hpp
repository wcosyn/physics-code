/*! \file DeuteronCross.hpp 
 * \brief Contains declaration of class DeuteronCross, computes deuteron exclusive cross sections
 * \author Wim Cosyn
 * \date 29/08/2012
 * 
 * \addtogroup DIS
 * @{
 */

#ifndef DEUTERONCROSS_HPP
#define DEUTERONCROSS_HPP

#include "DeuteronMomDistr.hpp"
#include <TKinematics2to2.h>
#include <DeuteronStructure.hpp>
#include <TElectronKinematics.h>

#include <string>

/*! \brief A class that computes semi-exclusive deuteron cross sections */
class DeuteronCross{
  
public:
    /*! Constructor
   * \param wfname Deuteron wave function name, see TDeuteron
   * \param proton photon interacts with proton [1] or neutron [0]
   * \param strucname which structure function parametrization ["CB"=Chrisy-Bosted, "SLAC", "Alekhin"]
   * \param sigmain [mb] total rescattering cross section in FSI
   * \param betain [GeV^-2] slope parameter
   * \param epsilonin real part of scattering amplitude
   * \param betaoffin [GeV^-2] off-shell beta parameter
   * \param lambdain [GeV^2] lambda cutoff off-shell parameter
   * \param offshellset  which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: no off-shell amplitude, fully suppressed
   * - 4: full off-shell amplitude, no suppression 
   */
  DeuteronCross(std::string wfname, bool proton, std::string strucname,
    double sigmain, double betain, double epsilonin, double betaoffin, double lambdain, int offshellset);
  ~DeuteronCross(); /*!<Destructor */
  /*! get the average cross section
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param electron has electron kinematics
   * \param pw plane-wave calculation [1] or not [0]
   * \param Einoff off-shell energy of the nucleon interacting with the photon
   * \return [nb/GeV^4] semi-inclusive cross section \f$ \frac{d\sigma}{dxdQ^2\frac{d^3p_s}{E_s}} \f$ integrated over phi
   */
  double getavgCross(TKinematics2to2 &kin, TElectronKinematics &electron, bool pw, double Einoff);
  /*! computes the results like they are presented in the Deeps data
   * \param Q2 [MeV^2] four-momentum transfer
   * \param W [MeV] invariant mass of X
   * \param Ein [MeV] beam energy
   * \param pr [MeV] spectator momentum
   * \param costhetar angle of spectator with q
   * \param proton DIS on proton (1) or neutron (0)
   * \param phi angle between hadron and electron plane
   * \param homedir SHARE dir with input
   * \param planewave plane wave result
   * \param fsi fsi result
   */
  void getDeepsresult(double Q2, double W, double Ein, double pr, double costhetar, bool proton, double phi, std::string homedir,
    double &planewave, double &fsi);
  /*! reads in the Deeps data set in an array >
   * \param pdeepsarray pointer to array with the deeps results, indices as follows <BR>
   * Q2 [2 values] <BR>
   * invariant mass W [5 values] <BR>
   * spectator momentum [5 values] <BR>
   * cos spectator with q [34 values] <BR>
   * cos angle, deeps value, stat error, sys error
   * \param dir SHAREDIR
   */
  static void readin_deeps(double ******pdeepsarray, std::string dir);
  /*! cleans up memory of array with deeps results
   * \param deepsarray pointer to array with deeps results
   */
  static void maint_deepsarray(double *****deepsarray);
private:
  double massi; /*!< mass of nucleon interacting with photon */
//   TElectronKinematics electron; /*!< electron kinematics */
  DeuteronMomDistr momdistr; /*!< object instance used to calculate the deuteron momentum distributions */
  DeuteronStructure structure; /*!< object used to calculate the structure functions needed in semi-inclusive deuteron scattering */
};
/*! @} */
#endif