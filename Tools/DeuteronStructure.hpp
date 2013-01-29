/*! \file NuclStructure.hpp 
 * \brief Contains declaration of class NuclStructure, contains functions to compute deuteron structure functions
 * \author Wim Cosyn
 * \date 29/08/2012
 * 
 * \addtogroup MePhys
 * @{
 */
#ifndef DEUTERONSTRUCTURE_HPP
#define DEUTERONSTRUCTURE_HPP

#include "TKinematics2to2.h"
#include "TElectronKinematics.h"
#include <string>

/*! \brief A class to compute a variety of deuteron structure functions */
class DeuteronStructure{

public:
  /*! constructor
   * \param el electron kinematics, see TElectronKinematics
   * \param proton photon interacts with proton [1] or neutron [0]
   * \param name string for the structure function parametrization. Possibilities: <BR>
   * "CB": Christy & Bosted parametrization (see F1F209.f file) <BR>
   * "SLAC": SLAC paramtetrization from Bodek <BR>
   * "Alekhin": leading twist parametrization by Alekhin [see PRD 68,014002], also see alekhin.f file <BR>
   */
  DeuteronStructure(TElectronKinematics &el, int proton, std::string name);
  /*! return all four structure functions in semi-inclusive kinematics
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param FL FL structure function
   * \param FT FT structure function
   * \param FTT FTT structure function
   * \param FTL FTL structure function
   * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  void getStructureFunctions(TKinematics2to2 &kin, double &FL, double &FT, double &FTT, double &FTL, double Einoff) const;
  /*! gives you the combination of all structure functions
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param phi [rad] angle between electron scattering plane and hadron reaction plane
   * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  double getStructure(TKinematics2to2 &kin, double phi, double Einoff) const;
  /*! combination of all structure functions averaged over phi (angle between hadron and electron plane)
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
  * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  double getavgStructure(TKinematics2to2 &kin, double Einoff) const;
  /*! Deuteron F2 structure functions, used in the inclusive calculations, defined as
   * \f$ F_2^D = F_L^D + \frac{Q^2}{2|q|^2}\frac{\nu}{m_N}F_T^D
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
  * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  double getInclStructure(TKinematics2to2 &kin, double Einoff) const;
  
private:
  TElectronKinematics electron; /*!< electron kinematics */
  int proton; /*!< photon interacts with proton [1] or neutron [0] */
  std::string name; /*!< string that has the structure function parametrization name */
  double massi; /*!< mass of nucleon interacting with photon */
};
/** @} */
#endif