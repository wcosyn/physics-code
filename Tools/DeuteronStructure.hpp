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
#include "LightConeKin2to2.hpp"

/*! \brief A class to compute a variety of deuteron structure functions */
class DeuteronStructure{

public:
  /*! constructor
   * \param proton photon interacts with proton [1] or neutron [0]
   * \param name string for the structure function parametrization. Possibilities: <BR>
   * "CB": Christy & Bosted parametrization (see F1F209.f file) <BR>
   * "SLAC": SLAC paramtetrization from Bodek <BR>
   * "Alekhin": leading twist parametrization by Alekhin [see PRD 68,014002], also see alekhin.f file <BR>
   */
  DeuteronStructure(const bool proton, std::string name);
  
  DeuteronStructure(); /*!< Default Constructor */
  DeuteronStructure(const DeuteronStructure&); /*!< Copy Constructor */
  DeuteronStructure& operator=(const DeuteronStructure&); /*!< assignment operator */

  
  /*! return all four structure functions in semi-inclusive kinematics <BR>
   * For details see Phys.Rev. C84 (2011) 014601 Eqs. (24-30)
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param[out] FL FL structure function
   * \param[out] FT FT structure function
   * \param[out] FTT FTT structure function
   * \param[out] FTL FTL structure function
   * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  void getStructureFunctions(TKinematics2to2 &kin, double &FL, double &FT, double &FTT, double &FTL, double Einoff) const;
  /*! return all four structure functions in semi-inclusive kinematics and F2 (to extract prefactor)<BR>
   * For details see Phys.Rev. C84 (2011) 014601 Eqs. (24-30)
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param[out] FL FL structure function
   * \param[out] FT FT structure function
   * \param[out] FTT FTT structure function
   * \param[out] FTL FTL structure function
   * \param[out] F2 F2N nucleon structure function
   * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  void getStructureFunctions(TKinematics2to2 &kin, double &FL, double &FT, double &FTT, double &FTL, double &F2, double Einoff) const;
  /*! return all four structure functions in semi-inclusive kinematics <BR>
   * For details see Phys.Rev. C84 (2011) 014601 Eqs. (24-30)<BR>
   * This one is for DIS where an off-shell X is produced, W is the mass and enters as input
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param[out] FL FL structure function
   * \param[out] FT FT structure function
   * \param[out] FTT FTT structure function
   * \param[out] FTL FTL structure function
   * \param Wsq [MeV^2] mass of produced X squared, differs from invariant p_X^2 due to off-shellness
   * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  void getStructureFunctions_off(TKinematics2to2 &kin, double &FL, double &FT, double &FTT, double &FTL, double Wsq, double Einoff) const;
  /*! gives you the combination of all structure functions
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param[in] el electron kinematics, see TElectronKinematics
   * \param phi [rad] angle between electron scattering plane and hadron reaction plane
   * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  double getStructure(TKinematics2to2 &kin, TElectronKinematics &el, double phi, double Einoff) const;
  /*! combination of all structure functions averaged over phi (angle between hadron and electron plane)
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
    * \param[in] el electron kinematics, see TElectronKinematics
     * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  double getavgStructure(TKinematics2to2 &kin, TElectronKinematics &el, double Einoff) const;
  /*! combination of all structure functions averaged over phi (angle between hadron and electron plane) for LC formalism
   * \param kin semi-exclusive gamma+D->X+N kinematics, see LightConeKin2to2 <BR>
   */
  double getavgStructureLC(LightConeKin2to2 &kin) const;
  /*! combination of all structure functions NOT averaged over any angle for LC formalism
   * \param kin semi-exclusive gamma+D->X+N kinematics, see LightConeKin2to2 <BR>
   */
  double getStructureLC(LightConeKin2to2 &kin) const;
  /*! prefactor needed in neutron structure function extraction,
   * obtained from combination of all structure functions averaged over phi (angle between hadron and electron plane) 
   * divided by F2N
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
    * \param[in] el electron kinematics, see TElectronKinematics
     * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  double getavgPrefactor(TKinematics2to2 &kin, TElectronKinematics &el, double Einoff) const;
  /*! combination of all structure functions averaged over phi (angle between hadron and electron plane)
   * Off-shell produced X case here!!!  mass W is known, but differs from invariant mass 
   * \param[in] kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
    * \param[in] el electron kinematics, see TElectronKinematics
    * \param Wsq [MeV^2] invariant mass squared of produced X
     * \param Einoff [MeV] off-shell energy of interacting nucleon
     * \result [] prefactor
   */
  double getavgStructure_off(TKinematics2to2 &kin, TElectronKinematics &el, double Wsq, double Einoff) const;
  /*! Deuteron F2 structure functions, used in the inclusive calculations, defined as
   * \f$ F_2^D = F_L^D + \frac{Q^2}{2|q|^2}\frac{\nu}{m_N}F_T^D \f$
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
  * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  double getInclStructure(TKinematics2to2 &kin, double Einoff) const;
  /*! Deuteron F2 structure functions, used in the inclusive calculations, defined as
   * \f$ F_2^D = F_L^D + \frac{Q^2}{2|q|^2}\frac{\nu}{m_N}F_T^D \f$ <BR>
   * Off-shell produced X case here!!!  mass W is known, but differs from invariant mass 
   * \param kin semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param Wsq [MeV^2] mass of produced X squared, differs from invariant p_X^2 due to off-shellness
  * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  double getInclStructure_off(TKinematics2to2 &kin, double Wsq, double Einoff) const;
  /*! Deuteron F2 response functions including electron factors, split in dependence on phi spectator
   * 
   * \param kin [in] semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param[in] el electron kinematics, see TElectronKinematics
   * \param[out] ResLandT longitudinal and transverse part
   * \param[out] ResTL TL part
   * \param [out] ResTT TT part
   * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  void getResponses(TKinematics2to2 &kin, TElectronKinematics &el, double &ResLandT, double &ResTT, double &ResTL, double Einoff) const;
  /*! Deuteron F2 response functions including electron factors, split in dependence on phi spectator
   * 
   * \param kin [in] semi-exclusive gamma+D->X+N kinematics, see TKinematics2to2 <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param[in] el electron kinematics, see TElectronKinematics
   * \param Wsq [MeV^2] mass of produced X squared, differs from invariant p_X^2 due to off-shellness
   * \param[out] ResLandT longitudinal and transverse part
   * \param[out] ResTL TL part
   * \param [out] ResTT TT part
   * \param Einoff [MeV] off-shell energy of interacting nucleon
   */
  void getResponses_off(TKinematics2to2 &kin, TElectronKinematics &el, double Wsq,
			double &ResLandT, double &ResTT, double &ResTL, double Einoff) const;
  const std::string & getName() const{return name;}
  
private:
//   TElectronKinematics electron; /*!< electron kinematics */
  bool proton; /*!< photon interacts with proton [1] or neutron [0] */
  std::string name; /*!< string that has the structure function parametrization name */
  double massi; /*!< mass of nucleon interacting with photon */
};
/** @} */
#endif