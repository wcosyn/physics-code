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

#include <DeuteronMomDistr.hpp>
#include <TKinematics2to2.h>
#include <DeuteronStructure.hpp>
#include <TElectronKinematics.h>
#include <LightConeKin2to2.hpp>

#include <string>

/*! \brief A class that computes tagged spectator deuteron cross sections */
class DeuteronCross{
  
public:
    /*! Constructor
   * \param wfname Deuteron wave function name, see TDeuteron constructor for possibilities ["AV18", "Paris", etc.]
   * \param proton photon interacts with proton [1] or neutron [0]
   * \param strucname which structure function parametrization ["CB"=Chrisy-Bosted, "SLAC", "Alekhin"], see NucleonStructure 
   * \param sigmain [mb] total rescattering cross section in FSI
   * \param betain [GeV^-2] slope parameter
   * \param epsilonin real part of scattering amplitude
   * \param betaoffin [GeV^-2] off-shell beta parameter when using offshellset 2
   * \param lambdain [GeV^2] lambda cutoff off-shell parameter when using offshellset 1
   * \param offshellset  which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: no off-shell amplitude, fully suppressed
   * - 4: full off-shell amplitude, no suppression 
   * \param looplimit max number of tries in loop to get prz pole [for FSI calculation]
   */
  DeuteronCross(std::string wfname, bool proton, std::string strucname,
    double sigmain, double betain, double epsilonin, double betaoffin, double lambdain, int offshellset, int looplimit);
  ~DeuteronCross(); /*!<Destructor */
  /*! get the average cross section in VNA formalism, only valid in lab frame
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param electron has electron kinematics
   * \param pw plane-wave calculation [1] or not [0]
   * \param Einoff off-shell energy of the nucleon interacting with the photon
   * \return [nb/GeV^4] semi-inclusive cross section \f$ \frac{d\sigma}{dxdQ^2\frac{d^3p_s}{E_s}} \f$ averaged over phi
   */
  double getavgVNALabCross(TKinematics2to2 &kin, TElectronKinematics &electron, bool pw, double Einoff);
  /*! get the average cross section in LC formalism, only valid for collinear LightConeKin2to2!!!
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * \param pw plane-wave calculation [1] or including FWI [0]
   * \return [nb/GeV^4] semi-inclusive cross section \f$ \frac{d\sigma}{dx_AdQ^2\frac{d^3p_s}{E_s}} \f$ 
   * averaged over phi (hadron-lepton plane angle) [see Wim's lightconeDIS.pdf note for formulas (eq. 42)]
   */
  double getavgLCCross(LightConeKin2to2 &kin, bool pw);
  /*! get the cross section in LC formalism
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * \param pw plane-wave calculation [1] or including FSI [0]
   * \return [nb/GeV^4] semi-inclusive cross section \f$ \frac{d\sigma}{dx_AdQ^2\frac{d^3p_s}{E_s}} \f$ 
   * not(!) averaged over phi (hadron-lepton plane angle) [see Wim's lightconeDIS.pdf note for formulas (eq. 40)]
   */
  double getLCCross(LightConeKin2to2 &kin, bool pw);
  /*! get the average cross section in VNA formalism, only valid for collinear LightConeKin2to2!!!
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * \param pw plane-wave calculation [1] or including FWI [0]
   * \return [nb/GeV^4] semi-inclusive cross section \f$ \frac{d\sigma}{dx_AdQ^2\frac{d^3p_s}{E_s}} \f$ 
   * averaged over phi (hadron-lepton plane angle) 
   * [see Wim's lightconeDIS.pdf note for formulas (eq. 42), only with VNA deuteron density now]
   */
  double getavgVNACross(LightConeKin2to2 &kin, bool pw);
  /*! get the cross section in VNA formalism
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * \param pw plane-wave calculation [1] or including FSI [0]
   * \return [nb/GeV^4] semi-inclusive cross section \f$ \frac{d\sigma}{dx_AdQ^2\frac{d^3p_s}{E_s}} \f$ 
   * not(!) averaged over phi (hadron-lepton plane angle) 
   * [see Wim's lightconeDIS.pdf note for formulas (eq. 40), only with VNA deuteron density now]
   */
  double getVNACross(LightConeKin2to2 &kin, bool pw);
  /*! computes the results like they are presented in the Deeps data
   * \param Q2 [MeV^2] four-momentum transfer
   * \param W [MeV] invariant mass of X
   * \param Ein [MeV] beam energy
   * \param pr [MeV] spectator momentum
   * \param costhetar angle of spectator with q
   * \param proton DIS on proton (1) or neutron (0)
   * \param[out] planewave [MeV^-3] plane wave result
   * \param[out]  fsi [MeV^-3]  fsi result
   */
  void setScatter(double sigmain, double betain, double epsin);
private:
  double massi; /*!< mass of nucleon interacting with photon */
  std::string strucname; /*!< structure function parametrization */
//   TElectronKinematics electron; /*!< electron kinematics */
  DeuteronMomDistr momdistr; /*!< object instance used to calculate the deuteron momentum distributions */
  DeuteronStructure structure; /*!< object used to calculate the structure functions needed in semi-inclusive deuteron scattering */

  /*! get the average cross section like it was generated in the BONUS MC
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param electron has electron kinematics
   * \param lc lightcone density to restore data or VNA one
   * \return [MeV^-6] semi-inclusive cross section \f$ \frac{d\sigma}{d\Omega dE'd^3p_s} \f$ averaged over phi
   */
  double getavgBonus(TKinematics2to2 &kin, TElectronKinematics &electron, bool lc);
  
  void getDeepsresult(double Q2, double W, double Ein, double pr, double costhetar, bool proton,
    double &planewave, double &fsi);
  /*! computes the results like they are presented in the Deeps data in the LC formalism
   * \param Q2 [MeV^2] four-momentum transfer
   * \param W [MeV] invariant mass of X
   * \param Ein [MeV] beam energy
   * \param pr [MeV] spectator momentum
   * \param costhetar angle of spectator with q
   * \param proton DIS on proton (1) or neutron (0)
   * \param[out] planewave [MeV^-3] plane wave result
   * \param[out]  fsi [MeV^-3]  fsi result
   */
  void getDeepsresultLC(double Q2, double W, double Ein, double pr, double costhetar, bool proton,
    double &planewave, double &fsi);
  /*! computes the avg cross section with the formula as used in the Bonus MC (TKachenko arXiv:1402.2477)
   * \param Q2 [MeV^2] four-momentum transfer
   * \param W [MeV] invariant mass of X
   * \param Ein [MeV] beam energy
   * \param pr [MeV] spectator momentum
   * \param costhetar angle of spectator with q
   * \param proton DIS on proton (1) or neutron (0)
   * \param pw [1] plane wave only
   * \param lc [1] lightcone or [0] VNA deuteron density to restore data
   * \param[out] MCresult [GeV^-6] plane wave result using formulas from the BONUS MC simulation.
   * Semi-inclusive cross section \f$ \frac{d\sigma}{d\Omega dE'd^3p_s} \f$ integrated over phi.
   * \param[out] modelresultpw [GeV^-6] plane wave result using our model.
   * Semi-inclusive cross section \f$ \frac{d\sigma}{d\Omega dE'd^3p_s} \f$ integrated over phi.
   * \param[out] modelresultfsi [GeV^-6] fsi result using our model.
   * Semi-inclusive cross section \f$ \frac{d\sigma}{d\Omega dE'd^3p_s} \f$ integrated over phi.
   * 
   */
  void getBonusMCresult(double &MCresult, double &modelresultpw, double &modelresultfsi, 
			double Q2, double W, double Ein, double pr, double costhetar, bool proton, bool pw, bool lc);
    
  /*! computes the avg cross section with the formula as used in the Bonus MC (TKachenko arXiv:1402.2477)
   * \param Q2 [MeV^2] four-momentum transfer
   * \param W [MeV] invariant mass of X
   * \param Ein [MeV] beam energy
   * \param pr [MeV] spectator momentum
   * \param costhetar angle of spectator with q
   * \param proton DIS on proton (1) or neutron (0)
   * \param lc [1] lightcone or [0] VNA deuteron density to restore data
   * \param xref bjorken x to which we want to extrapolate
   * \param norm normalization correction for the data from the fits we did earlier
   * \param Rdata R ratio from the bonusdata
   * \param error stat + sys error associated with R [data] above
   */
  void getBonusextrapolate(double Q2, double W, double Ein, double pr, double costhetar, bool proton, 
			   bool lc, double xref, double norm, double Rdata, double error);

  /*! get the average Azz observable or cross section ratio
   * \param kin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in TKinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param elec has electron kinematics
   * \param pw plane-wave calculation [1] or not [0]
   * \param Azz [0] cross section or [1] Azz
   * \param Einoff off-shell energy of the nucleon interacting with the photon
   * \return [nb/GeV^4] semi-inclusive cross section \f$ \frac{d\sigma}{dxdQ^2\frac{d^3p_s}{E_s}} \f$ averaged over phi
   * or Azz observable
   */
  double getavgAzz(TKinematics2to2 &kin,TElectronKinematics &elec, bool Azz, bool pw, double Einoff); 
  /*! computes the results like they are presented in the Deeps data
   * \param Q2 [MeV^2] four-momentum transfer
   * \param W [MeV] invariant mass of X
   * \param Ein [MeV] beam energy
   * \param pr [MeV] spectator momentum
   * \param costhetar angle of spectator with q
   * \param proton DIS on proton (1) or neutron (0)
   * \param[out] Azz [] Azz plane wave result
   * \param[out]  Azzfsi [] Azz fsi result
   * \param[out] planewave [MeV^-3] cross section ratio plane wave result
   * \param[out]  fsi [MeV^-3]  cross section ratio fsi result
   */
  void getDeepsAzz(double Q2, double W, double Ein, double pr, double costhetar, bool proton, 
		    double &Azz, double &Azzfsi, double &planewave, double &fsi); 

  
  
  /*! reads in the Deeps data set in an array >
   * \param[out] pdeepsarray pointer to array with the deeps results, indices as follows <BR>
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
  /*! Set scattering parameters for FSI
   * \param sigmain total cross section [mb]
   * \param betain slope parameter [GeV^-2]
   * \param epsin real part of amplitude
   */  
  
  
};
/*! @} */
#endif