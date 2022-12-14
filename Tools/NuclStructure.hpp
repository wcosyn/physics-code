/*! \defgroup MePhys libMePhys: Library that contains all classes for general medium physics purposes
 * \author Wim Cosyn, Pieter Vancraeyveld
 * \date 29/08/2012
 * \brief This code implements classes to do a lot of general medium physics stuff
 * 
 * \details 
 * 
 * I jacked and hacked most of the stuff included here from the strangecalc libraries.
 * Generally I did not doxygen-ize these classes.  The relevant API can be found at 
 * http://rprmodel.ugent.be/api/. <BR><BR>
 * 
 * I also added some new classes: <BR>
 * - A class that contains nucleon structure functions [NuclStructure] <BR>
 * - A class that contains deuteron structure functions [DeuteronStructure] <BR>
 * - A leptonkinematics class for neutrino scattering [TLeptonKinematics] <BR>
 * - A class for he3 wave functions [He3wf] <BR>
 * -A class that has lightcone kinematics for a two body scattering reaction  [LightConeKin2to2] <BR><BR>
 * -A class that has lightcone kinematics for a two to three body scattering reaction  [LightConeKin2to3] <BR><BR>
 * 
 * Utilfunctions.hpp contains a few helper functions (mostly math related)<BR>
 * --------------------------------------<BR>
 * 
 * Some changes I did in Pieter's files and classes: <BR>
  * -moved from kinematics2to3.h and.cpp to fourvector.h: <BR>
  * // Global functions for conversion <BR>
  * // from ROOT 4vectors to strangecalc 4vectors.<BR>
  * template<typename T> <BR>
  * FourVector<T> operator*(const TLorentzRotation&, const FourVector<T>&); <BR>
  * <BR>
  * FourVector<double> ToFourVector(const TLorentzVector&); <BR>
  * TLorentzVector ToLorentzVector(const FourVector<double>&); <BR>
  * <BR>
  * -moved definitions of ToFourVector and ToLorentzVector from kinematics2to3.cpp to a new file FourVector.cpp <BR>
  * <BR>
  * -commented out include "2to3kinematics.h" in TDeuteron.cpp <BR>
  * <BR>
  * -commented out all hash functions <BR>
  * -commented out numtoa.h includes <BR>
  * <BR>
  * -moved all the necessary wrapper classes in the dir with deuteron code, changed include <> to include "" <BR>
  * <BR>
  * -commented out all includes of Structures.h as only defines of constants were needed <BR>
  * <BR>
  * -included needed constants (STRANGEUFLOW, M_P, M_N) to my file constants.hpp (some  double definitions now of course) <BR>
  * <BR>
  * -replaced all includes of Structures.h with include "constants.hpp" (was not even needed in TYukawaPWF.cpp) <BR>
  * <BR>
  * -changed TYukawaPWF to also include parmetrization with non-equidistant Mass(i) (AV18 f.i.) <BR>
  *  -added double *fM member <BR>
  * -added fM(NULL) to initializer list of constructor<BR>
  * -added new constructor with array of m values<BR>
  * -check in constructor to delete fM when necessary<BR>
  * -modified Mass(int i) for this case<BR>
  * <BR>
  *-added AV18 and Av18b Deuteron wf to TDeuteron.h & TDeuteron.cpp <BR>
  *<BR>
  *-added coordinate state wavefunction functions to TDeuteron<BR>
  *-some more vector types to DeuteronRState and DeuteronPState<BR>
  *-added off-shell momentum wave function parts to TWaveFunctionImplementation and derived classes
  * (GetUpoff etc, DeuteronPStateOff functions...)<BR>
  *<BR>
  *-commented out hash functions in TKinematics.cpp<BR>

  *-added hypron kin vars in TKinematics.cpp: updated constructors, getters, updatekinematics<BR>
  *-pklab: ++energycount changed to ++anglecount<BR>
  *-l1250 added pklab kinematics<BR>
  *<BR>
  *-added tan2 to Electronkin<BR>
  *<BR>
  *-electronkinematics changed all 2to3 to 2to2<BR>
  *<BR>
  *-added TKinematics2to3 to Tools, commented out functions that I moved previously...<BR>
  *<BR>
  *-added masses to constructors to Tkinematics2to3(withlabangles)<BR>
  *changed setisospin function for masses functionality<BR>
 *
 */



/*! \file NuclStructure.hpp 
 * \brief Contains declaration of class NuclStructure, contains a variety of nucleon structure function parametrizations
 * \author Wim Cosyn
 * \date 29/08/2012
 * 
 * \addtogroup MePhys
 * @{
 */
#ifndef NUCLSTRUCTURE_HPP
#define NUCLSTRUCTURE_HPP

#include <string>
#include <cmath>
class c_mstwpdf; //forward declaration

/*! \brief A class for a variety of nucleon structure parametrizations.  Has Christy&Bosted / Alekhin (leading twist) / SLAC. */
class NuclStructure{
  
public:
  /*! Constructor
   * \param proton proton [1] or neutron [0]
   * \param var1 first variable, see below
   * \param var2 second variable, see below
   * \param switchvar Which input variables? Possibilities: <BR>
   * [0]: var1 is Q2 [MeV^2] and var2 is Bjorken x [] <BR>
   * [1]: var1 is Q2 [MeV^2] and var2 is W^2 [MeV^2] <BR>
   * [2]: var1 is Bjorken x [] and var2 is W^2 [MeV^2] <BR>
   * \param name Name of the parametriation. Possibilities: <BR>
   * "CB": Christy & Bosted parametrization (see F1F209.f file) if invariant mass is >5 GeV, auto switch to SLAC!!!<BR>
   * "SLAC": SLAC paramtetrization from Bodek <BR>
   * "Alekhin": leading twist parametrization by Alekhin [see PRD 68,014002], also see alekhin.f file <BR>
   * "CTEQ": F2 based on the pdf's from CTEQ (code from Misak, see cteq.f file) <BR>
   * "HMRS": F2 from pdfs from unknown source (Shunzo Kumano is the source, not enough info in his code file, see hmrs-b.f) <BR>
   * "MSTW": F2 from LO (or nlo etc) from MSTW pdfs (see mstwpdf.h & mstw.cpp) <BR>
   */
  NuclStructure(const bool proton, double var1, double var2, int switchvar, std::string name);
  /*! Constructor 
   * \param proton proton [1] or neutron [0]
   * \param Q2in [MeV^2] Q^2 of the virtual photon
   * \param xin [] Bjorken x
   * \param Wsqin [MeV^2] invariant mass squared W^2
   * \param name Name of the parametriation. Possibilities: <BR>
   * "CB": Christy & Bosted parametrization (see F1F209.f file) if invariant mass is >5 GeV, auto switch to SLAC!!!<BR>
   * "SLAC": SLAC paramtetrization from Bodek <BR>
   * "Alekhin": leading twist parametrization by Alekhin [see PRD 68,014002], also see alekhin.f file <BR>
   * "CTEQ": F2 based on the pdf's from CTEQ (code from Misak, see cteq.f file) <BR>
   * "HMRS": F2 from pdfs from unknown source (Shunzo Kumano is the source, not enough info in his code file, see hmrs-b.f) <BR>
   * "MSTW": F2 from LO (or nlo etc) from MSTW pdfs (see mstwpdf.h & mstw.cpp) <BR>
   */
  NuclStructure(const bool proton, double Q2in, double xin, double Wsqin, std::string name);

  NuclStructure(const NuclStructure&); /*!< Copy Constructor */
  NuclStructure& operator=(const NuclStructure&); /*!< assignment operator */
  ~NuclStructure(); /*!< Destructor */
  
  
  /*! Returns the F1 and F2 structure functions */
  void getF(double &F1, double &F2);
  double getF1(); /*!< returns the F1 structure function */
  double getF2(); /*!< returns the F2 structure function */
  const std::string getName() const{return name;} /*!< returns the name of the chosen parametrization */
  double getW(){return sqrt(Wsq);} /*!< [MeV] returns the invariant mass of X*/
  bool getProton() const{return proton;}
  
  /*! returns the g1 structure function from the grsv2000 parametrization (leading order,  hep-ph/0011215, Tools/grsv2000pdf_g1.f )
   * \param proton 1=proton 0=neutron
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return g1n structure functions
   */
  static double getG1_grsv2000(const bool proton, double x, double Q2);
  /*! returns the sigma_L/sigma_T ratio, from SLAC fits (code from Shunzo Kumano, see R199x.f)
   * This is the 1990 parametrization 0.1 < x < 0.9, 0.6 < Q^2 < 20 GeV^2
   * L. W. Withlow et al., Phys. Lett. B250 194(1990)
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return R=sigma_L/sigma_T (indicates scaling violations from Callan-Gross LT relation)
   */
  static double getr1990(double x, double Q2);
  /*! returns the sigma_L/sigma_T ratio, from SLAC fits (code from Shunzo Kumano, see R199x.f)
   * This is the 1998 parametrization 0.005 < x < 0.86, 0.5 < Q^2 < 130 GeV^2
   * (E143) K. Abe et al., Phys. Lett B452 194 (1999)
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return R=sigma_L/sigma_T (indicates scaling violations from Callan-Gross LT relation)
   */
  static double getr1998(double x, double Q2);


  /*! returns the g1+g2 structure function using the Wandura-Wilczek relation [g1+g2=int dy/y g1(y)] from the grsv2000 parametrization (leading order,  hep-ph/0011215, Tools/grsv2000pdf_g1.f )
   * This sum of structure functions is related to the transversity pdfs (and higher twist).
   * \param proton 1=proton 0=neutron
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return g1+g2 structure functions
   */
  static double getG1plusg2_grsv2000(const bool proton, double x, double Q2);
  
private:
  std::string name; /*!< name of the paramtetrization */
  bool proton; /*!< proton [1] or neutron [0] structure functions */
  double mass; /*!< [MeV] nucleon mass */
  std::string dir; /*!< share dir to read input files */
  double x; /*!< Bjorken x */
  double Q2; /*!< [MeV^2] Q^2 of the virtual photon */
  double Wsq; /*!< [MeV^2] invariant mass squared W^2 */
  c_mstwpdf *mstw;
  
  /*! proton F2 SLAC parametrization
   * \param massi [GeV] initial mass, off-shell
   * \param xp bjorken x
   * \param q2 [GeV^2] Q^2 of virtual photon
   * \param fm [GeV] invariant mass W
   * \return F2 of proton
   */
  double f2p_b(double massi, double xp, double q2, double fm) const;
  /*! neutron F2 SLAC parametrization
   * \param massi [GeV] initial mass, off-shell
   * \param xp bjorken x
   * \param q2 [GeV^2] Q^2 of virtual photon
   * \param fm [GeV] invariant mass W
   * \return F2 of neutron
   */
  double f2n_b(double massi, double xp, double q2, double fm) const;
  
  /*! SLAC parametrization helper function
   * \param wm [GeV] invariant mass
   * \param qsq [GeV^2] Q^2 of virtual photon
   * \return something :p
   */
  double bodek(double wm, double qsq) const;
  
  
  /*! returns F1 and F2 in the alekhin parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   */
  void getF_Alekhin(double &F1, double &F2);
  double getF1_Alekhin(); /*!< return F1 structure function in the Alekhin parametrization */
  double getF2_Alekhin(); /*!< return F2 structure function in the Alekhin parametrization */
  /*! returns F1 and F2 in the SLAC parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   */
  void getF_SLAC(double &F1, double &F2) const;
  double getF1_SLAC() const; /*!< return F1 structure function in the SLAC parametrization */
  double getF2_SLAC() const; /*!< return F2 structure function in the SLAC parametrization */
  /*! returns F1 and F2 in the Christy & Bosted parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   */
  void getF_CB(double &F1, double &F2) const;
  double getF1_CB() const;/*!< return F1 structure function in the Christy & Bosted parametrization */
  double getF2_CB() const;/*!< return F2 structure function in the Chrsity & Bosted parametrization */
  /*! returns F1 and F2 in the CTEQ parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   */
  void getF_CTEQ(double &F1, double &F2);
  double getF1_CTEQ();/*!< return F1 structure function in the CTEQ parametrization */
  double getF2_CTEQ();/*!< return F2 structure function in the CTEQ parametrization */
  
  /*! returns F1 and F2 in the HRMS-B parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   */
  void getF_HMRS(double &F1, double &F2);
  
 /*! returns F1 and F2 in the MSTW parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   */
  void getF_MSTW(double &F1, double &F2);
  
};
/*! @} */
#endif