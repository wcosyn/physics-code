/*! \file NucleonStructure.hpp 
 * \brief Contains declaration of class NucleonStructure, contains a variety of nucleon structure function parametrizations
 * \author Wim Cosyn
 * \date 12/09/2016
 * 
 * \addtogroup MePhys
 * @{
 */
#ifndef NUCLEONSTRUCTURE_HPP
#define NUCLEONSTRUCTURE_HPP

#include <string>
#include <cmath>
class c_mstwpdf; //forward declaration

/*! \brief A class for a variety of nucleon structure parametrizations.  Has Christy&Bosted / Alekhin (leading twist) / SLAC. */
class NucleonStructure{
  
public:
  NucleonStructure(); /*!< default constructor, SLAC parametrization taken */
  /*! Constructor
   * \param name Name of the parametriation. Possibilities: <BR>
   * "CB": Christy & Bosted parametrization (see F1F209.f file)<BR>
   * "SLAC": SLAC paramtetrization from Bodek <BR>
   * "Alekhin": leading twist parametrization by Alekhin [see PRD 68,014002], also see alekhin.f file <BR>
   * "CTEQ": F2 based on the pdf's from CTEQ (code from Misak, see cteq.f file) <BR>
   * "HMRS": F2 from pdfs from unknown source (Shunzo Kumano is the source, not enough info in his code file, see hmrs-b.f) <BR>
   * "MSTW": F2 from LO (or nlo etc) from MSTW pdfs (see mstwpdf.h & mstw.cpp) WARNING: for some reason class instances can only be made as pointers, dynamic memory)<BR>
   */
  NucleonStructure(const std::string name);
  NucleonStructure(const NucleonStructure&); /*!< Copy Constructor */
  NucleonStructure& operator=(const NucleonStructure&); /*!< assignment operator */
  ~NucleonStructure(); /*!< Destructor */
  
  
  /*! Returns the F1 and F2 structure functions 
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_xQ(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
  /*! Returns the F1 structure functions starting from bjorken x and Q^2
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 structure function
   */
  double getF1_xQ(const bool proton, const double x, const double Q2) const; 
  /*! Returns the F2 structure functions starting from bjorken x and Q^2
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 structure function
   */
  double getF2_xQ(const bool proton, const double x, const double Q2) const; 
  /*! Returns the F1 and F2 structure functions 
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param Wsq [MeV^2] invariant mass produced X
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_WQ(double &F1, double &F2, const bool proton, const double Wsq, const double Q2) const;
  /*! Returns the F1 structure functions starting from bjorken x and Q^2
   * \param proton 0 neutron, 1 proton
   * \param Wsq [MeV^2] invariant mass produced X
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 structure function
   */
  double getF1_WQ(const bool proton, const double Wsq, const double Q2) const; 
  /*! Returns the F2 structure functions starting from bjorken x and Q^2
   * \param proton 0 neutron, 1 proton
   * \param Wsq [MeV^2] invariant mass produced X
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 structure function
   */
  double getF2_WQ(const bool proton, const double Wsq, const double Q2) const; 
  const std::string getName() const{return name;} /*!< returns the name of the chosen parametrization */
  
  /*! returns the g1 structure function from the grsv2000 parametrization (leading order,  hep-ph/0011215, Tools/grsv2000pdf_g1.f )
   * \param proton 1=proton 0=neutron
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return g1n structure functions
   */
  static double getG1_grsv2000(const bool proton, const double x, const double Q2);
  /*! returns the sigma_L/sigma_T ratio, from SLAC fits (code from Shunzo Kumano, see R199x.f)
   * This is the 1990 parametrization 0.1 < x < 0.9, 0.6 < Q^2 < 20 GeV^2
   * L. W. Withlow et al., Phys. Lett. B250 194(1990)
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return R=sigma_L/sigma_T (indicates scaling violations from Callan-Gross LT relation)
   */
  static double getr1990(const double x, const double Q2);
  /*! returns the sigma_L/sigma_T ratio, from SLAC fits (code from Shunzo Kumano, see R199x.f)
   * This is the 1998 parametrization 0.005 < x < 0.86, 0.5 < Q^2 < 130 GeV^2
   * (E143) K. Abe et al., Phys. Lett B452 194 (1999)
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return R=sigma_L/sigma_T (indicates scaling violations from Callan-Gross LT relation)
   */
  static double getr1998(const double x, const double Q2);


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
  std::string dir; /*!< share dir to read input files */
  c_mstwpdf *mstw;
  
  /*! proton F2 SLAC parametrization
   * \param massi [GeV] initial mass, off-shell
   * \param xp bjorken x
   * \param q2 [GeV^2] Q^2 of virtual photon
   * \param fm [GeV] invariant mass W
   * \return F2 of proton
   */
  double f2p_b(const double massi, const double xp, const double q2, const double fm) const;
  /*! neutron F2 SLAC parametrization
   * \param massi [GeV] initial mass, off-shell
   * \param xp bjorken x
   * \param q2 [GeV^2] Q^2 of virtual photon
   * \param fm [GeV] invariant mass W
   * \return F2 of neutron
   */
  double f2n_b(const double massi, const double xp, const double q2, const double fm) const;
  
  /*! SLAC parametrization helper function
   * \param wm [GeV] invariant mass
   * \param qsq [GeV^2] Q^2 of virtual photon
   * \return something :p
   */
  double bodek(const double wm, const double qsq) const;
  
  
  /*! returns F1 and F2 in the alekhin parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_Alekhin(double &F1, double &F2, const bool proton, const double x, const double Q2)const;
  /* \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * return F1 in Alekhin parametriation
   */
  double getF1_Alekhin(const bool proton, const double x, const double Q2) const; 
  /* \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * return F2 in Alekhin parametriation
   */
  double getF2_Alekhin(const bool proton, const double x, const double Q2) const; 
  /*! returns F1 and F2 in the SLAC parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_SLAC(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
    /* \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * return F2 in SLAC parametriation
   */
  double getF2_SLAC(const bool proton, const double x, const double Q2) const; 
  /*! returns F1 and F2 in the Christy & Bosted parametrization
   * \param[out] F1 F1 structure function, determined through Callan-Gross relation incl scaling violation
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_CB(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
    /* \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * return F1 in Christy & Bosted parametriation
   */
  double getF1_CB(const bool proton, const double x, const double Q2) const;
    /* \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * return F2 in Christy & Bosted parametriation
   */
double getF2_CB(const bool proton, const double x, const double Q2) const;
  /*! returns F1 and F2 in the CTEQ parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_CTEQ(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
    /* \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * return F2 in CTEQ parametriation
   */
double getF2_CTEQ(const bool proton, const double x, const double Q2) const;/*!< return F2 structure function in the CTEQ parametrization */
  
  /*! returns F1 and F2 in the HRMS-B parametrization
   * \param[out] F1 F1 structure function, determined through Callan-Gross relation incl scaling violation
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_HMRS(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
  
 /*! returns F1 and F2 in the MSTW parametrization
   * \param[out] F1 F1 structure function, dermined through Callan-Gross relation incl scaling violations
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_MSTW(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
  
};
/*! @} */
#endif