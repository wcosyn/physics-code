/*! \file DeuteronCross.hpp 
 * \brief Contains declaration of class He3Cross, computes he3 exclusive DIS cross sections
 * \author Wim Cosyn
 * \date 26/09/2014
 * 
 * \addtogroup DIS
 * @{
 */

#ifndef HE3CROSS_HPP
#define HE3CROSS_HPP

#include <He3wf.hpp>
#include <DeuteronStructure.hpp>
#include <LightConeKin2to3.hpp>

#include <string>

/*! \brief A class that computes tagged spectator he3 cross sections */
class He3Cross{
public:
    /*! Constructor
   * \param wfname He3 wave function name, see He3wf class constructor for possibilities ["AV18", etc.]
   * \param inputdir share dir [full path please, else fortran code breaks for some reason]
   * \param strucname which structure function parametrization ["CB"=Chrisy-Bosted, "SLAC", "Alekhin"], see NucleonStructure
   * \param proton photon interacts with proton [1] or neutron [0]
   */
  He3Cross(std::string wfname, std::string inputdir, std::string strucname, bool proton);
  He3Cross(const He3Cross&); /*!< Copy Constructor */
  He3Cross& operator=(const He3Cross&); /*!< assignment operator */
  ~He3Cross(); /*!<Destructor */
  
  /*! Calculate He3 momentum distribution for a certain kinematics of the spectators
   * \param[in] kin LightconeKin2to3 object that has all the kinematics
   * \return [MeV-4] momentum distribution S(ps1,ps2) times E_sp1_lab*E_sp2_lab, normed so that \f$ \int d3p1 d3p2 S(ps1,ps2) = 1 \f$. 
   */  
  double getDensity(LightConeKin2to3 &kin);
    /*! get the average cross section
   * \param kin kinematics object containing the gamma+He3->X+N1+N2 kinematics <BR>
   * \return [nb/GeV^6] tagged spectators cross section \f$ \frac{d\sigma}{dx_A dQ^2 \frac{d^3p_{s_1}}{E_{s_1}}\frac{d^3p_{s_2}}{E_{s_2}}} \f$ 
   * NOT averaged over phi [see he3.pdf note for formula details]
   */
  double getCross (LightConeKin2to3 &kin);
  
private:
  He3wf wavefunction;
  DeuteronStructure strucfunc;
  bool proton;

};
/*! @} */
#endif  
