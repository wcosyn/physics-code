/*! \file He3wf.hpp 
 * \brief Contains declaration of class He3wf, has a variety of CD-Bonn and AV18 he3 wave functions
 * uses fortran code obtained from Misak Sargsian [getwavemod_tm.f,clebsch.f,spline.f].  This class is basically a wrapper for that code.
 * \author Wim Cosyn
 * \date 19/9/2014
 * 
 * \addtogroup MePhys
 * @{
 */
#ifndef HE3WF_HPP
#define HE3WF_HPP

#include <string>
#include <TVector3.h>
#include <complex>

/*! \brief A class that returns the momentum space he3 wave function for a variety of parametrizations
 */
class He3wf{

public:
  /*!Constructor
   * \param wf_name parametrization of the wave function.  Possibilities (Case sensitive!!):
   * "AV18", "AV18-UIX", "AV18-TML", "CDBonn", "CDBonn-TML"
   * \param dir dir where the data files are located
   */
  He3wf(std::string wf_name, std::string dir); 
  /*! Calculate the he3 wf for a certain momentum, spin and isospin configuration
   * \param p1 [MeV] momentum of particle 1
   * \param p2 [MeV] momentum of particle 2
   * \param ms1 [-1,1] twice the spin projection of particle 1
   * \param ms2 [-1,1] twice the spin projection of particle 2
   * \param mt1 [-1,1] twice the isospin projection of particle 1 (1 proton, -1 neutron)
   * \param ms2 [-1,1] twice the isospin projection of particle 2 (1 proton, -1 neutron)
   * \param mA [-1,1] twice the spin projection of the nucleus
   * \return [MeV^-3] returns the complex valued He3 wave function
   */
  std::complex<double> getWF(TVector3 &p1, TVector3 &p2, int ms1, int ms2, int mt1, int mt2, int mA);

  /*! Calculate the he3 wf for a certain momentum, spin and isospin configuration
   * \param p1 [MeV] momentum of particle 1
   * \param p2 [MeV] momentum of particle 2
   * \param ms1 [-1,1] twice the spin projection of particle 1
   * \param ms2 [-1,1] twice the spin projection of particle 2
   * \param mt1 [-1,1] twice the isospin projection of particle 1 (1 proton, -1 neutron)
   * \param ms2 [-1,1] twice the isospin projection of particle 2 (1 proton, -1 neutron)
   * \param mA [-1,1] twice the spin projection of the nucleus
   * \return [MeV^-3] returns the complex valued He3 wave function
   */
  std::complex<double> getWF(double p1[3], double p2[3], int ms1, int ms2, int mt1, int mt2, int mA);
  
private:
  int mA_set; /*!< [-1,1] twice (!) the spin mA of the he3 nucleus.  This needs to be set before calling the wave function */
};
/** @} */
#endif