/*! \file MeanFieldNucleusThick.hpp 
 * \brief Contains MeanFieldNucleusThick class declaration
 * \author Wim Cosyn
 * \date 16/08/2011
 * 
 */

#ifndef MEANFIELDNUCLEUSTHICK_H
#define MEANFIELDNUCLEUSTHICK_H

#include <string>

using namespace std;

#include <constants.hpp>
#include "MeanFieldNucleus.hpp"

/*! \brief A Class for a nucleus with mean-field wave functions for the individual nucleons and densities. 
 * 
 * Inherited from MeanFieldNucleus.
 * Contains all useful variables for the nucleus (A,Z,N,excitation energy of the levels, etc)
 * Provides all quantum numbers for individual nucleons.  Functions to compute the wave function for a certain
 * nucleon at a certain coordinate.  
 * Also provides proton, neutron and total densities.  
 * 
 * ALL THESE DENSITIES INCLUDE A FACTOR r^2 AND ARE NORMALIZED TO ONE!!!!!!!!!
 * ALL THESE DENSITIES INCLUDE A FACTOR r^2 AND ARE NORMALIZED TO ONE!!!!!!!!!
 * ALL THESE DENSITIES INCLUDE A FACTOR r^2 AND ARE NORMALIZED TO ONE!!!!!!!!!
 */
class MeanFieldNucleusThick : public MeanFieldNucleus{
public:
  /*! Constructor
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param dir string containing the dir were all input is located
   */
  MeanFieldNucleusThick(const int nucleus = 0, const string & dir = "."); 
  
  ~MeanFieldNucleusThick(); /*!< Destructor */
  //ALL DENSITIES ARE TIMES R^2 ALREADY!!!!!!!!!!!!
  //I REPEAT
  //ALL DENSITIES ARE TIMES R^2 ALREADY!!!!!!!!!!!!
  double getProtonDensity(const double r) const; /*!< get Proton density for r[fm] \param r [fm] coord (spherical) \return r*r*rho_p(r)*/
  double getNeutronDensity(const double r) const;  /*!< get Neutron density for r[fm] \param r [fm] coord (spherical) \return r*r*rho_n(r)*/
  double getTotalDensity(const double r) const; /*!< get Total density for r[fm] \param r [fm] coord (spherical)  \return r*r*rho(r)*/
  /*! get proton or neutron density for r[fm] 
   *\param r [fm] coord (spherical)
   *\param proton selects proton (1) or neutron(0) density
   *\return r*r*rho_(n or p)(r)
   */
  double getDensity(const double r, bool proton) const; 
  /*! get density for r[fm] 
   *\param r [fm] coord (spherical)
   *\param density density array: proton, neutron or total
   *\return r*r*density(r)
   */  
  double getDensity(const double r, const double *density) const; //get total density for r
  const double * getProtonDensity() const; /*!< returns the proton density array */
  const double * getNeutronDensity() const; /*!< returns the neutron density array */
  const double * getTotalDensity() const; /*!< returns the total density array */
  const double * getDensity(int proton) const; /*!< returns the proton (1) or neutron (0) density array \param proton selects proton (1) or neutron(0) density*/

private:
  double *protondensity; /*!< neutron density array*/
  double *neutrondensity; /*!< proton density array*/
  double *totaldensity; /*!< total density array*/
  void constructDensity(); /*!< Construct all density arrays, includes a factor r^2!!!!*/
  
  
};

#endif