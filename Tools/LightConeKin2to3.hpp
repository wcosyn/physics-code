/*! \file LightConeKin2to3.hpp 
 * \brief Contains declaration of class LightConeKin2to3, contains (possibly) all necessary lightcone variables for a two to three process
 * \author Wim Cosyn
 * \date 9/9/2014
 * 
 * \addtogroup MePhys
 * @{
 */
#ifndef LIGHTCONEKIN2TO3_HPP
#define LIGHTCONEKIN2TO3_HPP

#include <string>
#include <TVector3.h>
#include <FourVector.h>
#include <LightConeKin2to2.hpp>

/*! \brief A class that initializes all lightcone variables for a two body scattering process, 
 *written with DIS on he3 with two spectators in mind!
 */
class LightConeKin2to3: public LightConeKin2to2{

public:
  /*!Constructor for the most general case
   *\param massA [MeV] mass of hadron target
   *\param Q2 [MeV^2] minus four-momentum squared of virtual photon (in case this is used for hadrons, put -Mass^2 here)
   *\param massSp1 [MeV] mass of spectator 1
   *\param massSp2 [MeV] mass of spectator 1
   *\param[in] vecA [MeV] momentum vector of hadron target in the frame under consideration
   *\param[in] vecq [MeV] momentum vector of virtual photon
   *\param[in] vecSp1 [MeV] momentum vector of spectator1
   *\param[in] vecSp2 [MeV] momentum vector of spectator2
   *\param[in] vecBeam [MeV] momentum vector of beam (assumed massless electron)
   * \param[in] vecSpcm [MeV] momentum vector of cm of spectators
   */
  LightConeKin2to3(double massA, double Q2, double massSp1, double massSp2, 
		   TVector3 &vecA, TVector3 &vecq, TVector3 &vecSp1, TVector3 &vecSp2, TVector3 &vecBeam, TVector3 &vecSpcm); 
  /*!Constructor for the collinear case, where hadron and virtual photon are on the same axis, taken to be the z-axis!!!
    *\param massA [MeV] mass of hadron target
    *\param Q2 [MeV^2] minus four-momentum squared of virtual photon (in case this is used for hadrons, put -Mass^2 here)
    *\param massSp1 [MeV] mass of spectator 1
    *\param massSp2 [MeV] mass of spectator 1
    *\param pA [MeV] momentum of hadron target along the collinear axis (be wary of sign!)
    *\param vecq [MeV] momentum of virtual photon along the collinear axis (be wary of sign!)
    *\param[in] vecSp1 [MeV] momentum vector of spectator1
    *\param[in] vecSp2 [MeV] momentum vector of spectator2
    *\param[in] vecBeam [MeV] momentum vector of beam (assumed massless electron)
   * \param[in] vecSpcm [MeV] momentum vector of cm of spectators
   */
  LightConeKin2to3(double massA, double Q2, double massSp1, double massSp2, 
		   double pA, double vecq, TVector3 &vecSp1, TVector3 &vecSp2, TVector3 &vecBeam, TVector3 &vecSpcm); 
  /*! Copy Constructor */
  LightConeKin2to3(LightConeKin2to3 &rhs);
  
  const FourVector<double> &getSp1_mu() const{return Sp1_mu;} /*!< [MeV] fourvector spectator 1 */
  const FourVector<double> &getSp2_mu() const{return Sp2_mu;}/*!< [MeV] fourvector spectator 2 */
  
private:
  FourVector<double> Sp1_mu;/*!< [MeV] fourvector spectator 1 */
  FourVector<double> Sp2_mu;/*!< [MeV] fourvector spectator 2 */
};
/** @} */
#endif