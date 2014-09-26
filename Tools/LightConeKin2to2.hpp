/*! \file LightConeKin2to2.hpp 
 * \brief Contains declaration of class LightConeKin2to2, contains (possibly) all necessary lightcone variables for a two to two process
 * \author Wim Cosyn
 * \date 9/9/2014
 * 
 * \addtogroup MePhys
 * @{
 */
#ifndef LIGHTCONEKIN2TO2_HPP
#define LIGHTCONEKIN2TO2_HPP

#include <string>
#include <TVector3.h>
#include <FourVector.h>

/*! \brief A class that initializes all lightcone variables for a two body scattering process, 
 *written with DIS on a deuteron with spectator in mind!
 */
class LightConeKin2to2{

public:
  /*!Constructor for the most general case
   *\param massA [MeV] mass of hadron target
   *\param Q2 [MeV^2] minus four-momentum squared of virtual photon (in case this is used for hadrons, put -Mass^2 here)
   *\param massC [MeV] mass of spectator
   *\param[in] vecA [MeV] momentum vector of hadron target in the frame under consideration
   *\param[in] vecq [MeV] momentum vector of virtual photon
   *\param[in] vecC [MeV] momentum vector of spectator
   *\param[in] vecBeam [MeV] momentum vector of beam (assumed massless electron)
   */
  LightConeKin2to2(double massA, double Q2, double massC, TVector3 &vecA, TVector3 &vecq, TVector3 &vecC, TVector3 &vecBeam); 
  /*!Constructor for the collinear case, where hadron and virtual photon are on the same axis, taken to be the z-axis!!!
   *\param massA [MeV] mass of hadron target
   *\param Q2 [MeV^2] minus four-momentum squared of virtual photon (in case this is used for hadrons, put -Mass^2 here)
   *\param massC [MeV] mass of spectator
   *\param pA [MeV] momentum of hadron target along the collinear axis (be wary of sign!)
   *\param vecq [MeV] momentum of virtual photon along the collinear axis (be wary of sign!)
   *\param[in] vecC [MeV] momentum vector of spectator
   *\param[in] vecBeam [MeV] momentum vector of beam (assumed massless electron)
   */
  LightConeKin2to2(double massA, double Q2, double massC, double pA, double vecq, TVector3 &vecC, TVector3 &vecBeam); 

  /*! Copy constructor */
  LightConeKin2to2(LightConeKin2to2 &rhs);
  
  double getMassA() const{return massA;} /*!<[MeV] mass hadron target */
  double getQ2() const{return Q2;} /*!< [MeV^2] Q2 virtual photon */
  double getmassN () const{return massN;} /*!<[MeV] mass spectator */
  double getMassX() const{return massX;} /*!<[MeV] mass undetected particle, DIS X here */
  const FourVector<double> &getA_mu() const{return A_mu;} /*!< [MeV] fourvector hadron target */
  const FourVector<double> &getQ_mu() const{return q_mu;}/*!< [MeV] fourvector virtual photon */
  const FourVector<double> &getPs_mu() const{return ps_mu;} /*!< [MeV] fourvector spectator */
  const FourVector<double> &getX_mu() const{return X_mu;} /*!< [MeV] fourvector undetected DIS product X */
  const FourVector<double> &getBeam_mu() const{return Beam_mu;}/*!< [MeV] fourvector beam */
  const FourVector<double> &getPi_mu() const{return pi_mu;} /*!< [MeV] fourvector of struck nucleon */
  double getYA() const{return yA;} /*!< [] A_mu*q_mu/beam_mu*A_mu */	
  double getXA() const{return xA;} /*!< [] Q^2/(A_mu*q_mu) Bjorken x for target, careful with the factor here!!! */
  double getYN() const{return yN;} /*!< [] pi_mu*q_mu/beam_mu*pi_mu */
  double getXN() const{return xN;} /*!< [] Q^2/(2pi_mu*q_mu) Bjorken x for struck nucleon */
  double getPA_plus() const{return pA_plus;} /*!< LC plus component of target */
  const TVector3 &getPA_perp() const{return pA_perp;} /*!< [MeV] perp component of target momentum */
  const TVector3 &getPs_perp() const{return ps_perp;} /*!< [MeV] perp component of spectator momentum */
  const TVector3 &getQ_perp() const{return ps_perp;} /*!< [MeV] perp component of virtual photon */
  double getAlpha_s() const{return alpha_s;} /*!< LC plus fraction of spectator 2p_s^+/p_A^+ */
  double getAlpha_i() const{return alpha_i;} /*!< LC plus fraction of struck nucleon */
  double getEpsilon() const{return epsilon;} /*!< [] longitudinal to transverse polarization ratio of virtual photon */
  const TVector3 &getK_perp() const{return k_perp;} /*!< [MeV] perp component of LC rescaled deuteron momentum */
  double getK_z() const{return k_z;} /*!< [MeV] z-component of LC rescaled deuteron momentum */
  double getK() const{return k;} /*!< [MeV] norm of LC rescaled deuteron momentum */
  const TVector3 &getKvec() const{return kvec;} /*!< [MeV] threevector of LC rescaled deuteron momentum */
  double getEk() const{return Ek;} /*!< [MeV] on-shell energy of LC rescaled deuteron momentum */
  double getZs() const{return zs;} /*!< [MeV] (A_mu*ps_mu)/(A_mu*q_mu) */
  
  double getPs_norm2() const{return ps_mu[1]*ps_mu[1]+ps_mu[2]*ps_mu[2]+ps_mu[3]*ps_mu[3];}
  
  double getS() const{return (q_mu+A_mu)*(q_mu+A_mu);} /*!< [MeV^2] Mandelstam s */
  
  bool getIsCollinear() const{return isCollinear;} /*!< [0] general frame, [1] collinear frame*/


  
private:
  double massA; /*!<[MeV] mass hadron target */
  double Q2; /*!< [MeV^2] Q2 virtual photon */
  double massN; /*!<[MeV] mass spectator */
  double massX; /*!<[MeV] mass undetected particle, DIS X here */
  FourVector<double> A_mu; /*!< [MeV] fourvector hadron target */
  FourVector<double> q_mu;/*!< [MeV] fourvector virtual photon */
  FourVector<double> ps_mu; /*!< [MeV] fourvector spectator */
  FourVector<double> X_mu; /*!< [MeV] fourvector undetected DIS product X */
  FourVector<double> Beam_mu;/*!< [MeV] fourvector beam */
  FourVector<double> pi_mu; /*!< [MeV] fourvector of struck nucleon */
  double yA; /*!< [] A_mu*q_mu/beam_mu*A_mu */	
  double xA; /*!< [] Q^2/(A_mu*q_mu) Bjorken x for target, careful with the factor here!!! */
  double yN; /*!< [] pi_mu*q_mu/beam_mu*pi_mu */
  double xN; /*!< [] Q^2/(2pi_mu*q_mu) Bjorken x for struck nucleon */
  double pA_plus; /*!< [MeV] LC plus component of target */
  TVector3 pA_perp; /*!< [MeV] perp component of target momentum */
  TVector3 ps_perp; /*!< [MeV] perp component of spectator momentum */
  TVector3 q_perp; /*!< [MeV] perp component of virtual photon */
  double alpha_s; /*!< LC plus fraction of spectator 2p_s^+/p_A^+ */
  double alpha_i; /*!< LC plus fraction of struck nucleon */
  double epsilon; /*!< [] longitudinal to transverse polarization ratio of virtual photon */
  TVector3 k_perp; /*!< [MeV] perp component of LC rescaled deuteron momentum */
  double k_z; /*!< [MeV] z-component of LC rescaled deuteron momentum */
  double k; /*!< [MeV] norm of LC rescaled deuteron momentum */
  TVector3 kvec; /*!< [MeV] threevector of LC rescaled deuteron momentum */
  double Ek; /*!< [MeV] on-shell energy of LC rescaled deuteron momentum */
  double zs; /*!< [MeV] (A_mu*ps_mu)/(A_mu*q_mu) */
  
  bool isCollinear; /*!< [0] general frame, [1] collinear frame*/
};
/** @} */
#endif