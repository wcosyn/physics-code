/*! \file FastParticle.hpp 
 * \brief Contains FastParticle class declaration
 * \author Wim Cosyn
 * \date 16/08/2011
 * 
 */

#ifndef FASTPARTICLE_H
#define FASTPARTICLE_H

#include <string>
#include <TVector3.h>

#include "MeanFieldNucleus.hpp"


/*! \brief A Class for a fast particle ejected from a nucleus
 * 
 * Contains all useful kinematic variables as well as scattering parameters sigma, epsilon and beta^2 for 
 * scattering with a proton and neutron.  Also some coordinate transformations to obtain the (b,z) coordinates relative to 
 * this particle.
 */
class FastParticle{
public:
  /*! Constructor
   * \param particletype which particle [0-3] = [proton, neutron, pi+, pi-] can be extended of course
   * \param incoming is it the beam hadron particle?
   * \param momentum [MeV] momentum of the particle
   * \param theta [rad] spherical coord, theta angle of momentum with the z-axis
   * \param phi [rad] spherical coord, phi angle of momentum
   * \param hard_scale [GeV^2] for CT calculations, the hard scale associated with the particle (like Q^2 or |t|, etc)
   * \param Gamma [MeV] Decay width in rest frame (main branching assumed)
   * \param dir homedir where all input is located
   */
  FastParticle(const int particletype, const int incoming, const double momentum, const double theta, const double phi, 
	       const double hard_scale, const double Gamma, string dir);
  /*! Constructor
   * \param particletype which particle [0-3] = [proton, neutron, pi+, pi-] can be extended of course
   * \param incoming is it the beam hadron particle?
   * \param pvec [MeV] 3-vector with momentum
   * \param hard_scale [GeV^2] for CT calculations, the hard scale associated with the particle (like Q^2 or |t|, etc)
   * \param Gamma [MeV] Decay width in rest frame (main branching assumed)
   * \param dir homedir where all input is located
   */
  FastParticle(const int particletype, const int incoming, const TVector3 &pvec, 
	       const double hard_scale, const double Gamma, string dir);
  /*! Copy Constructor
   * \param Copy copy
   */
  FastParticle(const FastParticle &Copy);
  ~FastParticle(); /*!< Destructor */
  FastParticle& operator=(const FastParticle&);/*!< assignment overloading */

  int getParticletype() const; /*!< get type of particle */
  bool getIncoming() const;/*!<  get if the particle is a beam particle */
  double getP() const;/*!< [MeV] get momentum */
  double getTheta() const;/*!< [rad] get theta angle of momentum */
  double getCosTheta() const{return cos(theta);}/*!< [] get cos theta angle of momentum */
  double getPhi() const;/*!<  [rad] get phi angle of momentum */
  double getEx() const;/*!<  get x-component of momentum direction */
  double getEy() const;/*!< get y-component of momentum direction  */
  double getEz() const;/*!<  get z-component of momentum direction */
  double getHitz() const;/*!<  get z coord relative to this momentum of the hard interaction coordinate */
  double getHitbnorm() const;/*!< get norm of b coord relative to this momentum of the hard interaction coordinate */
  const double* getHitbvec() const;/*!<  get b vector in regular coord system for hard interaction coord*/
  /*! Set (b,z) coordinates of the hard interaction point relative to this particle
   * \param r [fm] r coord of hard interaction point, spherical
   * \param costheta cos(theta) coord of hard interaction point, spherical
   * \param sintheta sin(theta) coord of hard interaction point, spherical
   * \param cosphi cos(phi) coord of hard interaction point, spherical
   * \param sinphi cos(phi) coord of hard interaction point, spherical
   */
  void setHitcoord(double r, double costheta, double sintheta, double cosphi, double sinphi); 
  /*!  returns the z coord relative to this particle's momentum of a coordinate
   * \param r [fm] r coord , spherical
   * \param costheta cos(theta) coord , spherical
   * \param sintheta sin(theta) coord , spherical
   * \param cosphi cos(phi) coord , spherical
   * \param sinphi cos(phi) coord , spherical
   * \return z(r), relative to this momentum
   */
  double calcZ(double r, double costheta, double sintheta, double cosphi, double sinphi);
  double getSigmap() const; /*!< [fm^2] return sigma parameter for scattering with proton */
  double getBeta2p() const;  /*!< [fm^2] return beta^2 parameter for scattering with proton */
  double getEpsilonp() const;  /*!< [] return epsilon parameter for scattering with proton */
  double getSigman() const; /*!< [fm^2] return sigma parameter for scattering with neutron */
  double getBeta2n() const;  /*!< [fm^2] return beta^2 parameter for scattering with neutron */
  double getEpsilonn() const; /*!< [] return epsilon parameter for scattering with neutron */
  double getSigma(bool proton) const;/*!< [fm^2] return sigma parameter for scattering with proton(1) or neutron(0) \param proton selects proton or neutron */
  double getEpsilon(bool proton) const;/*!< [] return epsilon parameter for scattering with proton(1) or neutron(0) \param proton selects proton or neutron */
  double getBetasq(bool proton) const;/*!< [fm^2] return beta^2 parameter for scattering with proton(1) or neutron(0) \param proton selects proton or neutron */
  double getMass() const{ return mass;}
  double getE() const{return E;}
  double getDecay_dil() const{return decay_dil;}
  double getSigma_decay_p() const{ return sigma_decay_p;}
  double getSigma_decay_n() const{ return sigma_decay_n;}
  double getSigma_decay(bool proton) const{ return proton? sigma_decay_p : sigma_decay_n;}
  
  
  /*! [fm^2] return sigma parameter for scattering with nucleon from specified level
   * \param level selects level in nucleus
   * \param pnucleus pointer to a class instance of the nucleus 
   * \return [fm^2] sigma parameter
   */
  double getSigma(int level, MeanFieldNucleus *pnucleus) const;
  /*! [fm^2] return sigma parameter with decay products for scattering with nucleon from specified level
   * \param level selects level in nucleus
   * \param pnucleus pointer to a class instance of the nucleus 
   * \return [fm^2] sigma parameter
   */
  double getSigma_decay(int level, MeanFieldNucleus *pnucleus) const;
  /*! [] return epsilon parameter for scattering with nucleon from specified level
   * \param level selects level in nucleus
   * \param pnucleus pointer to a class instance of the nucleus 
   * \return [fm^2] sigma parameter
   */
  double getEpsilon(int level, MeanFieldNucleus *pnucleus) const;
  /*! [fm^2] return beta^2 parameter for scattering with nucleon from specified level
   * \param level selects level in nucleus
   * \param pnucleus pointer to a class instance of the nucleus 
   * \return [fm^2] sigma parameter
   */
  double getBetasq(int level, MeanFieldNucleus *pnucleus) const;
  /*! [fm^2] return scattering parameter front factor (\sigma(1-I*eps)/(4\pi\beta^2)
   * \param level selects level in nucleus
   * \param pnucleus pointer to a class instance of the nucleus 
   * \return [fm^2] frontfactor scattering
   */
  complex<double> getScatterfront(int level, MeanFieldNucleus *pnucleus) const;
  /*! [fm^2] return scattering parameter front factor (\sigma(1-I*eps)/(4\pi\beta^2)
   * \param proton scattering with proton (1) or neutron (0)
   * \return [fm^2] frontfactor scattering
   */
  complex<double> getScatterfront(bool proton) const;
  /*! for CT calc, return ratio between sigma_eff and regular sigma
   * \param zmom distance along path of outgoing particle
   * \return ratio
   */
  double getCTsigma(double zmom) const;

  /*! for CT calc including decay, return ratio between sigma_eff and regular sigma
   * \param zmom distance along path of outgoing particle
   * \param level selects level in nucleus
   * \param pnucleus pointer to a class instance of the nucleus 
   * \return ratio
   */
  double getCT_decay_sigma(double zmom, int level, MeanFieldNucleus *pnucleus) const;

  /*! for CT calc including decay, return ratio between sigma_eff and regular sigma
   * \param zmom distance along path of outgoing particle
   * \param proton scattering with proton (1) or neutron (0)
   * \return ratio
   */
  double getCT_decay_sigma(double zmom, bool proton) const;
  /*! for CT calc including decay, return ratio between sigma_eff and regular sigma
   * \param zmom distance along path of outgoing particle
   * \param level selects level in nucleus
   * \param pnucleus pointer to a class instance of the nucleus 
   * \return ratio
   */
  double getDecay_sigma(double zmom, int level, MeanFieldNucleus *pnucleus) const;
  /*! for CT calc including decay, return ratio between sigma_eff and regular sigma
   * \param zmom distance along path of outgoing particle
   * \param proton scattering with proton (1) or neutron (0)
   * \return ratio
   */
  double getDecay_sigma(double zmom, bool proton) const;

  double getHardScale() const; /*!< get hard scale associated with particle  */
  double getLc() const; /*!< get coherence length of particle (CT parameter)  */
  double getNkt_sq() const; /*!< [GeV^2] get <n_kt^2>, 0.35^2n^2 with n number of constituent quarks of particle  */
  void printParticle() const; /*!< Prints relevant information for the particle  */
  /*! used to calc (b-b')^2 for the glauberphase in the ref frame of this particle
   * \param r [fm] r coord , spherical
   * \param costheta cos(theta) coord , spherical
   * \param sintheta sin(theta) coord , spherical
   * \param cosphi cos(phi) coord , spherical
   * \param sinphi cos(phi) coord , spherical
   * \param zmom [fm] z coord of the coord relative to this momentum
   * \return (b-b')^2, relative to this momentum, with b' the hard interaction point b' coordinate
   */
  double getBdist(double r, double costheta, double sintheta, double cosphi, double sinphi, double zmom);
  
  void setSigma(double sigma){ sigman=sigmap=sigma/10.; return;}

    
private:
  int particletype; /*!< which particle [0-4] = [proton, neutron, pi+, pi-,rho0] can be extended of course */
  
  bool incoming;  /*!< incoming beam particle? */
  double p; /*!<[MeV] momentum */
  double theta; /*!< [rad]  theta angle of momentum*/
  double phi; /*!< [rad]  phi angle of momentum*/
  double ex;  /*!<  x-component of momentum direction*/
  double ey; /*!<  y-component of momentum direction*/
  double ez; /*!<  z-component of momentum direction*/
  double hitz; /*!<  z coord relative to this momentum of the hard interaction coordinate*/
  double hitbnorm; /*!<  norm of b coord relative to this momentum of the hard interaction coordinate*/
  double hitb[3]; /*!<  b vector in regular coord system for hard interaction coord*/
  double beta2p; /*!< [fm^2]  beta^2 parameter for scattering with proton*/
  double sigmap; /*!<  [fm^2]  sigma parameter for scattering with proton*/
  double epsilonp; /*!< []  epsilon parameter for scattering with proton*/
  double beta2n; /*!< [fm^2]  beta^2 parameter for scattering with neutron*/
  double sigman; /*!< [fm^2]  sigma parameter for scattering with neutron*/
  double epsilonn; /*!< []  epsilon parameter for scattering with neutron*/
  double hardscale; /*!<  hard scale associated with particle */
  double nkt_sq; /*!< [GeV^2]  <n_kt^2>, 0.35^2n^2 with n number of constituent quarks of particle*/
  double lc; /*!<  coherence length of particle (CT parameter) */
  double mass; /*!<  [MeV] mass */
  double E; /*!<  [MeV] Energy */
  double decay_dil; /*!<  [MeV] dilated decay width */
  double sigma_decay_p; /*!<  [fm^2]  sigma parameter for decay products scattering with proton*/
  double sigma_decay_n; /*!<  [fm^2]  sigma parameter for decay products scattering with neutron*/
  
  /*! Sets the glauber parameters for a nucleon fast particle
   * \param sigmap [fm^2]  sigma parameter for scattering with proton
   * \param beta2p [fm^2]  beta^2 parameter for scattering with proton
   * \param epsp []  epsilon parameter for scattering with proton
   * \param sigman [fm^2]  sigma parameter for scattering with neutron
   * \param beta2n [fm^2]  beta^2 parameter for scattering with neutron
   * \param epsn []  epsilon parameter for scattering with neutron
   */
  void setGlauberParameters(double &sigmap, double &beta2p, double &epsp, double &sigman, double &beta2n, double &epsn);
  /*! Sets the glauber parameters for a pion fast particle
   * \param sigmap [fm^2]  sigma parameter for scattering with proton
   * \param beta2p [fm^2]  beta^2 parameter for scattering with proton
   * \param epsp []  epsilon parameter for scattering with proton
   * \param sigman [fm^2]  sigma parameter for scattering with neutron
   * \param beta2n [fm^2]  beta^2 parameter for scattering with neutron
   * \param epsn []  epsilon parameter for scattering with neutron
   * \param dir string that contains the dir where all input is located
   */
  void setPionGlauberData(double &sigmap, double &beta2p, double &epsp, double &sigman, double &beta2n, double &epsn, string dir);
  
};
#endif