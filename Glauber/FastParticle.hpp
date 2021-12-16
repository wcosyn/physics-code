/*! \file FastParticle.hpp 
 * \brief Contains FastParticle class declaration
 * \author Wim Cosyn
 * \date 16/08/2011
 * \addtogroup Glauber
 * @{
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
   * \param type which particle [0-4] can be extended of course <BR>
   * 0: proton <BR>
   * 1: neutron <BR>
   * 2: pi+ <BR>
   * 3: pi- <BR>
   * 4: rho0 <BR>
   * 5: rho0 with more CT effects <BR>
   * 6: rho0 with less CT effects <BR>
   * 7: double pion along one momentum (rho decay simulation) <BR>
   * 8: proton for charge exchange reactions <BR>
   * 9: neutron for charge exchange reactions <BR>
   * \param inc is it the beam hadron particle?
   * \param momentum [MeV] momentum of the particle
   * \param ptheta [rad] spherical coord, theta angle of momentum with the z-axis
   * \param pphi [rad] spherical coord, phi angle of momentum
   * \param hard_scale [GeV^2] for CT calculations, the hard scale associated with the particle (like Q^2 or |t|, etc)
   * \param Gamma [MeV] Decay width in rest frame (main branching assumed)
   * \param lc_mod [] modification of coherence length in CT 
   * \param nkt_mod [] modification of initial CT sigma value
   * \param dir homedir where all input is located
   */
  FastParticle(const int type, const int inc, const double momentum,
			   const double ptheta, const double pphi, const double hard_scale, const double Gamma, 
         const double lc_mod, const double nkt_mod, const std::string dir);
  /*! Constructor
   * \param type which particle [0-4] can be extended of course <BR>
   * 0: proton <BR>
   * 1: neutron <BR>
   * 2: pi+ <BR>
   * 3: pi- <BR>
   * 4: rho0 <BR>
   * 5: rho0 with more CT effects <BR>
   * 6: rho0 with less CT effects <BR>
   * 7: double pion along one momentum (rho decay simulation) <BR>
   * 8: proton for charge exchange reactions <BR>
   * 9: neutron for charge exchange reactions <BR>
   * \param inc is it the beam hadron particle?
   * \param pvec [MeV] 3-vector with momentum
   * \param hard_scale [GeV^2] for CT calculations, the hard scale associated with the particle (like Q^2 or |t|, etc)
   * \param Gamma [MeV] Decay width in rest frame (main branching assumed)
   * \param lc_mod [] modification of coherence length in CT 
   * \param nkt_mod [] modification of initial CT sigma value
   * \param dir homedir where all input is located
   */
  FastParticle(const int type, const int inc, const TVector3 &pvec, 
			   const double hard_scale, const double Gamma, 
         const double lc_mod, const double nkt_mod, const std::string dir);
  /*! Copy Constructor
   * \param Copy copy
   */
  FastParticle(const FastParticle &Copy);
  ~FastParticle(); /*!< Destructor */
  FastParticle& operator=(const FastParticle&);/*!< assignment overloading */
  enum Particletype { P_EL=0    , N_EL=1, 
                      PI_PLUS=2 , PI_MIN=3, 
                      RHO0=4    , RHO0_MCT=5,
                      RHO0_LCT=6, DBL_PI=7,
                      P_SCX=8   ,N_SCX=9,
                      P_CLASS_SCX=10, N_CLASS_SCX=11 }; /*!< enum for particle types so that everythin a bit more readable */
  
  int getParticletype() const{ return particletype;} /*!< get type of particle */
  bool getIncoming() const {return incoming;}/*!<  get if the particle is a beam particle */
  double getP() const {return p;}/*!< [MeV] get momentum */
  double getTheta() const {return theta;}/*!< [rad] get theta angle of momentum */
  double getCosTheta() const{return costheta;}/*!< [] get cos theta angle of momentum */
  double getPhi() const {return phi;}/*!<  [rad] get phi angle of momentum */
  double getEx() const {return ex;}/*!<  get x-component of momentum direction */
  double getEy() const {return ey;}/*!< get y-component of momentum direction  */
  double getEz() const {return ez;}/*!<  get z-component of momentum direction */
  double getHitz() const {return hitz;}/*!<  get z coord relative to this momentum of the hard interaction coordinate */
  double getHitbnorm() const {return hitbnorm;}/*!< get norm of b coord relative to this momentum of the hard interaction coordinate */
  const double* getHitbvec() const {return hitb;}/*!<  get b vector in regular coord system for hard interaction coord*/
  TVector3 getPvec() const{ return p*TVector3(ex,ey,ez);} /*!<returns the momentumvector */
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
  double getSigmap() const {return sigmap;} /*!< [fm^2] return sigma parameter for scattering with proton */
  double getBeta2p() const {return beta2p;}  /*!< [fm^2] return beta^2 parameter for scattering with proton */
  double getEpsilonp() const {return epsilonp;}  /*!< [] return epsilon parameter for scattering with proton */
  double getSigman() const {return sigman;} /*!< [fm^2] return sigma parameter for scattering with neutron */
  double getBeta2n() const {return beta2n;}  /*!< [fm^2] return beta^2 parameter for scattering with neutron */
  double getEpsilonn() const {return epsilonn;} /*!< [] return epsilon parameter for scattering with neutron */
  /*! [fm^2] return sigma parameter for scattering with proton(1) or neutron(0) \param proton selects proton or neutron */
  double getSigma(bool proton) const{return proton? sigmap:sigman;}
  /*! [] return epsilon parameter for scattering with proton(1) or neutron(0) \param proton selects proton or neutron */
  double getEpsilon(bool proton) const {return proton? epsilonp:epsilonn;}
  /*! [fm^2] return beta^2 parameter for scattering with proton(1) or neutron(0) \param proton selects proton or neutron */
  double getBetasq(bool proton) const {return proton? beta2p:beta2n; }
  double getMass() const{ return mass;} /*!< [MeV] return mass */
  double getE() const{return E;} /*!< [MeV] on-shell energy*/
  double getDecay_dil() const{return decay_dil;} /*!< [MeV] return time dilatation factor */
  /*! [fm^2] return sigma parameter corrected for decaying particle for scattering with proton */
  double getSigma_decay_p() const{ return sigma_decay_p;}
  /*! [fm^2] return sigma parameter corrected for decaying particle for scattering with neutron */
  double getSigma_decay_n() const{ return sigma_decay_n;}
  /*! [fm^2] return sigma (corrected for decay) parameter for scattering with proton(1) or neutron(0) \param proton selects proton or neutron */
  double getSigma_decay(bool proton) const{ return proton? sigma_decay_p : sigma_decay_n;}
  double get_lc_mod(){return lc_mod;}
  double get_nkt_mod(){return nkt_mod;}
  
  
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
  /*! [fm^2] return scattering parameter front factor (\f$\frac{\sigma(1-I\epsilon)}{4\pi\beta^2}\f$)
   * \param level selects level in nucleus
   * \param pnucleus pointer to a class instance of the nucleus 
   * \return [fm^2] frontfactor scattering
   */
  std::complex<double> getScatterfront(int level, MeanFieldNucleus *pnucleus) const;
  /*! [fm^2] return scattering parameter front factor (\f$\frac{\sigma(1-I\epsilon)}{4\pi\beta^2}\f$)
   * \param proton scattering with proton (1) or neutron (0)
   * \return [fm^2] frontfactor scattering
   */
  std::complex<double> getScatterfront(bool proton) const;

 /*! gives you the Feynman scatter amplitude of the final-state interaction
   * \param t [MeV^2] momentum transfer squared
   * \param proton scattering with proton (1) or neutron (0)
   * \return [MeV^-2] \f$ \sigma_{tot} (I+\epsilon) e^{\beta t/2} \f$
   */
  std::complex<double> scatter(double t, bool proton) const{
    return getSigma(proton)*(I_UNIT+getEpsilon(proton))*exp(getBetasq(proton)*t/2.*INVHBARC*INVHBARC)*INVHBARC*INVHBARC; 
  }


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

  double getHardScale() const {return hardscale;} /*!< get hard scale associated with particle  */
  double getLc() const {return lc;} /*!< get coherence length of particle (CT parameter)  */
  double getNkt_sq() const {return nkt_sq;} /*!< [GeV^2] get <n_kt^2>, 0.35^2n^2 with n number of constituent quarks of particle  */
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
  
  /*! sets the scatter sigma_tot parameter
   * \param sigma [mb] tot cross section
   */
  void setSigma(double sigma){ sigman=sigmap=sigma/10.; userset=1; return;}
  /*! sets a screening for the sigma parameters
   * \param screening [%] amount of screening (can be anti-screening too of course!)
   */
  void setScreening(double screening){ sigman*=screening; sigmap*=screening; userset=1; return;}
  /*! sets the scatter parameters
   * \param sigma [mb] tot cross section
   * \param beta [GeV^2] slope param
   * \param eps [] real to imag ratio
   * \param p_dil [MeV] momentum used for decay time dilatation
   */
  void setScatter(double sigma, double beta, double eps, double p_dil);
  bool getUserset() const{return userset;} /*!< returns if user set the scattering parameters */
  /*! Sets the glauber parameters for a nucleon fast particle
   * \param mom [MeV]  particle momentum
   * \param sigmap [fm^2]  sigma parameter for scattering with proton
   * \param beta2p [fm^2]  beta^2 parameter for scattering with proton
   * \param epsp []  epsilon parameter for scattering with proton
   * \param sigman [fm^2]  sigma parameter for scattering with neutron
   * \param beta2n [fm^2]  beta^2 parameter for scattering with neutron
   * \param epsn []  epsilon parameter for scattering with neutron
   */
  static void setGlauberParameters(double mom, double &sigmap, double &beta2p, double &epsp, double &sigman, double &beta2n, double &epsn);
  /*! Sets the glauber parameters for a pion fast particle
   * \param mom [MeV]  particle momentum
   * \param sigmap [fm^2]  sigma parameter for scattering with proton
   * \param beta2p [fm^2]  beta^2 parameter for scattering with proton
   * \param epsp []  epsilon parameter for scattering with proton
   * \param sigman [fm^2]  sigma parameter for scattering with neutron
   * \param beta2n [fm^2]  beta^2 parameter for scattering with neutron
   * \param epsn []  epsilon parameter for scattering with neutron
   * \param dir std::string that contains the dir where all input is located
   */
  static void setPionGlauberData(double mom, double &sigmap, double &beta2p, 
				 double &epsp, double &sigman, double &beta2n, double &epsn, std::string dir);
  /*! Sets the glauber parameters for a pion fast particle through interpolation of static arrays
   * \param particletype which particle [0-3] = [proton, neutron, pi+, pi-] can be extended of course
   * \param mom [MeV]  particle momentum
   * \param sigmap [fm^2]  sigma parameter for scattering with proton
   * \param beta2p [fm^2]  beta^2 parameter for scattering with proton
   * \param epsp []  epsilon parameter for scattering with proton
   * \param sigman [fm^2]  sigma parameter for scattering with neutron
   * \param beta2n [fm^2]  beta^2 parameter for scattering with neutron
   * \param epsn []  epsilon parameter for scattering with neutron
   */
  static void interpPionGlauberData(int particletype, double mom, double &sigmap, 
				    double &beta2p, double &epsp, double &sigman, double &beta2n, double &epsn);
    
private:
  int particletype; /*!< which particle [0-4] = [proton, neutron, pi+, pi-,rho0] can be extended of course */
  
  bool incoming;  /*!< incoming beam particle? */
  double p; /*!<[MeV] momentum */
  double theta; /*!< [rad]  theta angle of momentum*/
  double costheta; /*!< []  cos(theta) angle of momentum*/
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
  bool userset; /*!<  if user has set scatt parameters */
  double Gamma; /*!< [MeV] decay width in rest frame */
  double lc_mod; /*!< [] modification factor of coherence length (playing around with parameters) */
  double nkt_mod; /*!< [] modification factor of initial CT sigma value (playing around with parameters) */
  static const double sigmap_array[];				      
  static const double beta2p_array[];
  static const double epsp_array[];
  static const double sigman_array[];
  static const double beta2n_array[];
  static const double epsn_array[];

  
};




/** @} */
#endif
