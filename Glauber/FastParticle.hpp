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
  complex<double> getScatterfront(int level, MeanFieldNucleus *pnucleus) const;
  /*! [fm^2] return scattering parameter front factor (\f$\frac{\sigma(1-I\epsilon)}{4\pi\beta^2}\f$)
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
   */
  void setScatter(double sigma, double beta, double eps);
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
   * \param dir string that contains the dir where all input is located
   */
  static void setPionGlauberData(double mom, double &sigmap, double &beta2p, 
				 double &epsp, double &sigman, double &beta2n, double &epsn, string dir);
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

  
};

static double sigmap_array[] = { 17.2781, 16.306, 15.2518, 14.1561, 13.0553, 11.9802, 10.9544, 9.99498, 9.11227, 8.31128,
				    7.59246, 6.95286, 6.3871, 5.8882, 5.44828, 5.05905, 4.71225, 4.39995, 4.1149, 3.85078,
				    3.60257, 3.36681, 3.14198, 2.92868, 2.72983, 2.55078, 2.3995, 2.28718, 2.21259, 2.07051,
				    1.93132, 1.80425, 1.69566, 1.60908, 1.54576, 1.50547, 1.48705, 1.48886, 1.50895, 1.54515,
				    1.59509, 1.65623, 1.72583, 1.80112, 1.87943, 1.95834, 2.03596, 2.11099, 2.18295, 2.2521,
				    2.31941, 2.38646, 2.4552, 2.52776, 2.60623, 2.69242, 2.7877, 2.89274, 3.00733, 3.1302,
				    3.25887, 3.38956, 3.51728, 3.63598, 3.73898, 3.81961, 3.87195, 3.89164, 3.8766, 3.82745,
				    3.74765, 3.64312, 3.52165, 3.39199, 3.26299, 3.14286, 3.03859, 2.95563, 2.89769, 2.86662,
				    2.86223, 2.88189, 2.91401, 2.94705, 2.97388, 2.99496, 3.01079, 3.02185, 3.02861, 3.03154,
				    3.03106, 3.02761, 3.02157, 3.0133, 3.00315, 2.99142, 2.9784, 2.96433, 2.94946, 2.93398,
				    2.91807, 2.90191, 2.88562, 2.86933, 2.85315, 2.83717, 2.82146, 2.80608, 2.7911, 2.77654,
				    2.76244, 2.74883, 2.73572, 2.72312, 2.71103, 2.69947, 2.68841, 2.67785, 2.66779, 2.6582,
				    2.64907, 2.64038, 2.63212, 2.62425, 2.61676, 2.60963, 2.60282, 2.59633, 2.59013, 2.58418,
				    2.57848, 2.573, 2.56771, 2.5626, 2.55764, 2.55282, 2.54812, 2.54353, 2.53902, 2.53457,
				    2.53019, 2.52585, 2.52154, 2.51725, 2.51297, 2.50869, 2.50442, 2.50013, 2.49582, 2.4915,
				    2.48716} ;

				    
static double beta2p_array[] =	{ -0.0485909, 0.143607, 0.286071, 0.387333, 0.454878, 0.495238, 0.51407, 0.516242, 0.505898, 0.486539,
0.46108, 0.431917, 0.400986, 0.369813, 0.339567, 0.311108, 0.285027, 0.261692, 0.241277, 0.223805,
0.209168, 0.197164, 0.187519, 0.179907, 0.173971, 0.169342, 0.165654, 0.162554, 0.159713, 0.15684,
0.153682, 0.150031, 0.145731, 0.140676, 0.134811, 0.128134, 0.120689, 0.11257, 0.10391, 0.0948794,
0.0856826, 0.076548, 0.0677235, 0.0594694, 0.0520508, 0.0457307, 0.0407623, 0.037382, 0.0358023, 0.036205,
0.0387351, 0.0434948, 0.0505389, 0.0598702, 0.0714361, 0.0851269, 0.100774, 0.118149, 0.136968, 0.156892,
0.177532, 0.198456, 0.219196, 0.239257, 0.258134, 0.275319, 0.290324, 0.302695, 0.312037, 0.318036,
0.320487, 0.319323, 0.314648, 0.306772, 0.296252, 0.28393, 0.259021, 0.251092, 0.245262, 0.241262,
0.23884, 0.237765, 0.237824, 0.238821, 0.240579, 0.242935, 0.245743, 0.248871, 0.252202, 0.255632,
0.259069, 0.262434, 0.26566, 0.268689, 0.271475, 0.27398, 0.276176, 0.278042, 0.279565, 0.28074,
0.281567, 0.282052, 0.282208, 0.28205, 0.281599, 0.28088, 0.27992, 0.278749, 0.2774, 0.275906,
0.274303, 0.272626, 0.270914, 0.269202, 0.267527, 0.265925, 0.264431, 0.263077, 0.261896, 0.260918,
0.260171, 0.25968, 0.259467, 0.259552, 0.259954, 0.260685, 0.261755, 0.263173, 0.264941, 0.267061,
0.269528, 0.272335, 0.275473, 0.278927, 0.282679, 0.28671, 0.290994, 0.295505, 0.300213, 0.305084,
0.310085, 0.315177, 0.320321, 0.325476, 0.3306, 0.33565, 0.340583, 0.345354, 0.349921, 0.354241,
0.358274} ;

static double epsp_array[] = { 40.2787, 27.7732, 18.5131, 11.7765, 6.97527, 3.63519, 1.37815, -0.0934038, -1.01015, -1.54774,
-1.83704, -1.973, -2.02214, -2.0289, -2.02098, -2.01371, -2.01368, -2.02155, -2.0344, -2.04739,
-2.05512, -2.05251, -2.0354, -2.00089, -1.94752, -1.87523, -1.78524, -1.67986, -1.5622, -1.43595,
-1.30504, -1.17342, -1.0448, -0.922449, -0.809057, -0.706621, -0.616385, -0.538837, -0.473743, -0.420227,
-0.376878, -0.341891, -0.313216, -0.288726, -0.266379, -0.244371, -0.221272, -0.196136, -0.168578, -0.138811,
-0.107641, -0.0764157, -0.0469264, -0.0212679, -0.00165783, 0.00977428, 0.0112181, 0.00141802, -0.0201046, -0.0528668,
-0.0953037, -0.144745, -0.197527, -0.249278, -0.295439, -0.332052, -0.350241, -0.371676, -0.38652, -0.393895,
-0.39396, -0.387572, -0.376013, -0.360774, -0.343385, -0.325295, -0.307782, -0.291901, -0.278453, -0.267978,
-0.260759, -0.256849, -0.256101, -0.258199, -0.262705, -0.269098, -0.276811, -0.285273, -0.29394, -0.302321,
-0.310004, -0.316671, -0.322104, -0.32619, -0.328917, -0.330365, -0.33069, -0.330105, -0.328861, -0.327221,
-0.325437, -0.327752, -0.325751, -0.323747, -0.32174, -0.319731, -0.31772, -0.315706, -0.313689, -0.31167,
-0.309648, -0.307624, -0.305597, -0.303568, -0.301536, -0.299502, -0.297465, -0.295426, -0.293384, -0.29134,
-0.289293, -0.287243, -0.285191, -0.283137, -0.28108, -0.27902, -0.276958, -0.274893, -0.272826, -0.270756,
-0.268684, -0.266609, -0.264532, -0.262452, -0.26037, -0.258285, -0.256197, -0.254107, -0.252015, -0.24992,
-0.247822, -0.245722, -0.24362, -0.241515, -0.239407, -0.237297, -0.235184, -0.233069, -0.230951, -0.228831,
-0.226708} ;

static double sigman_array[] = { 6.01242, 5.70712, 5.38087, 5.04634, 4.71471, 4.39525, 4.09517, 3.81963, 3.57196, 3.35396,
3.16623, 3.00846, 2.87971, 2.77865, 2.70364, 2.65289, 2.62443, 2.61612, 2.62557, 2.65014,
2.68683, 2.73236, 2.78323, 2.83602, 2.88772, 2.93638, 2.98186, 3.02666, 3.0771, 3.14463,
3.24794, 3.51793, 3.75167, 4.06071, 4.35431, 4.55696, 4.62355, 4.54971, 4.36772, 4.13054,
3.8935, 3.70243, 3.58992, 3.57727, 3.67775, 3.89777, 4.23356, 4.66272, 5.13272, 5.55444,
5.81446, 5.81624, 5.53835, 5.02546, 4.57165, 4.21632, 3.96029, 3.78954, 3.68587, 3.63125,
3.60932, 3.60609, 3.61021, 3.61323, 3.60977, 3.5973, 3.57588, 3.54748, 3.5153, 3.48311,
3.45458, 3.43288, 3.42036, 3.41836, 3.42719, 3.44609, 3.47332, 3.50628, 3.54167, 3.57582,
3.60495, 3.62559, 3.63497, 3.63136, 3.61435, 3.58487, 3.54516, 3.49847, 3.44867, 3.39982,
3.35566, 3.31922, 3.29238, 3.27558, 3.26757, 3.26535, 3.26429, 3.25874, 3.23093, 3.21291,
3.19527, 3.17801, 3.16112, 3.14459, 3.12842, 3.11259, 3.0971, 3.08194, 3.0671, 3.05258,
3.03836, 3.02445, 3.01083, 2.99749, 2.98444, 2.97166, 2.95915, 2.9469, 2.93491, 2.92316,
2.91167, 2.90041, 2.88939, 2.87859, 2.86802, 2.85767, 2.84753, 2.8376, 2.82788, 2.81836,
2.80903, 2.7999, 2.79095, 2.78219, 2.77361, 2.7652, 2.75696, 2.74889, 2.74099, 2.73325,
2.72566, 2.71823, 2.71095, 2.70382, 2.69683, 2.68998, 2.68328, 2.6767, 2.67026, 2.66395,
2.65777} ;

static double beta2n_array[] = { -1.87567, -0.90376, -0.318568, 0.0128662, 0.184313, 0.259972, 0.282194, 0.277594, 0.26181, 0.243137,
0.225239, 0.209102, 0.194393, 0.180349, 0.166306, 0.151965, 0.137469, 0.12336, 0.110459, 0.0997147,
0.0920482, 0.0882169, 0.0887069, 0.0936658, 0.102877, 0.115773, 0.131486, 0.148926, 0.166879, 0.184116,
0.199505, 0.212112, 0.221296, 0.226768, 0.228634, 0.227401, 0.223952, 0.219488, 0.21544, 0.213362,
0.214793, 0.22112, 0.233436, 0.252405, 0.278151, 0.310178, 0.347331, 0.387808, 0.429233, 0.468787,
0.503414, 0.530072, 0.546046, 0.54929, 0.538782, 0.514848, 0.479421, 0.436178, 0.390469, 0.348979,
0.330735, 0.30365, 0.292228, 0.293364, 0.302434, 0.314748, 0.326398, 0.334669, 0.338128, 0.336513,
0.330494, 0.321373, 0.310778, 0.300387, 0.291688, 0.285819, 0.283458, 0.28479, 0.289526, 0.296978,
0.306158, 0.31591, 0.325046, 0.33248, 0.337336, 0.339045, 0.337397, 0.332561, 0.325063, 0.315732,
0.305605, 0.29582, 0.287483, 0.281535, 0.278633, 0.27905, 0.282615, 0.288696, 0.296247, 0.303914,
0.3102, 0.313691, 0.313314, 0.30861, 0.299981, 0.288849, 0.277669, 0.269695, 0.268391, 0.27635,
0.293187, 0.292679, 0.292596, 0.292911, 0.293596, 0.294623, 0.295967, 0.297603, 0.299505, 0.301651,
0.304016, 0.306579, 0.309318, 0.312212, 0.315241, 0.318385, 0.321626, 0.324946, 0.328327, 0.331753,
0.335209, 0.338679, 0.342149, 0.345605, 0.349035, 0.352428, 0.35577, 0.359053, 0.362265, 0.365398,
0.368444, 0.371395, 0.374245, 0.376986, 0.379614, 0.382123, 0.384511, 0.386773, 0.388908, 0.390914,
0.392789} ;

static double epsn_array[] = { 0.220782, 0.39322, 0.489897, 0.528353, 0.524508, 0.492505, 0.44458, 0.390976, 0.339892, 0.297477,
0.267874, 0.253323, 0.254303, 0.269736, 0.297233, 0.333387, 0.374099, 0.414926, 0.45144, 0.479586,
0.496009, 0.498344, 0.485439, 0.457499, 0.416129, 0.364263, 0.305979, 0.246191, 0.190234, 0.143366,
0.110213, 0.0942141, 0.097116, 0.118584, 0.155994, 0.204469, 0.257203, 0.306086, 0.342616, 0.359026,
0.349499, 0.311294, 0.245565, 0.157647, 0.0566115, -0.0459978, -0.138182, -0.209941, -0.255599, -0.275228,
-0.274377, -0.261664, -0.244639, -0.22576, -0.174621, -0.179313, -0.159157, -0.135739, -0.118666, -0.11035,
-0.109356, -0.112648, -0.116993, -0.119757, -0.119255, -0.11482, -0.106675, -0.0957192, -0.0832652, -0.0707917,
-0.0597236, -0.0512651, -0.0462866, -0.0452688, -0.0482955, -0.0550886, -0.0650731, -0.0774616, -0.0913482, -0.105801,
-0.119946, -0.133034, -0.144492, -0.153949, -0.161242, -0.166403, -0.169629, -0.171237, -0.171621, -0.171194,
-0.170344, -0.169394, -0.168575, -0.168015, -0.167737, -0.16768, -0.167728, -0.167745, -0.167615, -0.167279,
-0.166756, -0.166148, -0.165627, -0.165385, -0.165576, -0.166225, -0.167154, -0.161663, -0.161403, -0.161117,
-0.160805, -0.160468, -0.160105, -0.159718, -0.159306, -0.158869, -0.158409, -0.157924, -0.157416, -0.156884,
-0.156328, -0.15575, -0.155149, -0.154526, -0.15388, -0.153212, -0.152522, -0.151811, -0.151078, -0.150325,
-0.14955, -0.148755, -0.14794, -0.147104, -0.146248, -0.145373, -0.144479, -0.143565, -0.142632, -0.141681,
-0.140711, -0.139723, -0.138717, -0.137694, -0.136653, -0.135594, -0.134519, -0.133427, -0.132319, -0.131194,
-0.130053} ;



/** @} */
#endif