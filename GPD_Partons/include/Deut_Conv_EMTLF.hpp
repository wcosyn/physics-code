#ifndef DEUT_CONV_EMTLF
#define DEUT_CONV_EMTLF


#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <TDeuteron.h>
#include <complex>
#include <numint/numint.hpp>
#include <TInterpolatingWavefunction.h>

/**
 * @brief Class implements deuteron LF EMT densities, see draft paper with Adam
 * 
 */
class Deut_Conv_EMTLF{

public:

/**
 * @brief constructor
 * 
 * @param wfname deuteron wave function parametrization name, see TDeuteron for possibilities
 */
Deut_Conv_EMTLF(const std::string &wfname);

/**
 * @brief Destructor
 * 
 */
~Deut_Conv_EMTLF();


/**
 * @brief returns the melosh rotated lf wf, melosh rotation only acting on the first (active) nucleon, since the helicities of the spectator are summed over
 * 
 * @param Ek lf rel momentum on-shell energy
 * @param k lf rel momentum
 * @param nonrelwf deuteron light-front wf without the Melosh rotations implemented.  Helicites are as follows [0] down down, [1] up down, [2] down up, [3] up up
 * @return std::vector< std::complex<double> >  deuteron light-front wf with the Melosh rotation for the active nucleon implemented.  Helicites are as follows [0] down down, [1] up down, [2] down up, [3] up up
 */

static std::vector< std::complex<double> > lf_deut(const double Ek, const TVector3& k, const std::vector< std::complex<double> > &nonrelwf);



/**
 * @brief Computes the deuteron helicity amplitudes with the convolution formula, xi=0 due to LF choice
 * 
 * @param t [Gev^2] momentum transfer sq
 * @return std::vector< std::complex<double> > helicity amplitudes  deuteron helicities are (final,initial) [0]++,[1]00,[2]0+,[3]+0,[4]-+ 
 */
std::vector< std::complex<double> > EMT_conv(const double t);

/**
 * @brief Computes the deuteron helicity amplitudes with the convolution formula, xi=0 due to LF choice
 * 
 * @param t [Gev^2] momentum transfer sq
 * @return std::vector< std::complex<double> > helicity amplitudes  deuteron helicities are (final,initial) [0]++,[1]00,[2]0+,[3]+0,[4]-+ 
 */
std::vector< double > EMT_conv_real(const double t);

TInterpolatingWavefunction * getWf(){ return &wf;} ///< return deuteron wf object


/**
 * @brief computes nucleon matrix elements for chiral even vector GPDs.  
 * (\bar{u}(p') \gamma^+ u(p))/2P^+ etc for vector GPDs
 * (\bar{u}(p') \gamma^+gamma^5 u(p))/2P^+ etc for axial GPDs  [Diehl conventions, but shouldn't matter]
 * 
 * @param sigma_in polarization incoming nucleon (spin times two!!)
 * @param sigma_out polarization outgoing nucleon (spin times two!!)
 * @param xi_n skewness nucleon
 * @param t [MeV^2] momentum transfer sq.
 * @param t0 [MeV^2] minimum momentum transfer sq
 * @param phi [azimuthal angle \xiP+Delta fourvector]
 * @param gpd_H GPD H or Htilde
 * @param gpd_E GPD E or Etilde
 * @param [1] vector [0] axial GPDS
 * @return std::complex<double> nucleon helicity amplitude
 */
// static std::complex<double> getGPD_even_nucl(const int sigma_in, const int sigma_out, const double xi_n, 
//                                     const double t, const double t0, const double phi,
//                                     const double gpd_H, const double gpd_E, const bool vector);              




private: 
/**
 * @brief structure needed to carry out the convolution integration, contains parameters and integrandum function
 * 
 */
struct Ftor_conv {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_conv &p = * (Ftor_conv *) param;
      p.f(ret,x[0],x[1],x[2], *p.emt,p.t,p.pold_in, p.pold_out,p.deltax, p.GFF_A, p.GFF_J, p.GFF_D);
    }
    Deut_Conv_EMTLF *emt;
    double t; // [MeV^2]
    double scale; ///< [GeV]
    int pold_in;
    int pold_out;
    double deltax; //[MeV]
    double GFF_A, GFF_J, GFF_D;
    
    void (*f)(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_EMTLF &emt,
              double t, int pold_in, int pold_out,  double deltax, double GFF_A, double GFF_J, double GFF_D);
      };


/**
 * @brief structure needed to carry out the convolution integration, contains parameters and integrandum function
 * 
 */
struct Ftor_conv_real{

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor_conv_real &p = * (Ftor_conv_real *) param;
      p.f(ret,x[0],x[1],x[2], *p.emt,p.t,p.pold_in, p.pold_out,p.deltax, p.GFF_A, p.GFF_J, p.GFF_D);
    }
    Deut_Conv_EMTLF *emt;
    double t; // [MeV^2]
    double scale; ///< [GeV]
    int pold_in;
    int pold_out;
    double deltax; //[MeV]
    double GFF_A, GFF_J, GFF_D;
    
    void (*f)(numint::vector_d & res, double alpha_1, double kperp, double kphi, Deut_Conv_EMTLF &emt,
              double t, int pold_in, int pold_out,  double deltax, double GFF_A, double GFF_J, double GFF_D);
      };



/**
 * @brief integrandum for the convolution integral, integration variables are alpha and k_perp of the initial nucleon
 * 
 * @param[out] res integral result 
 * @param alpha_1 [] lc momentum fraction of initial nucleon
 * @param kperp [MeV] perp lf momentum entering in initial deuteron wf
 * @param kphi [] azimuthal angle of lf momentum entering in initial deuteron wf
 * @param[in] gpd object that contains all necessary info on the gpds 
 * @param t [MeV^2] momentum transfer sq
 * @param pold_in polarization initial deuteron state (-1//0//+1)
 * @param pold_out polarization final deuteron state (-1//0//+1)
 * @param right [1] R matrix element [0] L matrix element
 * @param deltax [MeV!!!] perp component of the momentum transfer, along x-axis
 */
static void int_k3(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_EMTLF &emt,
              double t, int pold_in, int pold_out, double deltax, double GFF_A, double GFF_J, double GFF_D);


/**
 * @brief integrandum for the convolution integral, integration variables are alpha and k_perp of the initial nucleon
 * 
 * @param[out] res integral result 
 * @param alpha_1 [] lc momentum fraction of initial nucleon
 * @param kperp [MeV] perp lf momentum entering in initial deuteron wf
 * @param kphi [] azimuthal angle of lf momentum entering in initial deuteron wf
 * @param[in] gpd object that contains all necessary info on the gpds 
 * @param t [MeV^2] momentum transfer sq
 * @param pold_in polarization initial deuteron state (-1//0//+1)
 * @param pold_out polarization final deuteron state (-1//0//+1)
 * @param right [1] R matrix element [0] L matrix element
 * @param deltax [MeV!!!] perp component of the momentum transfer, along x-axis
 */
static void int_k3_real(numint::vector_d & res, double alpha_1, double kperp, double kphi, Deut_Conv_EMTLF &emt,
              double t, int pold_in, int pold_out, double deltax, double GFF_A, double GFF_J, double GFF_D);


/**
 * @brief compute the EMT nucleon matrix elements (matrix elements Eq 28 Freese Miller PRD103)
 * 
 * @param sigma_in 
 * @param sigma_out 
 * @param alpha_1 
 * @param p1R [MeV]
 * @param p1L [MeV]
 * @param t [MeV^2]
 * @param GFF_A
 * @param GFF_J
 * @param GFF_D
 * @return vector< complex<double> > 
 */
static std::vector< std::complex<double> > getGFF_nucl(const int sigma_in, const int sigma_out, 
                        const double alpha_1, const std::complex<double> p1R, const std::complex<double> p1L, const double t,
                        const double GFF_A, const double GFF_J, const double GFF_D);

TDeuteron::Wavefunction *wfref; /*!< contains instance of deuteron wave function*/
TInterpolatingWavefunction wf;


};

#endif 