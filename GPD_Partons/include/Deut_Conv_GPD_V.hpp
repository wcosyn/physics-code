#ifndef DEUT_CONV_GPD_T
#define DEUT_CONV_GPD_T


#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <TDeuteron.h>
#include <complex>
#include <numint/numint.hpp>
#include <TInterpolatingWavefunction.h>


#include <partons/Partons.h>
#include <partons/services/GPDService.h>

#include "Deut_GPD_V_set.hpp"

/**
 * @brief Class implements vector GPDs for the nucleon and deuteron (using convolution formalism) for now
 * 
 */
class Deut_Conv_GPD_V{

public:

/**
 * @brief constructor
 * 
 * @param pGPDService GPD Service from PARTONS
 * @param pGPDModel GPD Module from PARTONS
 * @param wfname deuteron wave function parametrization name, see TDeuteron for possibilities
 */
Deut_Conv_GPD_V(PARTONS::GPDService* pGPDService, PARTONS::GPDModule* pGPDModel, const std::string &wfname);

/**
 * @brief Destructor
 * 
 */
~Deut_Conv_GPD_V();


/**
 * @brief conversion from helicity amplitudes to GPDs for spin 1 chiral even vector quark GPDs.  We compute in a frame where phi=0
 * See Cano Pire EPJA App A
 * 
 * @param xi [] skewness
 * @param t [MeV^2] momentum transfer sq
 * @param helamps helicity amplitudes, deuteron helicities are (final,initial) [0]++,[1]00,[2]0+,[3]+0,[4]-+
 * @return std::vector<double>  gpds H_1 to H_5
 */
static std::vector< std::complex<double> > helamps_to_gpds_V(const double xi, const double t, const std::vector< std::complex<double> > & helamps);

/**
 * @brief conversion from GPDs to helicity amplitudes for spin 1 chiral even vector quark GPDs.  We compute in a frame where phi=0
 * See Cano Pire EPJA App A
 * 
 * @param xi [] skewness
 * @param t [MeV^2] momentum transfer sq
 * @param gpds gpds H_1 to H_5
 * @return std::vector<double> helicity amplitudes, deuteron helicities are (final,initial) [0]++,[1]00,[2]0+,[3]+0,[4]-+
 */
static std::vector< std::complex<double> > gpds_to_helamps_V(const double xi, const double t, const std::vector< std::complex<double> > & gpds);

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
 * @brief Computes the deuteron helicity amplitudes with the convolution formula
 * 
 * @param xi [] skewness
 * @param x [] parton average lf momentum fraction
 * @param t [Mev^2] momentum transfer sq
 * @return std::vector< std::complex<double> > helicity amplitudes  deuteron helicities are (final,initial) [0]++,[1]--,[2]00,[3]0+,[4]-0,[5]+0,[6]0-,[7]-+,[8]+-
 */
std::vector< std::complex<double> > gpd_conv(const double xi, const double x, const double t);

TInterpolatingWavefunction * getWf(){ return &wf;} ///< return deuteron wf object

/**
 * @brief Get a Deut_GPD_V_set object containing deuteron vector helicity amplitudes
 * 
 * @param x [] avg lf momentum fraction
 * @param xi [] skewness
 * @param t [GeV^2] mom transfer sq
 * @param ERBL grid for only ERBL region or not?
 * @return Deut_GPD_V_set 
 */
Deut_GPD_V_set getDeut_GPD_V_set(const double x, const double xi, const double t, const bool ERBL);

private: 
/**
 * @brief structure needed to carry out the convolution integration, contains parameters and integrandum function
 * 
 */
struct Ftor_conv {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_conv &p = * (Ftor_conv *) param;
      p.f(ret,x[0],x[1],x[2], *p.gpd,p.x,p.xi,p.t,p.pold_in, p.pold_out,p.deltax);
    }
    Deut_Conv_GPD_V *gpd;
    double x;
    double xi;
    double t;
    int pold_in;
    int pold_out;
    double deltax;
    
    void (*f)(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_GPD_V &gpd,
              double x, double xi, double t, int pold_in, int pold_out,  double deltax);
      };


/**
 * @brief integrandum for the convolution integral, integration variables are alpha and k_perp of the initial nucleon
 * 
 * @param[out] res integral result 
 * @param alpha_1 [] lc momentum fraction of initial nucleon
 * @param kperp [MeV] perp lf momentum entering in initial deuteron wf
 * @param kphi [] azimuthal angle of lf momentum entering in initial deuteron wf
 * @param[in] gpd object that contains all necessary info on the gpds 
 * @param x [] avg lc momentum fraction of struck quark
 * @param xi [] skewness
 * @param t [MeV^2] momentum transfer sq
 * @param pold_in polarization initial deuteron state (-1//0//+1)
 * @param pold_out polarization final deuteron state (-1//0//+1)
 * @param model diff implementations of KG parametrization, see TransGPD_set for details
 * @param right [1] R matrix element [0] L matrix element
 * @param deltax [MeV] perp component of the momentum transfer, along x-axis
 */
static void int_k3(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_GPD_V &gpd,
              double x, double xi, double t, int pold_in, int pold_out, double deltax);




/**
 * @brief computes nucleon matrix elements for chiral even vector GPDs
 * 
 * @param sigma_in polarization incoming nucleon (spin times two)
 * @param sigma_out polarization outgoing nucleon (spin times two)
 * @param xi_n skewness nucleon
 * @param t [MeV^2] momentum transfer sq.
 * @param t0 [MeV^2] minimum momentum transfer sq
 * @param phi [azimuthal angle \xiP+Delta fourvector]
 * @param gpd_H GPD H
 * @param gpd_E GPD E
 * @return std::complex<double> nucleon helicity amplitude
 */
std::complex<double> getGPD_even_nucl(const int sigma_in, const int sigma_out, const double xi_n, 
                                    const double t, const double t0, const double phi,
                                    const double gpd_H, const double gpd_E) const;              

TDeuteron::Wavefunction *wfref; /*!< contains instance of deuteron wave function*/
TInterpolatingWavefunction wf;
PARTONS::GPDService* pGPDService; ///< GPD Service from PARTONS
PARTONS::GPDModule* pGPDModel; ///< GPD Model from PARTONS


double t_grid; ///< [MeV^2] momentum transfer sq value the grid has
double xi_grid; ////< [] skewness value the grid has
bool grid_set; ///< is the grid with transv gpds set or not
bool ERBL_set; ///< is the grid only ERBL or not
Deut_GPD_V_set grid[201];

};

#endif 