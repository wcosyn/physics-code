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
#include "GPD_T_Nucl_grid.hpp"
#include "Deut_GPD_T_set.hpp"

class c_mstwpdf; //forward declaration
class TransGPD_set;
/**
 * @brief Class implements transversity GPDs for the nucleon and deuteron (using convolution formalism) for now
 * 
 */
class Deut_Conv_GPD_T{

public:

/**
 * @brief constructor
 * 
 * @param pdf_name pdf parametrization name used for the forward limit ("MSTW" is the only valid one for now)
 * @param wfname deuteron wave function parametrization name, see TDeuteron for possibilities
 */
Deut_Conv_GPD_T(const std::string &pdf_name, const std::string &wfname);

/**
 * @brief Destructor
 * 
 */
~Deut_Conv_GPD_T();

/**
 * @brief conversion from helicity amplitudes to GPDs for spin 1 chiral odd quark GPDs.  We compute in a frame where phi=0
 * See Cosyn, Pire PRD '18 App C
 * 
 * @param xi [] skewness
 * @param t [MeV^2] momentum transfer sq
 * @param helamps helicity amplitudes, deuteron helicities are (final,initial) [0]++,[1]--,[2]00,[3]0+,[4]-0,[5]+0,[6]0-,[7]-+,[8]+-
 * @return std::vector<double>  gpds H^T_1 to H^T_9
 */
static std::vector< std::complex<double> > helamps_to_gpds_T(const double xi, const double t, const std::vector< std::complex<double> > & helamps);

/**
 * @brief conversion from GPDs to helicity amplitudes for spin 1 chiral odd quark GPDs.  We compute in a frame where phi=0
 * See Cosyn, Pire PRD '18 App C
 * 
 * @param xi [] skewness
 * @param t [MeV^2] momentum transfer sq
 * @param gpds gpds H^T_1 to H^T_9
 * @return std::vector<double> helicity amplitudes, deuteron helicities are (final,initial) [0]++,[1]--,[2]00,[3]0+,[4]-0,[5]+0,[6]0-,[7]-+,[8]+-
 */
static std::vector< std::complex<double> > gpds_to_helamps_T(const double xi, const double t, const std::vector< std::complex<double> > & gpds);


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
 * @param model [0-2] different models for nucleon chiral odd GPDs
 * @return std::vector< std::complex<double> > helicity amplitudes  deuteron helicities are (final,initial) [0]++,[1]--,[2]00,[3]0+,[4]-0,[5]+0,[6]0-,[7]-+,[8]+-
 */
std::vector< std::complex<double> > gpd_conv(const double xi, const double x, const double t, const int model);

TInterpolatingWavefunction * getWf(){ return &wf;} ///< return deuteron wf object

/**
 * @brief Get a Deut_GPD_V_set object containing deuteron vector helicity amplitudes
 * 
 * @param x [] avg lf momentum fraction
 * @param xi [] skewness
 * @param t [GeV^2] mom transfer sq
 * @param ERBL grid for only ERBL region or not?
 * @param model [0-2] different models for nucleon chiral odd GPDs
 * @return Deut_GPD_V_set 
 */
Deut_GPD_T_set getDeut_GPD_T_set(const double x, const double xi, const double t, const bool ERBL, const int model);


/**
 * @brief computes nucleon matrix elements through helicity amplitudes for tensor current
 * 
 * @param sigma_in polarization incoming nucleon
 * @param sigma_out polarization outgoing nucleon
 * @param xi_n skewness nucleon
 * @param t [MeV^2] momentum transfer sq.
 * @param t0 [MeV^2] minimum momentum transfer sq
 * @param phi [azimuthal angle \xiP+Delta fourvector]
 * @param model different implementations of the nucleon transversity GPDs
 * @param right [1] right or [0] left matrix element (see notes Eq 112)
 * @param gpd_nucl constains a set of chiral odd quark nucleon gpds computed according to the GK model
 * @return std::complex<double> nucleon helicity amplitude
 */
static std::complex<double> getGPD_odd_nucl(const int sigma_in, const int sigma_out, const double xi_n, 
                                    const double t, const double t0, const double phi, const int model, const bool right,
                                    const TransGPD_set &gpd_nucl);              


private: 
/**
 * @brief structure needed to carry out the convolution integration, contains parameters and integrandum function
 * 
 */
struct Ftor_conv {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_conv &p = * (Ftor_conv *) param;
      p.f(ret,x[0],x[1],x[2], *p.gpd,p.x,p.xi,p.t,p.pold_in, p.pold_out,p.model,p.right,p.deltax);
    }
    Deut_Conv_GPD_T *gpd;
    double x;
    double xi;
    double t;
    int pold_in;
    int pold_out;
    int model;
    bool right;
    double deltax;
    
    void (*f)(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_GPD_T &gpd,
              double x, double xi, double t, int pold_in, int pold_out, int model, bool right, double deltax);
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
static void int_k3(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_GPD_T &gpd,
              double x, double xi, double t, int pold_in, int pold_out, int model, bool right, double deltax);
static void int_kprime3(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_GPD_T &gpd,
              double x, double xi, double t, int pold_in, int pold_out, int model, bool right, double deltax);


/**
 * @brief testing: do the convolution in minimal fashion: no D-wave, no wave function, only S-wave spin sums, 
 * symmetric kinematics retain all the symmetries of helicity amplitudes
 * 
 * @param x avg parton momentum fraction
 * @param xi skewness
 * @param t [MeV^2] momentum transfer squared
 * @param pold_in initial deuteron helicity
 * @param pold_out final deuteron helicity
 * @param model [0-2] type of model for chiral odd nucleon GPDs
 * @param right [0] left transverse index [1] right transverse index
 * @param deltax [MeV] perp component of momentum transfer
 */
std::complex<double >test(double x, double xi, double t, int pold_in, int pold_out, int model, bool right, double deltax);




TDeuteron::Wavefunction *wfref; /*!< contains instance of deuteron wave function*/
TInterpolatingWavefunction wf;
GPD_T_Nucl_grid chiralodd_grid; ///< object that allows to interpolate or get chiral odd nucleon GPDs

double t_grid; ///< [MeV^2] momentum transfer sq value the grid has
double xi_grid; ////< [] skewness value the grid has
bool grid_set; ///< is the grid with transv gpds set or not
bool ERBL_set; ///< is the grid only ERBL or not
Deut_GPD_T_set grid[201];

};

#endif 