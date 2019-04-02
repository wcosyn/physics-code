#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <TDeuteron.h>
#include <complex>
#include <numint/numint.hpp>
#include "TransGPD_set.hpp"
#include <TInterpolatingWavefunction.h>

class c_mstwpdf; //forward declaration
/**
 * @brief Class implements transversity GPDs for the nucleon and deuteron (using convolution formalism) for now
 * 
 */
class GPD{

public:

/**
 * @brief constructor
 * 
 * @param pdf_name pdf parametrization name used for the forward limit ("MSTW" is the only valid one for now)
 * @param wfname deuteron wave function parametrization name, see TDeuteron for possibilities
 */
GPD(const std::string &pdf_name, const std::string &wfname);

/**
 * @brief Destructor
 * 
 */
~GPD();

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
 * @brief Computes the set H_T, \bar{E}_T for up and down quarks according to the parametrisations of Goloskokov and Kroll EPJA47:112
 * Grid in x and xi is constructed and then interpolated, was a lot faster than calculating every value that appears in the convolution integral
 * 
 * 
 * @param x  [] parton lf momentum fraction
 * @param xi [] skewness
 * @param t [MeV^2] momentum transfer sq
 * @return TransGPD_set constains the H_T and \bar{E}_T gpds for up and down quarks.  See TransGPD_set.
 */
TransGPD_set getGK_param(const double x, const double xi, const double t);


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
    GPD *gpd;
    double x;
    double xi;
    double t;
    int pold_in;
    int pold_out;
    int model;
    bool right;
    double deltax;
    
    void (*f)(numint::vector_z & res, double alpha_1, double kperp, double kphi, GPD &gpd,
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
static void int_k3(numint::vector_z & res, double alpha_1, double kperp, double kphi, GPD &gpd,
              double x, double xi, double t, int pold_in, int pold_out, int model, bool right, double deltax);
static void int_kprime3(numint::vector_z & res, double alpha_1, double kperp, double kphi, GPD &gpd,
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

/**
 * @brief structure used for integration appearing in GK model GPDs
 * 
 */
struct Ftor_doubledistr {

    /*! integrandum function */
    static void exec(const numint::array<double,1> &rho, void *param, numint::vector_d &ret) {
      Ftor_doubledistr &p = * (Ftor_doubledistr *) param;
      p.f(ret,rho[0],p.x,p.xi,p.t, *p.gpdobject);
    }
    double t; ///< [MeV^2] momentum transfer
    double x; ///< avg parton momentum fraction
    double xi; ///< skewness
    GPD* gpdobject; ///< GPD object

    void (*f)(numint::vector_d & result, double rho, double x, double xi, double t, GPD &gpdobject);
      };

/**
 * @brief integral of double distribution to obtain transversity GPD according to GK prescription [EPJA 47:112 Eq 14]
 * 
 * @param res [0] H_T^d, [1] H_T^u, [2] \bar{E}_T^d, [3] \bar{E}_T^u
 * @param rho integral variable
 * @param x  avg parton lf momentum fraction
 * @param xi skewness
 * @param t [MeV^2] momentum transfer sq.
 */
static void DD_int_rho(numint::vector_d & res, double rho, double x, double xi, double t, GPD &gpdobject);


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
std::complex<double> getGPD_odd_nucl(const int sigma_in, const int sigma_out, const double xi_n, 
                                    const double t, const double t0, const double phi, const int model, const bool right,
                                    const TransGPD_set &gpd_nucl) const;              

/**
 * @brief Obtains the chiral odd nucleon gpds for a certain kinematics (interpolated from a grid)
 * 
 * @param x average lf momentum fraction quark
 * @param xi skewness nucleon
 * @param t [MeV^2] momentum transfer sq.
 * @return TransGPD_set set of GPDs H, Ebar
 */
TransGPD_set getTransGPDSet(const double x, const double xi, const double t);

c_mstwpdf *mstw; ///< object that stores MSTW pdf grids
TDeuteron::Wavefunction *wfref; /*!< contains instance of deuteron wave function*/
TInterpolatingWavefunction wf;
/**
 * @brief compute forward limit frontfactors for H_T in double distribution of GK model (Eq.24 in EPJA 47:112
 * 
 * @param x [] parton avg lightfront momentum fraction
 * @param[out] HTdfront H_T^d fronfactor
 * @param[out] HTufront H_T^u fronfactor
 */
void getHTfront(const double x, double &HTdfront, double &HTufront) const;
void getEbarTfront(const double x, double &EbarTdfront, double &EbarTufront) const;

double t_grid; ///< [MeV^2] momentum transfer sq value the grid has
bool grid_set; ///< is the grid with transv gpds set or not
TransGPD_set grid[201][101];


};