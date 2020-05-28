#ifndef DEUT_CONV_GPD_V
#define DEUT_CONV_GPD_V


#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <TDeuteron.h>
#include <complex>
#include <numint/numint.hpp>
#include <TInterpolatingWavefunction.h>

#include "GPD_V_Nucl_grid.hpp"

#include <partons/Partons.h>
#include <partons/services/GPDService.h>

#include "Deut_GPD_V_set.hpp"

#include <NucleonEMOperator.hpp>

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
 * @brief conversion from GPDs to helicity amplitudes for spin 1 chiral even vector quark FFs.  We compute in a frame where phi=0
 * See Cano Pire EPJA App A (restricted to ++,0+,-+ helamps) + mathematica sheet (gpd_helmatrix_newtensors)
 * 
 * @param xi [] skewness
 * @param t [GeV^2] momentum transfer sq
 * @param FFs gpds G_1 to G_3
 * @return std::vector<double> helicity amplitudes, deuteron helicities are (final,initial) [0]++,[1]0+,[2]-+
 */
static std::vector< std::complex<double> > FFs_to_helamps_V(const double xi, const double t, const std::vector< std::complex<double> > & FFs);

/**
 * @brief conversion from helicity amplitudes to vector FFs for spin 1.  We compute in a frame where phi=0
 * See mathematica sheet (gpd_helmatrix_newtensors)
 * 
 * @param xi [] skewness
 * @param t [GeV^2] momentum transfer sq
 * @param helamps helicity amplitudes, deuteron helicities are (final,initial) [0]++,[1]0+,[2]-+
 * @return std::vector<double>  FFs G_1 to G_3
 */
static std::vector< std::complex<double> > helamps_to_FFs_V(const double xi, const double t, const std::vector< std::complex<double> > & helamps);



/**
 * @brief conversion from helicity amplitudes to GPDs for spin 1 chiral even vector quark GPDs.  We compute in a frame where phi=0
 * See Cano Pire EPJA App A
 * 
 * @param xi [] skewness
 * @param t [GeV^2] momentum transfer sq
 * @param helamps helicity amplitudes, deuteron helicities are (final,initial) [0]++,[1]00,[2]0+,[3]+0,[4]-+
 * @return std::vector<double>  gpds H_1 to H_5
 */
static std::vector< std::complex<double> > helamps_to_gpds_V(const double xi, const double t, const std::vector< std::complex<double> > & helamps);


/**
 * @brief conversion from GPDs to helicity amplitudes for spin 1 chiral even vector quark GPDs.  We compute in a frame where phi=0
 * See Cano Pire EPJA App A
 * 
 * @param xi [] skewness
 * @param t [GeV^2] momentum transfer sq
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
 * @param t [Gev^2] momentum transfer sq
 * @return std::vector< std::complex<double> > helicity amplitudes  deuteron helicities are (final,initial) [0]++,[1]00,[2]0+,[3]+0,[4]-+ 
 */
std::vector< std::complex<double> > gpd_conv(const double xi, const double x, const double t, const double scale);

/**
 * @brief Computes the deuteron helicity amplitudes with the convolution formula for vector form factors
 * 
 * @param xi [] skewness
 * @param t [Gev^2] momentum transfer sq
 * @return std::vector< std::complex<double> > helicity amplitudes  deuteron helicities are (final,initial) [0]++,[1]00,[2]0+,[3]+0,[4]-+ 
 */
std::vector< std::complex<double> > FF_conv(const double xi, const double t);

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
Deut_GPD_V_set getDeut_GPD_V_set(const double x, const double xi, const double t, const double scale, const bool ERBL);

void setH(const bool H){incH=H;} ///< [1] include or [0] exclude isoscalar nucleon GPD H
void setE(const bool E){incE=E;} ///< [1] include or [0] exclude isoscalar nucleon GPD E


/**
 * @brief Calculates three deuteron vector form factors in a non-relativistic approximation
 * 
 * @param t [GeV^2] momentum transfer squared (negative!)
 * @return vector<double> [0] G_C [1] G_M [2] G_Q (See Gross, Gilman JPG '01 review Eq(15) for FF definitions )
 */
std::vector<double> calc_NR_ffs(double t);

/**
 * @brief computes nucleon matrix elements for chiral even vector GPDs
 * 
 * @param sigma_in polarization incoming nucleon (spin times two!!)
 * @param sigma_out polarization outgoing nucleon (spin times two!!)
 * @param xi_n skewness nucleon
 * @param t [MeV^2] momentum transfer sq.
 * @param t0 [MeV^2] minimum momentum transfer sq
 * @param phi [azimuthal angle \xiP+Delta fourvector]
 * @param gpd_H GPD H
 * @param gpd_E GPD E
 * @param [1] vector [0] axial GPDS
 * @return std::complex<double> nucleon helicity amplitude
 */
static std::complex<double> getGPD_even_nucl(const int sigma_in, const int sigma_out, const double xi_n, 
                                    const double t, const double t0, const double phi,
                                    const double gpd_H, const double gpd_E, const bool vector);              




private: 
/**
 * @brief structure needed to carry out the convolution integration, contains parameters and integrandum function
 * 
 */
struct Ftor_conv {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_conv &p = * (Ftor_conv *) param;
      p.f(ret,x[0],x[1],x[2], *p.gpd,p.x,p.xi,p.t,p.scale,p.pold_in, p.pold_out,p.deltax);
    }
    Deut_Conv_GPD_V *gpd;
    double x;
    double xi;
    double t; // [MeV^2]
    double scale; ///< [GeV]
    int pold_in;
    int pold_out;
    double deltax; //[MeV]
    
    void (*f)(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_GPD_V &gpd,
              double x, double xi, double t, double scale, int pold_in, int pold_out,  double deltax);
      };

/**
 * @brief structure needed to carry out the non relativistic evalution of deuteron form factors, contains parameters and integrandum function
 * 
 */
struct Ftor_conv_FF_nr {

    /*! integrandum function */
    static void exec(const numint::array<double,1> &x, void *param, numint::vector_d &ret) {
      Ftor_conv_FF_nr &p = * (Ftor_conv_FF_nr *) param;
      p.f(ret,x[0], *p.gpd,p.t);
    }
    Deut_Conv_GPD_V *gpd;
    double t;
    
    void (*f)(numint::vector_d & res, double r, Deut_Conv_GPD_V &gpd, double t);
};

/**
 * @brief integration function for NR evaluation of deuteron form factors, see Gross, Gilman JPG '01 review, around Eq (34)
 * 
 * @param[out] res returns [0] D_C [1] D_Q [1] D_M [2] D_E 
 * @param r integration variable
 * @param gpd object that contains all usefull things
 * @param t [GeV^2] momentum transfer squared (negative!)
 */
static void int_r(numint::vector_d & res, double r, Deut_Conv_GPD_V &gpd, double t);




/**
 * @brief structure needed to carry out the convolution integration for the form factors, contains parameters and integrandum function
 * 
 */
struct Ftor_conv_FF {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_conv_FF &p = * (Ftor_conv_FF *) param;
      p.f(ret,x[0],x[1],x[2], *p.gpd,p.xi,p.t,p.pold_in, p.pold_out,p.deltax,p.F1, p.F2);
    }
    Deut_Conv_GPD_V *gpd;
     double xi;
    double t; // [MeV^2]
    int pold_in;
    int pold_out;
    double deltax; //[MeV]
    double F1, F2;
    
    void (*f)(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_GPD_V &gpd,
               double xi, double t, int pold_in, int pold_out,  double deltax, double F1, double F2);
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
 * @param scale [GeV] factorization = renorm scale
 * @param pold_in polarization initial deuteron state (-1//0//+1)
 * @param pold_out polarization final deuteron state (-1//0//+1)
 * @param model diff implementations of KG parametrization, see TransGPD_set for details
 * @param right [1] R matrix element [0] L matrix element
 * @param deltax [MeV!!!] perp component of the momentum transfer, along x-axis
 */
static void int_k3(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_GPD_V &gpd,
              double x, double xi, double t, double scale, int pold_in, int pold_out, double deltax);

/**
 * @brief integrandum for the convolution integral, integration variables are alpha and k_perp of the initial nucleon
 * 
 * @param[out] res integral result 
 * @param alpha_1 [] lc momentum fraction of initial nucleon
 * @param kperp [MeV] perp lf momentum entering in initial deuteron wf
 * @param kphi [] azimuthal angle of lf momentum entering in initial deuteron wf
 * @param[in] gpd object that contains all necessary info on the gpds 
 * @param xi [] skewness
 * @param t [GeV^2] momentum transfer sq
 * @param pold_in polarization initial deuteron state (-1//0//+1)
 * @param pold_out polarization final deuteron state (-1//0//+1)
 * @param model diff implementations of KG parametrization, see TransGPD_set for details
 * @param right [1] R matrix element [0] L matrix element
 * @param deltax [MeV!!!] perp component of the momentum transfer, along x-axis
 * @param F1 isoscalar nucleon FF F1
 * @param F2 isoscalar nucleon FF F2
 */
static void int_k3_FF(numint::vector_z & res, double alpha_1, double kperp, double kphi, Deut_Conv_GPD_V &gpd,
               double xi, double t, int pold_in, int pold_out, double deltax, double F1, double F2);






TDeuteron::Wavefunction *wfref; /*!< contains instance of deuteron wave function*/
TInterpolatingWavefunction wf;
PARTONS::GPDService* pGPDService; ///< GPD Service from PARTONS
PARTONS::GPDModule* pGPDModel; ///< GPD Model from PARTONS

GPD_V_Nucl_grid chiraleven_grid;


double t_grid; ///< [GeV^2] momentum transfer sq value the grid has
double xi_grid; ////< [] skewness value the grid has
bool grid_set; ///< is the grid with transv gpds set or not
bool ERBL_set; ///< is the grid only ERBL or not
Deut_GPD_V_set grid[201];

bool incH; ///< include GPD H
bool incE; ///< include GPD E

};

#endif 