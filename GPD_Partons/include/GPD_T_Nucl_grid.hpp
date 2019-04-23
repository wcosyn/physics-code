#ifndef GPD_T_NUCL_GRID
#define GPD_T_NUCL_GRID

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <complex>
#include <numint/numint.hpp>
#include "TransGPD_set.hpp"

class c_mstwpdf; //forward declaration


/**
 * @brief contains a grid with chiral odd nucleon GPDs calculated with the GK model
 * 
 */
class GPD_T_Nucl_grid{

public:
/**
 * @brief Construct a new gpd_T nucl grid object
 * 
 * @param pdf_name pdf parametrization name used for the forward limit ("MSTW" is the only valid one for now)
 */
GPD_T_Nucl_grid(const std::string &pdf_name);

~GPD_T_Nucl_grid();


/**
 * @brief Computes the set H_T, \bar{E}_T for up and down quarks according to the parametrisations of Goloskokov and Kroll EPJA47:112
 * Using Grid in x and xi is constructed and then interpolated, was a lot faster than calculating every value that appears in the convolution integral
 * 
 * 
 * @param x  [] parton lf momentum fraction
 * @param xi [] skewness
 * @param t [MeV^2] momentum transfer sq
 * @param scale [GeV] factorization = renorm scale
 * @return TransGPD_set constains the H_T and \bar{E}_T gpds for up and down quarks.  See TransGPD_set.
 */
TransGPD_set getGK_param(const double x, const double xi, const double t, const double scale);

/**
 * @brief Obtains the chiral odd nucleon gpds for a certain kinematics (interpolated from a grid)
 * 
 * @param x average lf momentum fraction quark
 * @param xi skewness nucleon
 * @param t [MeV^2] momentum transfer sq.
 * @param scale [GeV] factorization = renorm scale
 * @return TransGPD_set set of GPDs H, Ebar
 */
TransGPD_set getTransGPDSet(const double x, const double xi, const double t, const double scale);


private:

/**
 * @brief structure used for integration appearing in GK model GPDs
 * 
 */
struct Ftor_doubledistr {

    /*! integrandum function */
    static void exec(const numint::array<double,1> &rho, void *param, numint::vector_d &ret) {
      Ftor_doubledistr &p = * (Ftor_doubledistr *) param;
      p.f(ret,rho[0],p.x,p.xi,p.t, p.scale, *p.gpdobject);
    }
    double t; ///< [MeV^2] momentum transfer
    double x; ///< avg parton momentum fraction
    double xi; ///< skewness
    double scale; ///< [GeV] factorization = renorm scale
    GPD_T_Nucl_grid* gpdobject; ///< GPD object

    void (*f)(numint::vector_d & result, double rho, double x, double xi, double t, double scale, GPD_T_Nucl_grid &gpdobject);
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
static void DD_int_rho(numint::vector_d & res, double rho, double x, double xi, double t, double scale, GPD_T_Nucl_grid &gpdobject);


c_mstwpdf *mstw; ///< object that stores MSTW pdf grids
/**
 * @brief compute forward limit frontfactors for H_T in double distribution of GK model (Eq.24 in EPJA 47:112
 * 
 * @param x [] parton avg lightfront momentum fraction
 * @param[out] HTdfront H_T^d fronfactor
 * @param[out] HTufront H_T^u fronfactor
 * @param scale [GeV] factorization = renorm scale
 */
void getHTfront(const double x, double &HTdfront, double &HTufront, const double scale) const;
void getEbarTfront(const double x, double &EbarTdfront, double &EbarTufront) const;

double t_grid; ///< [MeV^2] momentum transfer sq value the grid has
bool grid_set; ///< is the grid with transv gpds set or not
bool ERBL_set; ///< is the grid only ERBL or not
TransGPD_set grid[201][101];


};


#endif