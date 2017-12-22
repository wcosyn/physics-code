#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <TDeuteron.h>
#include <complex>
#include <numint/numint.hpp>


/**
 * @brief conversion from helicity amplitudes to GPDs for spin 1 chiral odd quark GPDs.  We compute in a frame where phi=0
 * 
 * @param xi [] skewness
 * @param x [] parton lf momentum fraction
 * @param t [MeV^2] momentum transfer sq
 * @param helamps 
 * @return std::vector<double> 
 */
std::vector< std::complex<double> > helamps_to_gpds(const double xi, const double x, const double t, const std::vector< std::complex<double> > & helamps);

/**
 * @brief 
 * 
 */
struct Ftor_conv {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_conv &p = * (Ftor_conv *) param;
      p.f(ret,x[0],x[1],x[2], p.wfref,p.x,p.xi,p.t,p.pold_in, p.pold_out);
    }
    TDeuteron::Wavefunction *wfref;
    double x;
    double xi;
    double t;
    int pold_in;
    int pold_out;
    
    void (*f)(numint::vector_z & res, double knorm, double kcosth, double kphi, TDeuteron::Wavefunction *wfref,
              double x, double xi, double t, int pold_in, int pold_out);
      };


/**
 * @brief 
 * 
 * @param res 
 * @param knorm 
 * @param kcosth 
 * @param kphi 
 * @param wfref 
 * @param x 
 * @param xi 
 * @param t 
 * @param pold_in polarization initial state (-1//0//+1)
 * @param pold_out polarization final state state (-1//0//+1)
 */
void int_k3(numint::vector_z & res, double knorm, double kcosth, double kphi, TDeuteron::Wavefunction *wfref,
              double x, double xi, double t, int pold_in, int pold_out);


std::complex<double> getGPD_odd_nucl(int sigma_in, int sigma_out, double x_n, double xi_n, double t, double t0);              

double getGPD_HT(double x,double xi, double t);
double getGPD_HtildeT(double x,double xi, double t);
double getGPD_ET(double x,double xi, double t);
double getGPD_EtildeT(double x,double xi, double t);