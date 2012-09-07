/**
 * @file quad_tr.cpp
 * @brief Trapezoidal quadrature rule (1D)
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt) (2010-)
 * @author k.vantournhout@gsi.de
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <iomanip>

#include <numint/macros.hpp>
#include <numint/typedef.hpp>
#include <numint/numint_tools.hpp>

OPEN_NUMINT_NAMESPACE

#define TRMAX 15

template<typename T>
struct tr_work {
  double a; // upper limit
  double b; // lower limit
  double h; // h
  T sum;    // intermediate sum
  T t;      // temporary
};

/* =========================================================================
   The trapezoidal rule
   ========================================================================= */

/**
   @brief a single step in the trapezoidal rule

   @param[in] f the function
   @param[in] a the lower intgral limit
   @param[in,out] s the refinement (in = older, out = new)
   @param[in] h the step size
   @param[in] n the total new steps
 */

template<typename T>
inline void _tr_step_(const function<T> &f, tr_work<T> &WORK, T &s, unsigned n) {
  
  T &sum = WORK.sum;
  T &t   = WORK.t;
  double a = WORK.a; double h = WORK.h;

  scal(0.5,s);
  FN_EVAL(f,a+h,sum);

  for (unsigned i = 1; i < n; ++i) {
    double x = a + h * (1+2*i);
    FN_EVAL(f,x,t);
    axpy(1.0,t,sum);
  }
  axpy(h,sum,s);
}

/**
   @brief trapezoidal integration

   @param[in] f the function
   @param[in] a the lower integral limit
   @param[in] b the upper integral limit
   @param[in] epsabs the absolute error
   @param[in] epsrel the relative error
   @param[out] result the integral
   @param neval number of evaluations
 */

template<typename T>
inline int _quad_tr_(const function<T> &f, double a, double b,
                     double epsabs, double epsrel, T &result, unsigned &neval,
                     int dirty, bool periodic,
                     typename Convergence<T>::type compare) {

  T s;   // the previous refinement value
  T r;   // temporary result;

  tr_work<T> WORK; WORK.a = a; WORK.b = b;
  double &h = WORK.h;
  h = b-a;

  std::cout << "===== START =====\n";

  // the first step
  if (periodic) {
    FN_EVAL(f,a,s);
    scal(h,s);
    r = s;
    neval += 1;
  } else {
    FN_EVAL(f,a,s); FN_EVAL(f,b,WORK.t);
    axpy(1.0,WORK.t,s);
    scal(0.5*h,s);
    r = s;
    neval += 2;
  }

  unsigned n = 1; // total new function calls

  bool conv = false;

  int RET = SUCCESS;

  for (unsigned j = 1; j < TRMAX && conv == false; ++j) {
    h *= 0.5;   // adjust the step size
    _tr_step_(f,WORK,s,n);
    neval += n; // total new function calls
    n *= 2;     // increase n for the next step;

    if (compare) {
      conv = (*compare)(r,s,epsabs,epsrel,dirty);
    } else {
      conv = convergence(r,s,epsabs,epsrel,dirty);
    }

    r = s;
    if (j == TRMAX-1 && conv == false) RET = EMAXITER;

  }
  std::cout << "\n===== STOP =====\n";
  numint::swap(result,r);
  
  return RET;

}

/* =========================================================================
   Determine the needed specializations
   ========================================================================= */
// single values

int quad_tr(const function<double> &f, double a, double b, double epsabs, double epsrel,
            double &result, unsigned &neval) {
  return _quad_tr_(f,a,b,epsabs,epsrel,result,neval,0,false,NULL);
}
int quad_tr(const function<dcomplex> &f, double a, double b, double epsabs, double epsrel,
            dcomplex &result, unsigned &neval) {
  return _quad_tr_(f,a,b,epsabs,epsrel,result,neval,0,false,NULL);
}
int quad_ptr(const function<double> &f, double a, double b, double epsabs, double epsrel,
             double &result, unsigned &neval) {
  return _quad_tr_(f,a,b,epsabs,epsrel,result,neval,0,true,NULL);
}
int quad_ptr(const function<dcomplex> &f, double a, double b, double epsabs, double epsrel,
             dcomplex &result, unsigned &neval) {
  return _quad_tr_(f,a,b,epsabs,epsrel,result,neval,0,true,NULL);
}

// vector values
int quad_tr(const function<vector_d> &f, double a, double b, double epsabs, double epsrel,
            vector_d &result, unsigned &neval, int dirty,
            typename Convergence<vector_d>::type compare) {
  return _quad_tr_(f,a,b,epsabs,epsrel,result,neval,dirty,false,compare);
}
int quad_tr(const function<vector_z> &f, double a, double b, double epsabs, double epsrel,
            vector_z &result, unsigned &neval, int dirty,
            typename Convergence<vector_z>::type compare) {
  return _quad_tr_(f,a,b,epsabs,epsrel,result,neval,dirty,false,compare);
}
int quad_ptr(const function<vector_d> &f, double a, double b, double epsabs, double epsrel,
             vector_d &result, unsigned &neval, int dirty,
             typename Convergence<vector_d>::type compare) {
  return _quad_tr_(f,a,b,epsabs,epsrel,result,neval,dirty,true,compare);
}
int quad_ptr(const function<vector_z> &f, double a, double b, double epsabs, double epsrel,
             vector_z &result, unsigned &neval, int dirty,
             typename Convergence<vector_z>::type compare) {
  return _quad_tr_(f,a,b,epsabs,epsrel,result,neval,dirty,true,compare);
}

CLOSE_NUMINT_NAMESPACE
