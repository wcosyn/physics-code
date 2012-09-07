/**
 * @file quad_romb.cpp
 * @brief Romberg + trapezoidal quadrature rule (1D)
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt) (2010-)
 * @author k.vantournhout@gsi.de
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iomanip>

#include <numint/macros.hpp>
#include <numint/typedef.hpp>
#include <numint/numint_tools.hpp>

OPEN_NUMINT_NAMESPACE

#define ROMBMAX 15
#define ROMBMIN 3

template<typename T>
struct romb_work {
  double a; // upper limit
  double b; // lower limit
  double h; // h
  T sum;    // intermediate sum
  T t;      // temporary

  T r1[ROMBMAX]; // The first Romberg line
  T r2[ROMBMAX]; // The second Romberg line
  T *R1;
  T *R2;
};

/* =========================================================================
   Romberg integration based on the trapezoidal rule
   ========================================================================= */

template<typename T>
inline void tr_step(const function<T> &f, double a,
                    T &s, double &h, unsigned &n, T &sum) {
  s *= 0.5;
  sum = FN_EVAL(f,a+h);
  for (unsigned i = 1; i < n; ++i) {
    double x = a + h * (1+2*i);
    sum += FN_EVAL(f,x);
  }
  sum *= h;
  s += sum;
}

template<typename T>
inline void _romb_step_(const function<T> &f, romb_work<T> &WORK, T &s,
                        unsigned &n, unsigned j) {
  
  // do a trapezoidal step
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


  // do Romberg integration
  T*& R1 = WORK.R1; T*& R2 = WORK.R2;
  T* R3; R3 = R1; R1 = R2; R2 = R3;

  R2[0] = s;

  // Romberg integration
  double fourm = 4.0;
  for (unsigned m = 1; m <= j; ++m) {
    //    R2[m] = (fourm * R2[m-1] - R1[m-1]) / (fourm-1.0);
    R2[m] = R2[m-1];
    scal(fourm,R2[m]);
    axpy(-1.0,R1[m-1],R2[m]);
    scal(1.0/(fourm-1.0),R2[m]);
    fourm *= 4.0;
  }
}

template<typename T>
inline int _quad_romb_(const function<T> &f, double a, double b,
                       double epsabs, double epsrel, T &result, unsigned &neval,
                       int dirty, bool periodic,
                       typename Convergence<T>::type compare) {

  T s; // the previous refinement value
  T rtr;   // temporary result (for trapezoidal convergence)
  T rromb; // temporary result (for romberg convergence)
 
  romb_work<T> WORK; WORK.a = a; WORK.b = b;
  WORK.R1 = WORK.r1; WORK.R2 = WORK.r2;
  double &h = WORK.h;
  h = b-a;

  // the first step
  if (periodic) {
    FN_EVAL(f,a,s);
    scal(h,s);
    rtr = rromb = s;
    neval += 1;
  } else {
    FN_EVAL(f,a,s); FN_EVAL(f,b,WORK.t);
    axpy(1.0,WORK.t,s);
    scal(0.5*h,s);
    rtr = rromb = s;
    neval += 2;
  }

  T*& R2 = WORK.R2;
  R2[0] = s;

  unsigned n = 1; // total new function calls

  bool conv_tr   = false;
  bool conv_romb = false;

  int RET = SUCCESS;

  for (unsigned j = 1; j < ROMBMAX && conv_tr == false && conv_romb == false; ++j) {
    h *= 0.5; // adjust the step size
    _romb_step_(f,WORK,s,n,j);
    neval += n; // total new function calls
    n *= 2;     // increase n for the next step;

    // we do two convergence checks, the romberg check and the trapezoidal check
    if (j > ROMBMIN) {
      conv_tr   = convergence(rtr,R2[0],epsabs,epsrel,dirty);   // trapezoidal check
      conv_romb = convergence(rromb,R2[j],epsabs,epsrel,dirty); // romberg check
    }

    rtr = R2[0];    rromb = R2[j];
    if (j == ROMBMAX-1 && conv_tr == false && conv_romb == false) RET = EMAXITER;
  }

  if (conv_romb) numint::swap(result,rromb);
  else numint::swap(result,rtr);

  return RET;

}

/* =========================================================================
   Determine the needed specializations
   ========================================================================= */
// single values

int quad_romb(const function<double> &f, double a, double b, double epsabs, double epsrel,
            double &result, unsigned &neval) {
  return _quad_romb_(f,a,b,epsabs,epsrel,result,neval,0,false,NULL);
}
int quad_romb(const function<dcomplex> &f, double a, double b, double epsabs, double epsrel,
            dcomplex &result, unsigned &neval) {
  return _quad_romb_(f,a,b,epsabs,epsrel,result,neval,0,false,NULL);
}
int quad_promb(const function<double> &f, double a, double b, double epsabs, double epsrel,
             double &result, unsigned &neval) {
  return _quad_romb_(f,a,b,epsabs,epsrel,result,neval,0,true,NULL);
}
int quad_promb(const function<dcomplex> &f, double a, double b, double epsabs, double epsrel,
             dcomplex &result, unsigned &neval) {
  return _quad_romb_(f,a,b,epsabs,epsrel,result,neval,0,true,NULL);
}

// vector values
int quad_romb(const function<vector_d> &f, double a, double b, double epsabs, double epsrel,
            vector_d &result, unsigned &neval, int dirty,
            typename Convergence<vector_d>::type compare) {
  return _quad_romb_(f,a,b,epsabs,epsrel,result,neval,dirty,false,compare);
}
int quad_romb(const function<vector_z> &f, double a, double b, double epsabs, double epsrel,
            vector_z &result, unsigned &neval, int dirty,
            typename Convergence<vector_z>::type compare) {
  return _quad_romb_(f,a,b,epsabs,epsrel,result,neval,dirty,false,compare);
}
int quad_promb(const function<vector_d> &f, double a, double b, double epsabs, double epsrel,
             vector_d &result, unsigned &neval, int dirty,
             typename Convergence<vector_d>::type compare) {
  return _quad_romb_(f,a,b,epsabs,epsrel,result,neval,dirty,true,compare);
}
int quad_promb(const function<vector_z> &f, double a, double b, double epsabs, double epsrel,
             vector_z &result, unsigned &neval, int dirty,
             typename Convergence<vector_z>::type compare) {
  return _quad_romb_(f,a,b,epsabs,epsrel,result,neval,dirty,true,compare);
}

CLOSE_NUMINT_NAMESPACE
