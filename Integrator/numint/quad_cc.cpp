/**
 * @file quad_cc.cpp
 * @brief Clenshaw-Curtis quadrature rule (1D)
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt) (2010-)
 * @author k.vantournhout@gsi.de
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <numint/macros.hpp>
#include <numint/typedef.hpp>
#include <numint/numint_tools.hpp>

OPEN_NUMINT_NAMESPACE

#define CCMAX   10    // maximum iterations
#define NMAX   513    // 2^(CCMAX-1) + 1

/* =========================================================================
   The Clenshaw-Curtis work
   ========================================================================= */

template<typename T>
struct cc_work {
  double x[NMAX];   // coordinates x in [-1,1]
  double w[NMAX];   // corresponding weight factors for [-1,1]
  T     fx[NMAX];   // function values

  double M,Q;       // transformation x' = Mx+Q with x[-1,1], x'[a,b]
  
};

/* =========================================================================
   The Clenshaw-Curtis rule
   ========================================================================= */

template<typename T>
inline void _cc_step_(const function<T> &f, cc_work<T> &WORK, T &s, unsigned SIZE) {

  double *x = WORK.x;
  double *w = WORK.w;
  T *fx = WORK.fx;

  double Q = WORK.Q;
  double M = WORK.M;

  unsigned N = SIZE/2;
  // move x and fx values
  for (unsigned i = N; i > 0; --i) {
    x[2*i] = x[i];
    numint::swap(fx[2*i],fx[i]);
  }
  // fill missing x and fx values
  for (unsigned i = 1; i < N; i += 2) {
    x[i] = -cos(i * 0.5 * PI / N);
    x[2*N-i] = -x[i];
    FN_EVAL(f,M*x[i]+Q,fx[i]);
    FN_EVAL(f,-M*x[i]+Q,fx[2*N-i]);
  }
  
  // creating the weights
  //
  // The weights can be evaluated through matrix multiplications
  // See http://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature
  //
  // w = Q D Q d / N
  // with Q = diag(0.5,1,..., 1,0.5)
  // D[k][n] = cos(k n PI/N)
  // d[k] = 2/(1-(2k)^2)
  
  double D[N+1][N+1], d[N+1];
  // construction of d
  double t = 0.0;
  for (unsigned int k = 0; k <= N; ++k) {
    d[k] = 2.0 / (1.0 - t*t) / N; t += 2.0;
  }
  // Q d
  d[0] *= 0.5;
  d[N] *= 0.5;
  
  // construction of D
  t = 1.0;
  for (unsigned int k = 0; k <= N; ++k) {
    D[0][k] = D[k][0] = 1.0;
    D[N][k] = D[k][N] = t;
    t = -t;
  }
  for (unsigned int k = 1; k < N; ++k)
    for (unsigned int l = k; l < N; ++l) {
      D[k][l] = cos(k*l*PI/N);
      D[l][k] = D[k][l];
    }
  // Calculate W and multiply with Q
  for (unsigned int k = 0; k <= N; ++k) {
    w[k] = 0.0;
    for (unsigned int l = 0; l <= N; ++l)
      w[k] += D[k][l] * d[l];
      w[2*N-k] = w[k];
  }
  // multiplication with Q
  // ABS[N].W * 0.5 * 2.0; << this is needed as 
  // I = sum(w_n (f(x_n) + f(-x_n)))
  // and thus 0 comes twice
  w[0] *= 0.5; w[2*N] *= 0.5;
  
  // doing the integration
  s = fx[0]; scal(w[0],s);
  for (unsigned i = 1; i < SIZE; ++i) axpy(w[i],fx[i],s);
}



template<typename T>
inline int _quad_cc_(const function<T> &f, double a, double b,
                     double epsabs, double epsrel, T &result, unsigned &neval,
                     int dirty, bool periodic,
                     typename Convergence<T>::type compare) {

  cc_work<T> WORK;
  WORK.M = 0.5 * (b-a); WORK.Q = 0.5 * (b+a); // x' = Mx+Q

  double *x = WORK.x;
  T *fx = WORK.fx;

  T s;   // the previous refinement value
  T r;   // temporary result;

  // the first three points (-1,0,1)
  x[0] = -1.0; x[1] = 0.0; x[2] = 1.0;
  FN_EVAL(f,a,fx[0]); FN_EVAL(f,WORK.Q,fx[1]);
  if (periodic) fx[2] = fx[0];
  else FN_EVAL(f,b,fx[2]);

  // r = (fx[0] + fx[2] + 4.0*fx[1])/3.0; // first order approximation
  r = fx[0]; axpy(1.0,fx[2],r); axpy(4.0,fx[1],r); scal(1.0/3.0,r);
  
  neval += periodic ? 2 : 3;

  // the second step
  unsigned n = 2; // total new function calls
  bool conv = false;

  int RET = SUCCESS;

  for (unsigned j = 2; j < CCMAX && conv == false; ++j) {
    unsigned SIZE = (0x1 << j) + 1;

    _cc_step_(f,WORK,s,SIZE);

    neval += n;
    n *= 2;

    if (compare) {
      conv = (*compare)(r,s,epsabs,epsrel,dirty);
    } else {
      conv = convergence(r,s,epsabs,epsrel,dirty);
    }
    r = s;
    if (j == CCMAX-1 && conv == false) RET = EMAXITER;
  }

  scal(WORK.M,r);
  numint::swap(result,r);
  return RET;
}

/* =========================================================================
   Determine the needed specializations
   ========================================================================= */
// single values

int quad_cc(const function<double> &f, double a, double b, double epsabs, double epsrel,
            double &result, unsigned &neval) {
  return _quad_cc_(f,a,b,epsabs,epsrel,result,neval,0,false,NULL);
}
int quad_cc(const function<dcomplex> &f, double a, double b, double epsabs, double epsrel,
            dcomplex &result, unsigned &neval) {
  return _quad_cc_(f,a,b,epsabs,epsrel,result,neval,0,false,NULL);
}
int quad_pcc(const function<double> &f, double a, double b, double epsabs, double epsrel,
             double &result, unsigned &neval) {
  return _quad_cc_(f,a,b,epsabs,epsrel,result,neval,0,true,NULL);
}
int quad_pcc(const function<dcomplex> &f, double a, double b, double epsabs, double epsrel,
             dcomplex &result, unsigned &neval) {
  return _quad_cc_(f,a,b,epsabs,epsrel,result,neval,0,true,NULL);
}

// vector values
int quad_cc(const function<vector_d> &f, double a, double b, double epsabs, double epsrel,
            vector_d &result, unsigned &neval, int dirty,
            typename Convergence<vector_d>::type compare) {
  return _quad_cc_(f,a,b,epsabs,epsrel,result,neval,dirty,false,compare);
}
int quad_cc(const function<vector_z> &f, double a, double b, double epsabs, double epsrel,
            vector_z &result, unsigned &neval, int dirty,
            typename Convergence<vector_z>::type compare) {
  return _quad_cc_(f,a,b,epsabs,epsrel,result,neval,dirty,false,compare);
}
int quad_pcc(const function<vector_d> &f, double a, double b, double epsabs, double epsrel,
             vector_d &result, unsigned &neval, int dirty,
             typename Convergence<vector_d>::type compare) {
  return _quad_cc_(f,a,b,epsabs,epsrel,result,neval,dirty,true,compare);
}
int quad_pcc(const function<vector_z> &f, double a, double b, double epsabs, double epsrel,
             vector_z &result, unsigned &neval, int dirty,
             typename Convergence<vector_z>::type compare) {
  return _quad_cc_(f,a,b,epsabs,epsrel,result,neval,dirty,true,compare);
}

CLOSE_NUMINT_NAMESPACE

