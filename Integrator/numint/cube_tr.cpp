/**
 * @file cube_tr.cpp
 * @brief trapezoidal cubature rule
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt)
 * @author k.vantournhout@gsi.de
 * @date creation: 11/04/2011
 * @date last update: 11/04/2011
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <numint/macros.hpp>
#include <numint/typedef.hpp>
#include <numint/numint_tools.hpp>

OPEN_NUMINT_NAMESPACE

#define TRMAX 15

template<typename T, unsigned N>
struct cube_tr_work {
  numint::array<double,N> a; // upper limit
  numint::array<double,N> b; // lower limit
  numint::array<double,N> h; // h
  T sum;    // intermediate sum
  T t;      // temporary
};

/**
   @brief a single step in the cubature trapezoidal rule

   @param[in] f the function
   @param[in] a the lower intgral limit
   @param[in,out] s the refinement (in = older, out = new)
   @param[in] h the step size
   @param[in] n the total points in one dimension
*/

template<typename T, unsigned N>
inline void _cube_tr_step(const mdfunction<T,N> &f, cube_tr_work<T,N> &WORK, T &s, unsigned n,
                          bool periodic) {
  T &sum = WORK.sum;
  T &t   = WORK.t;
  numint::array<double,N> &a = WORK.a;
  numint::array<double,N> &h = WORK.h;
  numint::array<double,N> x;
  static const int pow = (0x1 << N);
  
  double htot = 1.0;
  for (unsigned i = 0; i < N; ++i) htot *= h[i];
  
  scal(1.0/pow,s); // scale the previous sum
  
  // set sum to zero
  sum = s;
  axpy(-1.0,s,sum);

  // total number of points in the sum
  unsigned Ntot = 1;
  unsigned I[N];
  for (unsigned i = 0; i < N; ++i) Ntot *= 2*n+1; // n should be different from zero

  for (unsigned idx = 0; idx < Ntot; ++idx) { // run over all points
    unsigned itmp = idx;
    bool odd = false;
    for (unsigned i = 0; i < N; ++i) {
      I[i] = itmp % (2*n+1);
      itmp /= (2*n+1);
      odd = odd || (I[i] % 2 == 1);

    }
    if (odd) { // only calculate if there is an odd index
      double c = 1.0;
      for (unsigned i = 0; i < N; ++i) {
        if ((!periodic) && (I[i] == 0 || I[i] == 2*n)) c *= 0.5; // update if we are on the border
        x[i] = a[i] + h[i] * I[i];
      }

      FN_EVAL<T,N>(f,x,t);
      axpy(c,t,sum);
    }
  }

  axpy(htot,sum,s);
}

template<typename T, unsigned N>
inline int _cube_tr_(const mdfunction<T,N> &f, const numint::array<double,N> &a,
                     const numint::array<double,N> &b,
                     double epsabs, double epsrel, T &result, unsigned &neval,
                     int dirty, bool periodic, typename Convergence<T>::type compare) {

  T s;   // the previous refinement value
  T r;   // temporary result;

  cube_tr_work<T,N> WORK; WORK.a = a; WORK.b = b;
  numint::array<double,N> &h = WORK.h;
  double htot = 1.0;
  for (unsigned i = 0; i < N; ++i) { h[i] = b[i] - a[i]; htot *= h[i]; }
  
  // the first step
  if (periodic) {
    FN_EVAL<T,N>(f,a,s);
    scal(htot,s);
    r = s;
  } else {

    FN_EVAL<T,N>(f,a,s);

    numint::array<double,N> x;
    int Ntot = 0x1 << N;
    for (int idx = 1; idx < Ntot; ++idx) { // run over all points
      unsigned itmp = idx;
      for (unsigned i = 0; i < N; ++i) {
        x[i] = itmp % 2 ? b[i] : a[i];
        itmp /= 2;
      }

      FN_EVAL<T,N>(f,x,WORK.t);
      axpy(1.0,WORK.t,s);
    }
    scal(htot/Ntot,s);
    r = s;
  }

  unsigned n = 1; // total new function calls in 1 dimension

  bool conv = false;

  int RET = SUCCESS;

  for (unsigned j = 1; j < TRMAX && conv == false; ++j) {
    
    for (unsigned i = 0; i < N; ++i) h[i] *= 0.5;   // adjust the step size
    _cube_tr_step<T,N>(f,WORK,s,n,periodic);
    n *= 2;     // increase n for the next step;
    if (compare) {
      conv = (*compare)(r,s,epsabs,epsrel,dirty);
    } else {
      conv = convergence(r,s,epsabs,epsrel,dirty);
    }
    r = s;
    if (j == TRMAX-1 && conv == false) RET = EMAXITER;
  }

  unsigned Ntot = 1;
  for (unsigned i = 0; i < N; ++i) Ntot *= periodic ? n : n+1;
  neval += Ntot;

  numint::swap(result,r);
  
  return RET;

}


#define CUBATURE(N,TYPE)                                                \
  int cube_tr(const mdfunction<TYPE, N> &f, const numint::array<double,N> &a, \
              const numint::array<double,N> &b, double epsabs, double epsrel, \
              TYPE &result, unsigned &neval) {                          \
    return _cube_tr_<TYPE,N>(f,a,b,epsabs,epsrel,result,neval,0,false, NULL); \
  }                                                                     \
  int cube_ptr(const mdfunction<TYPE, N> &f, const numint::array<double,N> &a, \
               const numint::array<double,N> &b, double epsabs, double epsrel, \
               TYPE &result, unsigned &neval) {                         \
    return _cube_tr_<TYPE,N>(f,a,b,epsabs,epsrel,result,neval,0,true, NULL); \
  }
#define CUBATURE_LOOP(TYPE)                             \
  CUBATURE(1,TYPE); CUBATURE(2,TYPE); CUBATURE(3,TYPE); \
  CUBATURE(4,TYPE); CUBATURE(5,TYPE); CUBATURE(6,TYPE);

CUBATURE_LOOP(double);
CUBATURE_LOOP(dcomplex);
#undef CUBATURE_LOOP
#undef CUBATURE

#define CUBATURE(N,TYPE)                                                \
  int cube_tr(const mdfunction<TYPE, N> &f, const numint::array<double,N> &a, \
              const numint::array<double,N> &b, double epsabs, double epsrel, \
              TYPE &result, unsigned &neval, int dirty,                 \
              Convergence<TYPE>::type compare) {               \
    return _cube_tr_<TYPE,N>(f,a,b,epsabs,epsrel,result,neval,dirty,false,compare); \
  }                                                                     \
  int cube_ptr(const mdfunction<TYPE, N> &f, const numint::array<double,N> &a, \
               const numint::array<double,N> &b, double epsabs, double epsrel, \
               TYPE &result, unsigned &neval, int dirty,                \
               Convergence<TYPE>::type compare) {              \
    return _cube_tr_<TYPE,N>(f,a,b,epsabs,epsrel,result,neval,0,true,compare); \
  }

#define CUBATURE_LOOP(TYPE)                             \
  CUBATURE(1,TYPE); CUBATURE(2,TYPE); CUBATURE(3,TYPE); \
  CUBATURE(4,TYPE); CUBATURE(5,TYPE); CUBATURE(6,TYPE);

CUBATURE_LOOP(vector_d);
CUBATURE_LOOP(vector_z);

#undef CUBATURE_LOOP
#undef CUBATURE

CLOSE_NUMINT_NAMESPACE
