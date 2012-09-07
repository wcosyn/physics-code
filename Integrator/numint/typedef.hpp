/**
 * @file typedef.hpp
 * @brief general typedefs
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt) (2010-)
 * @author k.vantournhout@gsi.de
 */


#ifndef NUMINT_TYPEDEF_HPP
#define NUMINT_TYPEDEF_HPP 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <complex>
#include <vector>
#include <numint/macros.hpp>

#include <numint/array.hpp>

OPEN_NUMINT_NAMESPACE

typedef std::complex<float>      complex ; ///< single precision complex number
typedef std::complex<double>     dcomplex; ///< double precision complex number

typedef std::vector<bool>     vector_b; ///< bool dynamic vector
typedef std::vector<int>      vector_i; ///< 32bit dynamic integer vector
typedef std::vector<long>     vector_l; ///< 64bit dynamic integer vector
typedef std::vector<float>    vector_s; ///< single precision dynamic vector
typedef std::vector<double>   vector_d; ///< double precision dynamic vector
typedef std::vector<complex>  vector_c; ///< single precision dynamic complex vector
typedef std::vector<dcomplex> vector_z; ///< double precision dynamic complex vector



///   A functor used to evaluate a function with parameters
template<typename T>
struct function {
  /// a function with parameters
  /// @param x the evaluation point
  /// @param param pointer to parameters
  /// @param ret return value
  void (* func) (double x, void *param, T &ret);
  /// the parameters
  void * param; 
};

/// evaluates the functor
template<typename T>
inline void FN_EVAL(const function<T> &f, double x, T &ret) {
  (*(f.func))(x,f.param,ret);
}

///   A functor used to evaluate a multi-dimensional function with parameters
template<typename T, unsigned N>
struct mdfunction {
  /// a function with parameters
  /// @param x the evaluation point
  /// @param param pointer to parameters
  /// @param ret return value
  void (* func) (const numint::array<double,N> &x, void *param, T &ret);
  /// the parameters
  void * param;
};

/// evaluates the multidimensional functor
template<typename T, unsigned N>
inline void FN_EVAL(const mdfunction<T,N> &f, const numint::array<double,N> &x, T &ret) {
  (*(f.func))(x,f.param,ret);
}

/// convergence functor for type T
template<typename T>
struct Convergence {
  typedef bool (*type) (const T &r, const T &s, double epsabs, double epsrel, int dirty);
};

enum { 
  SUCCESS  = 0, 
  FAILURE  = -1,
  CONTINUE = -2,  /* iteration has not converged */
  EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  EFAULT   = 3,   /* invalid pointer */
  EINVAL   = 4,   /* invalid argument supplied by user */
  EFAILED  = 5,   /* generic failure */
  EFACTOR  = 6,   /* factorization failed */
  ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  ENOMEM   = 8,   /* malloc failed */
  EBADFUNC = 9,   /* problem with user-supplied function */
  ERUNAWAY = 10,  /* iterative process is out of control */
  EMAXITER = 11,  /* exceeded max number of iterations */
  EZERODIV = 12,  /* tried to divide by zero */
  EBADTOL  = 13,  /* user specified an invalid tolerance */
  ETOL     = 14,  /* failed to reach the specified tolerance */
  EUNDRFLW = 15,  /* underflow */
  EOVRFLW  = 16,  /* overflow  */
  ELOSS    = 17,  /* loss of accuracy */
  EROUND   = 18,  /* failed because of roundoff error */
  EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  ENOTSQR  = 20,  /* matrix not square */
  ESING    = 21,  /* apparent singularity detected */
  EDIVERGE = 22,  /* integral or series is divergent */
  EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  ECACHE   = 25,  /* cache limit exceeded */
  ETABLE   = 26,  /* table limit exceeded */
  ENOPROG  = 27,  /* iteration is not making progress towards solution */
  ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  EOF      = 32   /* end of file */
} ;


CLOSE_NUMINT_NAMESPACE

#endif
