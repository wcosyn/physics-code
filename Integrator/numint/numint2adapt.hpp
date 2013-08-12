/**
 * @file numint2adapt.hpp 

 * @brief This file is a wrapper that uses the same style of function
 * calling but implements the adaptive cubuture rule which is obtained
 * fromt he code that is found on :
 * http://ab-initio.mit.edu/wiki/index.php/Cubature
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt)
 * @author k.vantournhout@gsi.de
 */

#ifndef NUMINT_ADAPTIVE_CONVERSION
#define NUMINT_ADAPTIVE_CONVERSION 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <numint/macros.hpp>
#include <numint/typedef.hpp>
#include <numint/numint_tools.hpp>

#include <adaptive/cubature.h>


OPEN_NUMINT_NAMESPACE

template <typename T>
inline unsigned get_ada_size(const T &a) { return 0; }

template <>
inline unsigned get_ada_size(const double &a) { return 1; }

template <>
inline unsigned get_ada_size(const dcomplex &a) { return 2; }

template <>
inline unsigned get_ada_size(const vector_d &a) { return a.size(); }

template <>
inline unsigned get_ada_size(const vector_z &a) { return 2 * a.size(); }

template <typename T>
inline void ada2numint(unsigned fdim, double *res, T &a) {}

template<>
inline void ada2numint(unsigned fdim, double *res, double &a) { a = res[0]; }

template<>
inline void ada2numint(unsigned fdim, double *res, dcomplex &a) { a=dcomplex(res[0],res[1]); } // old standard lang: a.real() = res[0]; a.imag() = res[1];

template<>
inline void ada2numint(unsigned fdim, double *res, vector_d &a) {
  for (unsigned i = 0; i < fdim; ++i) a[i] = res[i];
}

template<>
inline void ada2numint(unsigned fdim, double *res, vector_z &a) {
  for (unsigned i = 0; i < fdim/2; ++i) {
    a[i] = dcomplex(res[2*i],res[2*i+1]);
    //a[i].real() = res[2*i]; old standard lang, functions have become pure getters
    //a[i].imag() = res[2*i+1]; idem
  }
}

template <typename T>
inline void numint2ada(unsigned fdim, double *res, T&a) {}

template <>
inline void numint2ada(unsigned fdim, double *res, double &a) {res[0] = a;}

template <>
inline void numint2ada(unsigned fdim, double *res, dcomplex &a) {res[0] = a.real(); res[1] = a.imag();}

template <>
inline void numint2ada(unsigned fdim, double *res, vector_d &a) {
  for (unsigned i = 0; i < fdim; ++i)
    res[i] = a[i];
}

template <>
inline void numint2ada(unsigned fdim, double *res, vector_z &a) {
  for (unsigned i = 0; i < fdim/2; ++i) {
    res[2*i  ] = a[i].real();
    res[2*i+1] = a[i].imag();
  }
}



template <typename T, unsigned N>
int cube_adaptive(const mdfunction<T,N> &f, const numint::array<double, N> &a,
                  const numint::array<double,N> &b, double epsabs, double epsrel, int maxEval,
                  T &result, unsigned &neval, int dirty,
                  typename Convergence<T>::type compare=NULL) {

  struct mdf2ada {

    static void exec(unsigned ndim, const double *x, void *fdata,
                  unsigned fdim, double *fval) {

      mdf2ada &p = * (mdf2ada *) fdata;

      numint::array<double,N> X;
      for (unsigned i = 0; i < N; ++i) X[i] = x[i];

      FN_EVAL<T,N>(p.f,X,p.ret);

      numint2ada(fdim,fval,p.ret);
      p.Ncall++;
    }

    mdfunction<T,N> f;
    unsigned Ncall;
    T ret;
  };

  mdf2ada F;
  F.f = f;
  F.Ncall = 0;

  //double xmin[N], xmax[N];

  unsigned fdim;

  FN_EVAL<T,N>(f,a,result); F.Ncall++;
  fdim = get_ada_size(result);

  double val[fdim], err[fdim];

  int code = adapt_integrate(fdim, mdf2ada::exec, &F,
                             N,a.begin(),b.begin(),4*pow(10,int(fdim)),
                             maxEval,epsabs,epsrel,
                             val,err);


  ada2numint(fdim,val,result);

  neval = F.Ncall;

  return code;

}

CLOSE_NUMINT_NAMESPACE

#endif
