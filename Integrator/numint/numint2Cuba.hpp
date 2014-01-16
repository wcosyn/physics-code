/**
 * @file numint2Cuba.hpp 
 *
 * @brief This file is a wrapper that uses the same style of function
 * calling but implements the Cuba library which is obtained
 * from the code that is found on :
 * http://www.feynarts.de/cuba/
 * The code of this wrapper is heavily inspired on the one 
 * Klaas Vantournhout wrote > numint2adapt.hpp
 * 
 * ! currently supported CUBA version is v3.3 (08/01/2014) !
 * 
 * @author Camille Colle
 * @author Ghent University
 */

#ifndef NUMINT_CUBA_CONVERSION
#define NUMINT_CUBA_CONVERSION 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <numint/macros.hpp>
#include <numint/typedef.hpp>
#include <numint/numint_tools.hpp>

#include <Cuba/cuba.h> // include the header with the actual integration declarations

OPEN_NUMINT_NAMESPACE

template <typename T>
inline unsigned get_cuba_size(const T &a) { return 0; }

template <>
inline unsigned get_cuba_size(const double &a) { return 1; }

template <>
inline unsigned get_cuba_size(const dcomplex &a) { return 2; }

template <>
inline unsigned get_cuba_size(const vector_d &a) { return a.size(); }

template <>
inline unsigned get_cuba_size(const vector_z &a) { return 2 * a.size(); }

template <typename T>
inline void numint2cuba(unsigned fdim, double *res, T&a) {}

template <>
inline void numint2cuba(unsigned fdim, double *res, double &a) {res[0] = a;}

template <>
inline void numint2cuba(unsigned fdim, double *res, dcomplex &a) {res[0] = a.real(); res[1] = a.imag();}

template <>
inline void numint2cuba(unsigned fdim, double *res, vector_d &a) {
  for (unsigned i = 0; i < fdim; ++i)
    res[i] = a[i];
}

template <>
inline void numint2cuba(unsigned fdim, double *res, vector_z &a) {
  for (unsigned i = 0; i < fdim/2; ++i) {
    res[2*i  ] = a[i].real();
    res[2*i+1] = a[i].imag();
  }
}

template <typename T>
inline void cuba2numint(unsigned fdim, double *res, T &a) {}

template<>
inline void cuba2numint(unsigned fdim, double *res, double &a) { a = res[0]; } // *res and res[0] are equivalent

template<>
inline void cuba2numint(unsigned fdim, double *res, dcomplex &a) { a=dcomplex(res[0],res[1]); }

template<>
inline void cuba2numint(unsigned fdim, double *res, vector_d &a) {
  for (unsigned i = 0; i < fdim; ++i) a[i] = res[i];
}

template<>
inline void cuba2numint(unsigned fdim, double *res, vector_z &a) {
  for (unsigned i = 0; i < fdim/2; ++i) {
    a[i] = dcomplex(res[2*i],res[2*i+1]);
  }
}

template<typename T,unsigned N>
struct mdf2cuba{
  static int exec(const int* ndim, const double* xx, const int* ncomp, double* ff, void *fdata)
  {
    mdf2cuba &p = *(mdf2cuba*) fdata;
    
    numint::array<double,N> x; // N should be equal to ndim!
    for (unsigned i=0; i<N; i++) x[i] = xx[i]*(p.b[i]-p.a[i])+p.a[i]; // remember cuba integrates in a unit hypercube
    FN_EVAL(p.f,x,p.ret); // evaluate function p.f @ coordinates x and store result in ret
    numint2cuba(*ncomp,ff,p.ret); // convert the <T> return value to something cuba integrator can use ff (array of doubles)
    for (unsigned i=0; i<N; i++) { // Jacobian of unit hypercube transformation
      for (int j=0; j<*ncomp; j++) {
          ff[j] *= (p.b[i]-p.a[i]);
      }
    }
    return 0; // No error checks, always "succes"
  }
  mdfunction<T,N> f;
  numint::array<double, N> a,b; // boundaries, needed to calculate transformation to and from hypercube
  T ret;
};

template<typename T, unsigned N>
void cuhre( const mdfunction<T,N> &f, const numint::array<double, N> &a,
           const numint::array<double,N> &b,int nvec, double epsrel, double epsabs,int flags,
           int minEval, int maxEval,int key, char* statefile,int &nregions, int &neval, int &fail, T &result,double &err, double &prob ){
  
  struct mdf2cuba<T,N> F;
  F.f = f; // set the integrand
  F.a = a; // set the lower boundary
  F.b = b; // set the upper boundary
  unsigned fdim;
  FN_EVAL<T,N>(f,a,result);
  fdim = get_cuba_size(result); // call function once here to find out the size of fdim, for example complex value gets treated as a two component function
  
  double *val = new double[fdim];
  
  Cuhre(N, fdim, mdf2cuba<T,N>::exec, &F, nvec,
        epsrel, epsabs, flags,
        minEval, maxEval, key, statefile,
        &nregions, &neval, &fail, val, &err, &prob);
  
  cuba2numint(fdim,val,result); // convert the multicomponent function back to normal variables like complex ones
  delete [] val; // after conversion from val to result we don't need this anymore
}


template<typename T, unsigned N>
void divonne( const mdfunction<T,N> &f, const numint::array<double, N> &a,
           const numint::array<double,N> &b,int nvec, double epsrel, double epsabs,int flags, int seed,
           int minEval, int maxEval,int key1, int key2, int key3, int maxpass, double border, double maxchisq,
           double mindeviation, int ngiven, int ldxgiven,double *xgiven, int nextra, peakfinder_t peakfinder,
           char* statefile,int &nregions, int &neval, int &fail, T &result,double &err, double &prob ){
  
  struct mdf2cuba<T,N> F;
  F.f = f; // set the integrand
  F.a = a; // set the lower boundary
  F.b = b; // set the upper boundary
  unsigned fdim;
  FN_EVAL<T,N>(f,a,result);
  fdim = get_cuba_size(result); // call function once here to find out the size of fdim, for example complex value gets treated as a two component function
  
  double *val = new double[fdim];
  
  Divonne(N, fdim, mdf2cuba<T,N>::exec, &F, nvec,
        epsrel, epsabs, flags,seed,
        minEval, maxEval, key1, key2, key3, maxpass,border,
        maxchisq,mindeviation,ngiven,ldxgiven,xgiven,nextra,peakfinder,
        statefile, &nregions, &neval, &fail, val, &err, &prob);
  
  cuba2numint(fdim,val,result); // convert the multicomponent function back to normal variables like complex ones
  delete [] val; // after conversion from val to result we don't need this anymore
}

template<typename T, unsigned N>
void suave( const mdfunction<T,N> &f, const numint::array<double, N> &a,
           const numint::array<double,N> &b,int nvec, double epsrel, double epsabs,int flags, int seed,
           int minEval, int maxEval,int nnew, double flatness, char* statefile, int &nregions, int &neval, int &fail, T &result,double &err, double &prob ){
  
  struct mdf2cuba<T,N> F;
  F.f = f; // set the integrand
  F.a = a; // set the lower boundary
  F.b = b; // set the upper boundary
  unsigned fdim;
  FN_EVAL<T,N>(f,a,result);
  fdim = get_cuba_size(result); // call function once here to find out the size of fdim, for example complex value gets treated as a two component function
  
  double *val = new double[fdim];
  
  Suave(N, fdim, mdf2cuba<T,N>::exec, &F, nvec,
        epsrel, epsabs, flags,seed,
        minEval, maxEval,nnew,flatness,
        statefile,&nregions, &neval, &fail, val, &err, &prob);
  
  cuba2numint(fdim,val,result); // convert the multicomponent function back to normal variables like complex ones
  delete [] val; // after conversion from val to result we don't need this anymore
}

template<typename T, unsigned N>
void vegas( const mdfunction<T,N> &f, const numint::array<double, N> &a,
           const numint::array<double,N> &b,int nvec, double epsrel, double epsabs,int flags, int seed,
           int minEval, int maxEval,int nstart, int nincrease, int nbatch, int gridno,char* statefile, int &neval, int &fail, T &result,double &err, double &prob ){
  
  struct mdf2cuba<T,N> F;
  F.f = f; // set the integrand
  F.a = a; // set the lower boundary
  F.b = b; // set the upper boundary
  unsigned fdim;
  FN_EVAL<T,N>(f,a,result);
  fdim = get_cuba_size(result); // call function once here to find out the size of fdim, for example complex value gets treated as a two component function
  
  double *val = new double[fdim];
  
  Vegas(N, fdim, mdf2cuba<T,N>::exec, &F, nvec,
        epsrel, epsabs, flags,seed,
        minEval, maxEval, nstart,nincrease,nbatch,gridno,
        statefile, &neval, &fail, val, &err, &prob);
  
  cuba2numint(fdim,val,result); // convert the multicomponent function back to normal variables like complex ones
  delete [] val; // after conversion from val to result we don't need this anymore
}

CLOSE_NUMINT_NAMESPACE

#endif // NUMINT_CUBA_CONVERSION