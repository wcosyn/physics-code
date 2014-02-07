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

/*template <typename T>
inline unsigned get_cuba_size(const T &a) { return 0; }*/

//template <>
inline unsigned get_cuba_size(const double &a) { return 1; }

//template <>
inline unsigned get_cuba_size(const dcomplex &a) { return 2; }

//template <>
inline unsigned get_cuba_size(const vector_d &a) { return a.size(); }

//template <>
inline unsigned get_cuba_size(const vector_z &a) { return 2 * a.size(); }

template <unsigned N>
inline unsigned get_cuba_size(const numint::array<double,N> &a) { return a.size(); }

/*template <typename T>
inline void numint2cuba(unsigned fdim, double *res, T&a) {}*/

//template <>
inline void numint2cuba(unsigned fdim, double *res, double &a) {res[0] = a;}

//template <>
inline void numint2cuba(unsigned fdim, double *res, dcomplex &a) {res[0] = a.real(); res[1] = a.imag();}

//template <>
inline void numint2cuba(unsigned fdim, double *res, vector_d &a) {
  for (unsigned i = 0; i < fdim; ++i)
    res[i] = a.at(i); // important to use the vector::at rather than [] acces operator because vector::at does range checking! e.g. This will throw errors if you pass incorrectly initialized vectors!
}

//template <>
inline void numint2cuba(unsigned fdim, double *res, vector_z &a) {
  for (unsigned i = 0; i < fdim/2; ++i) {
    res[2*i  ] = a.at(i).real();
    res[2*i+1] = a.at(i).imag(); // important to use the vector::at rather than [] acces operator because vector::at does range checking! e.g. This will throw errors if you pass incorrectly initialized vectors!
  }
}

template <unsigned N>
inline void numint2cuba(unsigned fdim,double *res, numint::array<double,N> &a){
  for (unsigned i=0; i<N; i++){
    res[i] = a[i];
  }
}

/*template <typename T>
inline void cuba2numint(unsigned fdim, double *res, T &a) {}*/

//template<>
inline void cuba2numint(unsigned fdim, double *res, double &a) { a = res[0]; } // *res and res[0] are equivalent

//template<>
inline void cuba2numint(unsigned fdim, double *res, dcomplex &a) { a=dcomplex(res[0],res[1]); }

//template<>
inline void cuba2numint(unsigned fdim, double *res, vector_d &a) {
  for (unsigned i = 0; i < fdim; ++i) { a.at(i) = res[i]; } // important to use the vector::at rather than [] acces operator because vector::at does range checking! e.g. This will throw errors if you pass incorrectly initialized vectors!
}

//template<>
inline void cuba2numint(unsigned fdim, double *res, vector_z &a) {
  for (unsigned i = 0; i < fdim/2; ++i) {
    a.at(i) = dcomplex(res[2*i],res[2*i+1]); // important to use the vector::at rather than [] acces operator because vector::at does range checking! e.g. This will throw errors if you pass incorrectly initialized vectors!
  }
}

template<unsigned N>
inline void cuba2numint(unsigned fdim, double *res, numint::array<double,N> &a){
    for (unsigned i=0; i<fdim; i++){
        a[i] = res[i];
    }
}

template<typename T,unsigned N>
struct mdf2cuba{
  static int exec(const int* ndim, const double* xx, const int* ncomp, double* ff, void *fdata)
  {
    mdf2cuba &p = *(mdf2cuba*) fdata;
    numint::array<double,N> x; // N should be equal to ndim!
    for (unsigned i=0; i<N; i++) x[i] = xx[i]*(p.b[i]-p.a[i])+p.a[i]; // remember cuba integrates in a unit hypercube
    FN_EVAL(p.f,x,p.ret); // evaluate mdfunction p.f @ coordinates x and store result in ret, eg: mdf.func(x,mdf.param,ret)
    numint2cuba(*ncomp,ff,p.ret); // convert the <T> return value to something cuba integrator can use ff (array of doubles)
    for (unsigned i=0; i<N; i++) { // Jacobian of unit hypercube transformation
      for (int j=0; j<*ncomp; j++) {
          ff[j] *= (p.b[i]-p.a[i]);
      }
    }
    
    return 0; // No error checks, always "succes"
  }
  void set_ncomp(int n) {} // here this method does nothing because numint::array covered by typename T is auto correctly initialized, but keep it here to comply with required interface
  mdfunction<T,N> f;
  numint::array<double, N> a,b; // boundaries, needed to calculate transformation to and from hypercube
  T ret;
};

/** specialisation mdf2cuba if T is a vector, because vector needs initialisation by resizing!
 *  it is quite embarrassing that this is necessary but to i did it to comply with the interface
 *  the core problem is that the number of dimensions and components are handled completely differently...
 * 
 */
template< typename T, unsigned N>
struct mdf2cuba<std::vector<T>,N >{
  static int exec(const int* ndim, const double* xx, const int* ncomp, double* ff, void *fdata)
  {
    mdf2cuba &p = *(mdf2cuba*) fdata;
    numint::array<double,N> x; // N should be equal to ndim!
    for (unsigned i=0; i<N; i++) x[i] = xx[i]*(p.b[i]-p.a[i])+p.a[i]; // remember cuba integrates in a unit hypercube
    FN_EVAL(p.f,x,p.ret); // evaluate mdfunction p.f @ coordinates x and store result in ret, eg: mdf.func(x,mdf.param,ret)
    numint2cuba(*ncomp,ff,p.ret); // convert the <T> return value to something cuba integrator can use ff (array of doubles)
    for (unsigned i=0; i<N; i++) { // Jacobian of unit hypercube transformation
      for (int j=0; j<*ncomp; j++) {
          ff[j] *= (p.b[i]-p.a[i]);
      }
    }
    
    return 0; // No error checks, always "succes"
  }
  mdfunction<std::vector<T>,N> f;
  numint::array<double, N> a,b; // boundaries, needed to calculate transformation to and from hypercube
  void set_ncomp(unsigned n){ ret.resize(n); } // because ret is a vector here we must explicitly resize it to the required size, resize method is garantueed to exist on type vector
  std::vector<T> ret;
};

template<typename T, unsigned N>
void cuhre( const mdfunction<T,N> &f, const numint::array<double, N> &a,
           const numint::array<double,N> &b,int nvec, double epsrel, double epsabs,int flags,
           int minEval, int maxEval,int key, char* statefile,int &nregions, int &neval, int &fail, T &result, T &err, T &prob ){
  
  struct mdf2cuba<T,N> F;
  F.f = f; // set the integrand
  F.a = a; // set the lower boundary
  F.b = b; // set the upper boundary
  
  FN_EVAL<T,N>(f,a,result);
  unsigned ncomp = get_cuba_size(result); // call function once here to find out the size of ncomp, for example complex value gets treated as a two component function
  F.set_ncomp(ncomp); // this only does something if you passed vector<...> as template T
  
  double *result_t = new double[ncomp];
  double *err_t    = new double[ncomp]; // error has same dimensions as result
  double *prob_t   = new double[ncomp]; // prob  has same dimensions as result
  
  Cuhre(N, ncomp, mdf2cuba<T,N>::exec, &F, nvec,
        epsrel, epsabs, flags,
        minEval, maxEval, key, statefile,
        &nregions, &neval, &fail, result_t, err_t, prob_t);
  
  cuba2numint(ncomp,result_t,result); // convert the multicomponent function back to normal variables like complex ones
  cuba2numint(ncomp,err_t,err);       // convert the multicomponent error back to the template variables
  cuba2numint(ncomp,prob_t,prob);     // convert the multicomponent prob back to the template variables
  delete [] result_t; // after conversion from the *_t arrays we don't need these anymore
  delete [] err_t;
  delete [] prob_t;
}


template<typename T, unsigned N>
void divonne( const mdfunction<T,N> &f, const numint::array<double, N> &a,
           const numint::array<double,N> &b,int nvec, double epsrel, double epsabs,int flags, int seed,
           int minEval, int maxEval,int key1, int key2, int key3, int maxpass, double border, double maxchisq,
           double mindeviation, int ngiven, int ldxgiven,double *xgiven, int nextra, peakfinder_t peakfinder,
           char* statefile,int &nregions, int &neval, int &fail, T &result,T &err,T &prob ){
  
  struct mdf2cuba<T,N> F;
  F.f = f; // set the integrand
  F.a = a; // set the lower boundary
  F.b = b; // set the upper boundary
  
  FN_EVAL<T,N>(f,a,result);
  unsigned ncomp = get_cuba_size(result); // call function once here to find out the size of ncomp, for example complex value gets treated as a two component function
  F.set_ncomp(ncomp); // this only does something if you passed vector<...> as template T
  
  double *result_t = new double[ncomp];
  double *err_t    = new double[ncomp]; // error has same dimensions as result
  double *prob_t   = new double[ncomp]; // prob  has same dimensions as result
  
  Divonne(N, ncomp, mdf2cuba<T,N>::exec, &F, nvec,
        epsrel, epsabs, flags,seed,
        minEval, maxEval, key1, key2, key3, maxpass,border,
        maxchisq,mindeviation,ngiven,ldxgiven,xgiven,nextra,peakfinder,
        statefile, &nregions, &neval, &fail, result_t, err_t, prob_t);
  
  cuba2numint(ncomp,result_t,result); // convert the multicomponent function back to normal variables like complex ones
  cuba2numint(ncomp,err_t,err);       // convert the multicomponent error back to the template variables
  cuba2numint(ncomp,prob_t,prob);     // convert the multicomponent prob back to the template variables
  delete [] result_t; // after conversion from the *_t arrays we don't need these anymore
  delete [] err_t;
  delete [] prob_t;
}

template<typename T, unsigned N>
void suave( const mdfunction<T,N> &f, const numint::array<double, N> &a,
           const numint::array<double,N> &b,int nvec, double epsrel, double epsabs,int flags, int seed,
           int minEval, int maxEval,int nnew, double flatness, char* statefile, int &nregions, int &neval, int &fail, T &result,T &err, T &prob ){
  
  struct mdf2cuba<T,N> F;
  F.f = f; // set the integrand
  F.a = a; // set the lower boundary
  F.b = b; // set the upper boundary
  
  FN_EVAL<T,N>(f,a,result);
  unsigned ncomp = get_cuba_size(result); // call function once here to find out the size of ncomp, for example complex value gets treated as a two component function
  F.set_ncomp(ncomp); // this only does something if you passed vector<...> as template T
  
  double *result_t = new double[ncomp];
  double *err_t    = new double[ncomp]; // error has same dimensions as result
  double *prob_t   = new double[ncomp]; // prob  has same dimensions as result
  
  Suave(N, ncomp, mdf2cuba<T,N>::exec, &F, nvec,
        epsrel, epsabs, flags,seed,
        minEval, maxEval,nnew,flatness,
        statefile,&nregions, &neval, &fail, result_t, err_t, prob_t);
  
  cuba2numint(ncomp,result_t,result); // convert the multicomponent function back to normal variables like complex ones
  cuba2numint(ncomp,err_t,err);       // convert the multicomponent error back to the template variables
  cuba2numint(ncomp,prob_t,prob);     // convert the multicomponent prob back to the template variables
  delete [] result_t; // after conversion from the *_t arrays we don't need these anymore
  delete [] err_t;
  delete [] prob_t;
}

template<typename T, unsigned N>
void vegas( const mdfunction<T,N> &f, const numint::array<double, N> &a,
           const numint::array<double,N> &b,int nvec, double epsrel, double epsabs,int flags, int seed,
           int minEval, int maxEval,int nstart, int nincrease, int nbatch, int gridno,char* statefile, int &neval, int &fail, T &result,T &err, T &prob ){
  
  struct mdf2cuba<T,N> F;
  F.f = f; // set the integrand
  F.a = a; // set the lower boundary
  F.b = b; // set the upper boundary
  
  FN_EVAL<T,N>(f,a,result);
  unsigned ncomp = get_cuba_size(result); // call function once here to find out the size of ncomp, for example complex value gets treated as a two component function
  F.set_ncomp(ncomp); // this only does something if you passed vector<...> as template T
  
  
  double *result_t = new double[ncomp];
  double *err_t    = new double[ncomp]; // error has same dimensions as result
  double *prob_t   = new double[ncomp]; // prob  has same dimensions as result
  
  Vegas(N, ncomp, mdf2cuba<T,N>::exec, &F, nvec,
        epsrel, epsabs, flags,seed,
        minEval, maxEval, nstart,nincrease,nbatch,gridno,
        statefile, &neval, &fail, result_t, err_t, prob_t);
  
  cuba2numint(ncomp,result_t,result); // convert the multicomponent function back to normal variables like complex ones
  cuba2numint(ncomp,err_t,err);       // convert the multicomponent error back to the template variables
  cuba2numint(ncomp,prob_t,prob);     // convert the multicomponent prob back to the template variables
  delete [] result_t; // after conversion from the *_t arrays we don't need these anymore
  delete [] err_t;
  delete [] prob_t;
}

CLOSE_NUMINT_NAMESPACE

#endif // NUMINT_CUBA_CONVERSION
