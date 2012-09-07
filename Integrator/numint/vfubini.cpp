/**
 * @file vfubini.cpp
 * @brief Fubini type multi-dimensional vector-integration
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt) (2010-)
 * @author k.vantournhout@gsi.de
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <numint/typedef.hpp>
#include <numint/quadrature.hpp>

OPEN_NUMINT_NAMESPACE

/**
   quad<T>::type is a pointer to a non-vector quadrature rule
*/
template<typename T>
struct quad {
  typedef int (* type)(const function<T> &f, double a, double b,
                       double epsabs, double epsrel, T &, unsigned &neval, int dirty,
                       typename Convergence<T>::type compare);
};
/**
   The set of integraiton parameters for multi-dimensional integration
*/
template<typename T, unsigned N>
struct IntegrationParameters {
  IntegrationParameters (const mdfunction<T,N> &F, const numint::array<double,N> &A,
                         const numint::array<double,N> &B, double EPSABS, double EPSREL,
                         unsigned NEVAL, int DIRTY, typename Convergence<T>::type COMPARE,
                         int ERROR, typename quad<T>::type QUAD) 
    : func(F), a(A), b(B), epsabs(EPSABS), epsrel(EPSREL),
      neval(NEVAL), dirty(DIRTY), compare(COMPARE), error(ERROR), quad(QUAD) {}
  
  const mdfunction<T,N> &func;      ///< the function to integrate
  
  const numint::array<double,N> &a; ///< the lower limits
  const numint::array<double,N> &b; ///< the upper limits
  double epsabs;                    ///< absolute error
  double epsrel;                    ///< relative error
  unsigned neval;                   ///< total function evaluations
  int dirty;                        ///< how dirty should the convergence be
  typename Convergence<T>::type compare; ///< use the following convergence function
  
  int error;                        ///< the error obtained during sub-integrations
  
  typename quad<T>::type quad;      ///< the used quadrature function
  
  numint::array<double,N> x;        ///< the current point in which func is evaluated
};


/**
   evaluates the integrandum of the Mth integral of the multidimensional integral
*/
template<typename T, unsigned M, unsigned N>
struct RecursiveFubini {
  static void exec(double xm, void *param, T &ret) {
    IntegrationParameters<T,N> *p = (IntegrationParameters<T,N> *) param;
    if (M < N) p->x[M] = xm;
    
    function<T> f;
    f.func = &RecursiveFubini<T,M-1,N>::exec;
    f. param = param;
    
    unsigned dummy = 0;
    int error;
    error = (*(p->quad))(f,p->a[M-1],p->b[M-1],p->epsabs,p->epsrel,ret,dummy,p->dirty,p->compare);
    if (error != SUCCESS && p->error == SUCCESS) p->error = error;
  }
};

/**
   specialization : evaluates the function value
*/
template<typename T, unsigned N>
struct RecursiveFubini<T,0,N> {
  static void exec(double x0, void *param, T &ret) {
    IntegrationParameters<T,N> *p = (IntegrationParameters<T,N> *) param;
    p->x[0] = x0;
    FN_EVAL<T,N>(p->func,p->x,ret);
    p->neval++;
  }
};

/**
   starts the general fubini rule
*/
template<typename T, unsigned N>
int _fubini_(const mdfunction<T,N> &f, const numint::array<double,N> &a,
             const numint::array<double,N> &b, double epsabs, double epsrel,
             T &result, unsigned &neval, int dirty,
             typename Convergence<T>::type compare, typename quad<T>::type quad) {
  
  int error = SUCCESS;
  
  IntegrationParameters<T,N> param(f,a,b,epsabs,epsrel,neval,dirty,compare,error,quad);
  
  RecursiveFubini<T,N,N>::exec(0.0,(void *) &param, result);
  neval = param.neval;
  return param.error;
}


#define FUBINI(CC,N,TYPE)                                               \
  int fubini_##CC                                                       \
  (const mdfunction<TYPE,N> &f, const numint::array<double,N> &a,       \
   const numint::array<double,N> &b, double epsabs, double epsrel,      \
   TYPE &result, unsigned &neval, int dirty,                            \
   typename Convergence<TYPE>::type compare) {                          \
    return _fubini_<TYPE,N>(f,a,b,epsabs,epsrel,                        \
                            result,neval,dirty,compare,quad_##CC);      \
  }                                                                     \
  
#define FUBINI_LOOP(CC,TYPE)                                    \
  FUBINI(CC,1,TYPE);  FUBINI(CC,2,TYPE);  FUBINI(CC,3,TYPE);    \
  FUBINI(CC,4,TYPE);  FUBINI(CC,5,TYPE);  FUBINI(CC,6,TYPE);

#define FUBINI_TYPE(CC)                                 \
  FUBINI_LOOP(CC,vector_d); FUBINI_LOOP(CC,vector_z);

FUBINI_TYPE(tr)
FUBINI_TYPE(romb)
FUBINI_TYPE(cc)

FUBINI_TYPE(ptr)
FUBINI_TYPE(promb)
FUBINI_TYPE(pcc)

#undef FUBINI_TYPE
#undef FUBINI_LOOP
#undef FUBINI

CLOSE_NUMINT_NAMESPACE
