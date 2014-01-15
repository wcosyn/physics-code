/**
 * @file lib/numint/numint_tools.hpp
 * @brief simple functions used in the numint library
 *
 * @author Klaas Vantournhout
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt)
 * @author k.vantournhout@gsi.de
 * @date creation: 11/04/2011
 * @date last update: 11/04/2011
 */

#ifndef NUMINT_TOOLS_HPP
#define NUMINT_TOOLS_HPP 1

#include <numint/macros.hpp>
#include <numint/typedef.hpp>
#include <iomanip>
#include <iostream>
#include <cstdlib>

OPEN_NUMINT_NAMESPACE

/** 
 * @brief convergence function
 * checks if two values are numerically comparable

 * @tparam T type of the values
 * @param r old value
 * @param s new value
 * @param epsabs absolute error
 * @param epsrel relative error
 * @param dirty not in use
 */
template<typename T>
inline bool convergence(const T& r, const T& s, double epsabs, double epsrel, int dirty) {
  double error = fabs(r-s);
  return error < epsabs || error < epsrel * fabs(s);
}

/** 
 * @brief convergence function
 * checks if two complex values are numerically comparable
 * We make use of the L1 norm to compare both values
 * @tparam T type of the complex values
 * @param r old value
 * @param s new value
 * @param epsabs absolute error
 * @param epsrel relative error
 * @param dirty not in use
 */
template<typename T>
inline bool convergence(const std::complex<T> &r, const std::complex<T> &s,
                        double epsabs, double epsrel, int dirty) {
  double error = fabs(real(r-s)) + fabs(imag(r-s));
  return error < epsabs || error < epsrel * (fabs(real(s))+fabs(imag(s)));
  
  // return convergence(r.real(),s.real(),epsabs,epsrel,dirty) &&
  //   convergence(r.imag(),s.imag(),epsabs,epsrel,dirty);
}

/** 
 * @brief convergence function
 * checks if two vectors of values are numerically comparable

 * @tparam T type of the vector values
 * @param r old vector
 * @param s new vector
 * @param epsabs absolute error
 * @param epsrel relative error
 * @param dirty compare the first dirty values, if dirty is negative check all.
 */

template<typename T>
inline bool convergence (const std::vector<T> &r, const std::vector<T> &s,
                         double epsabs, double epsrel, int dirty) {
  bool conv = true;
  if (dirty <= 0) dirty = s.size();
  if (s.size() == 0) return false;
  for (unsigned i = 0; i < (unsigned) dirty && conv; ++i)
    conv = conv && convergence(r[i],s[i],epsabs,epsrel,dirty);
  return conv;
}

/// simple axpy-like function
template<typename T>
inline void axpy(T alpha, const T &X, T &Y) { Y += alpha*X; }
/// simple axpy-like function
template<typename T>
inline void axpy(T alpha, const std::complex<T> &X, std::complex<T> &Y) { Y += alpha*X; }
/// simple axpy-like function
template<typename T>
inline void axpy(T alpha, const std::vector<T> &X, std::vector<T> &Y) {
  for (unsigned i = 0; i < X.size(); ++i)
    Y[i] += alpha * X[i];
}
template<typename T>
inline void axpy(T alpha, const std::vector<std::complex<T> > &X,
                 std::vector<std::complex<T> > &Y) {
  for (unsigned i = 0; i < X.size(); ++i)
    Y[i] += alpha * X[i];
}

/// simple scal-like function
template<typename T>
inline void scal(T alpha, T &X) { X *= alpha; }
/// simple scal-like function
template<typename T>
inline void scal(T alpha, std::complex<T> &X) { X *= alpha; }
/// simple scal-like function
template<typename T>
inline void scal(T alpha, std::vector<T> &X) {
  for (unsigned i = 0; i < X.size(); ++i)
    X[i] *= alpha;
}
/// simple scal-like function
template<typename T>
inline void scal(T alpha, std::vector<std::complex<T> > &X) {
  for (unsigned i = 0; i < X.size(); ++i)
    X[i] *= alpha;
}

/// simple swap
template<typename T>
inline void swap(T &x, T&y) { T t = x; x = y; y = t; }

/// simple swap
template<typename T>
inline void swap(std::vector<T> &x, std::vector<T> &y) { x.swap(y); }

/// simple print
template<typename T>
inline void print(T &x) { std::cout << std::setprecision(15) << x; }

template<typename T>
inline void print(std::vector<T> &x) { std::cout << "[";
  for (unsigned i = 0; i < x.size()-1; ++i)
    std::cout << std::setprecision(15) << x[i] << ",";
  std::cout << std::setprecision(15) << x[x.size()-1] << "]";
}

CLOSE_NUMINT_NAMESPACE

#endif // NUMINT_TOOLS_HPP
