/**
 * @file array.hpp
 * @brief c++ version of the c-style fixed-sized array
 *
 * The definition of a c++ version of the c-style fixed-sized array
 *
 * @author Klaas Vantournhout
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt)
 * @author k.vantournhout@gsi.de
 * @date 16/03/2011
 */


#ifndef LINALG_ARRAY_HPP
#define LINALG_ARRAY_HPP 1

#include <cstring>
#include <algorithm>

namespace numint {

/* =========================================================================
   The array class definition
   ========================================================================= */
/**
 * \class array
 * \brief a c++ version of the C-style fixed-sized array
 * \note
 * This class has no defined constructors and destructors.
 * This allows direct initialisation
 * \code
 * array<int,3> a = {{1,2,3}};
 * \endcode
 */

template<class T, unsigned N>
class array {
public:
  T data_[N];    ///< fixed-size array of elements of type T

private:
  typedef array<T,N> self_type; ///< self explanitory

public:

  /* -----------------------------------------------------------------------
     typedef
     ----------------------------------------------------------------------- */ 

  typedef unsigned size_type; ///< self explanitory
  typedef int difference_type;  ///< self explanitory

  typedef T value_type; ///< self explanitory

  typedef value_type *pointer; ///< self explanitory
  typedef value_type &reference; ///< self explanitory
  typedef const value_type *const_pointer; ///< self explanitory
  typedef const value_type &const_reference; ///< self explanitory

  typedef pointer       iterator; ///< self explanitory
  typedef const_pointer const_iterator; ///< self explanitory

  /* -----------------------------------------------------------------------
     data access
     ----------------------------------------------------------------------- */  
  /// data access C-style array
  const_reference operator [] (size_type i) const {
    return data_[i];
  }
  /// data access C-style array
  reference operator [] (size_type i) {
    return data_[i];
  }
  /// data access
  const_reference operator () (size_type i) const {
    return data_[i];
  }
  /// data access
  reference operator () (size_type i) {
    return data_[i];
  }

  /* -----------------------------------------------------------------------
     get size
     ----------------------------------------------------------------------- */  
  /// return size
  static size_type size(void) { return N; }
  /* -----------------------------------------------------------------------
     assignments
     ----------------------------------------------------------------------- */
  /**
   * \fn self_type &operator = (const self_type &a)
   * \brief copies one array into an other
   * \param a the array to be copied
   * \fn self_type &assign_temporary(self_type &a)
   * \brief swaps a temporary array
   * \param a the array to swap with
   * \fn void swap (self_type &y)
   * \brief swaps two arrays
   * \param y the array to swap with
   */
  self_type &operator = (const self_type &a) {
    if (this != &a)
      memcpy(data_,a.data_,N*sizeof(value_type));
    return *this; // concat
  }

  self_type &assign_temporary(self_type &a) {
    swap(a);
    return *this; // concat
  }

  void swap (array<T,N>& y) {
    for (size_type i = 0; i < N; ++i)
      std::swap(data_[i],y.data_[i]);
  }

  /* -----------------------------------------------------------------------
     check zero
     ----------------------------------------------------------------------- */
  /// set all values to binary zero
  void zero(void) {
    memset(data_,0,N*sizeof(T));
  }
  /// checks if all values are zero
  bool iszero(void) const {
    bool tmp(true);
    for (size_type i = 0; i < N && tmp; ++i)
      tmp = data_[i] == 0.0;
    return tmp;
  }
  /// initialize all values to a constant
  void init(T a) {
    for (size_type i = 0; i < N; ++i) data_[i] = a;
  }

  /* -----------------------------------------------------------------------
     conversion
     ----------------------------------------------------------------------- */
  /// conversion from array<T,N> into array<U,N>
  template<typename U>
  operator array<U,N> (void) const {
    array<U,N> t;
    for (size_type i = 0; i < N; ++i)
      t[i] = (U) (*this)[i];
    return t;
  }

  /* -----------------------------------------------------------------------
     iterators
     ----------------------------------------------------------------------- */
  /**
   * \fn iterator begin()
   * \brief returns an iterator pointing to the first element of an array
   * \fn const_iterator begin() const
   * \brief returns a const_iterator pointing to the first element of a const array
   * \fn iterator end()
   * \brief returns a iterator pointing to the last element of an array
   * \fn const_iterator end() const
   * \brief returns a const_iterator pointing to the last element of a const array
   */
  iterator begin() { return data_; }
  const_iterator begin() const { return data_; }
  iterator end() { return data_+N; }
  const_iterator end() const { return data_+N; }

};

/* =========================================================================
   comparison operators related to the array class
   ========================================================================= */

/**
   @fn bool linalg::operator== (const array& x, const array& y)
   @brief Equal-to operator
   @param[in] x an array
   @param[in] y an array
   @return Returns true if both arrays are identical when compared
   element-by-element, and otherwise returns false

   @fn bool operator< (const array& x, const array& y)
   @brief Less-than operator
   @param[in] x an array 
   @param[in] y an array 
   @return Returns true if array x is lexicographical less then array y,
   and otherwise returns false.
*/


template<class T, unsigned N>
inline bool operator== (const array<T,N>& x, const array<T,N>& y) {
  for (unsigned i = 0; i < N; ++i)
    if (x[i] != y[i]) return false;
  return true;
  //  return std::equal(x.begin(), x.end(), y.begin());
}
template<class T, unsigned N>
inline bool operator< (const array<T,N>& x, const array<T,N>& y) {
  for (unsigned i = 0; i < N; ++i) {
    if (x[i] < y[i]) return true;
    else if (y[i] < x[i]) return false;
  }
  return false;
  //  return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
}
template<class T, unsigned N>
inline bool operator!= (const array<T,N>& x, const array<T,N>& y) {
  return !(x==y);
}
template<class T, unsigned N>
inline bool operator> (const array<T,N>& x, const array<T,N>& y) {
  return y<x;
}
template<class T, unsigned N>
inline bool operator<= (const array<T,N>& x, const array<T,N>& y) {
  return !(y<x);
}
template<class T, unsigned N>
inline bool operator>= (const array<T,N>& x, const array<T,N>& y) {
  return !(x<y);
}

template<class T, unsigned N>
inline array<T,N> operator- (const array<T,N> &x) {
  array<T,N> t;
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t[n] = -x[n];
  return t;
}


/* =========================================================================
   arithmetic operators related to the array class
   ========================================================================= */

template<class T, class U, unsigned N>
inline array<T,N> operator+ (const array<T,N> &x, const array<U,N> &y) {
  array<T,N> t;
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t[n] = x[n] + y[n];
  return t;
}

template<class T, unsigned N>
inline array<T,N> operator- (const array<T,N> &x, const array<T,N> &y) {
  array<T,N> t;
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t[n] = x[n] - y[n];
  return t;
}

template<class T, unsigned N>
inline array<T,N> operator* (const array<T,N> &x, const array<T,N> &y) {
  array<T,N> t;
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t[n] = x[n] * y[n];
  return t;
}

template<class T, unsigned N>
inline array<T,N> operator/ (const array<T,N> &x, const array<T,N> &y) {
  array<T,N> t;
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t[n] = x[n] / y[n];
  return t;
}

template<class T, unsigned N>
inline array<T,N> operator* (const array<T,N> &x, const T &y) {
  array<T,N> t;
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t[n] = x[n] * y;
  return t;
}

template<class T, unsigned N>
inline array<T,N> operator* (const T &x, const array<T,N> &y) {
  array<T,N> t;
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t[n] = x * y[n];
  return t;
}
 
template<class T, unsigned N>
inline array<T,N> operator/ (const array<T,N> &x, const T &y) {
  array<T,N> t;
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t[n] = x[n] / y;
  return t;
}

template<class T, unsigned N>
inline array<T,N> operator/ (const T &x, const array<T,N> &y) {
  array<T,N> t;
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t[n] = x / y[n];
  return t;
}


template<class T, unsigned N>
inline T dot(const array<T,N> &x, const array<T,N> &y) {
  T t(0);
  for (typename array<T,N>::size_type n = 0; n < N; ++n)
    t += x[n] * y[n];
  return t;
}

template <class T, unsigned I, unsigned J>
inline array<T,I*J> KroneckerProduct(const array<T,I> &x, const array<T,J> &y) {
  array<T,I*J> t;
  for (typename array<T,I>::size_type i = 0; i < I; ++i)
    for (typename array<T,J>::size_type j = 0; j < J; ++j)
      t(i*J+j) = x(i) * y(j);
  return t;
}

template <class T, unsigned N>
inline void split(const array<T,2*N> &x, array<T,N> &y1, array<T,N> &y2) {
  for (typename array<T,N>::size_type i = 0; i < N; ++i) {
    y1[i] = x[i]; y2[i] = x[i+N];
  }
};

}

#endif
