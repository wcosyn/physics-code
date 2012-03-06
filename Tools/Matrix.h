/*
  Matrix.h
  Class representing a general (K x L) matrix
  The components are complex doubles

  author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)

*/

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <complex>
#include <cstdlib>

//-------------------------------------------------------------------------
//---- Certain classes and functions need to be pre-declared
//---- See C++ FAQ lite 35.16 as to why

template<int K,int L> class Matrix;  // pre-declare Matrix<K,L>
// pre-declare friend template functions of Matrix<K,L>
template<int K,int L> Matrix<K,L> operator*(double,const Matrix<K,L>&);
template<int K,int L> Matrix<K,L> operator*(const std::complex<double>&,
					    const Matrix<K,L>&);
template<int K,int L> std::ostream& operator<<(std::ostream&,const Matrix<K,L>&);

//-------------------------------------------------------------------------
//---- class Matrix<K,L> declaration

template< int K, int L >
class Matrix
{
 private:
  // Data member
  // -----------
  std::complex<double>* m;

 public:
  // Constructors
  // ------------
  Matrix();                                // Default Constructor
  //Matrix(std::complex<double>);               // Ambiguous
  Matrix(std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>); // Constructor for 2x2/4x1/1x4 matrix
  Matrix(std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>); // Constructor for 4x4/8x2/2x8/1x16/16x1 matrix
  Matrix(const Matrix&);                   // copy constructor

  // Destructor
  // ----------
  ~Matrix();

  // Assignment
  // ----------
  Matrix& operator=(const Matrix&);

  // Access individual components
  // ----------------------------
  std::complex<double>& operator()(int,int);
  const std::complex<double> operator()(int,int) const;

  // Addition
  // --------
  Matrix& operator+=(const Matrix&);
  Matrix& operator-=(const Matrix&);

  Matrix operator+(const Matrix&) const;
  Matrix operator-(const Matrix&) const;
  
  // Scalar multiplication
  // ---------------------
  Matrix& operator*=(double);
  Matrix& operator*=(const std::complex<double>&);

  friend Matrix operator*<>(double,const Matrix&);
  friend Matrix operator*<>(const std::complex<double>&,const Matrix&);
  
  Matrix operator*(double) const;
  Matrix operator*(const std::complex<double>&) const;

  // Addition with Scalar multiplication
  // -----------------------------------
  Matrix& multiAdd(double,const Matrix&);
  Matrix& multiAdd(const std::complex<double>&,const Matrix&);

  // Matrix multiplication
  // ---------------------
  template< int K1, int L1, int K2, int L2 >
    friend Matrix<K1,L2> operator*(const Matrix<K1,L1>& left,
				   const Matrix<K2,L2>& right);
  template< int L1 , int K2 >
    friend std::complex<double> operator*(const Matrix<1,L1>& left,
				     const Matrix<K2,1>& right);
  template < int K1, int L1 >
    friend Matrix<K1,L1>& operator*=(Matrix<K1,L1>& left,
				     const Matrix<L1,L1>& right);

  // Matrix manipulations
  // --------------------
  Matrix<L,K> T() const;                       // Transpose
  Matrix<L,K> H() const;                       // Hermitian conjugate
  friend Matrix<2,2> sqrt(const Matrix<2,2>&); // sqrt of 2x2 matrix
  template< int K1 >
    friend std::complex<double> Trace(const Matrix<K1,K1>&); // Trace

  // Output
  // ------
  void print(int, int) const;       // Print the (x,y) component
  void print() const;               // Print all components
  void printReal() const;           // Print the real part of all components 
  friend std::ostream& operator<< <>(std::ostream&, const Matrix&); // Overload <<
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//---- class Matrix<K,L> definitions
//---- We do not have a separate source file
//---- See C++ FAQ lite 35.12 for more info

// Constructors
// ------------

// Default constructor
template< int K, int L >
Matrix<K,L>::Matrix()
  : m(new std::complex<double>[K*L])
{
  static std::complex<double> zero(0,0);

  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      m[i*L+j] = zero;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Constructor for 1x1 matrix
/*
template< int K, int L >
Matrix<K,L>::Matrix(std::complex<double> number)
{
  if(K==1 && L==1) { // check dimensions of Matrix
      
    m = new std::complex<double>[1*1];
    m[0*1+0] = number; // Assign value
  }
  else {
    cerr << "Error in Matrix::Matrix(std::complex,std::complex,std::complex,std::complex): a ("
	 << K << "x" << L << ") matrix can not be initialized "
	 << "with 1 argument!" << std::endl;
    exit(1);
  }  
}*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Constructor for 2x2 matrix
template< int K, int L >
Matrix<K,L>::Matrix(std::complex<double> v11,std::complex<double> v12,
		    std::complex<double> v21,std::complex<double> v22)
{
  if( K*L == 4 ) { // check dimensions of Matrix
      
    m = new std::complex<double>[2*2];
    // Assign values
    m[0*2+0] = v11;
    m[0*2+1] = v12;
    m[1*2+0] = v21;
    m[1*2+1] = v22;
  }
  else {
    std::cerr << "Error in Matrix::Matrix(complex,complex,complex,complex): a ("
	 << K << "x" << L << ") matrix can not be initialized "
	 << "with 4 arguments!" << std::endl;
    exit(1);
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Constructor for 4x4 matrix
template< int K, int L >
Matrix<K,L>::Matrix(std::complex<double> v11,std::complex<double> v12,
		    std::complex<double> v13,std::complex<double> v14,
		    std::complex<double> v21,std::complex<double> v22,
		    std::complex<double> v23,std::complex<double> v24,
		    std::complex<double> v31,std::complex<double> v32,
		    std::complex<double> v33,std::complex<double> v34,
		    std::complex<double> v41,std::complex<double> v42,
		    std::complex<double> v43,std::complex<double> v44)
{
  if( K*L == 16 ) { // check dimensions of Matrix
    
    m = new std::complex<double>[4*4];
    // Assign values
    m[0*4+0] = v11;
    m[0*4+1] = v12;
    m[0*4+2] = v13;
    m[0*4+3] = v14;
    m[1*4+0] = v21;
    m[1*4+1] = v22;
    m[1*4+2] = v23;
    m[1*4+3] = v24;
    m[2*4+0] = v31;
    m[2*4+1] = v32;
    m[2*4+2] = v33;
    m[2*4+3] = v34;
    m[3*4+0] = v41;
    m[3*4+1] = v42;
    m[3*4+2] = v43;
    m[3*4+3] = v44;
  }
  else {
    std::cerr << "Error in Matrix::Matrix(complex(16times)): a ("
	 << K << "x" << L << ") matrix can not be initialized "
	 << "with 16 arguments!" << std::endl;
    exit(1);
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Copy constructor
template< int K, int L >
Matrix<K,L>::Matrix(const Matrix<K,L>& toCopy)
  : m(new std::complex<double>[K*L])
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      m[i*L+j] = toCopy.m[i*L+j];
}

//-------------------------------------------------------------------------

// Destructor
// ----------

template< int K, int L >
Matrix<K,L>::~Matrix()
{
  delete[] m;
}

//-------------------------------------------------------------------------

// Assignment
// ----------

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator=(const Matrix<K,L>& toCopy)
{
  if(this != &toCopy) // avoid self-assignment
    {
      for(int i=0; i<K; i++)
	for(int j=0; j<L; j++)
	  m[i*L+j] = toCopy.m[i*L+j];
    }

  return *this;
}

//-------------------------------------------------------------------------

// Access individual components
// ----------------------------

template< int K, int L >
std::complex<double>& Matrix<K,L>::operator()(int k, int l)
{
  if( !(k>=0 && k<K && l>=0 && l<L) ) {
    std::cerr << "Error in Matrix::operator()(int,int): "
	 << "Out of bounds!" << std::endl;
    exit(1);
  }
  return m[k*L+l];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
const std::complex<double> Matrix<K,L>::operator()(int k, int l) const
{
  if( !(k>=0 && k<K && l>=0 && l<L) ) {
    std::cerr << "Error in Matrix::operator()(int,int): "
	 << "Out of bounds!" << std::endl;
    exit(1);
  }
  return m[k*L+l];
}

//-------------------------------------------------------------------------

// Addition
// --------

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator+=(const Matrix<K,L>& that)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] += that.m[i*L+j];
     
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator-=(const Matrix<K,L>& that)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] -= that.m[i*L+j];
     
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L> Matrix<K,L>::operator+(const Matrix<K,L>& that) const
{
  Matrix<K,L> sum = *this;
  sum += that;
  return sum;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L> Matrix<K,L>::operator-(const Matrix<K,L>& that) const
{
  Matrix<K,L> sum = *this;
  sum -= that;
  return sum;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::multiAdd(double factor, const Matrix<K,L>& that)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] += factor * that.m[i*L+j];
  
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::multiAdd(const std::complex<double>& factor, 
				   const Matrix<K,L>& that)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] += factor * that.m[i*L+j];
  
  return *this;
}

//-------------------------------------------------------------------------

// Scalar multiplication 
// ---------------------

// (matrix*double)
template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator*=(double factor)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] *= factor;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (matrix*complex)
template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator*=(const std::complex<double>& factor)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] *= factor;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (matrix*double)
template< int K, int L >
Matrix<K,L> Matrix<K,L>::operator*(double factor) const
{
  Matrix<K,L> product = *this;
  product *= factor;
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (matrix*complex)
template< int K, int L >
Matrix<K,L> Matrix<K,L>::operator*(const std::complex<double>& factor) const
{
  Matrix<K,L> product = *this;
  product *= factor;
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (double*matrix)
template< int K, int L >
Matrix<K,L> operator*(double factor, const Matrix<K,L>& matrix)
{
  Matrix<K,L> product = matrix;
  product *= factor;
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (complex*matrix)
template< int K, int L >
Matrix<K,L> operator*(const std::complex<double>& factor,
		      const Matrix<K,L>& matrix)
{
  Matrix<K,L> product = matrix;
  product *= factor;
  
  return product;
}

//-------------------------------------------------------------------------

// Matrix multiplication
// ---------------------

template< int K1, int L1, int K2, int L2 >
Matrix<K1,L2> operator*(const Matrix<K1,L1>& left, const Matrix<K2,L2>& right)
{
  if(L1 != K2) { // Check whether multiplication is meaningful
    std::cerr << "Error in operator*(Matrix,Matrix): " 
	 << "Matrix sizes are not compatible for multiplication" << std::endl;
    exit(1);
  }
  
  Matrix<K1,L2> product;
  for(int i=0; i<K1; i++)
    for(int j=0; j<L2; j++)
      for(int x=0; x<K2; x++)
	product.m[i*L2+j] += (left.m[i*L1+x]*right.m[x*L2+j]);
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Matrix multiplication resulting in a scalar
template< int L1 , int K2 >
std::complex<double> operator*(const Matrix<1,L1>& left, const Matrix<K2,1>& right)
{
  if(L1 != K2) { // Check whether multiplication is meaningful
    std::cerr << "Error in operator*(Matrix,Matrix): " 
	      << "Matrix sizes are not compatible for multiplication" << std::endl;
    exit(1);
  }

  std::complex<double> scalar(0,0);
  for(int i=0; i<L1; i++)
    scalar += left.m[0*L1+i]*right.m[i*1+0];
  return scalar;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template < int K1, int L1 >
Matrix<K1,L1>& operator*=(Matrix<K1,L1>& left,const Matrix<L1,L1>& right)
{
  left = left * right;
  
  return left;
}

//-------------------------------------------------------------------------

// Transpose
// ---------
template< int K, int L >
Matrix<L,K> Matrix<K,L>::T() const
{
  Matrix<L,K> transposed;
  for(int i=0; i<L; i++)
    for(int j=0; j<K; j++)
      transposed(i,j) = m[j*L+i];

  return transposed;
}

//-------------------------------------------------------------------------

// Hermitian conjugate
// -------------------
template< int K, int L >
Matrix<L,K> Matrix<K,L>::H() const
{
  Matrix<L,K> conjugate;
  for(int i=0; i<L; i++)
    for(int j=0; j<K; j++)
      conjugate(i,j) = conj(m[j*L+i]);

  return conjugate;
}

//-------------------------------------------------------------------------

// Trace of square matrices
// ------------------------
template< int K>
std::complex<double> Trace(const Matrix<K,K>& m)
{
  std::complex<double> trace = 0.;
  
  for(int i=0; i<K; ++i) trace += m.m[i*K+i];

  return trace;
}

//-------------------------------------------------------------------------

// Print the (x,y) component
template< int K, int L >
void Matrix<K,L>::print(int x, int y) const
{
  std::cout << (*this)(x,y) << std::endl;
}

//-------------------------------------------------------------------------

// Print all components
template< int K, int L >
void Matrix<K,L>::print() const
{
  std::cout << "{";
  for(int i=0; i<K; i++)
    {
      std::cout << "{";
      for(int j=0; j<L; j++)
	{
	  std::cout << m[i*L+j]
	       << ( j<L-1 ? "," : "" );
	}
      std::cout << "}";
    }
  std::cout << "}" << std::endl;
}

//-------------------------------------------------------------------------

// Print the real part of all components
template< int K, int L >
void Matrix<K,L>::printReal() const
{
  std::cout << "{";
  for(int i=0; i<K; i++)
    {
      std::cout << "{";
      for(int j=0; j<L; j++)
	{
	  std::cout << real(m[i*L+j])
	       << ( j<L-1 ? "," : "" );
	}
      std::cout << "}";
    }
  std::cout << "}" << std::endl;
}

//-------------------------------------------------------------------------

// Enables output using << operator
template< int K, int L >
std::ostream& operator<<(std::ostream& output, const Matrix<K,L>& matrix)
{
  output << "{";
  for(int i=0; i<K; i++)
    {
      output << "{";
      for(int j=0; j<L; j++)
	{
	  output << matrix.m[i*L+j]
	       << ( j<L-1 ? "," : "" );
	}
      output << "}";
    }
  output << "}";

  return output;
}



#endif
