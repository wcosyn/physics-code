/* 
 * Matrix.cpp
 *
 * All declarations for class Matrix<K,L> are done inside
 * the header file (Matrix.h)
 *
 * Friend function that takes square root of 2x2 matrix is
 * declared here.
 *
 * author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * 
 */


#include "Matrix.h"
#include <iostream>
using std::cout; using std::cerr; using std::endl;
#include <complex>
using std::complex;
#include <cstdlib>
using std::exit;
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

/* Square root of a hermitian 2x2-matrix
 * This means taking square root of
 * eigenvalues and rotating this matrix
 * back
 */
Matrix<2,2> sqrt(const Matrix<2,2>& inputMatrix)
{
  // Test whether the matrix is hermitian
  if( inputMatrix.m[0*2+1] != conj(inputMatrix.m[1*2+0]) )
    {
      cerr << "ERROR in sqrt(const Matrix<2,2>&): "
	   << "Only sqrt of hermitian 2x2 matrix implemented." << endl;
      exit(1);
    }

  // Store the matrix in gsl format
  gsl_matrix_complex *matrix =  gsl_matrix_complex_alloc(2,2);
  for(int i=0; i<2; ++i) 
    for(int j=0; j<2; ++j)
      gsl_matrix_complex_set(matrix, i, j, 
			     gsl_complex_rect(inputMatrix.m[i*2+j].real(),
					      inputMatrix.m[i*2+j].imag()));

  //allocate workspace
  gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(2);
  //allocate eigenvalue vector
  gsl_vector *eval = gsl_vector_alloc(2);
  //allocate eigenvector matrix - mutually orthogonal, normalized to 1
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(2, 2);
  if ((gsl_eigen_hermv (matrix, eval, evec, w)) != 0) {
    //could not find eigenvalues
    cerr << "ERROR in sqrt(const Matrix<2,2>&): "
	 << "Failed at finding eigenvalues and eigenvectors." << endl;
    exit(1);
  }

  // sort the eigenvalues, take square root and put them into
  // a matrix with eigenvalues across the diagonal
  gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  gsl_matrix_complex *eval_matrix = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex_set_zero(eval_matrix);
  for (int i=0; i<2; i++) {
    gsl_matrix_complex_set(eval_matrix, i, i, 
			   gsl_complex_sqrt_real(gsl_vector_get(eval, i)));
  }
  
  // Rotate the solution backwards
  //perform V*D^1/2*inv(V)
  static gsl_complex zero = gsl_complex_rect(0.,0.);
  static gsl_complex one = gsl_complex_rect(1.,0.);
  gsl_matrix_complex *tmp = gsl_matrix_complex_alloc(2,2);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, evec, eval_matrix, zero, tmp);
  gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, one, tmp, evec, zero, matrix);

  // Store the solution in our format
  Matrix<2,2> root;
  for(int i=0; i<2; ++i) 
    for(int j=0; j<2; ++j)
      root.m[i*2+j] 
	= complex<double>(GSL_REAL(gsl_matrix_complex_get(matrix,i,j)),
			  GSL_IMAG(gsl_matrix_complex_get(matrix,i,j))); 
  
  //deallocate memory
  gsl_matrix_complex_free(eval_matrix);
  gsl_eigen_hermv_free(w);
  gsl_vector_free(eval);
  gsl_matrix_complex_free(evec);
  gsl_matrix_complex_free(tmp);
  gsl_matrix_complex_free(matrix);

  return root;
}
