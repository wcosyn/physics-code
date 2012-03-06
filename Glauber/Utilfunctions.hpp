/*! \file Utilfunctions.hpp
 * \brief Some useful utilities functions, mostly math related
 * \author Wim Cosyn
 * \date 16/08/2011
 */

#ifndef UTILFUNCTIONS_H
#define UTILFUNCTIONS_H

#include <cstdarg>
#include <sstream>
#include <complex>

/*! interpolation of an array
 * \param array array to interpolate
 * \param r interpolation value
 * \param rstep gridsize of the array
 * \param lines size of the array
 * \param offset units of rstep where the grid starts
 * \return interpolation value of the grid at r
 */
double interpolate(const double *array, double r, double rstep,int lines, int offset);
double power(double x,int n);/*!< returns x^n */
int power(int x,int n);/*!< returns x^n */
/*! returns the distance between points r1 and r2
 * \param r1 r value of point1
 * \param costheta1 cos(theta) value of point1
 * \param sintheta1 sin(theta) value of point1
 * \param cosphi1 cos(phi) value of point1
 * \param sinphi1 sin(phi) value of point1
 * \param r2 r value of point2
 * \param costheta2 cos(theta) value of point2
 * \param sintheta2 sin(theta) value of point2
 * \param cosphi2 cos(phi) value of point2
 * \param sinphi2 sin(phi) value of point2
 */
double normr(double r1, double costheta1, double sintheta1, double cosphi1, double sinphi1,
	     double r2, double costheta2, double sintheta2, double cosphi2, double sinphi2);
/*! solves matrix equation a*x=b with LU solve algorithm and store solution in b, from numerical recipes p46
 * \param a square matrix a
 * \param b column matrix b
 * \param order number of equations
 */
void LUsolve(double **a, double *b, int order);
/*! calculates sin and cos of value x in one go
 * \param x value to compute cos and sin of
 * \param sin pointer to sin value
 * \param cos pointer to cos value
 */
void sincos(const double x, double * sin, double * cos);
/*! performs a fit to data with a set of fitfunctions, also computes the chi squared of the fit
 * \param xarray data x array
 * \param yarray data y array
 * \param sigma data y errors
 * \param lower lower index of the data array for the fit
 * \param upper upper index of the data array for the fit
 * \param order order of the fit (in fitfunctions)
 * \param chi2 chi squard of the fit
 */
double* fitdata(double *xarray, double *yarray, double *sigma, int lower, int upper, int order, double & chi2);
/*!set of base functions used for the fit
 * \param x x value
 * \param i order index of the base function set
 * \return value of f_i(x)
 */
double fitfunction(double x, int i);
double **matrix(int n, int m); /*!< create a n*m matrix on the heap*/
void freematrix(double **matrix, int n); /*!< free a n*? matrix on the heap*/

/*! template function for romberg integration of N functions (slightly different preferably) T f(x)
 * \tparam T template return type of the function you're integrating
 * \param function function you want to integrate, must return N values (through array)
 * \param a lower limit
 * \param b upper limit
 * \param N number of functions you want to integrate
 * \param results array with integration results
 * \param acc relative accuracy that has to be reached
 * \param min min number of steps in the Romberg algorithm
 * \param max max number of steps in the Romberg algorithm
 * \param estimate pointer to a variable containing an upper limit on nested integration, to speed up things
 */
template <class T> void rombergerN(void (*function)(double, T*,  va_list), double a, double b, 
		       int N, T* results, double acc, int min, int max, double *estimate, ...)
{

  if(abs(a-b)<1e-06){for(int i=0;i<N;i++) results[i]=0.;return;}
  T **DN1= new T*[N];
  T **DN2= new T*[N];
  double x,h=b-a;
  T sum[N], value[N];
  double dev[N];
  
  va_list ap;
  va_start(ap,estimate);
  
  for(int i=0;i<N;i++) {
    DN1[i] = new T[1];
  }

  function(a,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++) DN1[i][0] = value[i]; 
  va_start(ap,estimate);
  function(b,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++){
    DN1[i][0] += value[i]; 
  DN1[i][0] *= 0.5*(b-a);
  }
  

  
  for (int n = 1; n <= max ; n++) {
    for(int i=0;i<N;i++) {
      sum[i] = 0.;
      //store new results
      DN2[i] = new T[n+1];
    }
    int ceiling = power(2,n);
    // Trapezium rule recursive rule
    for (int k = 1; k < ceiling; k += 2) {
      x = a + power(0.5,n)*k*h;
      va_start(ap,estimate);
      function(x,value,ap);
      va_end(ap);
      for(int i=0;i<N;i++) sum[i] += value[i];
    }
    for(int i=0;i<N;i++) DN2[i][0] = 0.5*DN1[i][0] + power(0.5,n)*h*sum[i];
    int p=4;
    for (int m = 1; m <= n; m++) {
      for(int i=0;i<N;i++) DN2[i][m] = (double(p)*DN2[i][m-1] - DN1[i][m-1])/(p - 1.);
      p*=4;
    }
    if (n >= min) {
      double deviation=0.;
      for(int i=0;i<N;i++) {
	dev[i] = (DN2[i][n] == complex<double>(0.,0.)) ? 0. : abs((DN2[i][n]-DN1[i][n-1]))/abs(DN2[i][n]);
	if(dev[i]>deviation) deviation=dev[i];
      }
     if (((deviation < acc ) ) || ((abs(DN2[0][n]-DN1[0][n-1]) <acc*1e-07 ))) {
       for(int i=0;i<N;i++) {
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];
	 if(abs(results[i])>(*estimate)) (*estimate)=abs(results[i]);
       }
       delete [] DN1; delete [] DN2;
       return;
     }
     double dummy = 0.;
     for(int i=0;i<N;i++) dummy +=abs(DN2[i][n]);
     if(dummy <(*estimate)*1e-05){
       for(int i=0;i<N;i++){
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];	
       }
       delete [] DN1; delete [] DN2;
       return;
     }
    }
    for(int i=0;i<N;i++){
      delete [] DN1[i];
      DN1[i] = DN2[i];
    }
  }
    
  for(int i=0;i<N;i++){
    results[i] = DN2[i][max];
    delete [] DN2[i];
  }
  delete [] DN1; delete [] DN2;
  va_end(ap);
  return;  
  
}


/*! template function for romberg integration of N member functions (slightly different preferably) T F->f(x)
 * \tparam T template return type of the function you're integrating
 * \tparam F template class of which the integration functions are member function
 * \param object pointer to class object that contains the memberfunction
 * \param function member function you want to integrate, must return N values (through array)
 * \param a lower limit
 * \param b upper limit
 * \param N number of functions you want to integrate
 * \param results array with integration results
 * \param acc relative accuracy that has to be reached
 * \param min min number of steps in the Romberg algorithm
 * \param max max number of steps in the Romberg algorithm
 * \param estimate pointer to a variable containing an upper limit on nested integration, to speed up things
 */
template <class T, class F> void rombergerN(F* object, void (F::*function)(double, T*,  va_list), double a, double b, 
		       int N, T* results, double acc, int min, int max, double *estimate, ...)
{

  if(abs(a-b)<1e-06){for(int i=0;i<N;i++) results[i]=0.;return;}
  T **DN1= new T*[N];
  T **DN2= new T*[N];
  double x,h=b-a;
  T sum[N], value[N];
  double dev[N];
  
  va_list ap;
  va_start(ap,estimate);
  
  for(int i=0;i<N;i++) {
    DN1[i] = new T[1];
  }

  (object->*function)(a,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++) DN1[i][0] = value[i]; 
  va_start(ap,estimate);
  (object->*function)(b,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++){
    DN1[i][0] += value[i]; 
  DN1[i][0] *= 0.5*(b-a);
  }
  

  
  for (int n = 1; n <= max ; n++) {
    for(int i=0;i<N;i++) {
      sum[i] = 0.;
      //store new results
      DN2[i] = new T[n+1];
    }
    int ceiling = power(2,n);
    // Trapezium rule recursive rule
    for (int k = 1; k < ceiling; k += 2) {
      x = a + power(0.5,n)*k*h;
      va_start(ap,estimate);
      (object->*function)(x,value,ap);
      va_end(ap);
      for(int i=0;i<N;i++) sum[i] += value[i];
    }
    for(int i=0;i<N;i++) DN2[i][0] = 0.5*DN1[i][0] + power(0.5,n)*h*sum[i];
    int p=4;
    for (int m = 1; m <= n; m++) {
      for(int i=0;i<N;i++) DN2[i][m] = (double(p)*DN2[i][m-1] - DN1[i][m-1])/(p - 1.);
      p*=4;
    }
    if (n >= min) {
      double deviation=0.;
      for(int i=0;i<N;i++) {
	dev[i] = (DN2[i][n] == complex<double>(0.,0.)) ? 0. : abs((DN2[i][n]-DN1[i][n-1]))/abs(DN2[i][n]);
	if(dev[i]>deviation) deviation=dev[i];
      }
     if (((deviation < acc ) ) || ((abs(DN2[0][n]-DN1[0][n-1]) <acc*1e-07 ))) {
       for(int i=0;i<N;i++) {
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];
	 if(abs(results[i])>(*estimate)) (*estimate)=abs(results[i]);
       }
       delete [] DN1; delete [] DN2;
       return;
     }
     double dummy = 0.;
     for(int i=0;i<N;i++) dummy +=abs(DN2[i][n]);
     if(dummy <(*estimate)*1e-05){
       for(int i=0;i<N;i++){
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];	
       }
       delete [] DN1; delete [] DN2;
       return;
     }
    }
    for(int i=0;i<N;i++){
      delete [] DN1[i];
      DN1[i] = DN2[i];
    }
  }
    
  for(int i=0;i<N;i++){
    results[i] = DN2[i][max];
    delete [] DN2[i];
  }
  delete [] DN1; delete [] DN2;
  va_end(ap);
  return;  
  
}
/*!convert a certain type to a string
 * \tparam T type to convert to string
 * \return returns a string
 */
template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}
#endif

