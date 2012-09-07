/**
 * @file quadrature.hpp
 * @brief header with quadrature rules
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt) (2010-)
 * @author k.vantournhout@gsi.de
 */

#ifndef NUMINT_QUADRATURE_HPP
#define NUMINT_QUADRATURE_HPP 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <numint/macros.hpp>
#include <numint/typedef.hpp>


OPEN_NUMINT_NAMESPACE

/* =========================================================================
   Quadrature rules
   ========================================================================= */
/**
   @name Trapezoidal quadrature
*/
/**@{
   Quadrature rules for functions returning double precision, complex
   double precision or vectors. P indicates that integration assumes a
   periodic function.

   @param f the functor to integrate
   @param a the lower limit
   @param b the upper limit
   @param epsabs the absolute error
   @param epsrel the relative error
   @param result the integral value
   @param neval  value incremented with total function evaluations
   @param dirty check convergence on the first dirty values
   @param compare use another function to determain convergence
*/
int quad_tr(const function<double> &f, double a, double b,
            double epsabs, double epsrel, double &result, unsigned &neval);
int quad_tr(const function<dcomplex> &f, double a, double b,
            double epsabs, double epsrel, dcomplex &result, unsigned &neval);
int quad_tr(const function<vector_d> &f, double a, double b,
            double epsabs, double epsrel, vector_d &result, unsigned &neval, int dirty,
            typename Convergence<vector_d>::type compare = NULL);
int quad_tr(const function<vector_z> &f, double a, double b,
            double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
            typename Convergence<vector_z>::type compare = NULL);

int quad_ptr(const function<double> &f, double a, double b,
             double epsabs, double epsrel, double &result, unsigned &neval);
int quad_ptr(const function<dcomplex> &f, double a, double b,
             double epsabs, double epsrel, dcomplex &result, unsigned &neval);
int quad_ptr(const function<vector_d> &f, double a, double b,
             double epsabs, double epsrel, vector_d &result, unsigned &neval, int dirty,
             typename Convergence<vector_d>::type compare = NULL);
int quad_ptr(const function<vector_z> &f, double a, double b,
             double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
             typename Convergence<vector_z>::type compare = NULL);
/**@}*/

/**
   @name Romberg quadrature
*/
/**@{
   Romberg quadrature rules for functions returning double precision, complex
   double precision or vectors. P indicates that integration assumes a
   periodic function.

   @param f the functor to integrate
   @param a the lower limit
   @param b the upper limit
   @param epsabs the absolute error
   @param epsrel the relative error
   @param result the integral value
   @param neval  value incremented with total function evaluations
   @param dirty check convergence on the first dirty values
   @param compare use another function to determain convergence
*/
int quad_romb(const function<double> &f, double a, double b,
              double epsabs, double epsrel, double &result, unsigned &neval);
int quad_romb(const function<dcomplex> &f, double a, double b,
              double epsabs, double epsrel, dcomplex &result, unsigned &neval);
int quad_romb(const function<vector_d> &f, double a, double b,
              double epsabs, double epsrel, vector_d &result, unsigned &neval, int dirty,
              typename Convergence<vector_d>::type compare = NULL);
int quad_romb(const function<vector_z> &f, double a, double b,
              double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
              typename Convergence<vector_z>::type compare = NULL);

int quad_promb(const function<double> &f, double a, double b,
               double epsabs, double epsrel, double &result, unsigned &neval);
int quad_promb(const function<dcomplex> &f, double a, double b,
               double epsabs, double epsrel, dcomplex &result, unsigned &neval);
int quad_promb(const function<vector_d> &f, double a, double b,
               double epsabs, double epsrel, vector_d &result, unsigned &neval, int dirty,
               typename Convergence<vector_d>::type compare = NULL);
int quad_promb(const function<vector_z> &f, double a, double b,
               double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
               typename Convergence<vector_z>::type compare = NULL);
/**@}*/

/**
   @name Clenshaw Curtis quadrature
*/
/**@{
   Clenshaw Curtis quadrature rules for functions returning double precision, complex
   double precision or vectors. P indicates that integration assumes a
   periodic function.

   @param f the functor to integrate
   @param a the lower limit
   @param b the upper limit
   @param epsabs the absolute error
   @param epsrel the relative error
   @param result the integral value
   @param neval  value incremented with total function evaluations
   @param dirty check convergence on the first dirty values
   @param compare use another function to determain convergence
*/
int quad_cc(const function<double> &f, double a, double b,
            double epsabs, double epsrel, double &result, unsigned &neval);
int quad_cc(const function<dcomplex> &f, double a, double b,
            double epsabs, double epsrel, dcomplex &result, unsigned &neval);
int quad_cc(const function<vector_d> &f, double a, double b,
            double epsabs, double epsrel, vector_d &result, unsigned &neval, int dirty,
            typename Convergence<vector_d>::type compare = NULL);
int quad_cc(const function<vector_z> &f, double a, double b,
            double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
            typename Convergence<vector_z>::type compare = NULL);

int quad_pcc(const function<double> &f, double a, double b,
             double epsabs, double epsrel, double &result, unsigned &neval);
int quad_pcc(const function<dcomplex> &f, double a, double b,
             double epsabs, double epsrel, dcomplex &result, unsigned &neval);
int quad_pcc(const function<vector_d> &f, double a, double b,
             double epsabs, double epsrel, vector_d &result, unsigned &neval, int dirty,
             typename Convergence<vector_d>::type compare = NULL);
int quad_pcc(const function<vector_z> &f, double a, double b,
             double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
             typename Convergence<vector_z>::type compare = NULL);

/**@}*/

CLOSE_NUMINT_NAMESPACE

#endif
