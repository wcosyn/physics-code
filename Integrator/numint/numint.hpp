/**
 * @file lib/numint/numint.hpp
 * @brief various numerical integrators
 *
 * This file contains various integrator rules, one-dimensional and
 * multi-dimensional for various types of functions.
 *
 * changelog :
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt) (2010-)
 * @author k.vantournhout@gsi.de
 *
 * @date last update : Thu Aug 16 08:57:18 CEST 2012
 */

/**
   @namespace numint
   blablabla
*/

#ifndef NUMINT_HPP
#define NUMINT_HPP 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <numint/macros.hpp>
#include <numint/typedef.hpp>
#include <numint/numint2adapt.hpp>
#include <numint/numint2Cuba.hpp>

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
             Convergence<vector_d>::type compare = NULL);
int quad_tr(const function<vector_z> &f, double a, double b,
            double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
             Convergence<vector_z>::type compare = NULL);

int quad_ptr(const function<double> &f, double a, double b,
             double epsabs, double epsrel, double &result, unsigned &neval);
int quad_ptr(const function<dcomplex> &f, double a, double b,
             double epsabs, double epsrel, dcomplex &result, unsigned &neval);
int quad_ptr(const function<vector_d> &f, double a, double b,
             double epsabs, double epsrel, vector_d &result, unsigned &neval, int dirty,
              Convergence<vector_d>::type compare = NULL);
int quad_ptr(const function<vector_z> &f, double a, double b,
             double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
              Convergence<vector_z>::type compare = NULL);
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
               Convergence<vector_d>::type compare = NULL);
int quad_romb(const function<vector_z> &f, double a, double b,
              double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);

int quad_promb(const function<double> &f, double a, double b,
               double epsabs, double epsrel, double &result, unsigned &neval);
int quad_promb(const function<dcomplex> &f, double a, double b,
               double epsabs, double epsrel, dcomplex &result, unsigned &neval);
int quad_promb(const function<vector_d> &f, double a, double b,
               double epsabs, double epsrel, vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int quad_promb(const function<vector_z> &f, double a, double b,
               double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
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
             Convergence<vector_d>::type compare = NULL);
int quad_cc(const function<vector_z> &f, double a, double b,
            double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
             Convergence<vector_z>::type compare = NULL);

int quad_pcc(const function<double> &f, double a, double b,
             double epsabs, double epsrel, double &result, unsigned &neval);
int quad_pcc(const function<dcomplex> &f, double a, double b,
             double epsabs, double epsrel, dcomplex &result, unsigned &neval);
int quad_pcc(const function<vector_d> &f, double a, double b,
             double epsabs, double epsrel, vector_d &result, unsigned &neval, int dirty,
              Convergence<vector_d>::type compare = NULL);
int quad_pcc(const function<vector_z> &f, double a, double b,
             double epsabs, double epsrel, vector_z &result, unsigned &neval, int dirty,
              Convergence<vector_z>::type compare = NULL);

/**@}*/

/* =========================================================================
   Fubini integration
   ========================================================================= */
/**
   @name Fubini integration using Trapezoidal quadrature
*/
/**@{
   Quadrature rules for 1,2 and 3 dimensional functions returning double
   precision, complex double precision or vectors. P indicates that
   integration assumes a periodic function.

   @param f the functor to integrate
   @param a the lower limits
   @param b the upper limits
   @param epsabs the absolute error
   @param epsrel the relative error
   @param result the integral value
   @param neval  value incremented with total function evaluations
   @param dirty check convergence on the first dirty values
   @param compare use another function to determain convergence
*/

int fubini_tr(const mdfunction<double,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int fubini_tr(const mdfunction<double,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int fubini_tr(const mdfunction<double,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);

int fubini_tr(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int fubini_tr(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int fubini_tr(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);

int fubini_tr(const mdfunction<vector_d,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int fubini_tr(const mdfunction<vector_d,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int fubini_tr(const mdfunction<vector_d,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);

int fubini_tr(const mdfunction<vector_z,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int fubini_tr(const mdfunction<vector_z,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int fubini_tr(const mdfunction<vector_z,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);

int fubini_ptr(const mdfunction<double,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int fubini_ptr(const mdfunction<double,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int fubini_ptr(const mdfunction<double,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);

int fubini_ptr(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int fubini_ptr(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int fubini_ptr(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);

int fubini_ptr(const mdfunction<vector_d,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int fubini_ptr(const mdfunction<vector_d,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int fubini_ptr(const mdfunction<vector_d,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);

int fubini_ptr(const mdfunction<vector_z,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int fubini_ptr(const mdfunction<vector_z,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int fubini_ptr(const mdfunction<vector_z,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);

int fubini_tr(const mdfunction<double,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int fubini_tr(const mdfunction<double,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int fubini_tr(const mdfunction<double,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);

int fubini_tr(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int fubini_tr(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int fubini_tr(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);

int fubini_tr(const mdfunction<vector_d,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int fubini_tr(const mdfunction<vector_d,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int fubini_tr(const mdfunction<vector_d,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);

int fubini_tr(const mdfunction<vector_z,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int fubini_tr(const mdfunction<vector_z,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int fubini_tr(const mdfunction<vector_z,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);

int fubini_ptr(const mdfunction<double,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int fubini_ptr(const mdfunction<double,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int fubini_ptr(const mdfunction<double,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);

int fubini_ptr(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int fubini_ptr(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int fubini_ptr(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);

int fubini_ptr(const mdfunction<vector_d,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int fubini_ptr(const mdfunction<vector_d,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int fubini_ptr(const mdfunction<vector_d,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);

int fubini_ptr(const mdfunction<vector_z,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int fubini_ptr(const mdfunction<vector_z,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int fubini_ptr(const mdfunction<vector_z,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);

/**@}*/

/**
   @name Fubini integration using Romberg quadrature
*/
/**@{
   Quadrature rules for 1,2 and 3 dimensional functions returning double
   precision, complex double precision or vectors. P indicates that
   integration assumes a periodic function.

   @param f the functor to integrate
   @param a the lower limits
   @param b the upper limits
   @param epsabs the absolute error
   @param epsrel the relative error
   @param result the integral value
   @param neval  value incremented with total function evaluations
   @param dirty check convergence on the first dirty values
   @param compare use another function to determain convergence
*/

int fubini_romb(const mdfunction<double,1> &f, const array<double,1> &a,
                const array<double,1> &b, double epsabs, double epsrel,
                double &result, unsigned &neval);
int fubini_romb(const mdfunction<double,2> &f, const array<double,2> &a,
                const array<double,2> &b, double epsabs, double epsrel,
                double &result, unsigned &neval);
int fubini_romb(const mdfunction<double,3> &f, const array<double,3> &a,
                const array<double,3> &b, double epsabs, double epsrel,
                double &result, unsigned &neval);

int fubini_romb(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
                const array<double,1> &b, double epsabs, double epsrel,
                dcomplex &result, unsigned &neval);
int fubini_romb(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
                const array<double,2> &b, double epsabs, double epsrel,
                dcomplex &result, unsigned &neval);
int fubini_romb(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
                const array<double,3> &b, double epsabs, double epsrel,
                dcomplex &result, unsigned &neval);

int fubini_romb(const mdfunction<vector_d,1> &f, const array<double,1> &a,
                const array<double,1> &b, double epsabs, double epsrel,
                vector_d &result, unsigned &neval, int dirty,
                 Convergence<vector_d>::type compare = NULL);
int fubini_romb(const mdfunction<vector_d,2> &f, const array<double,2> &a,
                const array<double,2> &b, double epsabs, double epsrel,
                vector_d &result, unsigned &neval, int dirty,
                 Convergence<vector_d>::type compare = NULL);
int fubini_romb(const mdfunction<vector_d,3> &f, const array<double,3> &a,
                const array<double,3> &b, double epsabs, double epsrel,
                vector_d &result, unsigned &neval, int dirty,
                 Convergence<vector_d>::type compare = NULL);

int fubini_romb(const mdfunction<vector_z,1> &f, const array<double,1> &a,
                const array<double,1> &b, double epsabs, double epsrel,
                vector_z &result, unsigned &neval, int dirty,
                 Convergence<vector_z>::type compare = NULL);
int fubini_romb(const mdfunction<vector_z,2> &f, const array<double,2> &a,
                const array<double,2> &b, double epsabs, double epsrel,
                vector_z &result, unsigned &neval, int dirty,
                 Convergence<vector_z>::type compare = NULL);
int fubini_romb(const mdfunction<vector_z,3> &f, const array<double,3> &a,
                const array<double,3> &b, double epsabs, double epsrel,
                vector_z &result, unsigned &neval, int dirty,
                 Convergence<vector_z>::type compare = NULL);

int fubini_promb(const mdfunction<double,1> &f, const array<double,1> &a,
                 const array<double,1> &b, double epsabs, double epsrel,
                 double &result, unsigned &neval);
int fubini_promb(const mdfunction<double,2> &f, const array<double,2> &a,
                 const array<double,2> &b, double epsabs, double epsrel,
                 double &result, unsigned &neval);
int fubini_promb(const mdfunction<double,3> &f, const array<double,3> &a,
                 const array<double,3> &b, double epsabs, double epsrel,
                 double &result, unsigned &neval);

int fubini_promb(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
                 const array<double,1> &b, double epsabs, double epsrel,
                 dcomplex &result, unsigned &neval);
int fubini_promb(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
                 const array<double,2> &b, double epsabs, double epsrel,
                 dcomplex &result, unsigned &neval);
int fubini_promb(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
                 const array<double,3> &b, double epsabs, double epsrel,
                 dcomplex &result, unsigned &neval);

int fubini_promb(const mdfunction<vector_d,1> &f, const array<double,1> &a,
                 const array<double,1> &b, double epsabs, double epsrel,
                 vector_d &result, unsigned &neval, int dirty,
                  Convergence<vector_d>::type compare = NULL);
int fubini_promb(const mdfunction<vector_d,2> &f, const array<double,2> &a,
                 const array<double,2> &b, double epsabs, double epsrel,
                 vector_d &result, unsigned &neval, int dirty,
                  Convergence<vector_d>::type compare = NULL);
int fubini_promb(const mdfunction<vector_d,3> &f, const array<double,3> &a,
                 const array<double,3> &b, double epsabs, double epsrel,
                 vector_d &result, unsigned &neval, int dirty,
                  Convergence<vector_d>::type compare = NULL);

int fubini_promb(const mdfunction<vector_z,1> &f, const array<double,1> &a,
                 const array<double,1> &b, double epsabs, double epsrel,
                 vector_z &result, unsigned &neval, int dirty,
                  Convergence<vector_z>::type compare = NULL);
int fubini_promb(const mdfunction<vector_z,2> &f, const array<double,2> &a,
                 const array<double,2> &b, double epsabs, double epsrel,
                 vector_z &result, unsigned &neval, int dirty,
                  Convergence<vector_z>::type compare = NULL);
int fubini_promb(const mdfunction<vector_z,3> &f, const array<double,3> &a,
                 const array<double,3> &b, double epsabs, double epsrel,
                 vector_z &result, unsigned &neval, int dirty,
                  Convergence<vector_z>::type compare = NULL);

int fubini_romb(const mdfunction<double,4> &f, const array<double,4> &a,
                const array<double,4> &b, double epsabs, double epsrel,
                double &result, unsigned &neval);
int fubini_romb(const mdfunction<double,5> &f, const array<double,5> &a,
                const array<double,5> &b, double epsabs, double epsrel,
                double &result, unsigned &neval);
int fubini_romb(const mdfunction<double,6> &f, const array<double,6> &a,
                const array<double,6> &b, double epsabs, double epsrel,
                double &result, unsigned &neval);

int fubini_romb(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
                const array<double,4> &b, double epsabs, double epsrel,
                dcomplex &result, unsigned &neval);
int fubini_romb(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
                const array<double,5> &b, double epsabs, double epsrel,
                dcomplex &result, unsigned &neval);
int fubini_romb(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
                const array<double,6> &b, double epsabs, double epsrel,
                dcomplex &result, unsigned &neval);

int fubini_romb(const mdfunction<vector_d,4> &f, const array<double,4> &a,
                const array<double,4> &b, double epsabs, double epsrel,
                vector_d &result, unsigned &neval, int dirty,
                 Convergence<vector_d>::type compare = NULL);
int fubini_romb(const mdfunction<vector_d,5> &f, const array<double,5> &a,
                const array<double,5> &b, double epsabs, double epsrel,
                vector_d &result, unsigned &neval, int dirty,
                 Convergence<vector_d>::type compare = NULL);
int fubini_romb(const mdfunction<vector_d,6> &f, const array<double,6> &a,
                const array<double,6> &b, double epsabs, double epsrel,
                vector_d &result, unsigned &neval, int dirty,
                 Convergence<vector_d>::type compare = NULL);

int fubini_romb(const mdfunction<vector_z,4> &f, const array<double,4> &a,
                const array<double,4> &b, double epsabs, double epsrel,
                vector_z &result, unsigned &neval, int dirty,
                 Convergence<vector_z>::type compare = NULL);
int fubini_romb(const mdfunction<vector_z,5> &f, const array<double,5> &a,
                const array<double,5> &b, double epsabs, double epsrel,
                vector_z &result, unsigned &neval, int dirty,
                 Convergence<vector_z>::type compare = NULL);
int fubini_romb(const mdfunction<vector_z,6> &f, const array<double,6> &a,
                const array<double,6> &b, double epsabs, double epsrel,
                vector_z &result, unsigned &neval, int dirty,
                 Convergence<vector_z>::type compare = NULL);

int fubini_promb(const mdfunction<double,4> &f, const array<double,4> &a,
                 const array<double,4> &b, double epsabs, double epsrel,
                 double &result, unsigned &neval);
int fubini_promb(const mdfunction<double,5> &f, const array<double,5> &a,
                 const array<double,5> &b, double epsabs, double epsrel,
                 double &result, unsigned &neval);
int fubini_promb(const mdfunction<double,6> &f, const array<double,6> &a,
                 const array<double,6> &b, double epsabs, double epsrel,
                 double &result, unsigned &neval);

int fubini_promb(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
                 const array<double,4> &b, double epsabs, double epsrel,
                 dcomplex &result, unsigned &neval);
int fubini_promb(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
                 const array<double,5> &b, double epsabs, double epsrel,
                 dcomplex &result, unsigned &neval);
int fubini_promb(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
                 const array<double,6> &b, double epsabs, double epsrel,
                 dcomplex &result, unsigned &neval);

int fubini_promb(const mdfunction<vector_d,4> &f, const array<double,4> &a,
                 const array<double,4> &b, double epsabs, double epsrel,
                 vector_d &result, unsigned &neval, int dirty,
                  Convergence<vector_d>::type compare = NULL);
int fubini_promb(const mdfunction<vector_d,5> &f, const array<double,5> &a,
                 const array<double,5> &b, double epsabs, double epsrel,
                 vector_d &result, unsigned &neval, int dirty,
                  Convergence<vector_d>::type compare = NULL);
int fubini_promb(const mdfunction<vector_d,6> &f, const array<double,6> &a,
                 const array<double,6> &b, double epsabs, double epsrel,
                 vector_d &result, unsigned &neval, int dirty,
                  Convergence<vector_d>::type compare = NULL);

int fubini_promb(const mdfunction<vector_z,4> &f, const array<double,4> &a,
                 const array<double,4> &b, double epsabs, double epsrel,
                 vector_z &result, unsigned &neval, int dirty,
                  Convergence<vector_z>::type compare = NULL);
int fubini_promb(const mdfunction<vector_z,5> &f, const array<double,5> &a,
                 const array<double,5> &b, double epsabs, double epsrel,
                 vector_z &result, unsigned &neval, int dirty,
                  Convergence<vector_z>::type compare = NULL);
int fubini_promb(const mdfunction<vector_z,6> &f, const array<double,6> &a,
                 const array<double,6> &b, double epsabs, double epsrel,
                 vector_z &result, unsigned &neval, int dirty,
                  Convergence<vector_z>::type compare = NULL);

/**@}*/

/**
   @name Fubini integration using Clenshaw-Curtis quadrature
*/
/**@{
   Quadrature rules for 1,2 and 3 dimensional functions returning double
   precision, complex double precision or vectors. P indicates that
   integration assumes a periodic function.

   @param f the functor to integrate
   @param a the lower limits
   @param b the upper limits
   @param epsabs the absolute error
   @param epsrel the relative error
   @param result the integral value
   @param neval  value incremented with total function evaluations
   @param dirty check convergence on the first dirty values
   @param compare use another function to determain convergence
*/

int fubini_cc(const mdfunction<double,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int fubini_cc(const mdfunction<double,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int fubini_cc(const mdfunction<double,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);

int fubini_cc(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int fubini_cc(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int fubini_cc(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);

int fubini_cc(const mdfunction<vector_d,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int fubini_cc(const mdfunction<vector_d,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int fubini_cc(const mdfunction<vector_d,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);

int fubini_cc(const mdfunction<vector_z,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int fubini_cc(const mdfunction<vector_z,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int fubini_cc(const mdfunction<vector_z,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);

int fubini_pcc(const mdfunction<double,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int fubini_pcc(const mdfunction<double,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int fubini_pcc(const mdfunction<double,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);

int fubini_pcc(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int fubini_pcc(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int fubini_pcc(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);

int fubini_pcc(const mdfunction<vector_d,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int fubini_pcc(const mdfunction<vector_d,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int fubini_pcc(const mdfunction<vector_d,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);

int fubini_pcc(const mdfunction<vector_z,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int fubini_pcc(const mdfunction<vector_z,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int fubini_pcc(const mdfunction<vector_z,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);

int fubini_cc(const mdfunction<double,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int fubini_cc(const mdfunction<double,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int fubini_cc(const mdfunction<double,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);

int fubini_cc(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int fubini_cc(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int fubini_cc(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);

int fubini_cc(const mdfunction<vector_d,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int fubini_cc(const mdfunction<vector_d,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int fubini_cc(const mdfunction<vector_d,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);

int fubini_cc(const mdfunction<vector_z,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int fubini_cc(const mdfunction<vector_z,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int fubini_cc(const mdfunction<vector_z,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);

int fubini_pcc(const mdfunction<double,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int fubini_pcc(const mdfunction<double,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int fubini_pcc(const mdfunction<double,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);

int fubini_pcc(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int fubini_pcc(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int fubini_pcc(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);

int fubini_pcc(const mdfunction<vector_d,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int fubini_pcc(const mdfunction<vector_d,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int fubini_pcc(const mdfunction<vector_d,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);

int fubini_pcc(const mdfunction<vector_z,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int fubini_pcc(const mdfunction<vector_z,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int fubini_pcc(const mdfunction<vector_z,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);


/**@}*/

/* =========================================================================
   Cubature integration
   ========================================================================= */
/**
   @name Cubature integration using Trapezoidal quadrature
*/
/**@{
   Quadrature rules for 1,2 and 3 dimensional functions returning double
   precision, complex double precision or vectors. P indicates that
   integration assumes a periodic function.

   @param f the functor to integrate
   @param a the lower limits
   @param b the upper limits
   @param epsabs the absolute error
   @param epsrel the relative error
   @param result the integral value
   @param neval  value incremented with total function evaluations
   @param dirty check convergence on the first dirty values
   @param compare use another function to determain convergence
*/

int cube_tr(const mdfunction<double,1> &f, const array<double,1> &a,
            const array<double,1> &b, double epsabs, double epsrel,
            double &result, unsigned &neval);
int cube_tr(const mdfunction<double,2> &f, const array<double,2> &a,
            const array<double,2> &b, double epsabs, double epsrel,
            double &result, unsigned &neval);
int cube_tr(const mdfunction<double,3> &f, const array<double,3> &a,
            const array<double,3> &b, double epsabs, double epsrel,
            double &result, unsigned &neval);
int cube_tr(const mdfunction<double,4> &f, const array<double,4> &a,
            const array<double,4> &b, double epsabs, double epsrel,
            double &result, unsigned &neval);
int cube_tr(const mdfunction<double,5> &f, const array<double,5> &a,
            const array<double,5> &b, double epsabs, double epsrel,
            double &result, unsigned &neval);
int cube_tr(const mdfunction<double,6> &f, const array<double,6> &a,
            const array<double,6> &b, double epsabs, double epsrel,
            double &result, unsigned &neval);

int cube_tr(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
            const array<double,1> &b, double epsabs, double epsrel,
            dcomplex &result, unsigned &neval);
int cube_tr(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
            const array<double,2> &b, double epsabs, double epsrel,
            dcomplex &result, unsigned &neval);
int cube_tr(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
            const array<double,3> &b, double epsabs, double epsrel,
            dcomplex &result, unsigned &neval);
int cube_tr(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
            const array<double,4> &b, double epsabs, double epsrel,
            dcomplex &result, unsigned &neval);
int cube_tr(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
            const array<double,5> &b, double epsabs, double epsrel,
            dcomplex &result, unsigned &neval);
int cube_tr(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
            const array<double,6> &b, double epsabs, double epsrel,
            dcomplex &result, unsigned &neval);

int cube_tr(const mdfunction<vector_d,1> &f, const array<double,1> &a,
            const array<double,1> &b, double epsabs, double epsrel,
            vector_d &result, unsigned &neval, int dirty,
             Convergence<vector_d>::type compare = NULL);
int cube_tr(const mdfunction<vector_d,2> &f, const array<double,2> &a,
            const array<double,2> &b, double epsabs, double epsrel,
            vector_d &result, unsigned &neval, int dirty,
             Convergence<vector_d>::type compare = NULL);
int cube_tr(const mdfunction<vector_d,3> &f, const array<double,3> &a,
            const array<double,3> &b, double epsabs, double epsrel,
            vector_d &result, unsigned &neval, int dirty,
             Convergence<vector_d>::type compare = NULL);
int cube_tr(const mdfunction<vector_d,4> &f, const array<double,4> &a,
            const array<double,4> &b, double epsabs, double epsrel,
            vector_d &result, unsigned &neval, int dirty,
             Convergence<vector_d>::type compare = NULL);
int cube_tr(const mdfunction<vector_d,5> &f, const array<double,5> &a,
            const array<double,5> &b, double epsabs, double epsrel,
            vector_d &result, unsigned &neval, int dirty,
             Convergence<vector_d>::type compare = NULL);
int cube_tr(const mdfunction<vector_d,6> &f, const array<double,6> &a,
            const array<double,6> &b, double epsabs, double epsrel,
            vector_d &result, unsigned &neval, int dirty,
             Convergence<vector_d>::type compare = NULL);

int cube_tr(const mdfunction<vector_z,1> &f, const array<double,1> &a,
            const array<double,1> &b, double epsabs, double epsrel,
            vector_z &result, unsigned &neval, int dirty,
             Convergence<vector_z>::type compare = NULL);
int cube_tr(const mdfunction<vector_z,2> &f, const array<double,2> &a,
            const array<double,2> &b, double epsabs, double epsrel,
            vector_z &result, unsigned &neval, int dirty,
             Convergence<vector_z>::type compare = NULL);
int cube_tr(const mdfunction<vector_z,3> &f, const array<double,3> &a,
            const array<double,3> &b, double epsabs, double epsrel,
            vector_z &result, unsigned &neval, int dirty,
             Convergence<vector_z>::type compare = NULL);
int cube_tr(const mdfunction<vector_z,4> &f, const array<double,4> &a,
            const array<double,4> &b, double epsabs, double epsrel,
            vector_z &result, unsigned &neval, int dirty,
             Convergence<vector_z>::type compare = NULL);
int cube_tr(const mdfunction<vector_z,5> &f, const array<double,5> &a,
            const array<double,5> &b, double epsabs, double epsrel,
            vector_z &result, unsigned &neval, int dirty,
             Convergence<vector_z>::type compare = NULL);
int cube_tr(const mdfunction<vector_z,6> &f, const array<double,6> &a,
            const array<double,6> &b, double epsabs, double epsrel,
            vector_z &result, unsigned &neval, int dirty,
             Convergence<vector_z>::type compare = NULL);

int cube_ptr(const mdfunction<double,1> &f, const array<double,1> &a,
             const array<double,1> &b, double epsabs, double epsrel,
             double &result, unsigned &neval);
int cube_ptr(const mdfunction<double,2> &f, const array<double,2> &a,
             const array<double,2> &b, double epsabs, double epsrel,
             double &result, unsigned &neval);
int cube_ptr(const mdfunction<double,3> &f, const array<double,3> &a,
             const array<double,3> &b, double epsabs, double epsrel,
             double &result, unsigned &neval);
int cube_ptr(const mdfunction<double,4> &f, const array<double,4> &a,
             const array<double,4> &b, double epsabs, double epsrel,
             double &result, unsigned &neval);
int cube_ptr(const mdfunction<double,5> &f, const array<double,5> &a,
             const array<double,5> &b, double epsabs, double epsrel,
             double &result, unsigned &neval);
int cube_ptr(const mdfunction<double,6> &f, const array<double,6> &a,
             const array<double,6> &b, double epsabs, double epsrel,
             double &result, unsigned &neval);

int cube_ptr(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
             const array<double,1> &b, double epsabs, double epsrel,
             dcomplex &result, unsigned &neval);
int cube_ptr(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
             const array<double,2> &b, double epsabs, double epsrel,
             dcomplex &result, unsigned &neval);
int cube_ptr(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
             const array<double,3> &b, double epsabs, double epsrel,
             dcomplex &result, unsigned &neval);
int cube_ptr(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
             const array<double,4> &b, double epsabs, double epsrel,
             dcomplex &result, unsigned &neval);
int cube_ptr(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
             const array<double,5> &b, double epsabs, double epsrel,
             dcomplex &result, unsigned &neval);
int cube_ptr(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
             const array<double,6> &b, double epsabs, double epsrel,
             dcomplex &result, unsigned &neval);

int cube_ptr(const mdfunction<vector_d,1> &f, const array<double,1> &a,
             const array<double,1> &b, double epsabs, double epsrel,
             vector_d &result, unsigned &neval, int dirty,
              Convergence<vector_d>::type compare = NULL);
int cube_ptr(const mdfunction<vector_d,2> &f, const array<double,2> &a,
             const array<double,2> &b, double epsabs, double epsrel,
             vector_d &result, unsigned &neval, int dirty,
              Convergence<vector_d>::type compare = NULL);
int cube_ptr(const mdfunction<vector_d,3> &f, const array<double,3> &a,
             const array<double,3> &b, double epsabs, double epsrel,
             vector_d &result, unsigned &neval, int dirty,
              Convergence<vector_d>::type compare = NULL);
int cube_ptr(const mdfunction<vector_d,4> &f, const array<double,4> &a,
             const array<double,4> &b, double epsabs, double epsrel,
             vector_d &result, unsigned &neval, int dirty,
              Convergence<vector_d>::type compare = NULL);
int cube_ptr(const mdfunction<vector_d,5> &f, const array<double,5> &a,
             const array<double,5> &b, double epsabs, double epsrel,
             vector_d &result, unsigned &neval, int dirty,
              Convergence<vector_d>::type compare = NULL);
int cube_ptr(const mdfunction<vector_d,6> &f, const array<double,6> &a,
             const array<double,6> &b, double epsabs, double epsrel,
             vector_d &result, unsigned &neval, int dirty,
              Convergence<vector_d>::type compare = NULL);

int cube_ptr(const mdfunction<vector_z,1> &f, const array<double,1> &a,
             const array<double,1> &b, double epsabs, double epsrel,
             vector_z &result, unsigned &neval, int dirty,
              Convergence<vector_z>::type compare = NULL);
int cube_ptr(const mdfunction<vector_z,2> &f, const array<double,2> &a,
             const array<double,2> &b, double epsabs, double epsrel,
             vector_z &result, unsigned &neval, int dirty,
              Convergence<vector_z>::type compare = NULL);
int cube_ptr(const mdfunction<vector_z,3> &f, const array<double,3> &a,
             const array<double,3> &b, double epsabs, double epsrel,
             vector_z &result, unsigned &neval, int dirty,
              Convergence<vector_z>::type compare = NULL);
int cube_ptr(const mdfunction<vector_z,4> &f, const array<double,4> &a,
             const array<double,4> &b, double epsabs, double epsrel,
             vector_z &result, unsigned &neval, int dirty,
              Convergence<vector_z>::type compare = NULL);
int cube_ptr(const mdfunction<vector_z,5> &f, const array<double,5> &a,
             const array<double,5> &b, double epsabs, double epsrel,
             vector_z &result, unsigned &neval, int dirty,
              Convergence<vector_z>::type compare = NULL);
int cube_ptr(const mdfunction<vector_z,6> &f, const array<double,6> &a,
             const array<double,6> &b, double epsabs, double epsrel,
             vector_z &result, unsigned &neval, int dirty,
              Convergence<vector_z>::type compare = NULL);
/**@}*/

/**
   @name Cubature integration with Romberg enhancement
*/
/**@{
   Quadrature rules for 1,2 and 3 dimensional functions returning double
   precision, complex double precision or vectors. P indicates that
   integration assumes a periodic function.

   @param f the functor to integrate
   @param a the lower limits
   @param b the upper limits
   @param epsabs the absolute error
   @param epsrel the relative error
   @param result the integral value
   @param neval  value incremented with total function evaluations
   @param dirty check convergence on the first dirty values
   @param compare use another function to determain convergence
*/

int cube_romb(const mdfunction<double,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int cube_romb(const mdfunction<double,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int cube_romb(const mdfunction<double,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int cube_romb(const mdfunction<double,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int cube_romb(const mdfunction<double,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);
int cube_romb(const mdfunction<double,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              double &result, unsigned &neval);

int cube_romb(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int cube_romb(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int cube_romb(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int cube_romb(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int cube_romb(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);
int cube_romb(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              dcomplex &result, unsigned &neval);

int cube_romb(const mdfunction<vector_d,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int cube_romb(const mdfunction<vector_d,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int cube_romb(const mdfunction<vector_d,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int cube_romb(const mdfunction<vector_d,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int cube_romb(const mdfunction<vector_d,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);
int cube_romb(const mdfunction<vector_d,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              vector_d &result, unsigned &neval, int dirty,
               Convergence<vector_d>::type compare = NULL);

int cube_romb(const mdfunction<vector_z,1> &f, const array<double,1> &a,
              const array<double,1> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int cube_romb(const mdfunction<vector_z,2> &f, const array<double,2> &a,
              const array<double,2> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int cube_romb(const mdfunction<vector_z,3> &f, const array<double,3> &a,
              const array<double,3> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int cube_romb(const mdfunction<vector_z,4> &f, const array<double,4> &a,
              const array<double,4> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int cube_romb(const mdfunction<vector_z,5> &f, const array<double,5> &a,
              const array<double,5> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);
int cube_romb(const mdfunction<vector_z,6> &f, const array<double,6> &a,
              const array<double,6> &b, double epsabs, double epsrel,
              vector_z &result, unsigned &neval, int dirty,
               Convergence<vector_z>::type compare = NULL);

int cube_promb(const mdfunction<double,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int cube_promb(const mdfunction<double,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int cube_promb(const mdfunction<double,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int cube_promb(const mdfunction<double,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int cube_promb(const mdfunction<double,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);
int cube_promb(const mdfunction<double,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               double &result, unsigned &neval);

int cube_promb(const mdfunction<dcomplex,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int cube_promb(const mdfunction<dcomplex,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int cube_promb(const mdfunction<dcomplex,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int cube_promb(const mdfunction<dcomplex,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int cube_promb(const mdfunction<dcomplex,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);
int cube_promb(const mdfunction<dcomplex,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               dcomplex &result, unsigned &neval);

int cube_promb(const mdfunction<vector_d,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int cube_promb(const mdfunction<vector_d,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int cube_promb(const mdfunction<vector_d,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int cube_promb(const mdfunction<vector_d,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int cube_promb(const mdfunction<vector_d,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);
int cube_promb(const mdfunction<vector_d,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               vector_d &result, unsigned &neval, int dirty,
                Convergence<vector_d>::type compare = NULL);

int cube_promb(const mdfunction<vector_z,1> &f, const array<double,1> &a,
               const array<double,1> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int cube_promb(const mdfunction<vector_z,2> &f, const array<double,2> &a,
               const array<double,2> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int cube_promb(const mdfunction<vector_z,3> &f, const array<double,3> &a,
               const array<double,3> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int cube_promb(const mdfunction<vector_z,4> &f, const array<double,4> &a,
               const array<double,4> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int cube_promb(const mdfunction<vector_z,5> &f, const array<double,5> &a,
               const array<double,5> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
int cube_promb(const mdfunction<vector_z,6> &f, const array<double,6> &a,
               const array<double,6> &b, double epsabs, double epsrel,
               vector_z &result, unsigned &neval, int dirty,
                Convergence<vector_z>::type compare = NULL);
/**@}*/


CLOSE_NUMINT_NAMESPACE

#endif
