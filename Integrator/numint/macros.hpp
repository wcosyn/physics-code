/**
 * @file lib/numint/macros.hpp
 * @brief preprocessor macros used in the numint namespace
 *
 * @author Klaas Vantournhout
 * @author Ghent University (2004-2009)
 * @author GSI Helmholtzzentrum für Schwerionenforschung GmbH (Darmstadt)
 * @author k.vantournhout@gsi.de
 * @date creation: 11/04/2011
 * @date last update: 11/04/2011
 */


#ifndef NUMINT_MACROS_HPP
#define NUMINT_MACROS_HPP 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// namespaces
#define OPEN_NUMINT_NAMESPACE namespace numint { //!< opens numint namespace
#define CLOSE_NUMINT_NAMESPACE }                 //!< closes numint namespace

#define OPEN_EMPTY_NAMESPACE namespace {         //!< opens an empty namespace
#define CLOSE_EMPTY_NAMESPACE }                  //!< closes an empty namespace

#ifdef TOTAL_DIM
#define CUBATURE_DIM TOTAL_DIM                   //!< the cubature is evaluated for CUBATURE_DIM dimensions
#else
#define CUBATURE_DIM 3                           //!< the cubature is evaluated for 3 dimensions
#endif

#define PI      3.1415926535897932385
#define TWOPI   6.2831853071795864769
#define SQRTPI  1.7724538509055160273
#define SQRT2PI 2.5066282746310005024
#define PI2     9.8696044010893586188
#define INVSQRTPI 0.564189583547756285
#define SQRT2   1.4142135623730950488
#define SQRT3   1.7320508075688772935
#define EXP1    2.7182818284590452354

#endif
