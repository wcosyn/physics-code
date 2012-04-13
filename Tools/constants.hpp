/*! \file constants.hpp 
 * \brief defines some useful constants
 * \author Wim Cosyn
 * \date 16/08/2011
 */

#ifndef CONST_H
#define CONST_H

#define PI 3.14159265 /*!< \def defines constant pi */
#define PREC 1e-03 /*!< \def defines precision of the integrations */
#define MASSN 939.56536 /*!< \def defines proton mass [MeV] */
#define MASSP 938.27231 /*!< \def defines neutron mass [MeV] */
#define MASSPI 139.57 /*!< \def defines charged pion mass [MeV] */
#define MASSRHO 775.49 /*!< \def defines rho0 mass [MeV] */
#define HBARC 197.327  /*!< \def defines hbar*c [MeV*fm] */
#define ALPHA 0.00729735253 /*!< \def defines finestructure constant [] */
#define I complex<double>(0.,1.) /*!< \def defines complex number I */
#define INVHBARC 0.00506770453255 /*!< \def defines (hbar*c)^-1 [(MeV*fm)^-1] */
#define MUP 2.79 /*!< \def defines proton anomalous magnetic moment */
#define DEGRTORAD 0.0174532925 /*!< \def defines Pi/180 */
#define RADTODEGR 57.2957795 /*!< \def defines 180/Pi */
#define MASSD 1875.6
/*! underflow in the code:
   * double smaller than underflow
   * should be treated as being zero */
#define STRANGEUFLOW 1.0e-9 
#define M_P 938.272       // proton
#define M_N 939.565       // neutron
#define M_KP 493.677      // kaon^+
#define M_K0 497.672      // kaon^0
#define M_L 1115.683      // lambda
#define M_S0 1192.642     // sigma^0
#define M_SP 1189.37      // sigma^+
#define M_SM 1197.449     // sigma^-
#define M_RHO 775.5       // rho
#define M_OMEGA 782.65    // omega
#define M_PHI 1019.46     // phi
#define M_PIP 139.57      // pi^+
#define M_PIM 139.57      // pi^-
#define M_PI0 134.977     // pi^0

#define HOMEDIR "/home/wim/Code/share"


#endif
