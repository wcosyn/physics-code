/*
 * TDeuteron.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: February,5 2009
 *
 */



#include "TDeuteron.h"
#include "TYukawaPWF.h"
#include "TGrossPWF.h"
#include "TCstPWF.h"
#include "TInterpolatingWavefunction.h"
// #include "TStrangePolarization.h"
//#include "TKinematics2to3.h"
//#include <numtoa.h>
#include <cstring>
#include <iostream>
#include <cassert>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include "constants.hpp"
#include <Math/SpecFuncMathMore.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>



using std::cout; using std::cerr; using std::strcmp;
using TMath::Sqrt; using TMath::Exp; using TMath::Pi;
using TMath::Power; using TMath::Abs;
using std::complex;

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TDeuteron                                                             //
//                                                                       //
// Namespace for deuteron related functionality                          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TDeuteron)

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TDeuteron::Wavefunction                                               //
//                                                                       //
// This class represents a generic deuteron wavefunction in strangecalc. //
//                                                                       //
// Different wavefunctions have different implementations. Wavefunction //
// deals with this by owning a concrete instance of the                  //
// TWavefunctionImplementation interface.                                //
//                                                                       //
// For safety reasons, users are not allowed to create a wavefunction    //
// of their own 'in code'. Therefore all meaningful constructors are     //
// hidden. The named constructors are the only way of creating a         //
// Wavefunction.                                                        //
//                                                                       //
// The deuteron with fourvector                                          
// BEGIN_LATEX
// p_{D} = (M_{D}, #vec{p}_{D} )
// END_LATEX
// is composed of two nucleons, labeled 1 and 2.
// We define
// BEGIN_LATEX
// p = #frac{p_{1}-p_{2}}{2}
// p_{D} = p_{1}+p_{2}
// END_LATEX
// and
// BEGIN_LATEX
// p_{1} = #frac{p_{D}}{2} + p
// p_{2} = #frac{p_{D}}{2} - p
// END_LATEX
//
// Unless explicitely specified, we evaluate the deuteron state in its
// restframe.
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TDeuteron::Wavefunction)

//_____________________________________________________________________
const double TDeuteron::Wavefunction::kMd = 1877.05;

//_____________________________________________________________________
TDeuteron::Wavefunction::Wavefunction(const TString& name,
			     const TString& title)
  : TNamed(name,title), fRelativistic(0), fImplementation(0), yukawa(0)
{
  // Default constructor is private
  // User are only allowed to use named constructors
}

//_____________________________________________________________________
TDeuteron::Wavefunction::Wavefunction(TRootIOCtor*)
  : TNamed(), fRelativistic(0), fImplementation(0), yukawa(0)
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
TDeuteron::Wavefunction::Wavefunction(const Wavefunction& rhs)
  : TNamed(rhs), fRelativistic(rhs.fRelativistic), fImplementation(0), yukawa(rhs.yukawa)
{
  // Copy constructor

  fImplementation = rhs.fImplementation->Clone();
}

//_____________________________________________________________________
TDeuteron::Wavefunction& TDeuteron::Wavefunction::operator=(const Wavefunction& rhs)
{
  // Assignment

  if( this!=&rhs ) { // avoid self-assignment
    TNamed::operator=(rhs);
    fRelativistic = rhs.fRelativistic;
    yukawa=rhs.yukawa;
    delete fImplementation;
    fImplementation = rhs.fImplementation->Clone();
  }

  return *this;
}

//_____________________________________________________________________
TDeuteron::Wavefunction::~Wavefunction()
{
  // Destructor

  delete fImplementation;
}

//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateCDBonnWavefunction()
{
  // Non-relativistic deuteron wavefunction generated with the CD-Bonn potential.
  // For more info see PRC63(2001)024001.

  Wavefunction *wf = new Wavefunction("CDBonn");

  wf->fRelativistic = false;
  wf->yukawa = true;

  static double c[11] = { +.88472985e+0,
			  -.26408759e+0,
			  -.44114404e-1,
			  -.14397512e+2,
			  +.85591256e+2,
			  -.31876761e+3,
			  +.70336701e+3,
			  -.90049586e+3,
			  +.66145441e+3,
			  -.25958894e+3,
			  0. };
  static double d[11] = { +.22623762e-1,
			  -.50471056e+0,
			  +.56278897e+0,
			  -.16079764e+2,
			  +.11126803e+3,
			  -.44667490e+3,
			  +.10985907e+4,
			  -.16114995e+4,
			  0., 0., 0.};
  
  wf->fImplementation = new TYukawaPWF(11,.2315380,.9,c,d);

  return wf;  
}

//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateCDBonn87Wavefunction()
{
//  R. Machleidt, K. Holinde, C. Elster, 
//  Phys. Rept. 149 (1987) 1-89. See pages 64-68.

  Wavefunction *wf = new Wavefunction("CDBonn87");

  wf->fRelativistic = false;
  wf->yukawa = true;

  static double c[11] = { +0.90457337E+00,
			  -0.35058661E+00,
			  -0.17635927E+00,
			  -0.10418261E+02,
			  +0.45089439E+02,
			  -0.14861947E+03,
			  +0.31779642E+03,
			  -0.37496518E+03,
			  +0.22560032E+03,
			  -0.54858290E+02,
			  0. };
  static double d[11] = { +0.24133026E-01,
			  -0.64430531E+00,
			  + 0.51093352E+00,
			  -0.54419065E+01,
			  +0.15872034E+02,
			  -0.14742981E+02,
			  +0.44956539E+01,
			  -0.71152863E-01,
			  0., 0., 0.};
  
  wf->fImplementation = new TYukawaPWF(11,0.231609,.9,c,d);

  return wf;  
}
//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateParisWavefunction()
{
  // Non-relativistic deuteron wavefunction generated with the Paris potential.
  // For more info see PLB101(1981)139 and PRC21(1980)861.

  Wavefunction *wf = new Wavefunction("Paris");
  wf->fRelativistic = false;
  wf->yukawa = true;
  
  static double c[13] = { +.88688076e-0,
			  -.34717093e+0,
			  -.30502380e+1,
			  +.56207766e+2,
			  -.74957334e+3,
			  +.53365279e+4,
			  -.22706863e+5,
			  +.60434469e+5,
			  -.10292058e+6,
			  +.11223357e+6,
			  -.75925226e+5,
			  +.29059715e+5,
			  0. };
  static double d[13] = { +.23135193e-1,
			  -.85604572e+0,
			  +.56068193e+1,
			  -.69462922e+2,
			  +.41631118e+3,
			  -.12546621e+4,
			  +.12387830e+4,
			  +.33739172e+4,
			  -.13041151e+5,
			  +.19512524e+5,
			  0., 0., 0.};
  
  wf->fImplementation = new TYukawaPWF(13,.23162461,1.,c,d);

  return wf;  
}
//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateAV18Wavefunction()
{
  // Non-relativistic deuteron wavefunction generated with the AV18 potential.

  Wavefunction *wf = new Wavefunction("AV18");

  wf->fRelativistic = false;
  wf->yukawa = true;
  
  static double c[12] = { 0.706699E+00/Sqrt(2./Pi()),
			  -0.169743E+00/Sqrt(2./Pi()),
			  0.112368E+01/Sqrt(2./Pi()),
			  -0.852995E+01/Sqrt(2./Pi()),
			  0.195033E+02/Sqrt(2./Pi()),
			  -0.757831E+02/Sqrt(2./Pi()),
			  0.283739E+03/Sqrt(2./Pi()),
			  -0.694734E+03/Sqrt(2./Pi()),
			  0.885257E+03/Sqrt(2./Pi()),
			  -0.720739E+03/Sqrt(2./Pi()),
			  0.412969E+03/Sqrt(2./Pi()),
			  -0.103336E+03/Sqrt(2./Pi())};
  static double d[12] = { 0.176655E-01/Sqrt(2./Pi()),
			  -0.124551E+00/Sqrt(2./Pi()),
			  -0.108815E+01/Sqrt(2./Pi()),
			  0.384848E+01/Sqrt(2./Pi()),
			  -0.852442E+01/Sqrt(2./Pi()),
			  0.209435E+02/Sqrt(2./Pi()),
			  -0.490728E+02/Sqrt(2./Pi()),
			  0.577382E+02/Sqrt(2./Pi()),
			  -0.127114E+01/Sqrt(2./Pi()),
			  -0.628361E+02/Sqrt(2./Pi()),
			  0.581016E+02/Sqrt(2./Pi()),
			  -0.177062E+02/Sqrt(2./Pi())};
			  
  static double m[12] = { 0.2316,
			  1.0,
			  1.5,
			  2.0,
			  2.5,
			  3.5,
			  4.5,
			  5.5,
			  6.5,
			  8.0,
			  9.5,
			  11.0};
  
  wf->fImplementation = new TYukawaPWF(12,c,d,m);

  return wf;  
}

//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateAV18bWavefunction()
{
  // Non-relativistic deuteron wavefunction generated with the AV1b8 potential.

  Wavefunction *wf = new Wavefunction("AV18b");

  wf->fRelativistic = false;
  wf->yukawa = true;
  
  static double c[12] = { 0.105252223e+02/(4.*Pi()),
			  0.124352529e+02/(4.*Pi()),
			  -0.687541641e+02/(4.*Pi()),
			  0.239111042e+03/(4.*Pi()),
			  -0.441014422e+03/(4.*Pi()),
			  0.300140328e+03/(4.*Pi()),
			  -0.230639939e+03/(4.*Pi()),
			  0.409671540e+03/(4.*Pi()),
			  -0.733453611e+03/(4.*Pi()),
			  0.123506081e+04/(4.*Pi()),
			  -0.120520606e+04/(4.*Pi()),0.};
  static double d[12] = { 0.280995496e+00/(4.*Pi()),
			  0.334117629e-01/(4.*Pi()),
			  -0.727192237e+00/(4.*Pi()),
			  -0.302809607e+01/(4.*Pi()),
			  -0.903824982e+01/(4.*Pi()),
			  0.496045967e+01/(4.*Pi()),
			  -0.271985613e+02/(4.*Pi()),
			  0.125334598e+03/(4.*Pi()),
			  -0.346742235e+03/(4.*Pi()),0.,0.,0.};
			  
  static double m[12] = {0.232500e+00,
			 0.500000e+00,
			 0.800000e+00,
			 0.120000e+01,
			 0.160000e+01,
			 0.200000e+01,
			 0.400000e+01,
			 0.600000e+01,
			 0.100000e+02,
			 0.140000e+02,
			 0.180000e+02,
			 0.220000e+02};
  
  wf->fImplementation = new TYukawaPWF(12,c,d,m);

  return wf;  
}
//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateGrossIIbWavefunction()
{
  // Relativistic deuteron wavefunction generated with OBE model.
  // Model IIb from Gross et al., PRC 45(1992)2094.
  // Parameters from J.W. Van Orden (email February,5 2009)
  
  Wavefunction *wf = new Wavefunction("GrossIIb");

  wf->fRelativistic = true;

  static double alpha[3] = { 45.702,   // L=0
			     138.0,    // L=1
			     45.702 }; // L=2

  static double par[4*12] = { +1.234739094966355e+01, // L=0
			      -9.791825483944739e-01,
			      -1.331496538056728e+01,
			      +2.035989847774199e+01,
			      -5.498836869256259e+01,
			      -3.352082107736260e+01,
			      +1.428164210651713e+01,
			      +4.586138938279015e+01,
			      +5.365661937998652e+01,
			      +4.425711958893110e+01,
			      +2.710010081195388e+01,
			      +1.035634235843613e+01,
			      +3.044454986909504e-01, // L=2
			      -2.283891432756719e+00,
			      -1.255358797238707e+01,
			      +4.754848431159639e+00,
			      +3.447858734346326e+00,
			      +1.063375115847049e+00,
			      -7.154490446866921e-03,
			      -2.600002172348260e-01,
			      -2.119574479663605e-01,
			      -1.117025885337745e-01,
			      -4.100604363231088e-02,
			      -8.544248703450209e-03,
			      +3.586491947965324e-02, // L=1 triplet
			      -1.119415185418387e+00,
			      +1.241002812912811e+00,
			      +3.010547062642895e+01,
			      -5.392128643947384e+01,
			      -2.968692230392706e+01,
			      +2.391384855293073e+01,
			      +5.341402580624263e+01,
			      +5.342806472289331e+01,
			      +3.719349628437146e+01,
			      +1.845155663025572e+01,
			      +5.341935033039075e+00,
			      +9.015742962943808e-03, // L=1 singlet
			      -5.622985139281625e-01,
			      +5.907329424032518e+00,
			      -1.884995361121991e+01,
			      +7.362525765548103e+00,
			      +1.352404338899928e+01,
			      +7.268503404936552e+00,
			      +5.589065626656743e-01,
			      -2.629918666125826e+00,
			      -2.818935581112417e+00,
			      -1.678734932366940e+00,
			      -5.386768666258961e-01 };

  wf->fImplementation = new TGrossPWF(12,alpha,138.,
				      par,par+12,par+(2*12),par+(3*12));

  return wf;
}

//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateGrossL4Wavefunction()
{
  // Relativistic deuteron wavefunction generated with OBE model.
  // Model with #lambda=.4 from Buck & Gross, PRD 20(1979)2361.
    
  Wavefunction *wf = new Wavefunction("GrossL4");

  wf->fRelativistic = true;
  wf->yukawa = false;
  
  static double alpha[3] = { 45.702,   // L=0
			     138.0,    // L=1
			     45.702 }; // L=2

  static double par[4*8] = { +.40627e+0*Sqrt(1e3), // L=0
			     -.48817e-1*Sqrt(1e3),
			     -.14180e+0*Sqrt(1e3),
			     -.24136e+1*Sqrt(1e3),
			     +.59324e+1*Sqrt(1e3),
			     -.13507e+2*Sqrt(1e3),
			     +.19513e+2*Sqrt(1e3),
			     -.12367e+2*Sqrt(1e3),
			     +.10589e-1*Sqrt(1e3), // L=2
			     -.70546e-1*Sqrt(1e3),
			     -.84650e+0*Sqrt(1e3),
			     +.35455e+1*Sqrt(1e3),
			     -.13965e+2*Sqrt(1e3),
			     +.29006e+2*Sqrt(1e3),
			     -.20446e+2*Sqrt(1e3),
			     -.21302e+2*Sqrt(1e3),
			     +.75280e-3*Sqrt(1e3), // L=1 triplet
			     +.10957e+0*Sqrt(1e3),
			     -.55089e+0*Sqrt(1e3),
			     +.41630e+1*Sqrt(1e3),
			     -.16092e+2*Sqrt(1e3),
			     +.34651e+2*Sqrt(1e3),
			     -.45555e+2*Sqrt(1e3),
			     +.34221e+2*Sqrt(1e3),
			     +.41566e-3*Sqrt(1e3), // L=1 singlet
			     +.61962e-1*Sqrt(1e3),
			     -.60043e+0*Sqrt(1e3),
			     +.48688e+1*Sqrt(1e3),
			     -.23789e+2*Sqrt(1e3),
			     +.69093e+2*Sqrt(1e3),
			     -.12215e+3*Sqrt(1e3),
			     +.12828e+3*Sqrt(1e3) };

  wf->fImplementation = new TGrossPWF(8,alpha,138.,
				      par,par+8,par+(2*8),par+(3*8));

  return wf;
}

//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateGrossWJC1Wavefunction()
{
  // Relativistic deuteron wavefunction generated with OBE model.
  // Model WJC-1 from Gross and Stadler, arXiv:1007.0778.
  
  // References for the parameters used in this function are to arXiv:1007.0778.
  Wavefunction *wf = new Wavefunction("GrossWJC1");

  wf->fRelativistic = true;
  wf->yukawa = false;
    
  // Table VII
  static int order[4] = { 17,   // u
			  16,   // w
			  12,   // v_t
			  12 }; // v_s
  
  // (3.44) for u and w,  m_{s0} in Table V for v_s and v_t
  static double leadingMass[4] = { 45.7159, // u
				   45.7159, // w
				   488.,    // v_t
				   633. };  // v_s

  // Table VII
  static double stepMass[4] = { 75.00,    // u
				80.00,    // w
				109.09,   // v_t
				109.09 }; // v_s

  // Table VII
  static double tailMass[4] = { 2000.,   // u
				2000.,   // w
				2000.,   // v_t
				2000. }; // v_s
  
  // Text below (3.45)
  static int exponent[4] = { 2,   // u
			     3,   // w
			     2,   // v_t
			     2 }; // v_s

  // Table VII
  static double par[57] = { 134.963 *sqrt(1.e-9),    // u
			    52.871 *sqrt(1.e-9),
			    -217.709 *sqrt(1.e-9),
			    1876.699 *sqrt(1.e-9),
			    -11369.449 *sqrt(1.e-9),
			    49427.176 *sqrt(1.e-9),
			    -156695.247 *sqrt(1.e-9),
			    369322.468 *sqrt(1.e-9),
			    -655367.189 *sqrt(1.e-9),
			    879178.453 *sqrt(1.e-9),
			    -887103.150 *sqrt(1.e-9),
			    662786.733 *sqrt(1.e-9),
			    -355720.583 *sqrt(1.e-9),
			    129738.381 *sqrt(1.e-9),
			    -28805.739 *sqrt(1.e-9),
			    2939.797 *sqrt(1.e-9),
			    -125.e-6 *sqrt(1.e-9),
			    23.813 *sqrt(1.e-9),     // w
			    32.709 *sqrt(1.e-9),
			    -111.381 *sqrt(1.e-9),
			    844.861 *sqrt(1.e-9),
			    -4376.965 *sqrt(1.e-9),
			    16489.297 *sqrt(1.e-9),
			    -45122.587 *sqrt(1.e-9),
			    91061.493 *sqrt(1.e-9),
			    -136672.748 *sqrt(1.e-9),
			    152411.033 *sqrt(1.e-9),
			    -124633.670 *sqrt(1.e-9),
			    72617.712 *sqrt(1.e-9),
			    -28549.974 *sqrt(1.e-9),
			    6790.953 *sqrt(1.e-9),
			    -738.436 *sqrt(1.e-9),
			    -100.e-6 *sqrt(1.e-9),
			    -27.120 *sqrt(1.e-9),    // v_t
			    233.155 *sqrt(1.e-9),
			    -998.893 *sqrt(1.e-9),
			    2679.013 *sqrt(1.e-9),
			    -4894.348 *sqrt(1.e-9),
			    6277.553 *sqrt(1.e-9),
			    -5674.481 *sqrt(1.e-9),
			    3550.952 *sqrt(1.e-9),
			    -1467.429 *sqrt(1.e-9),
			    360.796 *sqrt(1.e-9),
			    -40.010 *sqrt(1.e-9),
			    498.e-7 *sqrt(1.e-9),
			    -78.763 *sqrt(1.e-9),    // v_s
			    688.895 *sqrt(1.e-9),
			    -2910.298 *sqrt(1.e-9),
			    7690.084 *sqrt(1.e-9),
			    -13861.596 *sqrt(1.e-9),
			    17571.825 *sqrt(1.e-9),
			    -15720.926 *sqrt(1.e-9),
			    9747.293 *sqrt(1.e-9),
			    -3994.747 *sqrt(1.e-9),
			    975.086 *sqrt(1.e-9),
			    -107.489 *sqrt(1.e-9),
			    -394.e-7 *sqrt(1.e-9) };
  
  wf->fImplementation = new TCstPWF(order,leadingMass,stepMass,tailMass,exponent,
				    par,par+17,par+33,par+45);

  return wf;
}

//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateGrossWJC2Wavefunction()
{
  // Relativistic deuteron wavefunction generated with OBE model.
  // Model WJC-2 from Gross and Stadler, arXiv:1007.0778.
  
  // References for the parameters used in this function are to arXiv:1007.0778.
  Wavefunction *wf = new Wavefunction("GrossWJC2");

  wf->fRelativistic = true;
  wf->yukawa = false;
  
  // Table VIII
  static int order[4] = { 30,   // u
			  18,   // w
			  21,   // v_t
			  18 }; // v_s
  
  // (3.44) for u and w,  m_{s0} in Table VI for v_s and v_t
  static double leadingMass[4] = { 45.7159, // u
				   45.7159, // w
				   551.,    // v_t
				   977. };  // v_s

  // Table VIII
  static double stepMass[4] = { 344.83,  // u
				70.59,   // w
				60.00,   // v_t
				70.59 }; // v_s

  // Table VIII
  static double tailMass[4] = { 11000.,  // u
				2000.,   // w
				2000.,   // v_t
				2000. }; // v_s
  
  // Text below (3.45)
  static int exponent[4] = { 2,   // u
			     3,   // w
			     2,   // v_t
			     2 }; // v_s

  // Table VIII
  static double par[87] = { 181.569 *sqrt(1.e-9),  // u
			    -0.817 *sqrt(1.e-9),
			    11.779 *sqrt(1.e-9),
			    -64.782 *sqrt(1.e-9),
			    247.517 *sqrt(1.e-9),
			    -772.915 *sqrt(1.e-9),
			    1948.950 *sqrt(1.e-9),
			    -3813.276 *sqrt(1.e-9),
			    5397.722 *sqrt(1.e-9),
			    -4637.491 *sqrt(1.e-9),
			    533.832 *sqrt(1.e-9),
			    3806.187 *sqrt(1.e-9),
			    -3119.875 *sqrt(1.e-9),
			    -2058.321 *sqrt(1.e-9),
			    3773.554 *sqrt(1.e-9),
			    1092.507 *sqrt(1.e-9),
			    -3919.342 *sqrt(1.e-9),
			    -736.912 *sqrt(1.e-9),
			    4018.336 *sqrt(1.e-9),
			    612.635 *sqrt(1.e-9),
			    -4213.343 *sqrt(1.e-9),
			    -260.753 *sqrt(1.e-9),
			    4581.603 *sqrt(1.e-9),
			    -1171.113 *sqrt(1.e-9),
			    -4559.802 *sqrt(1.e-9),
			    5566.714 *sqrt(1.e-9),
			    -2984.961 *sqrt(1.e-9),
			    814.407 *sqrt(1.e-9),
			    -92.509 *sqrt(1.e-9),
			    -118.e-11 *sqrt(1.e-9),
			    19.228 *sqrt(1.e-9),      // w
			    30.600 *sqrt(1.e-9),
			    -107.606 *sqrt(1.e-9),
			    897.969 *sqrt(1.e-9),
			    -5434.481 *sqrt(1.e-9),
			    25662.551 *sqrt(1.e-9),
			    -93771.705 *sqrt(1.e-9),
			    265664.025 *sqrt(1.e-9),
			    -582644.240 *sqrt(1.e-9),
			    986160.259 *sqrt(1.e-9),
			    -1280768.137 *sqrt(1.e-9),
			    1262888.436 *sqrt(1.e-9),
			    -927821.391 *sqrt(1.e-9),
			    491712.394 *sqrt(1.e-9),
			    -177560.697 *sqrt(1.e-9),
			    39095.262 *sqrt(1.e-9),
			    -3959.947 *sqrt(1.e-9),
			    -172.e-6 *sqrt(1.e-9),
			    33002.947 *sqrt(1.e-9),     // v_t
			    -580810.178 *sqrt(1.e-9),
			    4678729.666 *sqrt(1.e-9),
			    -22669425.846 *sqrt(1.e-9),
			    72755203.760 *sqrt(1.e-9),
			    -159377601.760 *sqrt(1.e-9),
			    232305154.041 *sqrt(1.e-9),
			    -193963168.786 *sqrt(1.e-9),
			    15154785.470 *sqrt(1.e-9),
			    160238368.984 *sqrt(1.e-9),
			    -139567635.494 *sqrt(1.e-9),
			    -52643827.848 *sqrt(1.e-9),
			    170318699.570 *sqrt(1.e-9),
			    -69098137.722 *sqrt(1.e-9),
			    -117233319.092 *sqrt(1.e-9),
			    192985378.690 *sqrt(1.e-9),
			    -138452165.069 *sqrt(1.e-9),
			    56855870.564 *sqrt(1.e-9),
			    -13049893.573 *sqrt(1.e-9),
			    1310791.818 *sqrt(1.e-9),
			    514.e-7 *sqrt(1.e-9),
			    -74492.403 *sqrt(1.e-9),    // v_s
			    1001357.141 *sqrt(1.e-9),
			    -6147239.139 *sqrt(1.e-9),
			    22510592.546 *sqrt(1.e-9),
			    -53617329.732 *sqrt(1.e-9),
			    83750225.006 *sqrt(1.e-9),
			    -77842609.919 *sqrt(1.e-9),
			    20071103.929 *sqrt(1.e-9),
			    48069105.480 *sqrt(1.e-9),
			    -61037252.628 *sqrt(1.e-9),
			    8837981.631 *sqrt(1.e-9),
			    50533251.186 *sqrt(1.e-9),
			    -65739754.705 *sqrt(1.e-9),
			    43130588.287 *sqrt(1.e-9),
			    -16799889.802 *sqrt(1.e-9),
			    3717624.361 *sqrt(1.e-9),
			    -363261.256 *sqrt(1.e-9),
			    -110.e-6 *sqrt(1.e-9) };
  
  wf->fImplementation = new TCstPWF(order,leadingMass,stepMass,tailMass,exponent,
				    par,par+30,par+48,par+69);

  return wf;
}

//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateTestWavefunction()
{
  // Dummy wavefunction:
  // * r-space: exp(r)*r
  // * p-space: gaussian

  Wavefunction *wf = new Wavefunction("Test");

  wf->fRelativistic = false;
  wf->yukawa = false;
    
  TInterpolatingWavefunction *implementation = new TInterpolatingWavefunction;

  double x,q;
  for(int i=0; i<101; ++i) {
    q = i*.05;
    x = Exp(-q)*q;
    implementation->AddUr(q,x);
    implementation->AddWr(q,x);
  }
  for(int i=0; i<101; ++i) {
    q = i*5.;
    x = TMath::Gaus(q,0,200.);
    implementation->AddUp(q,x);
    implementation->AddWp(q,x);
  }

  wf->fImplementation = implementation;
  
  return wf;
}

//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateWavefunction(const TString& name)
{
  // known types are * CDBonn           (cdbonn)
  //                 * Paris            (paris)
  //                 * Nijmegen III     (nijm3)
  //                 * Gross IIb        (grossIIb)
  //                 * Gross #lambda=.4 (grossL4)
  //                 * Gross WJC-1      (grossWJC1)
  //                 * Gross WJC-2      (grossWJC2)

  Wavefunction *wf=0;

  if(!name.CompareTo("cdbonn",TString::kIgnoreCase))
    wf = CreateCDBonnWavefunction();
  else if(!name.CompareTo("cdbonn87",TString::kIgnoreCase))
    wf = CreateCDBonn87Wavefunction();
  else if(!name.CompareTo("paris",TString::kIgnoreCase))
    wf = CreateParisWavefunction();
  else if(!name.CompareTo("av18",TString::kIgnoreCase))
    wf = CreateAV18Wavefunction();
  else if(!name.CompareTo("av18b",TString::kIgnoreCase))
    wf = CreateAV18bWavefunction();
  else if(!name.CompareTo("nijm3",TString::kIgnoreCase))
    wf = CreateNijmegenWavefunction();
  else if(!name.CompareTo("grossiib",TString::kIgnoreCase))
    wf = CreateGrossIIbWavefunction();
  else if(!name.CompareTo("grossl4",TString::kIgnoreCase))
    wf = CreateGrossL4Wavefunction();
  else if(!name.CompareTo("grosswjc1",TString::kIgnoreCase))
    wf = CreateGrossWJC1Wavefunction();
  else if(!name.CompareTo("grosswjc2",TString::kIgnoreCase))
    wf = CreateGrossWJC2Wavefunction();
  else if(!name.CompareTo("test",TString::kIgnoreCase))
    wf = CreateTestWavefunction();
  else
    cerr << "WARNING in TDeuteron::Wavefunction::CreateWavefunction(const TString&): "
	 << "unknown wavefunction.\n";
  
  return wf;
}

//_____________________________________________________________________
FourVector<GammaStructure> TDeuteron::Wavefunction::Vertex(const FourVector<double>& p2,
							   const FourVector<double>& pD) const
{
  // The covariant Blankenbecler and Cook vertex with particle 2 (4momentum p2) on-mass shell:
  // BEGIN_LATEX
  // #Gamma = F_{1} #gamma - #frac{F_{2}}{M_{N}} p - #frac{M_{N}-#gamma.p_{1}}{M_{N}} #left( F_{3} #gamma - #frac{F_{4}}{M_{N}} p #right).
  // END_LATEX
  // with BEGIN_LATEX p=#frac{p_{D}}{2}-p_{2} END_LATEX.
  //
  // The scalar form factors of the vertex are linked to the radial wave functions. See Eq.(46) in PRD 20(1979)2361.

  static double mass = (M_P + M_N)/2.; // average nucleon mass

  // momentum 4vector of particle 1
  FourVector<double> p1 = pD - p2;

  // relative momentum 4vector
  FourVector<double> p4vect = pD*(1./2.) - p2;

  // To find the scalar form factors we need the relative momentum and energy
  // in the rest frame of the deuteron.
  double ep;
  double p;
  if( (Abs(pD[1])<STRANGEUFLOW) && (Abs(pD[2])<STRANGEUFLOW) && (Abs(pD[3])<STRANGEUFLOW) ) {
    // We are in the rest frame
    ep = p2[0];
    p = Sqrt( p2[1]*p2[1] + p2[2]*p2[2] + p2[3]*p2[3] );
  } else {
    TLorentzRotation boost(-pD[1]/pD[0],-pD[2]/pD[0],-pD[3]/pD[0]);
    FourVector<double> p2CM = boost*p2;
    ep = p2CM[0];
    p = Sqrt( p2CM[1]*p2CM[1] + p2CM[2]*p2CM[2] + p2CM[3]*p2CM[3] );
  }

  // vertex form factors: (46) in PRD 20(1979)2361.
  double f1 = Pi()*Sqrt(2.*kMd)*(2*ep-kMd)
    * ( GetUp(p) - GetWp(p)/Sqrt(2.) + mass/p*Sqrt(3./2.)*GetVTp(p) );

  double f2 = Pi()*Sqrt(2.*kMd)*(2*ep-kMd)
    * ( mass/(ep+mass)*GetUp(p) + mass*(2*ep+mass)/Sqrt(2.)/p/p*GetWp(p)
  	+ mass/p*Sqrt(3./2.)*GetVTp(p) );

  double f3 = Pi()*Sqrt(2.*kMd)*ep*mass/p*Sqrt(3./2.)*GetVTp(p);

  double f4 = -Pi()*Sqrt(2.*kMd)*mass*mass/kMd
    * ( (2*ep-kMd)*(GetUp(p)/(ep+mass) - (ep+2.*mass)/Sqrt(2.)/p/p*GetWp(p))
  	+ kMd*Sqrt(3.)/p*GetVSp(p) );

  static GammaStructure unity(1.); // unit 4x4 matrix

  // vertex: 
  return ( f1*GMU
  	   - (f2/mass*p4vect)*unity
  	   - (unity - (p1*(1./mass))*GMU) * (f3*GMU - (f4/mass*p4vect)*unity) );
}

//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::SphericalHarmonic(unsigned int l, int m,
						 double theta, double phi)
{
  // Compute the (normalized) spherical harmonic:
  // BEGIN_LATEX
  // Y_{lm}(#theta,#phi) = #sqrt{#frac{(2l+1)}{4#pi}#frac{(l-m)!}{(l+m)!}} P_{l}^{m}(cos#theta)e^{im#phi}
  // END_LATEX
  // with l>=0 and |m|<=l
  //
  // The Condon-Shortley phase convention is used.
 
  // first calculate Y_{l|m}}(theta,phi=0)
  double Ylm = ROOT::Math::sph_legendre(l,TMath::Abs(m),theta);
  
  // (-1)^m when m is negative
  if( m<0 ) Ylm *= (1+2*(m%2));
      
  // Add the phi-dependent phase before returning
  return Ylm * exp(complex<double>(0,m*phi));
}

complex<double> TDeuteron::Wavefunction::SphericalHarmonicCos(unsigned int l, int m,
						 double costheta, double phi)
{
  // Compute the (normalized) spherical harmonic:
  // BEGIN_LATEX
  // Y_{lm}(#theta,#phi) = #sqrt{#frac{(2l+1)}{4#pi}#frac{(l-m)!}{(l+m)!}} P_{l}^{m}(cos#theta)e^{im#phi}
  // END_LATEX
  // with l>=0 and |m|<=l
  //
  // The Condon-Shortley phase convention is used.
 
  // first calculate Y_{l|m}}(theta,phi=0)
  double Ylm = gsl_sf_legendre_sphPlm(l,abs(m),costheta);
  
  // (-1)^m when m is negative
  if( m<0 ) Ylm *= (1+2*(m%2));
      
  // Add the phi-dependent phase before returning
  return Ylm * exp(complex<double>(0,m*phi));
}

//_____________________________________________________________________
double TDeuteron::Wavefunction::ClebschGordan(int two_j1,int two_m1,int two_j2, 
					      int two_m2,int two_j, int two_m)
{
  // Compute the Clebsch-Gordan coefficient: 
  // BEGIN_LATEX
  // < j_{1}m_{1} j_{2}m_{2} | jm >
  // END_LATEX
  // The Condon-Shortley phase convention is used.
  //
  // The arguments are the j_i and m_i multiplied by two.
  
  return gsl_sf_coupling_3j(two_j1,two_j2,two_j,two_m1,two_m2,-two_m)
    * Power(-1.,(two_j1-two_j2+two_m)/2.)
    * Sqrt(two_j+1.);
}

//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronPState(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const TVector3& p) const
{
  // Calculate the state of a deuteron in momentum space. 
  // deuteronPol is the deuteron polarization (times two)
  // nucleon1pol fixes the polarization of the first nucleon (the one that interacts with the photon). (times two)
  //nucleon2pol fixes spectator nucleon polarization (times two)
  
  // The Dnp vertex matrix element is
  // BEGIN_LATEX
  // #frac{1}{2E_{1}} #bar{u}(p_{2},m_{2})u(#dot{p}_{1},m_{1})#Gamma_{#mu}#phi^{#mu}_{D}(m_{D})
  // = #sum_{L} u_{L}(|#vec{p}_{1}|) #sum_{m_{L},m_{S}} <Lm_{L}Sm_{S}|1m_{D}><#frac{1}{2}m_{1}#frac{1}{2} m_{2}|1m_{S}> Y_{Lm_{L}}(#hat{p}_{1})
  // END_LATEX
  // with BEGIN_LATEX #dot{p}_{1} = ( #sqrt{#vec{p}_{1}^{2}+M_{1}^{2}},#vec{p}_{1}) END_LATEX the on-shell 4-momentum of particle 1.
  // The labelling particle 1 and 2 is arbitrary in the non-relativistic case (You can do the math).
//   if(IsRelativistic()) {
//     cerr << "ERROR in TDeuteron::Wavefunction::DeuteronPState(TSpinor::Polarization&,"
// 	 << " TVector3&&): "
// 	 << "has no meaning for relativistic wave functions.\n";
//     exit(1);
//   }

  
  complex<double> state=0.;
// L=2, mL=mD-mS, S=1, mS=m1+m2
  for(int mS=-1; mS<=1; ++mS)
    state -=
      ClebschGordan(4,deuteronPol-mS*2,2,mS*2,2,deuteronPol)
      *ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,mS*2)
      *SphericalHarmonicCos(2,deuteronPol/2-mS,p.CosTheta(),p.Phi());
  state*=Radial_p(2,p.Mag());
  // L=0, mL=0, S=1, mS=mD
  state += 
    ClebschGordan(0,0,2,deuteronPol,2,deuteronPol)
    *ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,deuteronPol)
    *SphericalHarmonicCos(0,0,p.CosTheta(),p.Phi())
    *Radial_p(0,p.Mag());

  return state;
}
//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronPState(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const FourVector<double>& p1) const
{
  // See TDeuteron::Wavefunction::DeuteronPState for more info.
  return DeuteronPState(deuteronPol,nucleon2Pol,nucleon1Pol,TVector3(p1[1],p1[2],p1[3]));    
}
//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronPState(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const TLorentzVector& p1) const
{
  // See TDeuteron::Wavefunction::DeuteronPState for more info.
  return DeuteronPState(deuteronPol,nucleon2Pol,nucleon1Pol,TVector3(p1.X(),p1.Y(),p1.Z()));    
}
//_____________________________________________________________________
// complex<double> TDeuteron::Wavefunction::DeuteronPState(const TStrangePolarization& pol,
// 						       const TSpinor::Polarization& pol1, 
// 						       const TVector3& p) const
// {
//   // Calculate the state of a deuteron in momentum space. 
//   // The polarization state of reaction is given by pol. pol1 fixes the polarization
//   // of the first nucleon (the one that interacts with the photon).
// 
//   return DeuteronPState(pol.Deuteron(),pol.Nucleon(),pol1,p);
// }
// 
// //_____________________________________________________________________
// complex<double> TDeuteron::Wavefunction::DeuteronPState(const TStrangePolarization& pol,
// 						       const TSpinor::Polarization& pol1, 
// 						       const FourVector<double>& p1) const
// {
//   // See TDeuteron::Wavefunction::DeuteronPState for more info.
//   return DeuteronPState(pol.Deuteron(),pol.Nucleon(),pol1,TVector3(p1[1],p1[2],p1[3]));  
// }
// //_____________________________________________________________________
// complex<double> TDeuteron::Wavefunction::DeuteronPState(const TStrangePolarization& pol,
// 						       const TSpinor::Polarization& pol1, 
// 						       const TLorentzVector& p1) const
// {
//   // See TDeuteron::Wavefunction::DeuteronPState for more info.
//   return DeuteronPState(pol.Deuteron(),pol.Nucleon(),pol1,TVector3(p1.X(),p1.Y(),p1.Z()));
// }

//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronPStateOff(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const TVector3& p) const
{
  // Calculate the state of a deuteron in momentum space. 
  // deuteronPol is the deuteron polarization (times two)
  // nucleon1pol fixes the polarization of the first nucleon (the one that interacts with the photon). (times two)
  //nucleon2pol fixes spectator nucleon polarization (times two)
  
  // The Dnp vertex matrix element is
  // BEGIN_LATEX
  // #frac{1}{2E_{1}} #bar{u}(p_{2},m_{2})u(#dot{p}_{1},m_{1})#Gamma_{#mu}#phi^{#mu}_{D}(m_{D})
  // = #sum_{L} u_{L}(|#vec{p}_{1}|) #sum_{m_{L},m_{S}} <Lm_{L}Sm_{S}|1m_{D}><#frac{1}{2}m_{1}#frac{1}{2} m_{2}|1m_{S}> Y_{Lm_{L}}(#hat{p}_{1})
  // END_LATEX
  // with BEGIN_LATEX #dot{p}_{1} = ( #sqrt{#vec{p}_{1}^{2}+M_{1}^{2}},#vec{p}_{1}) END_LATEX the on-shell 4-momentum of particle 1.
  // The labelling particle 1 and 2 is arbitrary in the non-relativistic case (You can do the math).
  if(IsRelativistic()) {
    cerr << "ERROR in TDeuteron::Wavefunction::DeuteronPState(TSpinor::Polarization&,"
	 << " TVector3&&): "
	 << "has no meaning for relativistic wave functions.\n";
    assert(1==0);
    exit(1);
  }
  if(abs(p.Z())<1.E-09) return 0.;
  
  complex<double> state=0.;

  // L=0, mL=0, S=1, mS=mD
  state += 
    ClebschGordan(0,0,2,deuteronPol,2,deuteronPol)
    *ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,deuteronPol)
    *SphericalHarmonicCos(0,0,p.CosTheta(),p.Phi())
    *GetUpoff(p);

  double Wpoff2=  GetWpoff2(p.Perp2());
  double Wpoff1=  GetWpoff1(p);  
  // L=2, mL=mD-mS, S=1, mS=m1+m2
  for(int mS=-1; mS<=1; ++mS)
    state -=
      ClebschGordan(4,deuteronPol-mS*2,2,mS*2,2,deuteronPol)
      *ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,mS*2)
      *(SphericalHarmonicCos(2,deuteronPol/2-mS,p.CosTheta(),p.Phi())
      *(Wpoff1+p.Perp2()/(p.Z()*p.Z())*Wpoff2)
      -SphericalHarmonicCos(2,deuteronPol/2-mS,0.,p.Phi())
      *p.Perp2()/(p.Z()*p.Z())*Wpoff2);
      
  //cout << GetUpoff(p)*p.Z()/HBARC*sqrt(1.E09) << " " << GetWpoff1(p)*p.Z()/HBARC*sqrt(1.E09)  << " " << GetWpoff1(p)*p.Z()/HBARC*sqrt(1.E09)
      
  return -I_UNIT*p.Z()/HBARC*state;
}
//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronPStateOff(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const FourVector<double>& p1) const
{
  // See TDeuteron::Wavefunction::DeuteronPState for more info.
  return DeuteronPStateOff(deuteronPol,nucleon2Pol,nucleon1Pol,TVector3(p1[1],p1[2],p1[3]));    
}
//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronPStateOff(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const TLorentzVector& p1) const
{
  // See TDeuteron::Wavefunction::DeuteronPState for more info.
  return DeuteronPState(deuteronPol,nucleon2Pol,nucleon1Pol,TVector3(p1.X(),p1.Y(),p1.Z()));    
}

//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronRState(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const TVector3& r) const
{
  // Calculate the state of a deuteron in coordinate space. 
  // deuteronPol is the deuteron polarization (times two)
  // nucleon1pol fixes the polarization of the first nucleon (the one that interacts with the photon). (times two)
  //nucleon2pol fixes spectator nucleon polarization (times two)
  
  if(IsRelativistic()) {
    cerr << "ERROR in TDeuteron::Wavefunction::DeuteronRState(TSpinor::Polarization&,"
	 << " TVector3&&): "
	 << "has no meaning for relativistic wave functions.\n";
    exit(1);
  }

  
  complex<double> state=0.;

  // L=0, mL=0, S=1, mS=mD
  state += 
    ClebschGordan(0,0,2,deuteronPol,2,deuteronPol)
    *ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,deuteronPol)
    *SphericalHarmonic(0,0,r.Theta(),r.Phi())
    *Radial_r(0,r.Mag());

  // L=2, mL=mD-mS, S=1, mS=m1+m2
//   for(int mS=-1; mS<=1; ++mS)
    state -=
      ClebschGordan(4,deuteronPol-nucleon1Pol-nucleon2Pol,2,nucleon1Pol+nucleon2Pol,2,deuteronPol)
      *ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,nucleon1Pol+nucleon2Pol)
      *SphericalHarmonic(2,(deuteronPol-nucleon1Pol-nucleon2Pol)/2,r.Theta(),r.Phi())
      *Radial_r(2,r.Mag());
  
  return state;
}
//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronRState(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const FourVector<double>& r1) const
{
  // See TDeuteron::Wavefunction::DeuteronRState for more info.
  return DeuteronRState(deuteronPol,nucleon2Pol,nucleon1Pol,TVector3(r1[1],r1[2],r1[3]));    
}
//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronRState(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const TLorentzVector& r1) const
{
  // See TDeuteron::Wavefunction::DeuteronRState for more info.
  return DeuteronRState(deuteronPol,nucleon2Pol,nucleon1Pol,TVector3(r1.X(),r1.Y(),r1.Z()));    
}
//_____________________________________________________________________
/*complex<double> TDeuteron::Wavefunction::DeuteronRState(const TStrangePolarization& pol,
						       const TSpinor::Polarization& pol1, 
						       const TVector3& r) const
{
  // Calculate the state of a deuteron in momentum space. 
  // The polarization state of reaction is given by pol. pol1 fixes the polarization
  // of the first nucleon (the one that interacts with the photon).

  return DeuteronRState(pol.Deuteron(),pol.Nucleon(),pol1,r);
}

//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronRState(const TStrangePolarization& pol,
						       const TSpinor::Polarization& pol1, 
						       const FourVector<double>& r1) const
{
  // See TDeuteron::Wavefunction::DeuteronRState for more info.
  return DeuteronRState(pol.Deuteron(),pol.Nucleon(),pol1,TVector3(r1[1],r1[2],r1[3]));  
}
//_____________________________________________________________________
complex<double> TDeuteron::Wavefunction::DeuteronRState(const TStrangePolarization& pol,
						       const TSpinor::Polarization& pol1, 
						       const TLorentzVector& r1) const
{
  // See TDeuteron::Wavefunction::DeuteronRState for more info.
  return DeuteronRState(pol.Deuteron(),pol.Nucleon(),pol1,TVector3(r1.X(),r1.Y(),r1.Z()));
}*/



double TDeuteron::Wavefunction::getResidu() const{
  
  if(yukawa) return dynamic_cast<TYukawaPWF *>(fImplementation)->getResidu();
  else {
    cerr << "not a valid wf form for pole residu extraction" << std::endl;
    assert(1==0);
  }
}



///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TDeuteron::Polarization                                               //
//                                                                       //
// Polarization state of a deuteron                                      //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TDeuteron::Polarization)

//______________________________________________________________________________
TDeuteron::Polarization::Polarization(double theta, double phi, int polState)
  : fTheta(theta), fPhi(phi), fPolarizationState(polState)
{
  // Constructor
  SetTheta(theta);
  SetPhi(phi);
}

//______________________________________________________________________________
TDeuteron::Polarization::Polarization(TRootIOCtor *rio)
  : TObject(), fTheta(0.), fPhi(0.), fPolarizationState(-5)
{
  // ROOT I/O Constructor
}

//______________________________________________________________________________
TDeuteron::Polarization::Polarization(const Polarization& rhs)
  : fTheta(rhs.fTheta), fPhi(rhs.fPhi), fPolarizationState(rhs.fPolarizationState)
{
  // Copy constructor
}

//______________________________________________________________________________
TDeuteron::Polarization::~Polarization()
{
  // Destructor
}

//______________________________________________________________________________
TDeuteron::Polarization& TDeuteron::Polarization::operator=(const Polarization& rhs)
{
  // Assignment

  if( this!=&rhs) { // avoid self-assignment
    TObject::operator=(rhs);
    fTheta = rhs.fTheta;
    fPhi = rhs.fPhi;
    fPolarizationState = rhs.fPolarizationState;
  }

  return *this;
}

//______________________________________________________________________________
bool TDeuteron::Polarization::operator==(const Polarization& rhs)
{
  // Comparison

  return 
    fabs(fTheta-rhs.fTheta)<STRANGEUFLOW &&
    fabs(fPhi-rhs.fPhi)<STRANGEUFLOW &&
    fPolarizationState==rhs.fPolarizationState;
}

//______________________________________________________________________________
void TDeuteron::Polarization::SetTheta(double theta)
{
  // Set the polar angle of the spin quantization axis (in rad)

  if( theta<0. || theta>TMath::Pi() ) {
    cerr << "ERROR in TDeuteron::Polarization::SetTheta(double): "
	 << "polar angle not in range [0,pi].\n";
    exit(1);
  }

  fTheta = theta;
}

//______________________________________________________________________________
void TDeuteron::Polarization::SetPhi(double phi)
{
  // Set the azimuthal angle of the spin quantization axis (in rad)
  fPhi = phi;
}

//______________________________________________________________________________
void TDeuteron::Polarization::SetState(int polState)
{
  // Set the polarization state = -1, 0 or +1
  if( polState>1 || polState<-1 ) {
    cerr << "ERROR in TDeuteron::Polarization::SetState(int): "
	 << "polarization state should be -1, 0 or +1.\n";
    exit(1);
  }

  fPolarizationState = polState;
}

//______________________________________________________________________________
TDeuteron::Polarization& TDeuteron::Polarization::operator++()
{
  // Move the polarization state one-up like this ++spinorPolarization
  if( fPolarizationState!=1 ) ++fPolarizationState;

  return *this;
}

//______________________________________________________________________________
TDeuteron::Polarization& TDeuteron::Polarization::operator--()
{
  // Move the polarization state one-up like this ++spinorPolarization
  if( fPolarizationState!=-1 ) --fPolarizationState;

  return *this;
}

//______________________________________________________________________________
TDeuteron::Polarization::operator int() const
{
  // Cast polarization state to int: 
  // returns magnetic quantum number multiplied by two
  return fPolarizationState * 2;
}

//______________________________________________________________________________
FourVector<complex<double> > TDeuteron::Polarization::GetPolarizationVector() const
{
  // Return the polarization 4-vector for the current polarization state
  // in the rest frame of the deuteron.
  
  // Pick the correct polarization vector
  if( fPolarizationState==0 ) 
    return FourVector<complex<double> >(0.,-sin(fTheta),0.,cos(fTheta));
  
  else if( fPolarizationState==1 ) 
    return FourVector<complex<double> >(0.,
 complex<double>(-cos(fPhi)*cos(fTheta)/sqrt(2.),
		 sin(fPhi)*cos(fTheta)/sqrt(2.)),
 complex<double>(-sin(fPhi)/sqrt(2.),
		 -cos(fPhi)/sqrt(2.)),
 complex<double>(-cos(fPhi)*sin(fTheta)/sqrt(2.),
		 sin(fPhi)*sin(fTheta)/sqrt(2.)));
  
  else if( fPolarizationState==-1 ) 
    return FourVector<complex<double> >(0.,
 complex<double>(cos(fPhi)*cos(fTheta)/sqrt(2.),
		 sin(fPhi)*cos(fTheta)/sqrt(2.)),
 complex<double>(sin(fPhi)/sqrt(2.),
		 -cos(fPhi)/sqrt(2.)),
 complex<double>(cos(fPhi)*sin(fTheta)/sqrt(2.),
		 sin(fPhi)*sin(fTheta)/sqrt(2.)));
  
  else {
    cerr << "ERROR in TDeuteron::Polarization::GetPolarizationVector(): "
	 << "unknown polarization state.\n";
    assert(1==0);
  }

  
  return FourVector<complex<double> >();
}

/*//______________________________________________________________________________
TString TDeuteron::Polarization::HashName() const
{
  // Returns a unique string that specifies the current state of
  // the deuteron's polarization.

  // the polarization is fixed by the angles of the quantization axis
  // and the state.
  double key[3] = { fTheta, fPhi, (double)fPolarizationState };
  
  // We create a hash key using the function in the cachetools folder
  char hash[sizeof(double)*3*2 + 1];
  doublestochar(3, key, hash);

  return hash;
}
*/
//_____________________________________________________________________
TDeuteron::Wavefunction* TDeuteron::Wavefunction::CreateNijmegenWavefunction()
{
  // Non-relativistic deuteron wavefunction generated with the Nijm93 potential.
  // 
  // This potential is described in PRC49(1994)2950.
  // The numerical values of the wave function have been downloaded from
  // http://nn-online.org.

  Wavefunction *wf = new Wavefunction("Nijm93");

  wf->fRelativistic = false;
  wf->yukawa = false;
  TInterpolatingWavefunction *implementation = new TInterpolatingWavefunction;

  static double r[839] =
    { .00000000,.00100000,.00199936,.00299872,.00399808,.00499744,.00599680,
      .00699616,.00799552,.00899488,.00999424,.01099360,.01199296,.01299232,
      .01399168,.01499105,.01599041,.01698977,.01798913,.01898849,.01998785,
      .02198657,.02398529,.02598401,.02798273,.02998145,.03198017,.03397889,
      .03597761,.03797633,.03997505,
      .04197378,.04397250,.04597122,.04796994,.04996866,.05196738,.05396610,
      .05596482,.05796354,.05996226,.06196098,.06395970,.06595842,.06795715,
      .06995587,.07195459,.07395331,.07595203,.07795075,.07994947,.08194819,
      .08394691,.08594563,.08794435,.08994307,.09194179,.09394052,.09593924,
      .09793796,.09993668,.10393412,.10793156,.11192900,.11592644,.11992389,
      .12392133,.12791877,.13191621,.13591365,.13991109,.14390853,.14790598,
      .15190342,.15590086,.15989830,.16389574,.16789318,.17189063,.17588807,
      .17988551,.18388295,.18788039,.19187783,.19587527,.19987272,.20387016,
      .20786760,.21186504,.21586248,.21985992,.22385737,.22785481,.23185225,
      .23584969,.23984713,.24384457,.24784201,.25183946,.25583690,.25983434,
      .26383178,.26782922,.27182666,.27582411,.27982155,.28381899,.28781643,
      .29181387,.29581131,.29980875,.30780364,.31579852,.32379340,.33178829,
      .33978317,.34777805,.35577294,.36376782,.37176270,.37975759,.38775247,
      .39574735,.40374223,.41173712,.41973200,.42772688,.43572177,.44371665,
      .45171153,.45970642,.46770130,.47569618,.48369107,.49168595,.49968083,
      .50767571,.51567060,.52366548,.53166036,.53965525,.54765013,.55564501,
      .56363990,.57163478,.57962966,.58762455,.59561943,.60361431,.61160919,
      .61960408,.62759896,.63559384,.64358873,.65158361,.65957849,.66757338,
      .67556826,.68356314,.69155803,.69955291,.70754779,.71554267,.72353756,
      .73153244,.73952732,.74752221,.75551709,.76351197,.77150686,.77950174,
      .78749662,.79549151,.80348639,.81148127,.81947615,.82747104,.83546592,
      .84346080,.85145569,.85945057,.86744545,.87544034,.88343522,.89143010,
      .89942499,.90741987,.91541475,.92340963,.93140452,.93939940,.94739428,
      .95538917,.96338405,.97137893,.97937382,.98736870,.99536358,1.00735591,
      1.02334567,1.03933544,1.05532520,1.07131497,1.08730474,1.10329450,
      1.11928427,1.13527404,1.15126380,1.16725357,1.18324333,1.19923310,
      1.21522287,1.23121263,1.24720240,1.26319216,1.27918193,1.29517170,
      1.31116146,1.32715123,1.34314100,1.35913076,1.37512053,1.39111029,
      1.40710006,1.42308983,1.43907959,1.45506936,1.47105912,1.48704889,
      1.50303866,1.51902842,1.53501819,1.55100796,1.56699772,1.58298749,
      1.59897725,1.61496702,1.63095679,1.64694655,1.66293632,1.67892608,
      1.69491585,1.71090562,1.72689538,1.74288515,1.75887492,1.77486468,
      1.79085445,1.80684421,1.82283398,1.83882375,1.85481351,1.87080328,
      1.88679304,1.90278281,1.91877258,1.93476234,1.95075211,1.96674187,
      1.98273164,1.99872141,2.01471117,2.03070094,2.04669071,2.06268047,
      2.07867024,2.09466000,2.11064977,2.12663954,2.14262930,2.15861907,
      2.17460883,2.19059860,2.20658837,2.22257813,2.23856790,2.25455767,
      2.27054743,2.28653720,2.30252696,2.31851673,2.33450650,2.35049626,
      2.36648603,2.38247579,2.39846556,2.41445533,2.43044509,2.44643486,
      2.46242463,2.47841439,2.49440416,2.51039392,2.52638369,2.54237346,
      2.55836322,2.57435299,2.59034275,2.60633252,2.62232229,2.63831205,
      2.65430182,2.67029159,2.68628135,2.70227112,2.71826088,2.73425065,
      2.75024042,2.76623018,2.78221995,2.79820971,2.81419948,2.83018925,
      2.84617901,2.86216878,2.87815854,2.89414831,2.91013808,2.92612784,
      2.94211761,2.95810738,2.97409714,2.99008691,3.00607667,3.02206644,
      3.03805621,3.05404597,3.07003574,3.08602550,3.10201527,3.11800504,
      3.13399480,3.14998457,3.16597434,3.18196410,3.19795387,3.21394363,
      3.22993340,3.24592317,3.26191293,3.27790270,3.29389246,3.30988223,
      3.32587200,3.34186176,3.35785153,3.37384130,3.38983106,3.40582083,
      3.42181059,3.43780036,3.45379013,3.46977989,3.48576966,3.50175942,
      3.51774919,3.53373896,3.54972872,3.56571849,3.58170826,3.59769802,
      3.61368779,3.62967755,3.64566732,3.66165709,3.67764685,3.69363662,
      3.70962638,3.72561615,3.74160592,3.75759568,3.77358545,3.78957521,
      3.80556498,3.82155475,3.83754451,3.85353428,3.86952405,3.88551381,
      3.90150358,3.91749334,3.93348311,3.94947288,3.96546264,3.98145241,
      3.99744217,4.01343194,4.02942171,4.04541147,4.06140124,4.07739101,
      4.09338077,4.10937054,4.12536030,4.14135007,4.15733984,4.17332960,
      4.18931937,4.20530913,4.22129890,4.23728867,4.25327843,4.26926820,
      4.28525797,4.30124773,4.31723750,4.33322726,4.34921703,4.36520680,
      4.38119656,4.39718633,4.41317609,4.42916586,4.44515563,4.46114539,
      4.47713516,4.49312493,4.50911469,4.52510446,4.54109422,4.55708399,
      4.57307376,4.58906352,4.60505329,4.62104305,4.63703282,4.65302259,
      4.66901235,4.68500212,4.70099188,4.71698165,4.73297142,4.74896118,
      4.76495095,4.78094072,4.79693048,4.81292025,4.82891001,4.84489978,
      4.86088955,4.87687931,4.89286908,4.90885884,4.92484861,4.94083838,
      4.95682814,4.97281791,4.98880768,5.01279232,5.04477186,5.07675139,
      5.10873092,5.14071045,5.17268999,5.20466952,5.23664905,5.26862858,
      5.30060812,5.33258765,5.36456718,5.39654671,5.42852624,5.46050578,
      5.49248531,5.52446484,5.55644437,5.58842391,5.62040344,5.65238297,
      5.68436250,5.71634203,5.74832157,5.78030110,5.81228063,5.84426016,
      5.87623970,5.90821923,5.94019876,5.97217829,6.00415783,6.03613736,
      6.06811689,6.10009642,6.13207595,6.16405549,6.19603502,6.22801455,
      6.25999408,6.29197362,6.32395315,6.35593268,6.38791221,6.41989174,
      6.45187128,6.48385081,6.51583034,6.54780987,6.57978941,6.61176894,
      6.64374847,6.67572800,6.70770754,6.73968707,6.77166660,6.80364613,
      6.83562566,6.86760520,6.89958473,6.93156426,6.96354379,6.99552333,
      7.02750286,7.05948239,7.09146192,7.12344146,7.15542099,7.18740052,
      7.21938005,7.25135958,7.28333912,7.31531865,7.34729818,7.37927771,
      7.41125725,7.44323678,7.47521631,7.50719584,7.53917537,7.57115491,
      7.60313444,7.63511397,7.66709350,7.69907304,7.73105257,7.76303210,
      7.79501163,7.82699117,7.85897070,7.89095023,7.92292976,7.95490929,
      7.98688883,8.01886836,8.05084789,8.08282742,8.11480696,8.14678649,
      8.17876602,8.21074555,8.24272508,8.27470462,8.30668415,8.33866368,
      8.37064321,8.40262275,8.43460228,8.46658181,8.49856134,8.53054088,
      8.56252041,8.59449994,8.62647947,8.65845900,8.69043854,8.72241807,
      8.75439760,8.78637713,8.81835667,8.85033620,8.88231573,8.91429526,
      8.94627480,8.97825433,9.01023386,9.04221339,9.07419292,9.10617246,
      9.13815199,9.17013152,9.20211105,9.23409059,9.26607012,9.29804965,
      9.33002918,9.36200871,9.39398825,9.42596778,9.45794731,9.48992684,
      9.52190638,9.55388591,9.58586544,9.61784497,9.64982451,9.68180404,
      9.71378357,9.74576310,9.77774263,9.80972217,9.84170170,9.87368123,
      9.90566076,9.93764030,9.96961983,10.00159936,10.06555842,10.12951749,
      10.19347655,10.25743562,10.32139468,10.38535375,10.44931281,10.51327188,
      10.57723094,10.64119001,10.70514907,10.76910814,10.83306720,10.89702626,
      10.96098533,11.02494439,11.08890346,11.15286252,11.21682159,11.28078065,
      11.34473972,11.40869878,11.47265785,11.53661691,11.60057597,11.66453504,
      11.72849410,11.79245317,11.85641223,11.92037130,11.98433036,12.04828943,
      12.11224849,12.17620756,12.24016662,12.30412568,12.36808475,12.43204381,
      12.49600288,12.55996194,12.62392101,12.68788007,12.75183914,12.81579820,
      12.87975727,12.94371633,13.00767539,13.07163446,13.13559352,13.19955259,
      13.26351165,13.32747072,13.39142978,13.45538885,13.51934791,13.58330698,
      13.64726604,13.71122510,13.77518417,13.83914323,13.90310230,13.96706136,
      14.03102043,14.09497949,14.15893856,14.22289762,14.28685669,14.35081575,
      14.41477482,14.47873388,14.54269294,14.60665201,14.67061107,14.73457014,
      14.79852920,14.86248827,14.92644733,14.99040640,15.05436546,15.11832453,
      15.18228359,15.24624265,15.31020172,15.37416078,15.43811985,15.50207891,
      15.56603798,15.62999704,15.69395611,15.75791517,15.82187424,15.88583330,
      15.94979236,16.01375143,16.07771049,16.14166956,16.20562862,16.26958769,
      16.33354675,16.39750582,16.46146488,16.52542395,16.58938301,16.65334207,
      16.71730114,16.78126020,16.84521927,16.90917833,16.97313740,17.03709646,
      17.10105553,17.16501459,17.22897366,17.29293272,17.35689178,17.42085085,
      17.48480991,17.54876898,17.61272804,17.67668711,17.74064617,17.80460524,
      17.86856430,17.93252337,17.99648243,18.06044150,18.12440056,18.18835962,
      18.25231869,18.31627775,18.38023682,18.44419588,18.50815495,18.57211401,
      18.63607308,18.70003214,18.76399121,18.82795027,18.89190933,18.95586840,
      19.01982746,19.08378653,19.14774559,19.21170466,19.27566372,19.33962279,
      19.40358185,19.46754092,19.53149998,19.59545904,19.65941811,19.72337717,
      19.78733624,19.85129530,19.91525437,19.97921343,20.04317250,20.10713156,
      20.17109063,20.23504969,20.29900875,20.36296782,20.42692688,20.49088595,
      20.55484501,20.61880408,20.68276314,20.74672221,20.81068127,20.87464034,
      20.93859940,21.00255846,21.06651753,21.13047659,21.19443566,21.25839472,
      21.32235379,21.38631285,21.45027192,21.51423098,21.57819005,21.64214911,
      21.70610818,21.77006724,21.83402630,21.89798537,21.96194443,22.02590350,
      22.08986256,22.15382163,22.21778069,22.28173976,22.34569882,22.40965789,
      22.47361695,22.53757601,22.60153508,22.66549414,22.72945321,22.79341227,
      22.85737134,22.92133040,22.98528947,23.04924853,23.11320760,23.17716666,
      23.24112572,23.30508479,23.36904385,23.43300292,23.49696198,23.56092105,
      23.62488011,23.68883918,23.75279824,23.81675731,23.88071637,23.94467543,
      24.00863450,24.07259356,24.13655263,24.20051169,24.26447076,24.32842982,
      24.39238889,24.45634795,24.52030702,24.58426608,24.64822514,24.71218421,
      24.77614327,24.84010234,24.90406140,24.96802047 };

  static double ur[839] =
    { .000000000e+00,.157907255e-03,.315715669e-03,.473530214e-03,
      .631354727e-03,.789194063e-03,.947052302e-03,.110493311e-02,
      .126284026e-02,.142077751e-02,.157874861e-02,.173675734e-02,
      .189480745e-02,.205290271e-02,.221104689e-02,.236924373e-02,
      .252749701e-02,.268581050e-02,.284418794e-02,.300263312e-02,
      .316114979e-02,.347841268e-02,.379600672e-02,.411396205e-02,
      .443230880e-02,.475107711e-02,.507029712e-02,.538999901e-02,
      .571021291e-02,.603096902e-02,.635229751e-02,.667422856e-02,
      .699679237e-02,.732001916e-02,.764393913e-02,.796858251e-02,
      .829397954e-02,.862016047e-02,.894715555e-02,.927499505e-02,
      .960370925e-02,.993332843e-02,.102638829e-01,.105954029e-01,
      .109279189e-01,.112614611e-01,.115960599e-01,.119317456e-01,
      .122685487e-01,.126064993e-01,.129456281e-01,.132859652e-01,
      .136275412e-01,.139703864e-01,.143145312e-01,.146600061e-01,
      .150068415e-01,.153550678e-01,.157047154e-01,.160558148e-01,
      .164083965e-01,.171181282e-01,.178341543e-01,.185567182e-01,
      .192860635e-01,.200224337e-01,.207660721e-01,.215172220e-01,
      .222761262e-01,.230430276e-01,.238181685e-01,.246017909e-01,
      .253941364e-01,.261954462e-01,.270059606e-01,.278259196e-01,
      .286555621e-01,.294951265e-01,.303448502e-01,.312049695e-01,
      .320757197e-01,.329573349e-01,.338500479e-01,.347540902e-01,
      .356696915e-01,.365970803e-01,.375364831e-01,.384881246e-01,
      .394522276e-01,.404290127e-01,.414186984e-01,.424215008e-01,
      .434376335e-01,.444673075e-01,.455107310e-01,.465681095e-01,
      .476396452e-01,.487255373e-01,.498259817e-01,.509411708e-01,
      .520712933e-01,.532165342e-01,.543770748e-01,.555530920e-01,
      .567447587e-01,.579522435e-01,.591757104e-01,.604153187e-01,
      .616712230e-01,.629435729e-01,.642325130e-01,.668607155e-01,
      .695568791e-01,.723219621e-01,.751568273e-01,.780622381e-01,
      .810388544e-01,.840872285e-01,.872078018e-01,.904009006e-01,
      .936667338e-01,.970053895e-01,.100416833e+00,.103900903e+00,
      .107457313e+00,.111085645e+00,.114785353e+00,.118555761e+00,
      .122396059e+00,.126305308e+00,.130282439e+00,.134326252e+00,
      .138435418e+00,.142608482e+00,.146843864e+00,.151139861e+00,
      .155494650e+00,.159906293e+00,.164372737e+00,.168891820e+00,
      .173461275e+00,.178078736e+00,.182741739e+00,.187447730e+00,
      .192194071e+00,.196978042e+00,.201796851e+00,.206647637e+00,
      .211527479e+00,.216433398e+00,.221362369e+00,.226311324e+00,
      .231277161e+00,.236256749e+00,.241246935e+00,.246244553e+00,
      .251246430e+00,.256249391e+00,.261250270e+00,.266245913e+00,
      .271233188e+00,.276208987e+00,.281170239e+00,.286113913e+00,
      .291037021e+00,.295936631e+00,.300809869e+00,.305653923e+00,
      .310466051e+00,.315243585e+00,.319983938e+00,.324684603e+00,
      .329343164e+00,.333957293e+00,.338524761e+00,.343043432e+00,
      .347511275e+00,.351926359e+00,.356286860e+00,.360591059e+00,
      .364837345e+00,.369024217e+00,.373150281e+00,.377214256e+00,
      .381214966e+00,.385151346e+00,.389022438e+00,.392827393e+00,
      .396565463e+00,.400236008e+00,.403838487e+00,.407372457e+00,
      .410837577e+00,.414233594e+00,.417560351e+00,.420817777e+00,
      .424005887e+00,.427124777e+00,.431673731e+00,.437498990e+00,
      .443052855e+00,.448339428e+00,.453363689e+00,.458131380e+00,
      .462648891e+00,.466923152e+00,.470961528e+00,.474771717e+00,
      .478361667e+00,.481739486e+00,.484913367e+00,.487891524e+00,
      .490682126e+00,.493293251e+00,.495732835e+00,.498008637e+00,
      .500128206e+00,.502098857e+00,.503927649e+00,.505621370e+00,
      .507186527e+00,.508629339e+00,.509955734e+00,.511171349e+00,
      .512281534e+00,.513291351e+00,.514205591e+00,.515028770e+00,
      .515765148e+00,.516418735e+00,.516993301e+00,.517492390e+00,
      .517919330e+00,.518277245e+00,.518569069e+00,.518797554e+00,
      .518965284e+00,.519074684e+00,.519128034e+00,.519127477e+00,
      .519075030e+00,.518972591e+00,.518821952e+00,.518624804e+00,
      .518382747e+00,.518097297e+00,.517769892e+00,.517401900e+00,
      .516994622e+00,.516549301e+00,.516067127e+00,.515549238e+00,
      .514996727e+00,.514410647e+00,.513792011e+00,.513141798e+00,
      .512460956e+00,.511750401e+00,.511011024e+00,.510243691e+00,
      .509449243e+00,.508628501e+00,.507782266e+00,.506911319e+00,
      .506016423e+00,.505098325e+00,.504157756e+00,.503195431e+00,
      .502212052e+00,.501208303e+00,.500184858e+00,.499142376e+00,
      .498081503e+00,.497002873e+00,.495907106e+00,.494794810e+00,
      .493666582e+00,.492523005e+00,.491364650e+00,.490192078e+00,
      .489005837e+00,.487806462e+00,.486594478e+00,.485370397e+00,
      .484134720e+00,.482887937e+00,.481630525e+00,.480362951e+00,
      .479085670e+00,.477799126e+00,.476503752e+00,.475199968e+00,
      .473888186e+00,.472568806e+00,.471242216e+00,.469908795e+00,
      .468568911e+00,.467222922e+00,.465871175e+00,.464514008e+00,
      .463151748e+00,.461784714e+00,.460413214e+00,.459037548e+00,
      .457658005e+00,.456274867e+00,.454888407e+00,.453498887e+00,
      .452106564e+00,.450711685e+00,.449314490e+00,.447915209e+00,
      .446514067e+00,.445111280e+00,.443707059e+00,.442301604e+00,
      .440895112e+00,.439487771e+00,.438079764e+00,.436671268e+00,
      .435262451e+00,.433853478e+00,.432444509e+00,.431035694e+00,
      .429627183e+00,.428219117e+00,.426811633e+00,.425404864e+00,
      .423998937e+00,.422593975e+00,.421190096e+00,.419787416e+00,
      .418386042e+00,.416986082e+00,.415587638e+00,.414190807e+00,
      .412795685e+00,.411402363e+00,.410010928e+00,.408621464e+00,
      .407234053e+00,.405848773e+00,.404465699e+00,.403084903e+00,
      .401706456e+00,.400330423e+00,.398956869e+00,.397585856e+00,
      .396217443e+00,.394851687e+00,.393488643e+00,.392128364e+00,
      .390770899e+00,.389416298e+00,.388064607e+00,.386715870e+00,
      .385370131e+00,.384027430e+00,.382687807e+00,.381351299e+00,
      .380017942e+00,.378687772e+00,.377360820e+00,.376037119e+00,
      .374716699e+00,.373399589e+00,.372085816e+00,.370775407e+00,
      .369468386e+00,.368164779e+00,.366864607e+00,.365567893e+00,
      .364274656e+00,.362984918e+00,.361698696e+00,.360416009e+00,
      .359136873e+00,.357861304e+00,.356589316e+00,.355320926e+00,
      .354056145e+00,.352794986e+00,.351537462e+00,.350283584e+00,
      .349033361e+00,.347786805e+00,.346543923e+00,.345304725e+00,
      .344069217e+00,.342837409e+00,.341609306e+00,.340384914e+00,
      .339164239e+00,.337947286e+00,.336734059e+00,.335524564e+00,
      .334318802e+00,.333116778e+00,.331918494e+00,.330723952e+00,
      .329533154e+00,.328346102e+00,.327162796e+00,.325983238e+00,
      .324807427e+00,.323635363e+00,.322467047e+00,.321302477e+00,
      .320141652e+00,.318984571e+00,.317831232e+00,.316681634e+00,
      .315535773e+00,.314393648e+00,.313255255e+00,.312120592e+00,
      .310989655e+00,.309862441e+00,.308738946e+00,.307619166e+00,
      .306503097e+00,.305390734e+00,.304282073e+00,.303177108e+00,
      .302075835e+00,.300978249e+00,.299884345e+00,.298794115e+00,
      .297707556e+00,.296624661e+00,.295545424e+00,.294469838e+00,
      .293397898e+00,.292329597e+00,.291264928e+00,.290203885e+00,
      .289146461e+00,.288092648e+00,.287042440e+00,.285995829e+00,
      .284952808e+00,.283913369e+00,.282877506e+00,.281845209e+00,
      .280816472e+00,.279791287e+00,.278769646e+00,.277751540e+00,
      .276230993e+00,.274215895e+00,.272214792e+00,.270227614e+00,
      .268254294e+00,.266294761e+00,.264348945e+00,.262416775e+00,
      .260498177e+00,.258593078e+00,.256701407e+00,.254823087e+00,
      .252958044e+00,.251106204e+00,.249267491e+00,.247441829e+00,
      .245629142e+00,.243829353e+00,.242042385e+00,.240268163e+00,
      .238506608e+00,.236757644e+00,.235021193e+00,.233297178e+00,
      .231585522e+00,.229886147e+00,.228198976e+00,.226523932e+00,
      .224860937e+00,.223209914e+00,.221570786e+00,.219943476e+00,
      .218327908e+00,.216724004e+00,.215131688e+00,.213550885e+00,
      .211981516e+00,.210423508e+00,.208876784e+00,.207341268e+00,
      .205816887e+00,.204303563e+00,.202801224e+00,.201309795e+00,
      .199829201e+00,.198359369e+00,.196900225e+00,.195451697e+00,
      .194013711e+00,.192586195e+00,.191169077e+00,.189762285e+00,
      .188365747e+00,.186979392e+00,.185603150e+00,.184236951e+00,
      .182880723e+00,.181534398e+00,.180197905e+00,.178871177e+00,
      .177554144e+00,.176246739e+00,.174948893e+00,.173660539e+00,
      .172381611e+00,.171112040e+00,.169851762e+00,.168600709e+00,
      .167358818e+00,.166126021e+00,.164902256e+00,.163687457e+00,
      .162481560e+00,.161284503e+00,.160096221e+00,.158916652e+00,
      .157745733e+00,.156583403e+00,.155429600e+00,.154284263e+00,
      .153147330e+00,.152018742e+00,.150898439e+00,.149786360e+00,
      .148682447e+00,.147586640e+00,.146498881e+00,.145419113e+00,
      .144347276e+00,.143283315e+00,.142227172e+00,.141178789e+00,
      .140138112e+00,.139105085e+00,.138079651e+00,.137061756e+00,
      .136051345e+00,.135048363e+00,.134052758e+00,.133064474e+00,
      .132083460e+00,.131109661e+00,.130143026e+00,.129183503e+00,
      .128231039e+00,.127285583e+00,.126347084e+00,.125415492e+00,
      .124490756e+00,.123572826e+00,.122661652e+00,.121757185e+00,
      .120859376e+00,.119968177e+00,.119083539e+00,.118205414e+00,
      .117333755e+00,.116468514e+00,.115609644e+00,.114757099e+00,
      .113910833e+00,.113070800e+00,.112236953e+00,.111409248e+00,
      .110587639e+00,.109772083e+00,.108962534e+00,.108158949e+00,
      .107361284e+00,.106569495e+00,.105783540e+00,.105003376e+00,
      .104228959e+00,.103460249e+00,.102697203e+00,.101939780e+00,
      .101187938e+00,.100441637e+00,.997008350e-01,.989654926e-01,
      .982355695e-01,.975110259e-01,.967918223e-01,.960779195e-01,
      .953692784e-01,.946658606e-01,.939676275e-01,.932745411e-01,
      .925865635e-01,.919036573e-01,.912257850e-01,.905529099e-01,
      .898849949e-01,.892220038e-01,.885639003e-01,.879106485e-01,
      .872622126e-01,.859796474e-01,.847159244e-01,.834707672e-01,
      .822439038e-01,.810350659e-01,.798439890e-01,.786704126e-01,
      .775140801e-01,.763747384e-01,.752521383e-01,.741460340e-01,
      .730561835e-01,.719823483e-01,.709242932e-01,.698817866e-01,
      .688546002e-01,.678425092e-01,.668452919e-01,.658627300e-01,
      .648946081e-01,.639407144e-01,.630008397e-01,.620747783e-01,
      .611623272e-01,.602632866e-01,.593774595e-01,.585046518e-01,
      .576446722e-01,.567973323e-01,.559624465e-01,.551398317e-01,
      .543293076e-01,.535306968e-01,.527438240e-01,.519685169e-01,
      .512046055e-01,.504519223e-01,.497103025e-01,.489795834e-01,
      .482596048e-01,.475502090e-01,.468512404e-01,.461625458e-01,
      .454839742e-01,.448153768e-01,.441566071e-01,.435075206e-01,
      .428679751e-01,.422378302e-01,.416169480e-01,.410051921e-01,
      .404024286e-01,.398085252e-01,.392233517e-01,.386467799e-01,
      .380786832e-01,.375189372e-01,.369674191e-01,.364240079e-01,
      .358885846e-01,.353610316e-01,.348412334e-01,.343290760e-01,
      .338244469e-01,.333272357e-01,.328373332e-01,.323546321e-01,
      .318790264e-01,.314104120e-01,.309486859e-01,.304937471e-01,
      .300454956e-01,.296038333e-01,.291686633e-01,.287398900e-01,
      .283174196e-01,.279011594e-01,.274910180e-01,.270869056e-01,
      .266887335e-01,.262964144e-01,.259098622e-01,.255289923e-01,
      .251537211e-01,.247839662e-01,.244196466e-01,.240606824e-01,
      .237069949e-01,.233585065e-01,.230151407e-01,.226768224e-01,
      .223434773e-01,.220150322e-01,.216914152e-01,.213725553e-01,
      .210583826e-01,.207488281e-01,.204438241e-01,.201433035e-01,
      .198472005e-01,.195554501e-01,.192679884e-01,.189847523e-01,
      .187056798e-01,.184307095e-01,.181597813e-01,.178928356e-01,
      .176298140e-01,.173706588e-01,.171153131e-01,.168637209e-01,
      .166158270e-01,.163715772e-01,.161309178e-01,.158937960e-01,
      .156601599e-01,.154299582e-01,.152031404e-01,.149796568e-01,
      .147594584e-01,.145424968e-01,.143287245e-01,.141180947e-01,
      .139105610e-01,.137060781e-01,.135046010e-01,.133060856e-01,
      .131104884e-01,.129177664e-01,.127278773e-01,.125407796e-01,
      .123564323e-01,.121747947e-01,.119958273e-01,.118194906e-01,
      .116457460e-01,.114745555e-01,.113058814e-01,.111396868e-01,
      .109759352e-01,.108145908e-01,.106556181e-01,.104989822e-01,
      .103446489e-01,.101925843e-01,.100427550e-01,.989512811e-02,
      .974967135e-02,.960635278e-02,.946514096e-02,.932600494e-02,
      .918891419e-02,.905383865e-02,.892074870e-02,.878961515e-02,
      .866040925e-02,.853310265e-02,.840766743e-02,.828407609e-02,
      .816230152e-02,.804231702e-02,.792409627e-02,.780761335e-02,
      .769284271e-02,.757975917e-02,.746833795e-02,.735855459e-02,
      .725038504e-02,.714380556e-02,.703879278e-02,.693532367e-02,
      .683337555e-02,.673292604e-02,.663395312e-02,.653643509e-02,
      .644035056e-02,.634567845e-02,.625239801e-02,.616048878e-02,
      .606993059e-02,.598070360e-02,.589278823e-02,.580616520e-02,
      .572081552e-02,563672046e-02,.555386159e-02,.547222072e-02,
      .539177997e-02,.531252168e-02,.523442847e-02,.515748323e-02,
      .508166906e-02,.500696935e-02,.493336771e-02,.486084801e-02,
      .478939433e-02,.471899101e-02,.464962261e-02,.458127391e-02,
      .451392993e-02,.444757590e-02,.438219725e-02,.431777966e-02,
      .425430900e-02,.419177135e-02,.413015299e-02,.406944042e-02,
      .400962030e-02,.395067953e-02,.389260519e-02,.383538452e-02,
      .377900499e-02,.372345423e-02,.366872006e-02,.361479047e-02,
      .356165364e-02,.350929791e-02,.345771180e-02,.340688400e-02,
      .335680336e-02,.330745889e-02,.325883978e-02,.321093537e-02,
      .316373514e-02,.311722875e-02,.307140599e-02,.302625682e-02,
      .298177134e-02,.293793978e-02,.289475254e-02,.285220015e-02,
      .281027327e-02,.276896270e-02,.272825940e-02 };

  static double wr[839] = 
    { .000000000e+00,.458198607e-09,.366205604e-08,.123554909e-07,
      .292824853e-07,.571867520e-07,.988121003e-07,.156902217e-06,
      .234200703e-06,.333451047e-06,.457396614e-06,.608780627e-06,
      .790346151e-06,.100483608e-05,.125499312e-05,.154355976e-05,
      .187327830e-05,.224689076e-05,.266713894e-05,.313676435e-05,
      .365850823e-05,.486931481e-05,.632148219e-05,.803692827e-05,
      .100375647e-04,.123452966e-04,.149820215e-04,.179696294e-04,
      .213300017e-04,.250850108e-04,.292565195e-04,.338663805e-04,
      .389364355e-04,.444885150e-04,.505444371e-04,.571260074e-04,
      .642550179e-04,.719532466e-04,.802424566e-04,.891443957e-04,
      .986807951e-04,.108873369e-03,.119743815e-03,.131313810e-03,
      .143605013e-03,.156639063e-03,.170437577e-03,.185022151e-03,
      .200414358e-03,.216635746e-03,.233707842e-03,.251652143e-03,
      .270490123e-03,.290243226e-03,.310932869e-03,.332580441e-03,
      .355207296e-03,.378834760e-03,.403484125e-03,.429176650e-03,
      .455933556e-03,.512725228e-03,.574028176e-03,.640010800e-03,
      .710840816e-03,.786685214e-03,.867710215e-03,.954081221e-03,
      .104596277e-02,.114351847e-02,.124691097e-02,.135630189e-02,
      .147185176e-02,.159371999e-02,.172206476e-02,.185704299e-02,
      .199881028e-02,.214752082e-02,.230332736e-02,.246638109e-02,
      .263683159e-02,.281482678e-02,.300051280e-02,.319403396e-02,
      .339553266e-02,.360514931e-02,.382302221e-02,.404928754e-02,
      .428407921e-02,.452752877e-02,.477976540e-02,.504091573e-02,
      .531110378e-02,.559045091e-02,.587907566e-02,.617709371e-02,
      .648461775e-02,.680175741e-02,.712861916e-02,.746530620e-02,
      .781191838e-02,.816855212e-02,.853530026e-02,.891225204e-02,
      .929949292e-02,.969710459e-02,.101051648e-01,.105237472e-01,
      .109529215e-01,.113927530e-01,.118433031e-01,.127767811e-01,
      .137537557e-01,.147745535e-01,.158394242e-01,.169485393e-01,
      .181019891e-01,.192997813e-01,.205418391e-01,.218279996e-01,
      .231580128e-01,.245315407e-01,.259481566e-01,.274073447e-01,
      .289085001e-01,.304509293e-01,.320338510e-01,.336563963e-01,
      .353176113e-01,.370164576e-01,.387518151e-01,.405224842e-01,
      .423271883e-01,.441645770e-01,.460332295e-01,.479316579e-01,
      .498583116e-01,.518115811e-01,.537898025e-01,.557912619e-01,
      .578142008e-01,.598568203e-01,.619172870e-01,.639937373e-01,
      .660842838e-01,.681870197e-01,.703000246e-01,.724213699e-01,
      .745491243e-01,.766813585e-01,.788161511e-01,.809515934e-01,
      .830857943e-01,.852168850e-01,.873430241e-01,.894624015e-01,
      .915732428e-01,.936738134e-01,.957624220e-01,.978374241e-01,
      .998972253e-01,.101940284e+00,.103965115e+00,.105970289e+00,
      .107954438e+00,.109916257e+00,.111854501e+00,.113767992e+00,
      .115655616e+00,.117516324e+00,.119349136e+00,.121153137e+00,
      .122927477e+00,.124671373e+00,.126384108e+00,.128065029e+00,
      .129713546e+00,.131329135e+00,.132911328e+00,.134459722e+00,
      .135973968e+00,.137453779e+00,.138898917e+00,.140309202e+00,
      .141684503e+00,.143024738e+00,.144329874e+00,.145599922e+00,
      .146834935e+00,.148035010e+00,.149200279e+00,.150330915e+00,
      .151427122e+00,.152489139e+00,.153517234e+00,.154511704e+00,
      .155472872e+00,.156401088e+00,.157732440e+00,.159395871e+00,
      .160934785e+00,.162352839e+00,.163653872e+00,.164841868e+00,
      .165920909e+00,.166895140e+00,.167768734e+00,.168545863e+00,
      .169230674e+00,.169827262e+00,.170339654e+00,.170771795e+00,
      .171127527e+00,.171410589e+00,.171624598e+00,.171773049e+00,
      .171859309e+00,.171886612e+00,.171858061e+00,.171776623e+00,
      .171645135e+00,.171466302e+00,.171242700e+00,.170976780e+00,
      .170670873e+00,.170327188e+00,.169947824e+00,.169534768e+00,
      .169089904e+00,.168615015e+00,.168111791e+00,.167581830e+00,
      .167026642e+00,.166447661e+00,.165846239e+00,.165223660e+00,
      .164581136e+00,.163919818e+00,.163240796e+00,.162545102e+00,
      .161833716e+00,.161107569e+00,.160367543e+00,.159614479e+00,
      .158849175e+00,.158072392e+00,.157284855e+00,.156487254e+00,
      .155680248e+00,.154864468e+00,.154040516e+00,.153208967e+00,
      .152370373e+00,.151525261e+00,.150674139e+00,.149817492e+00,
      .148955786e+00,.148089471e+00,.147218975e+00,.146344713e+00,
      .145467084e+00,.144586469e+00,.143703239e+00,.142817748e+00,
      .141930337e+00,.141041337e+00,.140151064e+00,.139259823e+00,
      .138367908e+00,.137475603e+00,.136583180e+00,.135690902e+00,
      .134799021e+00,.133907780e+00,.133017414e+00,.132128146e+00,
      .131240193e+00,.130353762e+00,.129469053e+00,.128586256e+00,
      .127705555e+00,.126827126e+00,.125951137e+00,.125077749e+00,
      .124207116e+00,.123339385e+00,.122474697e+00,.121613186e+00,
      .120754979e+00,.119900199e+00,.119048961e+00,.118201374e+00,
      .117357544e+00,.116517570e+00,.115681544e+00,.114849555e+00,
      .114021687e+00,.113198018e+00,.112378623e+00,.111563571e+00,
      .110752927e+00,.109946752e+00,.109145103e+00,.108348033e+00,
      .107555591e+00,.106767822e+00,.105984768e+00,.105206467e+00,
      .104432956e+00,.103664265e+00,.102900424e+00,.102141458e+00,
      .101387392e+00,.100638244e+00,.998940345e-01,.991547772e-01,
      .984204854e-01,.976911698e-01,.969668391e-01,.962474994e-01,
      .955331551e-01,.948238087e-01,.941194608e-01,.934201100e-01,
      .927257537e-01,.920363871e-01,.913520044e-01,.906725979e-01,
      .899981589e-01,.893286772e-01,.886641411e-01,.880045381e-01,
      .873498544e-01,.867000749e-01,.860551839e-01,.854151643e-01,
      .847799984e-01,.841496675e-01,.835241519e-01,.829034314e-01,
      .822874850e-01,.816762909e-01,.810698267e-01,.804680694e-01,
      .798709954e-01,.792785806e-01,.786908003e-01,.781076294e-01,
      .775290425e-01,.769550134e-01,.763855158e-01,.758205231e-01,
      .752600083e-01,.747039439e-01,.741523023e-01,.736050558e-01,
      .730621761e-01,.725236350e-01,.719894039e-01,.714594542e-01,
      .709337571e-01,.704122835e-01,.698950044e-01,.693818905e-01,
      .688729127e-01,.683680416e-01,.678672477e-01,.673705016e-01,
      .668777738e-01,.663890348e-01,.659042552e-01,.654234053e-01,
      .649464558e-01,.644733771e-01,.640041398e-01,.635387145e-01,
      .630770720e-01,.626191828e-01,.621650178e-01,.617145478e-01,
      .612677439e-01,.608245769e-01,.603850181e-01,.599490388e-01,
      .595166101e-01,.590877036e-01,.586622908e-01,.582403435e-01,
      .578218334e-01,.574067325e-01,.569950129e-01,.565866467e-01,
      .561816064e-01,.557798644e-01,.553813934e-01,.549861662e-01,
      .545941558e-01,.542053352e-01,.538196778e-01,.534371570e-01,
      .530577463e-01,.526814196e-01,.523081508e-01,.519379139e-01,
      .515706833e-01,.512064333e-01,.508451386e-01,.504867740e-01,
      .501313144e-01,.497787350e-01,.494290111e-01,.490821181e-01,
      .487380317e-01,.483967277e-01,.480581823e-01,.477223715e-01,
      .473892717e-01,.470588596e-01,.467311118e-01,.464060052e-01,
      .460835169e-01,.457636243e-01,.454463046e-01,.451315357e-01,
      .448192952e-01,.445095612e-01,.442023117e-01,.438975252e-01,
      .435951802e-01,.432952552e-01,.429977293e-01,.427025814e-01,
      .424097907e-01,.421193367e-01,.418311988e-01,.415453568e-01,
      .412617906e-01,.409804802e-01,.407014060e-01,.404245483e-01,
      .401498876e-01,.398774047e-01,.396070806e-01,.393388963e-01,
      .390728330e-01,.388088721e-01,.385469952e-01,.382871840e-01,
      .379013009e-01,.373938523e-01,.368943456e-01,.364026421e-01,
      .359186056e-01,.354421025e-01,.349730013e-01,.345111729e-01,
      .340564908e-01,.336088304e-01,.331680696e-01,.327340885e-01,
      .323067693e-01,.318859964e-01,.314716562e-01,.310636374e-01,
      .306618305e-01,.302661281e-01,.298764249e-01,.294926173e-01,
      .291146038e-01,.287422847e-01,.283755621e-01,.280143400e-01,
      .276585240e-01,.273080217e-01,.269627422e-01,.266225964e-01,
      .262874968e-01,.259573575e-01,.256320942e-01,.253116242e-01,
      .249958664e-01,.246847409e-01,.243781696e-01,.240760757e-01,
      .237783839e-01,.234850202e-01,.231959120e-01,.229109882e-01,
      .226301786e-01,.223534148e-01,.220806293e-01,.218117559e-01,
      .215467299e-01,.212854874e-01,.210279659e-01,.207741041e-01,
      .205238415e-01,.202771191e-01,.200338788e-01,.197940636e-01,
      .195576175e-01,.193244855e-01,.190946136e-01,.188679490e-01,
      .186444396e-01,.184240344e-01,.182066831e-01,.179923367e-01,
      .177809467e-01,.175724658e-01,.173668473e-01,.171640455e-01,
      .169640154e-01,.167667129e-01,.165720947e-01,.163801182e-01,
      .161907417e-01,.160039240e-01,.158196250e-01,.156378050e-01,
      .154584251e-01,.152814471e-01,.151068336e-01,.149345476e-01,
      .147645530e-01,.145968142e-01,.144312963e-01,.142679650e-01,
      .141067866e-01,.139477279e-01,.137907564e-01,.136358401e-01,
      .134829477e-01,.133320482e-01,.131831114e-01,.130361075e-01,
      .128910071e-01,.127477815e-01,.126064024e-01,.124668421e-01,
      .123290731e-01,.121930688e-01,.120588026e-01,.119262487e-01,
      .117953816e-01,.116661763e-01,.115386080e-01,.114126526e-01,
      .112882864e-01,.111654859e-01,.110442281e-01,.109244904e-01,
      .108062506e-01,.106894868e-01,.105741776e-01,.104603018e-01,
      .103478387e-01,.102367678e-01,.101270690e-01,.100187226e-01,
      .991170913e-02,.980600952e-02,.970160497e-02,.959847699e-02,
      .949660742e-02,.939597836e-02,.929657225e-02,.919837177e-02,
      .910135992e-02,.900551995e-02,.891083540e-02,.881729007e-02,
      .872486803e-02,.863355359e-02,.854333132e-02,.845418606e-02,
      .836610285e-02,.827906701e-02,.819306407e-02,.810807980e-02,
      .802410021e-02,.794111150e-02,.785910011e-02,.777805270e-02,
      .769795613e-02,.761879746e-02,.754056397e-02,.746324313e-02,
      .738682261e-02,.731129027e-02,.723663416e-02,.716284251e-02,
      .708990375e-02,.701780647e-02,.694653945e-02,.687609163e-02,
      .680645213e-02,.673761023e-02,.666955539e-02,.660227722e-02,
      .653576548e-02,.647001010e-02,.640500115e-02,.634072888e-02,
      .627718365e-02,.615223653e-02,.603008570e-02,.591065923e-02,
      .579388734e-02,.567970230e-02,.556803840e-02,.545883183e-02,
      .535202066e-02,.524754473e-02,.514534567e-02,.504536675e-02,
      .494755288e-02,.485185056e-02,.475820778e-02,.466657405e-02,
      .457690027e-02,.448913873e-02,.440324306e-02,.431916818e-02,
      .423687027e-02,.415630670e-02,.407743604e-02,.400021799e-02,
      .392461333e-02,.385058395e-02,.377809273e-02,.370710358e-02,
      .363758136e-02,.356949190e-02,.350280192e-02,.343747902e-02,
      .337349168e-02,.331080919e-02,.324940165e-02,.318923997e-02,
      .313029577e-02,.307254146e-02,.301595012e-02,.296049553e-02,
      .290615218e-02,.285289517e-02,.280070026e-02,.274954381e-02,
      .269940279e-02,.265025473e-02,.260207776e-02,.255485053e-02,
      .250855222e-02,.246316255e-02,.241866171e-02,.237503040e-02,
      .233224980e-02,.229030154e-02,.224916769e-02,.220883077e-02,
      .216927372e-02,.213047989e-02,.209243303e-02,.205511729e-02,
      .201851718e-02,.198261759e-02,.194740376e-02,.191286131e-02,
      .187897615e-02,.184573457e-02,.181312314e-02,.178112878e-02,
      .174973868e-02,.171894035e-02,.168872159e-02,.165907046e-02,
      .162997533e-02,.160142479e-02,.157340774e-02,.154591329e-02,
      .151893083e-02,.149244997e-02,.146646054e-02,.144095266e-02,
      .141591660e-02,.139134289e-02,.136722226e-02,.134354566e-02,
      .132030420e-02,.129748923e-02,.127509226e-02,.125310503e-02,
      .123151941e-02,.121032746e-02,.118952145e-02,.116909377e-02,
      .114903699e-02,.112934389e-02,.111000732e-02,.109102033e-02,
      .107237612e-02,.105406805e-02,.103608958e-02,.101843435e-02,
      .100109610e-02,.984068725e-03,.967346250e-03,.950922810e-03,
      .934792677e-03,.918950246e-03,.903390002e-03,.888106595e-03,
      .873094729e-03,.858349254e-03,.843865126e-03,.829637404e-03,
      .815661224e-03,.801931851e-03,.788444635e-03,.775195021e-03,
      .762178553e-03,.749390859e-03,.736827660e-03,.724484770e-03,
      .712358067e-03,.700443525e-03,.688737223e-03,.677235278e-03,
      .665933895e-03,.654829387e-03,.643918101e-03,.633196480e-03,
      .622661026e-03,.612308316e-03,.602134999e-03,.592137781e-03,
      .582313435e-03,.572658808e-03,.563170793e-03,.553846368e-03,
      .544682549e-03,.535676408e-03,.526825092e-03,.518125765e-03,
      .509575707e-03,.501172192e-03,.492912592e-03,.484794323e-03,
      .476814826e-03,.468971592e-03,.461262190e-03,.453684225e-03,
      .446235321e-03,.438913184e-03,.431715574e-03,.424640199e-03,
      .417684920e-03,.410847622e-03,.404126153e-03,.397518512e-03,
      .391022628e-03,.384636546e-03,.378358315e-03,.372186034e-03,
      .366117836e-03,.360151911e-03,.354286384e-03,.348519517e-03,
      .342849565e-03,.337274873e-03,.331793740e-03,.326404511e-03,
      .321105607e-03,.315895428e-03,.310772457e-03,.305735112e-03,
      .300781962e-03,.295911497e-03,.291122312e-03,.286412945e-03,
      .281782063e-03,.277228270e-03,.272750278e-03,.268346756e-03,
      .264016369e-03,.259757905e-03,.255570068e-03,.251451654e-03,
      .247401513e-03,.243418470e-03,.239501302e-03,.235648942e-03,
      .231860214e-03,.228134117e-03,.224469519e-03,.220865364e-03,
      .217320641e-03,.213834284e-03,.210405372e-03,.207032896e-03,
      .203715899e-03,.200453371e-03,.197244484e-03,.194088300e-03,
      .190983864e-03,.187930274e-03,.184926797e-03,.181972511e-03,
      .179066591e-03,.176208209e-03,.173396556e-03,.170630842e-03,
      .167910368e-03,.165234252e-03,.162601776e-03,.160012280e-03,
      .157464935e-03,.154959117e-03,.152494072e-03,.150069121e-03,
      .147683673e-03,.145336940e-03,.143028425e-03,.140757356e-03,
      .138523108e-03,.136325129e-03,.134162886e-03,.132035552e-03,
      .129942686e-03,.127883754e-03,.125858188e-03,.123865324e-03,
      .121904796e-03,.119975813e-03,.118078062e-03,.116210845e-03,
      .114373886e-03,.112566524e-03,.110788316e-03 };

  static double p[500] =
    { .00,5.00,10.00,15.00,20.00,25.00,30.00,35.00,40.00,45.00,50.00,
      55.00,60.00,65.00,70.00,75.00,80.00,85.00,90.00,95.00,100.00,105.00,
      110.00,115.00,120.00,125.00,130.00,135.00,140.00,145.00,150.00,155.00,
      160.00,165.00,170.00,175.00,180.00,185.00,190.00,195.00,200.00,205.00,
      210.00,215.00,220.00,225.00,230.00,235.00,240.00,245.00,250.00,255.00,
      260.00,265.00,270.00,275.00,280.00,285.00,290.00,295.00,300.00,305.00,
      310.00,315.00,320.00,325.00,330.00,335.00,340.00,345.00,350.00,355.00,
      360.00,365.00,370.00,375.00,380.00,385.00,390.00,395.00,400.00,405.00,
      410.00,415.00,420.00,425.00,430.00,435.00,440.00,445.00,450.00,455.00,
      460.00,465.00,470.00,475.00,480.00,485.00,490.00,495.00,500.00,505.00,
      510.00,515.00,520.00,525.00,530.00,535.00,540.00,545.00,550.00,555.00,
      560.00,565.00,570.00,575.00,580.00,585.00,590.00,595.00,600.00,605.00,
      610.00,615.00,620.00,625.00,630.00,635.00,640.00,645.00,650.00,655.00,
      660.00,665.00,670.00,675.00,680.00,685.00,690.00,695.00,700.00,705.00,
      710.00,715.00,720.00,725.00,730.00,735.00,740.00,745.00,750.00,755.00,
      760.00,765.00,770.00,775.00,780.00,785.00,790.00,795.00,800.00,805.00,
      810.00,815.00,820.00,825.00,830.00,835.00,840.00,845.00,850.00,855.00,
      860.00,865.00,870.00,875.00,880.00,885.00,890.00,895.00,900.00,905.00,
      910.00,915.00,920.00,925.00,930.00,935.00,940.00,945.00,950.00,955.00,
      960.00,965.00,970.00,975.00,980.00,985.00,990.00,995.00,1000.00,1005.00,
      1010.00,1015.00,1020.00,1025.00,1030.00,1035.00,1040.00,1045.00,1050.00,1055.00,
      1060.00,1065.00,1070.00,1075.00,1080.00,1085.00,1090.00,1095.00,1100.00,1105.00,
      1110.00,1115.00,1120.00,1125.00,1130.00,1135.00,1140.00,1145.00,1150.00,1155.00,
      1160.00,1165.00,1170.00,1175.00,1180.00,1185.00,1190.00,1195.00,1200.00,1205.00,
      1210.00,1215.00,1220.00,1225.00,1230.00,1235.00,1240.00,1245.00,1250.00,1255.00,
      1260.00,1265.00,1270.00,1275.00,1280.00,1285.00,1290.00,1295.00,1300.00,1305.00,
      1310.00,1315.00,1320.00,1325.00,1330.00,1335.00,1340.00,1345.00,1350.00,1355.00,
      1360.00,1365.00,1370.00,1375.00,1380.00,1385.00,1390.00,1395.00,1400.00,1405.00,
      1410.00,1415.00,1420.00,1425.00,1430.00,1435.00,1440.00,1445.00,1450.00,1455.00,
      1460.00,1465.00,1470.00,1475.00,1480.00,1485.00,1490.00,1495.00,1500.00,1505.00,
      1510.00,1515.00,1520.00,1525.00,1530.00,1535.00,1540.00,1545.00,1550.00,1555.00,
      1560.00,1565.00,1570.00,1575.00,1580.00,1585.00,1590.00,1595.00,1600.00,1605.00,
      1610.00,1615.00,1620.00,1625.00,1630.00,1635.00,1640.00,1645.00,1650.00,1655.00,
      1660.00,1665.00,1670.00,1675.00,1680.00,1685.00,1690.00,1695.00,1700.00,1705.00,
      1710.00,1715.00,1720.00,1725.00,1730.00,1735.00,1740.00,1745.00,1750.00,1755.00,
      1760.00,1765.00,1770.00,1775.00,1780.00,1785.00,1790.00,1795.00,1800.00,1805.00,
      1810.00,1815.00,1820.00,1825.00,1830.00,1835.00,1840.00,1845.00,1850.00,1855.00,
      1860.00,1865.00,1870.00,1875.00,1880.00,1885.00,1890.00,1895.00,1900.00,1905.00,
      1910.00,1915.00,1920.00,1925.00,1930.00,1935.00,1940.00,1945.00,1950.00,1955.00,
      1960.00,1965.00,1970.00,1975.00,1980.00,1985.00,1990.00,1995.00,2000.00,2005.00,
      2010.00,2015.00,2020.00,2025.00,2030.00,2035.00,2040.00,2045.00,2050.00,2055.00,
      2060.00,2065.00,2070.00,2075.00,2080.00,2085.00,2090.00,2095.00,2100.00,2105.00,
      2110.00,2115.00,2120.00,2125.00,2130.00,2135.00,2140.00,2145.00,2150.00,2155.00,
      2160.00,2165.00,2170.00,2175.00,2180.00,2185.00,2190.00,2195.00,2200.00,2205.00,
      2210.00,2215.00,2220.00,2225.00,2230.00,2235.00,2240.00,2245.00,2250.00,2255.00,
      2260.00,2265.00,2270.00,2275.00,2280.00,2285.00,2290.00,2295.00,2300.00,2305.00,
      2310.00,2315.00,2320.00,2325.00,2330.00,2335.00,2340.00,2345.00,2350.00,2355.00,
      2360.00,2365.00,2370.00,2375.00,2380.00,2385.00,2390.00,2395.00,2400.00,2405.00,
      2410.00,2415.00,2420.00,2425.00,2430.00,2435.00,2440.00,2445.00,2450.00,2455.00,
      2460.00,2465.00,2470.00,2475.00,2480.00,2485.00,2490.00,2495.00};
 
  static double up[500] =
    { 12.7566650000,12.6010335622,12.1554909056,11.4771906262,
      10.6421873651,9.7270685394,8.7959746677,7.8951494666,
      7.0533352354,6.2852005618,5.5954865046,4.9826713081,
      4.4417188436,3.9659376537,3.5481395535,3.1813075374,
      2.8589402194,2.5751985935,2.3249427470,2.1037000539,
      1.9076117395,1.7333607699,1.5781031624,1.4394034886,
      1.3151736468,1.2036218234,1.1032091192,1.0126099540,
      .9306802413,.8564296361,.7889992997,.7276421533,
      .6717065066,.6206225086,.5738903851,.5310707181,
      .4917762279,.4556646766,.4224327918,.3918113741,
      .3635607578,.3374671676,.3133394867,.2910064332,
      .2703142358,.2511245873,.2333127773,.2167662495,
      .2013832596,.1870715304,.1737472075,.1613341524,
      .1497631038,.1389708335,.1288995350,.1194963914,
      .1107130947,.1025052827,.0948321774,.0876563974,
      .0809435195,.0746617303,.0687817999,.0632767276,
      .0581214409,.0532928865,.0487696759,.0445318889,
      .0405611894,.0368404523,.0333537849,.0300864662,
      .0270246760,.0241556772,.0214675052,.0189490317,
      .0165898955,.0143803515,.0123113970,.0103745020,
      .0085617971,.0068658206,.0052796599,.0037967852,
      .0024111132,.0011169204,-.0000911583,-.0012181529,
      -.0022687972,-.0032475220,-.0041585068,-.0050056653,
      -.0057926901,-.0065230462,-.0071999925,-.0078266052,
      -.0084057635,-.0089402059,-.0094324737,-.0098850080,
      -.0103000557,-.0106797925,-.0110262127,-.0113412466,
      -.0116266781,-.0118842113,-.0121154519,-.0123218991,
      -.0125050002,-.0126660804,-.0128064322,-.0129272452,
      -.0130296581,-.0131147522,-.0131835266,-.0132369584,
      -.0132759394,-.0133013322,-.0133139521,-.0133145540,
      -.0133038781,-.0132826004,-.0132513718,-.0132108154,
      -.0131615012,-.0131039902,-.0130388022,-.0129664255,
      -.0128873395,-.0128019797,-.0127107672,-.0126141081,
      -.0125123712,-.0124059189,-.0122950945,-.0121802126,
      -.0120615848,-.0119395017,-.0118142329,-.0116860454,
      -.0115551865,-.0114218867,-.0112863748,-.0111488620,
      -.0110095459,-.0108686228,-.0107262739,-.0105826678,
      -.0104379722,-.0102923435,-.0101459257,-.0099988624,
      -.0098512889,-.0097033282,-.0095551019,-.0094067276,
      -.0092583107,-.0091099541,-.0089617596,-.0088138182,
      -.0086662165,-.0085190411,-.0083723718,-.0082262813,
      -.0080808424,-.0079361250,-.0077921909,-.0076491001,
      -.0075069124,-.0073656817,-.0072254573,-.0070862886,
      -.0069482231,-.0068113016,-.0066755640,-.0065410498,
      -.0064077948,-.0062758307,-.0061451887,-.0060158997,
      -.0058879900,-.0057614832,-.0056364034,-.0055127731,
      -.0053906111,-.0052699344,-.0051507604,-.0050331051,
      -.0049169804,-.0048023977,-.0046893685,-.0045779028,
      -.0044680072,-.0043596881,-.0042529517,-.0041478031,
      -.0040442443,-.0039422773,-.0038419035,-.0037431239,
      -.0036459369,-.0035503403,-.0034563318,-.0033639090,
      -.0032730667,-.0031837993,-.0030961009,-.0030099657,
      -.0029253861,-.0028423538,-.0027608602,-.0026808966,
      -.0026024538,-.0025255207,-.0024500861,-.0023761389,
      -.0023036676,-.0022326597,-.0021631022,-.0020949819,
      -.0020282857,.0019630000,-.0018991102,-.0018366013,
      -.0017754587,-.0017156676,-.0016572124,-.0016000772,
      -.0015442461,-.0014897030,-.0014364319,-.0013844162,
      -.0013336392,-.0012840840,-.0012357339,-.0011885720,
      -.0011425810,-.0010977437,-.0010540428,-.0010114609,
      -.0009699808,-.0009295850,-.0008902561,-.0008519767,
      -.0008147292,-.0007784964,-.0007432608,-.0007090051,
      -.0006757120,-.0006433642,-.0006119443,-.0005814354,
      -.0005518205,-.0005230827,-.0004952051,-.0004681707,
      -.0004419630,-.0004165656,-.0003919622,-.0003681364,
      -.0003450720,-.0003227527,-.0003011631,-.0002802874,
      -.0002601103,-.0002406163,-.0002217897,-.0002036156,
      -.0001860795,-.0001691667,-.0001528627,-.0001371529,
      -.0001220229,-.0001074589,-.0000934474,-.0000799749,
      -.0000670278,-.0000545927,-.0000426564,-.0000312061,
      -.0000202297,-.0000097148,.0000003511,.0000099802,
      .0000191846,.0000279757,.0000363646,.0000443625,
      .0000519805,.0000592299,.0000661215,.0000726655,
      .0000788718,.0000847505,.0000903114,.0000955647,
      .0001005199,.0001051864,.0001095728,.0001136880,
      .0001175406,.0001211395,.0001244934,.0001276104,
      .0001304984,.0001331648,.0001356173,.0001378634,
      .0001399109,.0001417670,.0001434387,.0001449322,
      .0001462542,.0001474112,.0001484096,.0001492560,
      .0001499564,.0001505165,.0001509418,.0001512376,
      .0001514094,.0001514627,.0001514030,.0001512353,
      .0001509642,.0001505940,.0001501295,.0001495750,
      .0001489352,.0001482144,.0001474169,.0001465462,
      .0001456062,.0001446003,.0001435323,.0001424059,
      .0001412246,.0001399919,.0001387107,.0001373839,
      .0001360144,.0001346051,.0001331589,.0001316786,
      .0001301671,.0001286268,.0001270597,.0001254683,
      .0001238547,.0001222210,.0001205697,.0001189027,
      .0001172222,.0001155299,.0001138273,.0001121161,
      .0001103979,.0001086745,.0001069474,.0001052184,
      .0001034886,.0001017594,.0001000319,.0000983072,
      .0000965865,.0000948710,.0000931619,.0000914603,
      .0000897671,.0000880830,.0000864089,.0000847456,
      .0000830937,.0000814542,.0000798277,.0000782151,
      .0000766168,.0000750334,.0000734654,.0000719132,
      .0000703772,.0000688580,.0000673560,.0000658716,
      .0000644052,.0000629571,.0000615275,.0000601166,
      .0000587246,.0000573517,.0000559982,.0000546643,
      .0000533501,.0000520558,.0000507814,.0000495270,
      .0000482926,.0000470782,.0000458838,.0000447096,
      .0000435555,.0000424216,.0000413078,.0000402140,
      .0000391401,.0000380860,.0000370517,.0000360370,
      .0000350419,.0000340662,.0000331100,.0000321729,
      .0000312549,.0000303558,.0000294753,.0000286134,
      .0000277699,.0000269445,.0000261372,.0000253477,
      .0000245758,.0000238214,.0000230841,.0000223639,
      .0000216604,.0000209734,.0000203028,.0000196483,
      .0000190097,.0000183867,.0000177792,.0000171869,
      .0000166095,.0000160469,.0000154987,.0000149648,
      .0000144449,.0000139388,.0000134462,.0000129670,
      .0000125008,.0000120474,.0000116066,.0000111782,
      .0000107619,.0000103575,.0000099647,.0000095834,
      .0000092133,.0000088541,.0000085056,.0000081677,
      .0000078401,.0000075226,.0000072149,.0000069169,
      .0000066283,.0000063489,.0000060786,.0000058169,
      .0000055639,.0000053193,.0000050829,.0000048544,
      .0000046338,.0000044208,.0000042152,.0000040169,
      .0000038256,.0000036411,.0000034633,.0000032920,
      .0000031271,.0000029683,.0000028155,.0000026686,
      .0000025274,.0000023917,.0000022614,.0000021362,
      .0000020161,.0000019009,.0000017904,.0000016846,
      .0000015833,.0000014863,.0000013936,.0000013050,
      .0000012203,.0000011395,.0000010623,.0000009888,
      .0000009187,.0000008519,.0000007884,.0000007280,
      .0000006707,.0000006163,.0000005648,.0000005160 };

  static double wp[500] =
    { .0000000000,.0037707224,.0145447075,.0308822592,.0508683032,
      .0725761561,.0943949724,.1151672199,.1341785338,.1510710367,
      .1657384927,.1782336221,.1886985646,.1973177105,.2042883544,
      .2098035881,.2140435297,.2171714184,.2193325597,.2206550564,
      .2212510406,.2212184068,.2206424835,.2195976322,.2181487244,
      .2163524225,.2142582767,.2119096886,.2093447608,.2065969760,
      .2036957820,.2006671416,.1975339203,.1943163266,.1910321686,
      .1876971869,.1843252699,.1809286761,.1775182011,.1741033936,
      .1706926060,.1672932007,.1639116267,.1605535152,.1572237819,
      .1539266909,.1506659284,.1474446440,.1442655662,.1411309912,
      .1380428197,.1350026552,.1320118066,.1290713085,.1261819597,
      .1233443572,.1205589024,.1178258432,.1151452883,.1125171928,
      .1099414089,.1074176986,.1049457255,.1025250724,.1001552586,
      .0978357452,.0955659410,.0933452103,.0911728773,.0890482390,
      .0869705672,.0849391063,.0829530843,.0810117185,.0791142157,
      .0772597707,.0754475742,.0736768152,.0719466826,.0702563666,
      .0686050615,.0669919639,.0654162746,.0638772019,.0623739642,
      .0609057877,.0594719064,.0580715658,.0567040214,.0553685404,
      .0540644004,.0527908907,.0515473104,.0503329739,.0491472095,
      .0479893558,.0468587626,.0457547915,.0446768210,.0436242381,
      .0425964461,.0415928552,.0406128946,.0396560005,.0387216251,
      .0378092291,.0369182866,.0360482854,.0351987205,.0343691042,
      .0335589536,.0327678022,.0319951917,.0312406738,.0305038147,
      .0297841846,.0290813705,.0283949650,.0277245709,.0270698043,
      .0264302846,.0258056460,.0251955285,.0245995798,.0240174615,
      .0234488376,.0228933841,.0223507855,.0218207295,.0213029173,
      .0207970547,.0203028528,.0198200356,.0193483294,.0188874677,
      .0184371943,.0179972544,.0175674021,.0171474006,.0167370138,
      .0163360146,.0159441830,.0155612998,.0151871549,.0148215452,
      .0144642680,.0141151288,.0137739393,.0134405115,.0131146648,
      .0127962254,.0124850202,.0121808810,.0118836471,.0115931584,
      .0113092578,.0110317969,.0107606286,.0104956071,.0102365938,
      .0099834535,.0097360508,.0094942558,.0092579441,.0090269915,
      .0088012757,.0085806814,.0083650944,.0081544009,.0079484924,
      .0077472641,.0075506115,.0073584318,.0071706278,.0069871036,
      .0068077634,.0066325151,.0064612707,.0062939424,.0061304429,
      .0059706895,.0058146019,.0056620996,.0055131035,.0053675386,
      .0052253317,.0050864093,.0049506992,.0048181337,.0046886454,
      .0045621672,.0044386335,.0043179823,.0042001522,.0040850819,
      .0039727112,.0038629836,.0037558429,.0036512329,.0035490986,
      .0034493880,.0033520499,.0032570332,.0031642871,.0030737638,
      .0029854162,.0028991974,.0028150608,.0027329617,.0026528577,
      .0025747060,.0024984639,.0024240901,.0023515452,.0022807900,
      .0022117856,.0021444931,.0020788761,.0020148990,.0019525264,
      .0018917227,.0018324534,.0017746859,.0017183879,.0016635268,
      .0016100705,.0015579882,.0015072505,.0014578278,.0014096908,
      .0013628103,.0013171589,.0012727097,.0012294362,.0011873115,
      .0011463094,.0011064052,.0010675746,.0010297934,.0009930373,
      .0009572825,.0009225068,.0008886882,.0008558049,.0008238348,
      .0007927565,.0007625497,.0007331947,.0007046717,.0006769608,
      .0006500426,.0006238986,.0005985111,.0005738623,.0005499343,
      .0005267096,.0005041711,.0004823030,.0004610893,.0004405141,
      .0004205614,.0004012157,.0003824624,.0003642875,.0003466766,
      .0003296156,.0003130903,.0002970873,.0002815939,.0002665976,
      .0002520856,.0002380454,.0002244645,.0002113314,.0001986348,
      .0001863637,.0001745068,.0001630530,.0001519914,.0001413120,
      .0001310049,.0001210604,.0001114686,.0001022197,.0000933042,
      .0000847134,.0000764387,.0000684715,.0000608032,.0000534253,
      .0000463293,.0000395077,.0000329531,.0000266579,.0000206147,
      .0000148161,.0000092547,.0000039238,-.0000011830,-.0000060721,
      -.0000107500,-.0000152231,-.0000194977,-.0000235799,-.0000274752,
      -.0000311891,-.0000347271,-.0000380946,-.0000412972,-.0000443403,
      -.0000472286,-.0000499669,-.0000525597,-.0000550118,-.0000573278,
      -.0000595124,-.0000615700,-.0000635045,-.0000653198,-.0000670200,
      -.0000686088,-.0000700903,-.0000714682,-.0000727463,-.0000739279,
      -.0000750161,-.0000760142,-.0000769255,-.0000777532,-.0000785007,
      -.0000791708,-.0000797663,-.0000802898,-.0000807441,-.0000811318,
      -.0000814556,-.0000817182,-.0000819220,-.0000820694,-.0000821623,
      -.0000822032,-.0000821940,-.0000821370,-.0000820344,-.0000818882,
      -.0000817003,-.0000814725,-.0000812065,-.0000809039,-.0000805667,
      -.0000801964,-.0000797948,-.0000793636,-.0000789040,-.0000784175,
      -.0000779055,-.0000773692,-.0000768101,-.0000762295,-.0000756288,
      -.0000750091,-.0000743715,-.0000737171,-.0000730469,-.0000723619,
      -.0000716633,-.0000709520,-.0000702291,-.0000694955,-.0000687520,
      -.0000679993,-.0000672383,-.0000664697,-.0000656943,-.0000649130,
      -.0000641265,-.0000633354,-.0000625404,-.0000617420,-.0000609408,
      -.0000601373,-.0000593322,-.0000585261,-.0000577194,-.0000569127,
      -.0000561065,-.0000553011,-.0000544969,-.0000536944,-.0000528939,
      -.0000520958,-.0000513006,-.0000505085,-.0000497201,-.0000489354,
      -.0000481548,-.0000473785,-.0000466068,-.0000458399,-.0000450782,
      -.0000443218,-.0000435711,-.0000428262,-.0000420873,-.0000413544,
      -.0000406278,-.0000399077,-.0000391941,-.0000384873,-.0000377874,
      -.0000370945,-.0000364087,-.0000357302,-.0000350589,-.0000343950,
      -.0000337385,-.0000330895,-.0000324481,-.0000318143,-.0000311883,
      -.0000305700,-.0000299595,-.0000293568,-.0000287619,-.0000281748,
      -.0000275956,-.0000270242,-.0000264607,-.0000259050,-.0000253573,
      -.0000248174,-.0000242853,-.0000237611,-.0000232446,-.0000227359,
      -.0000222350,-.0000217417,-.0000212562,-.0000207782,-.0000203079,
      -.0000198452,-.0000193900,-.0000189422,-.0000185019,-.0000180689,
      -.0000176432,-.0000172247,-.0000168135,-.0000164093,-.0000160123,
      -.0000156223,-.0000152392,-.0000148630,-.0000144936,-.0000141310,
      -.0000137750,-.0000134256,-.0000130828,-.0000127464,-.0000124164,
      -.0000120927,-.0000117753,-.0000114641,-.0000111589,-.0000108598,
      -.0000105666,-.0000102793,-.0000099977,-.0000097219,-.0000094518,
      -.0000091872,-.0000089280,-.0000086743,-.0000084260,-.0000081829,
      -.0000079450,-.0000077122,-.0000074844,-.0000072616,-.0000070436,
      -.0000068305,-.0000066221,-.0000064184,-.0000062192,-.0000060246,
      -.0000058344,-.0000056486,-.0000054671,-.0000052898,-.0000051166,
      -.0000049475,-.0000047825,-.0000046214,-.0000044641,-.0000043107,
      -.0000041610,-.0000040150,-.0000038726,-.0000037338,-.0000035984 };

  for(int i=0; i<839; ++i) {
    implementation->AddUr(r[i],ur[i]);
    implementation->AddWr(r[i],wr[i]);
  }
  for(int i=0; i<500; ++i) {
    implementation->AddUp(p[i],
			  up[i]/Power(TWavefunctionImplementation::fgHbarc,3./2.));
    implementation->AddWp(p[i],
			  wp[i]/Power(TWavefunctionImplementation::fgHbarc,3./2.));
  }

  wf->fImplementation = implementation;
  
  return wf;
}



