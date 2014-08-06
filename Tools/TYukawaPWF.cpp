/*
 * TWavefunctionImplementation.cpp
 * 
 * Author: Pieter Vancraeyveld
 *         pieter.vancraeyveld@ugent.be
 *
 * Worked started February,5 2009
 */

////////////////////////////////////////////////////////////////////////
//                                                                    //
// TYukawaPWF                                                         //
//                                                                    //
// TYukawaPWF implements the abstract interface                       //
// TWavefunctionImplementation.                                       //
//                                                                    //
// This class parametrizes non-relativistic wavefunctions as a        //
// discrete superposition of Yukawa-type terms.                       //
// This parametrization was first introduced in PLB 101(1981)139.     //
//                                                                    //
// TWavefunctionImplementation fixes the following conventions        //
// concerning the units:                                              //
// * position in [fm]                                                 //
// * momentum in [MeV]                                                //
// * r-space wavefunctions in [fm^-1/2]                               //
// * p-space wavefunctions in [MeV^-3/2]                              //
//                                                                    //
// BEGIN_LATEX                                                        
// u(r) = #sum^{N}_{i=1} C_{i} e^{-m_{i}r}
// w(r) = #sum^{N}_{i=1} D_{i} e^{-m_{i}r} #left(1+#frac{3}{m_{i}r}+#frac{3}{m^{2}_{i}r^{2}}#right)
// v_{t}(r) = v_{s}(r) = 0
// END_LATEX                                                          //
// with                                                               //
// BEGIN_LATEX
// m_{i} = #alpha + (i-1) m_{0}
// END_LATEX
// In momentum space (after Fourier transformation) the parametrization is
// BEGIN_LATEX
// u(p) = #left(#frac{2}{#pi}#right)^{1/2} #sum^{N}_{i=1} #frac{C_{i}}{p^{2}+m_{i}^{2}}
// w(p) = #left(#frac{2}{#pi}#right)^{1/2} #sum^{N}_{i=1} #frac{-D_{i}}{p^{2}+m_{i}^{2}} 
// v_{t}(p) = v_{s}(p) = 0
// END_LATEX
// Imposing boundary conditions on the wavefunction in configuration  //
// space, leads to the following constraints for the parameters:      //
// BEGIN_LATEX
// #sum^{N}_{i} C_{i} = #sum^{N}_{i} D_{i} = #sum^{N}_{i} D_{i} m_{i}^{2} = #sum^{N}_{i} #frac{D_{i}}{m_{i}^{2}} = 0
// END_LATEX 
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TYukawaPWF.h"
// #include "Structures.h"
#include <iostream>
#include <cstdlib>
#include <TMath.h>
#include <TDecompLU.h>
#include <TVectorT.h>
#include <TMatrixT.h>

using std::cout; using std::cerr; using std::endl;

ClassImp(TYukawaPWF)

//_____________________________________________________________________
TYukawaPWF::TYukawaPWF(int order, double alpha, double m0,
		       double *c, double *d)
: TWavefunctionImplementation(0), fOrder(order), fAlpha(alpha), fM0(m0),
  fC(new double[order]), fD(new double[order]), fM(NULL)
{
  // Constructor
  //
  // order = number of terms N in the parametrization
  // alpha = Parameter BEGIN_LATEX #alpha END_LATEX in BEGIN_LATEX m_{i} END_LATEX formula [fm^-1]
  // m0 = Parameter BEGIN_LATEX m_{0} END_LATEX in BEGIN_LATEX m_{i} END_LATEX formula [fm^-1]
  // c,d = arrays containing parameters. [fm^1/2]

  if(!c || !d) {
    cerr << "Error in TYukawaPWF::TYukawaPWF(int,double*,double*): "
	 << "double arrays are NULL!\n";
    exit(1);
  }

  if(fOrder<1) {
    cerr << "Error in TYukawaPWF::TYukawaPWF(int,double*,double*): "
	 << "order should be >=1!\n";
    exit(1);
  }

  for(int i=0; i<fOrder; ++i) {
    fC[i] = c[i];
    fD[i] = d[i];
  }

  // C_n, D_n, D_(n-1) and D_(n-2) are fixed,
  // to satify the boundary conditions
  // Finding C_n is trivial
  fC[fOrder-1] = 0.;
  for(int i=0; i<fOrder-1; ++i) fC[fOrder-1] -= fC[i];

  // To find D_n, D_(n-1) and D_(n-2) we solve a 3x3
  // set of equations.
  // We use ROOT's LU decomposition solver method.
  double sumD=0., sumDm=0., sumD_m=0.;
  for(int i=0; i<fOrder-3; ++i) {
    sumD += fD[i];
    sumDm += fD[i]*Mass(i+1)*Mass(i+1);
    sumD_m += fD[i]/Mass(i+1)/Mass(i+1);
  }
   
  double elements[9] = {1/Mass(fOrder)/Mass(fOrder),
			1/Mass(fOrder-1)/Mass(fOrder-1),
			1/Mass(fOrder-2)/Mass(fOrder-2),
			1.,1.,1.,
			Mass(fOrder)*Mass(fOrder),
			Mass(fOrder-1)*Mass(fOrder-1),
			Mass(fOrder-2)*Mass(fOrder-2),
  };
  TMatrixD extMatrix(3,3,elements);
  elements[0]=-1.*sumD_m;
  elements[1]=-1.*sumD;
  elements[2]=-1.*sumDm;
  TVectorD sol(3,elements);
  TDecompLU lu(extMatrix);
  lu.SetTol(1e-16);
  lu.Solve(sol);

  fD[fOrder-3] = sol[2];
  fD[fOrder-2] = sol[1];
  fD[fOrder-1] = sol[0];
}

TYukawaPWF::TYukawaPWF(int order, double *c, double *d, double *m)
: TWavefunctionImplementation(0), fOrder(order),
  fC(new double[order]), fD(new double[order]), fM(new double[order])
{
  // Constructor
  //
  // order = number of terms N in the parametrization
  // c,d,m = arrays containing parameters. [fm^1/2]

  if(!c || !d || !m) {
    cerr << "Error in TYukawaPWF::TYukawaPWF(int,double*,double*): "
	 << "double arrays are NULL!\n";
    exit(1);
  }

  if(fOrder<1) {
    cerr << "Error in TYukawaPWF::TYukawaPWF(int,double*,double*): "
	 << "order should be >=1!\n";
    exit(1);
  }

  for(int i=0; i<fOrder; ++i) {
    fC[i] = c[i];
    fD[i] = d[i];
    fM[i] = m[i];
  }

  // C_n, D_n, D_(n-1) and D_(n-2) are fixed,
  // to satify the boundary conditions
  // Finding C_n is trivial
  fC[fOrder-1] = 0.;
  for(int i=0; i<fOrder-1; ++i) fC[fOrder-1] -= fC[i];

  // To find D_n, D_(n-1) and D_(n-2) we solve a 3x3
  // set of equations.
  // We use ROOT's LU decomposition solver method.
  double sumD=0., sumDm=0., sumD_m=0.;
  for(int i=0; i<fOrder-3; ++i) {
    sumD += fD[i];
    sumDm += fD[i]*Mass(i+1)*Mass(i+1);
    sumD_m += fD[i]/Mass(i+1)/Mass(i+1);
  }
   
  double elements[9] = {1/Mass(fOrder)/Mass(fOrder),
			1/Mass(fOrder-1)/Mass(fOrder-1),
			1/Mass(fOrder-2)/Mass(fOrder-2),
			1.,1.,1.,
			Mass(fOrder)*Mass(fOrder),
			Mass(fOrder-1)*Mass(fOrder-1),
			Mass(fOrder-2)*Mass(fOrder-2),
  };
  TMatrixD extMatrix(3,3,elements);
  elements[0]=-1.*sumD_m;
  elements[1]=-1.*sumD;
  elements[2]=-1.*sumDm;
  TVectorD sol(3,elements);
  TDecompLU lu(extMatrix);
  lu.SetTol(1e-16);
  lu.Solve(sol);

  fD[fOrder-3] = sol[2];
  fD[fOrder-2] = sol[1];
  fD[fOrder-1] = sol[0];
}

//_____________________________________________________________________
TYukawaPWF::TYukawaPWF(TRootIOCtor*)
  : TWavefunctionImplementation(0), fOrder(0), fAlpha(0), fM0(0), fC(0), fD(0)
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
TYukawaPWF::TYukawaPWF(const TYukawaPWF& rhs)
  : TWavefunctionImplementation(0),
    fOrder(rhs.fOrder), fAlpha(rhs.fAlpha), fM0(rhs.fM0),
    fC(new double[rhs.fOrder]), fD(new double[rhs.fOrder])
{
  // Copy constructor

  for(int i=0; i<fOrder; ++i) {
    fC[i] = rhs.fC[i];
    fD[i] = rhs.fD[i];
  }
}

//_____________________________________________________________________
TYukawaPWF& TYukawaPWF::operator=(const TYukawaPWF& rhs)
{
  // Assignment

  if(this!=&rhs) { // avoid self-assignment
    if( fOrder!=rhs.fOrder ) {
      delete[] fC; delete[] fD;
      fOrder = rhs.fOrder;
      fC = new double[fOrder];
      fD = new double[fOrder];
    }
    for(int i=0; i<fOrder; ++i) {
      fC[i] = rhs.fC[i];
      fD[i] = rhs.fD[i];
    }
    fAlpha = rhs.fAlpha;
    fM0 = rhs.fM0;
  }

  return *this;
}

//_____________________________________________________________________
TYukawaPWF::~TYukawaPWF()
{
  // Destructor

  delete[] fC;
  delete[] fD;
  if(fM) delete[] fM;
}

//_____________________________________________________________________
double TYukawaPWF::Mass(int i) const
{
  // Mass formula:
  // BEGIN_LATEX
  // m_{i} = #alpha + (i-1) m_{0}
  // END_LATEX
  // Return value in units [fm^-1]

  return (fM? fM[i-1]: fAlpha + (i-1) * fM0);
}

//_____________________________________________________________________
double TYukawaPWF::GetUp(double p) const
{
  // L=0 wave in momentum space.
  // p = momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double radial=0.;

  for(int i=0; i<fOrder; ++i)
      radial += fC[i] / ( p*p/fgHbarc/fgHbarc + Mass(i+1)*Mass(i+1) );

  return radial * TMath::Sqrt(2./TMath::Pi()) / pow(fgHbarc,3/2.);
}

//_____________________________________________________________________
double TYukawaPWF::GetWp(double p) const
{
  // L=2 wave in momentum space.
  // p = momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double radial=0.;

  for(int i=0; i<fOrder; ++i)
      radial -= fD[i] / ( p*p/fgHbarc/fgHbarc + Mass(i+1)*Mass(i+1) );

  return radial * TMath::Sqrt(2./TMath::Pi()) / pow(fgHbarc,3/2.);
}

double TYukawaPWF::GetUpoff(const TVector3& p) const{

  double radial=0.;
  
  for(int i=0; i<fOrder; ++i)
      radial += fC[i] / ( p.Mag()*p.Mag()/fgHbarc/fgHbarc + Mass(i+1)*Mass(i+1) )
		  /sqrt(p.Perp2()/fgHbarc/fgHbarc+Mass(i+1)*Mass(i+1));

  return radial * TMath::Sqrt(2./TMath::Pi()) / pow(fgHbarc,3/2.);
  
}

double TYukawaPWF::GetWpoff1(const TVector3& p) const{

  double radial=0.;
  
  for(int i=0; i<fOrder; ++i)
      radial -= fD[i] / ( p.Mag()*p.Mag()/fgHbarc/fgHbarc + Mass(i+1)*Mass(i+1) )
		  /sqrt(p.Perp2()/fgHbarc/fgHbarc+Mass(i+1)*Mass(i+1));

  return radial * TMath::Sqrt(2./TMath::Pi()) / pow(fgHbarc,3/2.);
  
}

double TYukawaPWF::GetWpoff2(double pperp2) const{
  double radial=0.;
  
  for(int i=0; i<fOrder; ++i)
      radial -= fD[i] / ( Mass(i+1)*Mass(i+1) )
		  /sqrt(pperp2/fgHbarc/fgHbarc+Mass(i+1)*Mass(i+1));

  return radial * TMath::Sqrt(2./TMath::Pi()) / pow(fgHbarc,3/2.);
  
}


//_____________________________________________________________________
double TYukawaPWF::GetVTp(double p) const
{
  // TYukawaPWF is used for non-relativistic wavefunctions.
  // L=1 wavefunctions have a relativistic origin.
  // GetVTp always returns 0.

  return 0;
}

//_____________________________________________________________________
double TYukawaPWF::GetVSp(double p) const
{
  // TYukawaPWF is used for non-relativistic wavefunctions.
  // L=1 wavefunctions have a relativistic origin.
  // GetVSp always returns 0.

  return 0;
}

//_____________________________________________________________________
double TYukawaPWF::GetUr(double r) const
{
  // L=0 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double radial=0.;

  for(int i=0; i<fOrder; ++i)
    radial += fC[i]*TMath::Exp(-1.*Mass(i+1)*r);

  return radial;
}

//_____________________________________________________________________
double TYukawaPWF::GetWr(double r) const
{
  // L=2 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double radial=0.;

  // convergence for r->0 is really bad. We have verified that
  // for r<10^-5 the wavefunction can be considered zero.
  if( r<1e-5 )
    return 0;

  for(int i=0; i<fOrder; ++i)
    radial += fD[i]*TMath::Exp(-1.*Mass(i+1)*r)
      * (1.+ 3./Mass(i+1)/r + 3/TMath::Power(Mass(i+1)*r,2.));

  return radial;
}

//_____________________________________________________________________
double TYukawaPWF::GetVTr(double r) const
{
  // TYukawaPWF is used for non-relativistic wavefunctions.
  // L=1 wavefunctions have a relativistic origin.
  // GetVTr always returns 0.

  return 0;
}

//_____________________________________________________________________
double TYukawaPWF::GetVSr(double r) const
{
  // TYukawaPWF is used for non-relativistic wavefunctions.
  // L=1 wavefunctions have a relativistic origin.
  // GetVSr always returns 0.

  return 0;
}
double TYukawaPWF::getResidu() const{
  return fC[0]*sqrt(2.*fgHbarc)/TMath::Pi(); //MeV^1/2, coefficient that goes with p^2+m[0]^2, times 2 because of t/(p^2+m[0]^2) jacobian
}