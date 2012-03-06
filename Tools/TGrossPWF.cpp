/*
 * TGrossPWF.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on February,12 2009
 *
 */

////////////////////////////////////////////////////////////////////////
//                                                                    //
// TGrossPWF                                                          //
//                                                                    //
// TGrossPWF implements the abstract interface                        //
// TWavefunctionImplementation.                                       //
//                                                                    //
// This class parametrizes relativistic wavefunctions.                //
// This parametrization was introduced by Buck and Gross              //
// in PRD 20(1979)2361 (section IB).                                  //
//                                                                    //
// TWavefunctionImplementation fixes the following conventions        //
// concerning the units:                                              //
// * position in [fm]                                                 //
// * momentum in [MeV]                                                //
// * r-space wavefunctions in [fm^-1/2]                               //
// * p-space wavefunctions in [MeV^-3/2]                              //
//                                                                    //
// BEGIN_LATEX
// z_{L}(q) = #sum_{n=1}^{N} b_{L,n} g_{L,n}(q)
// END_LATEX
// with                                                               //
// BEGIN_LATEX
// g_{L,n}(q) = f_{L,n} - #sum_{i=1}^{L+2} K^{n}_{L,i} f_{L,N+i}(q)
// END_LATEX
// The functions f in configuration space are given by                //
// BEGIN_LATEX
// * f_{0,i}(r) = e^{-m_{L,i}r}
// * f_{1,i}(r) = e^{-m_{L,i}r} #left(1+#frac{1}{m_{L,i}r}#right)
// * f_{2,i}(r) = e^{-m_{L,i}r} #left(1+#frac{3}{m_{L,i}r}+#frac{3}{m^{2}_{L,i}r^{2}}#right)
// END_LATEX
// In momentum space they are defined as                              //
// BEGIN_LATEX
// * f_{0,i}(p) = #left(#frac{2}{#pi}#right)^{1/2} #frac{1}{p^{2}+m_{L,i}^{2}}
// * f_{1,i}(p) = #left(#frac{2}{#pi}#right)^{1/2} #frac{1}{p^{2}+m_{L,i}^{2}} #left(#frac{p}{m_{L,i}}#right)
// * f_{2,i}(p) = #left(#frac{2}{#pi}#right)^{1/2} #frac{1}{p^{2}+m_{L,i}^{2}} #left(#frac{p}{m_{L,i}}#right)^{2}
// END_LATEX
// Furthermore we have                                                //
// BEGIN_LATEX
// m_{L,i} = #alpha_{L} + (i-1) m_{0}
// K^{n}_{L,i} = #left(#frac{m_{N+i}}{m_{n}}#right)^{L}#prod_{j#neq i}^{L+2} #frac{m_{N+j}^{2}-m_{n}^{2}}{m_{N+j}^{2}-m_{N+i}^{2}}
// END_LATEX
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TGrossPWF.h"
#include <iostream>
#include <TMath.h>
#include "constants.hpp"

using std::cout; using std::cerr; using std::endl;
using TMath::Sqrt; using TMath::Power; using TMath::Exp;

ClassImp(TGrossPWF)

//_____________________________________________________________________
TGrossPWF::TGrossPWF(int order, double *alpha, double m0, 
		     double *parU, double *parW, double *parVT, double *parVS)
: TWavefunctionImplementation(0),
  fOrder(order), fAlpha(new double[3]), fM0(m0), fPar(new double[4*order])
{
  // Constructor
  //
  // order = number of terms N in the parametrization
  // alpha = Parameters BEGIN_LATEX #alpha_{L} END_LATEX in BEGIN_LATEX m_{L,i} END_LATEX formula [MeV]
  // m0 = Parameter BEGIN_LATEX m_{0} END_LATEX in BEGIN_LATEX m_{L,i} END_LATEX formula [MeV]
  // parU,parW,parVT,parVS = arrays containing parameters BEGIN_LATEX b_{L,n} END_LATEX. [MeV^1/2]

  if(!alpha || !parU || !parW || !parVT || !parVS) {
    cerr << "ERROR in TGrossPWF::TGrossPWF(int,double*,double,double* x4): "
	 << "double arrays are NULL.\n";
    exit(1);
  }

  if(fOrder<1) {
    cerr << "ERROR in TGrossPWF::TGrossPWF(int,double*,double,double* x4): "
	 << "order should be >=1!\n";
    exit(1);
  }
  
  for(int i=0; i<3; ++i) fAlpha[i] = alpha[i];

  // all parameters are stored in one continuous array.
  // storage rule: [EWave*fOrder+term]
  for(int i=0; i<fOrder; ++i) fPar[kU*fOrder+i] = parU[i];
  for(int i=0; i<fOrder; ++i) fPar[kW*fOrder+i] = parW[i];
  for(int i=0; i<fOrder; ++i) fPar[kVT*fOrder+i] = parVT[i];
  for(int i=0; i<fOrder; ++i) fPar[kVS*fOrder+i] = parVS[i];
}

//_____________________________________________________________________
TGrossPWF::TGrossPWF(TRootIOCtor*)
  : TWavefunctionImplementation(0), fOrder(0), fAlpha(0), fM0(0), fPar(0)
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
TGrossPWF::TGrossPWF(const TGrossPWF& rhs)
  : TWavefunctionImplementation(rhs),
    fOrder(rhs.fOrder), fAlpha(new double[rhs.fOrder]), 
    fM0(rhs.fM0), fPar(new double[4*rhs.fOrder])
{
  // Copy constructor

  for(int i=0; i<3; ++i) fAlpha[i] = rhs.fAlpha[i];
  for(int i=0; i<4*fOrder; ++i) fPar[i] = rhs.fPar[i];
}

//_____________________________________________________________________
TGrossPWF& TGrossPWF::operator=(const TGrossPWF& rhs)
{
  // Assignment

  if(this!=&rhs) { // avoid self-assignment
    
    TWavefunctionImplementation::operator=(rhs);
    if(fOrder!=rhs.fOrder) delete[] fPar;
    fOrder = rhs.fOrder;
    fM0 = rhs.fM0;
    for(int i=0; i<3; ++i) fAlpha[i] = rhs.fAlpha[i];
    fPar = new double[4*fOrder];
    for(int i=0; i<4*fOrder; ++i) fPar[i] = rhs.fPar[i];
  }

  return *this;
}

//_____________________________________________________________________
TGrossPWF::~TGrossPWF()
{
  // Destructor

  delete[] fAlpha;
  delete[] fPar;
}

//_____________________________________________________________________
double TGrossPWF::GetUr(double r) const
{
  // L=0 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double ur = 0.;

  for(int i=0; i<fOrder; ++i)
    ur += fPar[kU*fOrder+i] * GetG(kR,0,i+1,r);
  
  return ur / Sqrt(fgHbarc);
}

//_____________________________________________________________________
double TGrossPWF::GetWr(double r) const
{
  // L=2 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double wr = 0.;

  for(int i=0; i<fOrder; ++i)
    wr += fPar[kW*fOrder+i] * GetG(kR,2,i+1,r);
  
  return wr / Sqrt(fgHbarc);
}

//_____________________________________________________________________
double TGrossPWF::GetVTr(double r) const
{
  // L=1 triplet wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double vtr = 0.;

  for(int i=0; i<fOrder; ++i)
    vtr += fPar[kVT*fOrder+i] * GetG(kR,1,i+1,r);
  
  return vtr / Sqrt(fgHbarc);
}

//_____________________________________________________________________
double TGrossPWF::GetVSr(double r) const
{
  // L=1 singlet wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double vsr = 0.;

  for(int i=0; i<fOrder; ++i)
    vsr += fPar[kVS*fOrder+i] * GetG(kR,1,i+1,r);
  
  return vsr / Sqrt(fgHbarc);
}

//_____________________________________________________________________
double TGrossPWF::GetUp(double p) const
{
  // L=0 wave in momentum space
  // p = relative momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double up = 0.;

  for(int i=0; i<fOrder; ++i)
    up += fPar[kU*fOrder+i] * GetG(kP,0,i+1,p);
  
  return up;
}

//_____________________________________________________________________
double TGrossPWF::GetWp(double p) const
{
  // L=2 wave in momentum space
  // p = relative momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double wp = 0.;

  for(int i=0; i<fOrder; ++i)
    wp += fPar[kW*fOrder+i] * GetG(kP,2,i+1,p);
  
  return wp;
}

//_____________________________________________________________________
double TGrossPWF::GetVTp(double p) const
{
  // L=1 triplet wave in momentum space
  // p = relative momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double vtp = 0.;

  for(int i=0; i<fOrder; ++i)
    vtp += fPar[kVT*fOrder+i] * GetG(kP,1,i+1,p);
  
  return vtp;
}

//_____________________________________________________________________
double TGrossPWF::GetVSp(double p) const
{
  // L=1 singlet wave in momentum space
  // p = relative momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double vsp = 0.;

  for(int i=0; i<fOrder; ++i)
    vsp += fPar[kVS*fOrder+i] * GetG(kP,1,i+1,p);
  
  return vsp;
}

double TGrossPWF::GetUpoff(const TVector3& p) const{
 
  cerr << "Off-shell wave function not implemented in TGrossPWF wave function class" << endl;
  return 0.;
  
}

double TGrossPWF::GetWpoff1(const TVector3& p) const{
 
  cerr << "Off-shell wave function not implemented in TGrossPWF wave function class" << endl;
  return 0.;
  
}

double TGrossPWF::GetWpoff2(double pperp2) const{
 
  cerr << "Off-shell wave function not implemented in TGrossPWF wave function class" << endl;
  return 0.;
  
}




//_____________________________________________________________________
double TGrossPWF::Mass(int L, int i) const
{
  // 'Mass' formula
  // BEGIN_LATEX
  // m_{L,i} = #alpha_{L} + (i-1) m_{0}
  // END_LATEX
  // Return value in units [MeV]

  return fAlpha[L] + (i-1)*fM0;
}

//_____________________________________________________________________
double TGrossPWF::GetK(int L, int term, int index) const
{
  // BEGIN_LATEX
  // K^{n}_{L,i} = #left(#frac{m_{N+i}}{m_{n}}#right)^{L}#prod_{j#neq i}^{L+2} #frac{m_{N+j}^{2}-m_{n}^{2}}{m_{N+j}^{2}-m_{N+i}^{2}}
  // END_LATEX
  //
  // L = orbital momentum
  // term = n
  // index = i
  // N = fOrder
  // 
  // Return value is dimensionless
  
  double K = Power(Mass(L,fOrder+index)/Mass(L,term),L);

  for(int j=1; j<=(L+2); ++j)
    if(j!=index)
      K *= (Mass(L,fOrder+j)*Mass(L,fOrder+j)-Mass(L,term)*Mass(L,term)) 
	/ (Mass(L,fOrder+j)*Mass(L,fOrder+j)-Mass(L,fOrder+index)*Mass(L,fOrder+index));

  return K;
}

//_____________________________________________________________________
double TGrossPWF::GetG(ESpace space, int L, int term, double q) const
{
  // BEGIN_LATEX
  // g_{L,n}(q) = f_{L,n} - #sum_{i=1}^{L+2} K^{n}_{L,i} f_{L,N+i}(q)
  // END_LATEX
  //
  // space = r-space or p-space
  // q in [fm] or [MeV]
  // L = orbital momentum
  // term = n
  // 
  // return value is dimensionless (r-space) or [MeV^-2] (p-space)

  double g = GetF(space,L,term,q);

  for(int i=1; i<=(L+2); ++i)
    g -= GetK(L,term,i) * GetF(space,L,fOrder+i,q);

  return g;
}

//_____________________________________________________________________
double TGrossPWF::GetF(ESpace space, int L, int index, double q) const
{
  // In r-space:
  // BEGIN_LATEX
  // * f_{0,i}(r) = e^{-m_{L,i}r}
  // * f_{1,i}(r) = e^{-m_{L,i}r} #left(1+#frac{1}{m_{L,i}r}#right)
  // * f_{2,i}(r) = e^{-m_{L,i}r} #left(1+#frac{3}{m_{L,i}r}+#frac{3}{m^{2}_{L,i}r^{2}}#right)
  // END_LATEX
  // In p-space:
  // BEGIN_LATEX
  // * f_{L,i}(p) = #left(#frac{2}{#pi}#right)^{1/2} #frac{1}{p^{2}+m_{L,i}^{2}} #left(#frac{p}{m_{L,i}}#right)^{L}
  // END_LATEX
  // 
  // space = r-space or p-space
  // q in [fm] or [MeV]
  // L = orbital momentum
  // index = i
  //
  // return value is dimensionless (r-space) or [MeV^-2] (p-space)

  double f = 0.;

  if(space==kR) { // r-space

    if(q<STRANGEUFLOW) // f(r) does not converge nicely
      return 0.;

    q /= fgHbarc; // convert position [fm] -> [MeV^-1]
    
    switch(L) {
    case 0:
      f = Exp(-1.*Mass(L,index)*q);
      break;
    case 1:
      f =  Exp(-1.*Mass(L,index)*q) *(1.+1./Mass(L,index)/q);
      break;
    case 2:
      f =  Exp(-1.*Mass(L,index)*q) 
	*(1.+3./Mass(L,index)/q+3./Power(Mass(L,index)*q,2.));
      break;
    default:
      cerr << "ERROR in TGrossPWF::GetF(ESpace,int,int,double): "
	   << "L=" << L << " wave unknown.\n";
      return 0;
    }
  } // end r-space
  else if(space==kP) { // p-space

    f = Sqrt(2./TMath::Pi())/(q*q + Mass(L,index)*Mass(L,index));
    switch(L) {
    case 0:
      break;
    case 1:
      f *= q/Mass(L,index);
      break;
    case 2:
      f *= Power(q/Mass(L,index),2.);
      break;
    default:
      cerr << "ERROR in TGrossPWF::GetF(ESpace,int,int,double): "
	   << "L=" << L << " wave unknown.\n";
      return 0;
    }
   } // end p-space

  return f;
}
