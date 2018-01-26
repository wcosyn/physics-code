/*
 * TCstPWF.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on July,23 2010
 *
 */

////////////////////////////////////////////////////////////////////////
//                                                                    //
// TCstPWF                                                            //
//                                                                    //
// TCstPWF implements the abstract interface 
// TWavefunctionImplementation.
//
// This class parametrizes relativistic wavefunctions.
// This parametrization was introduced by Gross and Stadler
// in arXiv:1007.0778 (section IIIG).
//
// TWavefunctionImplementation fixes the following conventions
// concerning the units:
// * position in [fm]
// * momentum in [MeV]
// * r-space wavefunctions in [fm^-1/2]
// * p-space wavefunctions in [MeV^-3/2]
//
// BEGIN_LATEX
// z_{L}(q) = #sum_{i=1}^{N} b_{z,i} G^{z}_{L,i}(q)
// END_LATEX
// with q either the relative momentum p or relative position r.
//
// The functions G in momentum space are given by
// BEGIN_LATEX
// * G^{z}_{L,i}(p) = #sqrt{#frac{2}{#pi}} #frac{p^{L}m_{z,i}^{2}M_{z,i}^{2n_{z}-L}}{#left(m_{z,i}^{2}+p^{2}#right)#left(M_{z,i}^{2}+p^{2}#right)^{n_{z}}}
// * G^{z}_{L,N}(p) = #sqrt{#frac{2}{#pi}} #frac{p^{L}M_{z,N}^{2n_{z}+1-L}}{#left(M_{z,N}^{2}+p^{2}#right)^{n_{z}+1/2}}
// END_LATEX
// In configuration space they are defined as
// BEGIN_LATEX
// * G^{z}_{0,i}(r) = A_{z,i} #left[e^{-z_{z,i}}-e^{-Z_{z,i}}#left(1+#frac{Z_{z,i}}{2}#left(1-R_{z,i}^{2}#right)#right)#right]
// * G^{z}_{1,i}(r) = A_{z,i} #left[R_{z,i}e^{-z_{z,i}}#left(1+#frac{1}{z_{z,i}}#right)-e^{-Z_{z,i}}#left(1+#frac{1}{z_{z,i}}+#frac{Z_{z,i}}{2}#left(1-R_{z,i}^{2}#right)#right)#right]
// * G^{z}_{2,i}(r) = B_{z,i} #left[R_{z,i}^{2}e^{-z_{z,i}}#left(1+#frac{3}{z_{z,i}}+#frac{3}{z_{z,i}^{2}}#right)-e^{-Z_{z,i}}#left(1+#frac{3}{z_{z,i}}+#frac{3}{z_{z,i}^{2}}+#frac{1+Z_{z,i}}{2}#left(1-R_{z,i}^{2}#right)+#frac{Z_{z,i}^{2}}{8}#left(1-R_{z,i}^{2}#right)^{2}#right)#right]
// END_LATEX
// and 
// BEGIN_LATEX
// * G^{z}_{0,N} = #frac{2}{3#pi} M_{z,N}^{2} Z_{z,N}^{2} K_{1}(Z_{z,n})
// * G^{z}_{1,N} = #frac{2}{3#pi} M_{z,N}^{2} Z_{z,N}^{2} K_{0}(Z_{z,n})
// * G^{z}_{2,N} = #frac{2}{15#pi} M_{z,N}^{2} Z_{z,N}^{3} K_{0}(Z_{z,n})
// END_LATEX
// with BEGIN_LATEX K_{i}(x) END_LATEX the modified Bessel functions of the 2nd kind.
// In the previous expressions, we have defined
// BEGIN_LATEX
// m_{z,i} = #alpha_{z} + (i-1) m^{z}_{x}
// M_{z,i} = m_{z,i+1}
// z_{z,i} = m_{z,i} r
// Z_{z,i} = M_{z,i} r
// R_{z,i} = #frac{m_{z,i}}{M_{z,i}}
// A_{z,i} = #frac{m_{z,i}^{2}M_{z,i}^{4}}{#left(M_{z,i}^{2}-m_{z,i}^{2}#right)^{2}}
// B_{z,i} = A_{z,i} #frac{M_{z,i}^{2}}{M_{z,i}^{2}-m_{z,i}^{2}}
// END_LATEX
//
//
// The parameters for each wave function component are:
// BEGIN_LATEX N, b_{z,i}, n_{z}, #alpha_{z}, m^{z}_{x} and M_{z,N} END_LATEX
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TCstPWF.h"
#include <iostream>
#include <TMath.h>
#include <gsl/gsl_sf_bessel.h>
#include "constants.hpp"
#include <cassert>
#include <cmath>
using namespace std;

ClassImp(TCstPWF)

//_____________________________________________________________________
TCstPWF::TCstPWF(int *order, double *leadingMass, double *stepMass, double *tailMass,
		 int *exponent, double *parU, double *parW, 
		 double *parVT, double *parVS)
: TWavefunctionImplementation(0), fOrder(new int[4]), fCoeff(new double*[4]),
  fAlpha(new double[4]), fMx(new double[4]), fMn(new double[4]), fExp(new int[4])
{
  // Constructor
  //
  // Users should provide the parameters of the parametrization. First parameters as arrays of size 4 (u,w,vt,vs): 
  // * BEGIN_LATEX N END_LATEX
  //   order = number of terms in the parametrization
  // * BEGIN_LATEX n_{z} END_LATEX
  //   exponent = parameter in the expansion functions G in momentum space
  // * BEGIN_LATEX #alpha_{z} END_LATEX
  //   leadingMass [MeV] = asymptotic behavior at large distances
  // * BEGIN_LATEX m^{z}_{x} END_LATEX
  //   stepMass [MeV]
  // * BEGIN_LATEX M_{z,N} END_LATEX
  //   tailMass [MeV]
  //
  // Finally 4 arrays of a size compatible with the order parameter given supra containing the expansion parameters:
  // * BEGIN_LATEX b_{z,i} END_LATEX
  //   parU,parW,parVT,parVS [MeV^-3/2] = arrays containing expansion parameters

  // Check memory allocation
  assert(fOrder);
  assert(fCoeff);
  assert(fAlpha);
  assert(fMx);
  assert(fMn);
  assert(fExp);

  // Check the arguments
  if(!order || !parU || !parW || !parVT || !parVS || 
     !leadingMass || !stepMass || !tailMass || !exponent) {
    cerr << "ERROR in TCstPWF::TCstPWF(int*,double* x8): "
	 << "arrays are NULL.\n";
    exit(1);
  }
    
  for(int i=0; i<4; ++i) {
    // Set the number of terms in the parametrization
    fOrder[i] = order[i];
    if(fOrder[i]<2) {
      cerr << "ERROR in TCstPWF::TCstPWF(int*,double* x8): "
	   << "order should be >=2!\n";
      exit(1);
    }
    // Set the masses and exponents
    fAlpha[i] = leadingMass[i];
    fMx[i] = stepMass[i];
    fMn[i] = tailMass[i];
    fExp[i] = exponent[i];
    // allocate memory for coefficients
    fCoeff[i] = new double[fOrder[i]];
    assert(fCoeff[i]);
  }

  // store the expansion coefficients  
  for(int i=0; i<fOrder[kU]; ++i) fCoeff[kU][i] = parU[i];
  for(int i=0; i<fOrder[kW]; ++i) fCoeff[kW][i] = parW[i];
  for(int i=0; i<fOrder[kVT]; ++i) fCoeff[kVT][i] = parVT[i];
  for(int i=0; i<fOrder[kVS]; ++i) fCoeff[kVS][i] = parVS[i];
}

//_____________________________________________________________________
TCstPWF::TCstPWF(TRootIOCtor*)
  : TWavefunctionImplementation(0), fOrder(0), fCoeff(0), fAlpha(0), fMx(0), fMn(0),
    fExp(0)
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
TCstPWF::TCstPWF(const TCstPWF& rhs)
  : TWavefunctionImplementation(rhs), fOrder(new int[4]), fCoeff(new double*[4]),
    fAlpha(new double[4]), fMx(new double[4]), fMn(new double[4]), 
    fExp(new int[4])    
{
  // Copy constructor
  
  // Check memory allocation
  assert(fOrder);
  assert(fCoeff);
  assert(fAlpha);
  assert(fMx);
  assert(fMn);
  assert(fExp);

  for(int i=0; i<4; ++i) {
    fOrder[i] = rhs.fOrder[i];
    fAlpha[i] = rhs.fAlpha[i];
    fMx[i] = rhs.fMx[i];
    fMn[i] = rhs.fMn[i];
    fExp[i] = rhs.fExp[i];
    fCoeff[i] = new double[fOrder[i]];
    assert(fCoeff[i]);
    for(int j=0; j<fOrder[i]; ++j) fCoeff[i][j] = rhs.fCoeff[i][j];
  }
}

//_____________________________________________________________________
TCstPWF& TCstPWF::operator=(const TCstPWF& rhs)
{
  // Assignment

  if(this!=&rhs) { // avoid self-assignment
    
    TWavefunctionImplementation::operator=(rhs);
    for(int i=0; i<4; ++i) {
      if( fOrder[i] != rhs.fOrder[i] ) {
	delete[] fCoeff[i];
	fCoeff[i] = new double[rhs.fOrder[i]];
	assert(fCoeff[i]);
      }
      fOrder[i] = rhs.fOrder[i];
      fAlpha[i] = rhs.fAlpha[i];
      fMx[i] = rhs.fMx[i];
      fMn[i] = rhs.fMn[i];
      fExp[i] = rhs.fExp[i];
      for(int j=0; j<fOrder[i]; ++j) fCoeff[i][j] = rhs.fCoeff[i][j];
    }
  }

  return *this;
}

//_____________________________________________________________________
TCstPWF::~TCstPWF()
{
  // Destructor

  delete[] fExp;
  delete[] fMn;
  delete[] fMx;
  delete[] fAlpha;
  for(int i=4; i>0; --i) delete[] fCoeff[i-1];
  delete[] fCoeff;
  delete[] fOrder;
}

//_____________________________________________________________________
double TCstPWF::GetUr(double r) const
{
  // L=0 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double ur = 0.;

  for(int i=1; i<=fOrder[kU]; ++i)
    ur += fCoeff[kU][i-1] * GetG(kR,kU,i,r);

  return ur / sqrt(fgHbarc);
}

//_____________________________________________________________________
double TCstPWF::GetWr(double r) const
{
  // L=2 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double wr = 0.;

  for(int i=1; i<=fOrder[kW]; ++i)
    wr += fCoeff[kW][i-1] * GetG(kR,kW,i,r);
  
  return wr / sqrt(fgHbarc);
}

//_____________________________________________________________________
double TCstPWF::GetVTr(double r) const
{
  // L=1 triplet wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double vtr = 0.;

  for(int i=1; i<=fOrder[kVT]; ++i)
    vtr += fCoeff[kVT][i-1] * GetG(kR,kVT,i,r);
  
  return vtr / sqrt(fgHbarc);
}

//_____________________________________________________________________
double TCstPWF::GetVSr(double r) const
{
  // L=1 singlet wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  double vsr = 0.;

  for(int i=1; i<=fOrder[kVS]; ++i)
    vsr += fCoeff[kVS][i-1] * GetG(kR,kVS,i,r);
    
  return vsr / sqrt(fgHbarc);
}

//_____________________________________________________________________
double TCstPWF::GetUp(double p) const
{
  // L=0 wave in momentum space
  // p = relative momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double up = 0.;

  for(int i=1; i<=fOrder[kU]; ++i)
    up += fCoeff[kU][i-1] * GetG(kP,kU,i,p);
  
  return up;
}

//_____________________________________________________________________
double TCstPWF::GetWp(double p) const
{
  // L=2 wave in momentum space
  // p = relative momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double wp = 0.;

  for(int i=1; i<=fOrder[kW]; ++i)
    wp += fCoeff[kW][i-1] * GetG(kP,kW,i,p);
  
  return wp;
}

//_____________________________________________________________________
double TCstPWF::GetVTp(double p) const
{
  // L=1 triplet wave in momentum space
  // p = relative momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double vtp = 0.;

  for(int i=1; i<=fOrder[kVT]; ++i)
    vtp += fCoeff[kVT][i-1] * GetG(kP,kVT,i,p);
  
  return vtp;
}

//_____________________________________________________________________
double TCstPWF::GetVSp(double p) const
{
  // L=1 singlet wave in momentum space
  // p = relative momentum in [MeV]
  // Return value in units [MeV^-3/2]

  double vsp = 0.;

  for(int i=1; i<=fOrder[kVS]; ++i)
    vsp += fCoeff[kVS][i-1] * GetG(kP,kVS,i,p);
  
  return vsp;
}

double TCstPWF::GetUpoff(const TVector3& p) const{
 
  cerr << "Off-shell wave function not implemented in TCstPWF wave function class" << endl;
  return 0.;
  
}

double TCstPWF::GetWpoff1(const TVector3& p) const{
 
  cerr << "Off-shell wave function not implemented in TCstPWF wave function class" << endl;
  return 0.;
  
}

double TCstPWF::GetWpoff2(double pperp2) const{
 
  cerr << "Off-shell wave function not implemented in TCstPWF wave function class" << endl;
  return 0.;
  
}




//_____________________________________________________________________
double TCstPWF::Mass(EWave wave, int i) const
{
  // 'Mass' formula
  // BEGIN_LATEX
  // m_{z,i} = #alpha_{z} + (i-1) m^{z}_{x}
  // END_LATEX
  // Return value in units [MeV]

  return fAlpha[wave] + (i-1)*fMx[wave];
}

//_____________________________________________________________________
int TCstPWF::AngularMomentum(EWave wave)
{
  // Returns the angular momentum of each wave:
  //   wave |  L
  //  ------|------
  //    u   |  0
  //    v_t |  1
  //    v_s |  1
  //    w   |  2

  switch(wave) {
  case kU: return 0; break;
  case kW: return 2; break;
  default: return 1; break;
  }
  return -1;
}

//_____________________________________________________________________
double TCstPWF::GetG(ESpace space, EWave wave, int term, double q) const
{
  // In momentum space:
  // BEGIN_LATEX
  // * G^{z}_{L,i}(p) = #sqrt{#frac{2}{#pi}} #frac{p^{L}m_{z,i}^{2}M_{z,i}^{2n_{z}-L}}{#left(m_{z,i}^{2}+p^{2}#right)#left(M_{z,i}^{2}+p^{2}#right)^{n_{z}}}
  // * G^{z}_{L,N}(p) = #sqrt{#frac{2}{#pi}} #frac{p^{L}M_{z,N}^{2n_{z}+1-L}}{#left(M_{z,N}^{2}+p^{2}#right)^{n_{z}+1/2}}
  // END_LATEX
  // In configuration space:
  // BEGIN_LATEX
  // * G^{z}_{0,i}(r) = A_{z,i} #left[e^{-z_{z,i}}-e^{-Z_{z,i}}#left(1+#frac{Z_{z,i}}{2}#left(1-R_{z,i}^{2}#right)#right)#right]
  // * G^{z}_{1,i}(r) = A_{z,i} #left[R_{z,i}e^{-z_{z,i}}#left(1+#frac{1}{z_{z,i}}#right)-e^{-Z_{z,i}}#left(1+#frac{1}{z_{z,i}}+#frac{Z_{z,i}}{2}#left(1-R_{z,i}^{2}#right)#right)#right]
  // * G^{z}_{2,i}(r) = B_{z,i} #left[R_{z,i}^{2}e^{-z_{z,i}}#left(1+#frac{3}{z_{z,i}}+#frac{3}{z_{z,i}^{2}}#right)-e^{-Z_{z,i}}#left(1+#frac{3}{z_{z,i}}+#frac{3}{z_{z,i}^{2}}+#frac{1+Z_{z,i}}{2}#left(1-R_{z,i}^{2}#right)+#frac{Z_{z,i}^{2}}{8}#left(1-R_{z,i}^{2}#right)^{2}#right)#right]
  // END_LATEX
  // and 
  // BEGIN_LATEX
  // * G^{z}_{0,N} = #frac{2}{3#pi} M_{z,N}^{2} Z_{z,N}^{2} K_{1}(Z_{z,n})
  // * G^{z}_{1,N} = #frac{2}{3#pi} M_{z,N}^{2} Z_{z,N}^{2} K_{0}(Z_{z,n})
  // * G^{z}_{2,N} = #frac{2}{15#pi} M_{z,N}^{2} Z_{z,N}^{3} K_{0}(Z_{z,n})
  // END_LATEX
  // with BEGIN_LATEX K_{i}(x) END_LATEX the modified Bessel functions of the 2nd kind
  //
  // space = r-space or p-space
  // wave = u,w,vt,vs
  // term = n
  // q in [fm] or [MeV]
  // 
  // return value is dimensionless (p-space) or [MeV^2] (r-space)

  double g = 0.;

  if( (term>0) && (term<fOrder[wave]) ) {
    if( space==kP ) {
      g = sqrt(2./TMath::Pi())
	* pow(q,AngularMomentum(wave))
	* Mass(wave,term)*Mass(wave,term)
	* pow(Mass(wave,term+1),2*fExp[wave]-AngularMomentum(wave))
	/ (Mass(wave,term)*Mass(wave,term) + q*q)
	/ pow( Mass(wave,term+1)*Mass(wave,term+1) + q*q, fExp[wave]);
    } 

    else {
      if(q<STRANGEUFLOW) // G(r) does not converge nicely
	return 0.;
      
      q /= fgHbarc; // convert position [fm] -> [MeV^-1]
      
      double A = Mass(wave,term)*Mass(wave,term)
	*Mass(wave,term+1)*Mass(wave,term+1)*Mass(wave,term+1)*Mass(wave,term+1)
	/ pow(Mass(wave,term+1)*Mass(wave,term+1)-Mass(wave,term)*Mass(wave,term),2.);
      double B = A
	*Mass(wave,term+1)*Mass(wave,term+1)
	/(Mass(wave,term+1)*Mass(wave,term+1)-Mass(wave,term)*Mass(wave,term));
      double z = Mass(wave,term) * q;
      double Z = Mass(wave,term+1) * q;
      double R = Mass(wave,term) / Mass(wave,term+1);

      if( wave==kU )
	g = A * ( exp(-z) - exp(-Z)*(1.+Z/2.*(1.-R*R)) );
      else if( wave==kW )
	g = B * 
	  ( R*R*exp(-z)*(1.+3./z+3./z/z) 
	    -exp(-Z)*(1.+3./Z+3./Z/Z+(1.+Z)/2.*(1.-R*R)+Z*Z/8.*pow(1.-R*R,2.)));
      else
	g = A * ( R*exp(-z)*(1.+1./z) - exp(-Z)*(1.+1./Z+Z/2.*(1.-R*R)) );
    } // else space==kR      
  }
  else if( term == fOrder[wave] ) {
    if( space==kP ) {
      g = sqrt(2./TMath::Pi())
	* pow(q,AngularMomentum(wave))
	* pow(fMn[wave],2*fExp[wave]+1-AngularMomentum(wave))
	/ pow( fMn[wave]*fMn[wave] + q*q, fExp[wave]+.5);
    }

    else {
      if(q<STRANGEUFLOW) // G(r) does not converge nicely
	return 0.;
      
      q /= fgHbarc; // convert position [fm] -> [MeV^-1]
      
      double Z = fMn[wave] * q;

      if( wave==kU )
	g = 2./3./TMath::Pi() * fMn[wave]*fMn[wave] *Z*Z *gsl_sf_bessel_k1_scaled(Z)*exp(-Z);
      else if( wave==kW )
	g = 2./15./TMath::Pi() * fMn[wave]*fMn[wave]*Z*Z*Z*gsl_sf_bessel_k0_scaled(Z)*exp(-Z);
      else
	g = 2./3./TMath::Pi() * fMn[wave]*fMn[wave] *Z*Z *gsl_sf_bessel_k0_scaled(Z)*exp(-Z);
    } // else space==kR
  }
  else {
    cerr << "ERROR in TCstPWF::GetG(ESpace,EWave,int,double): "
	 << "invalid index." << endl;
    exit(1);
  }

  return g;
}

//_____________________________________________________________________
void TCstPWF::Streamer(TBuffer &R__b){}
// {
//    // Stream an object of class TCstPWF.

//    UInt_t R__s, R__c;
//    if (R__b.IsReading()) {
//       Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
//       TWavefunctionImplementation::Streamer(R__b);
//       delete [] fOrder;
//       fOrder = new int[4];
//       R__b.ReadFastArray(fOrder,4);
//       delete [] fCoeff;
//       fCoeff = new double*[4];
//       for(int i=0; i<4; ++i) {
// 	fCoeff[i] = new double[fOrder[i]];
// 	R__b.ReadFastArray(fCoeff[i],fOrder[i]);
//       }
//       delete [] fAlpha;
//       fAlpha = new double[4];
//       R__b.ReadFastArray(fAlpha,4);
//       delete [] fMx;
//       fMx = new double[4];
//       R__b.ReadFastArray(fMx,4);
//       delete [] fMn;
//       fMn = new double[4];
//       R__b.ReadFastArray(fMn,4);
//       delete [] fExp;
//       fExp = new int[4];
//       R__b.ReadFastArray(fExp,4);
//       R__b.CheckByteCount(R__s, R__c, TCstPWF::IsA());
//    } else {
//       R__c = R__b.WriteVersion(TCstPWF::IsA(), kTRUE);
//       TWavefunctionImplementation::Streamer(R__b);
//       R__b.WriteFastArray(fOrder,4);
//       for(int i=0; i<4; ++i) R__b.WriteFastArray(fCoeff[i],fOrder[i]);
//       R__b.WriteFastArray(fAlpha,4);
//       R__b.WriteFastArray(fMx,4);
//       R__b.WriteFastArray(fMn,4);
//       R__b.WriteFastArray(fExp,4);
//       R__b.SetByteCount(R__c, kTRUE);
//    }
// }
