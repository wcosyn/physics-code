/*
 * TSpinor.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: August,10 2009
 *
 */

#include "TSpinor.h"
#include "TLorentzQuaternion.h"
#include "constants.hpp"
//#include <numtoa.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TRotation.h>
#include "FourVector.h"
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <climits>
#include <iostream>
#include <complex>

using TMath::Sign;
using namespace std;

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TSpinor                                                               //
//                                                                       //
// Dirac spinor for on-mass shell spin-1/2 particles                     //
//
//BEGIN_LATEX
//u(#vec{p},s) = #sqrt{E_{p}+m}#left( #splitline{1}{#frac{#vec{#sigma}#upoint#vec{p}}{E_{p}+m}} #right) #chi^{s},
//END_LATEX
// with BEGIN_LATEX E_{p} = #sqrt{#vec{p}^{2}+m^{2}} END_LATEX.
// See documentation of TSpinor::Polarization::GetPauliSpinor for the 
// definition of the two-component spinors BEGIN_LATEX #chi^{s} END_LATEX.
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TSpinor)

//_____________________________________________________________________
const FourVector< Matrix<2,2> > TSpinor::kSigmaPauli
(Matrix<2,2>(1.0,0.0,0.0,1.0),        // t-component
 Matrix<2,2>(0.0,1.0,1.0,0.0),        // x-component
 Matrix<2,2>(0.0,
	     complex<double>(0,-1.0),
	     complex<double>(0,1.0),
	     0.0),                    // y-component
 Matrix<2,2>(1.0,0.0,0.0,-1.0));      // z-component

//_____________________________________________________________________
TSpinor::TSpinor(double px,double py,double pz,double mass, const Polarization& pol,
		 Normalization norm, int energyState)
  : TObject(), fComponent(new Matrix<4,1>)
{
  // Protected Constructor
  // Allows to initialize the spinor for negative energy states. Needed by 
  // TVSpinor.
  //
  // (px,py,pz) is the momentum of the particle with mass 'mass'
  // It's polarization state is given by 'pol'
  // 'norm' defines the normalization convention of the spinors
  // 'energyState' specifies whether you want positive or negative energy solutions.
  assert(fComponent); // check allocation
  InitializeSpinor(px,py,pz,mass,pol,norm,energyState);
}

//_____________________________________________________________________
TSpinor::TSpinor(const TLorentzVector& v, double mass, const Polarization& pol, 
		 Normalization norm)
  : TObject(), fComponent(new Matrix<4,1>)
{
  // Constructor
  //
  // 'v' is the 4momentum of the particle with mass 'mass'. The particle
  // needn't be on-mass shell.
  // It's polarization state is given by 'pol'
  // 'norm' defines the normalization convention of the spinors
  assert(fComponent); // check allocation
  InitializeSpinor(v.X(),v.Y(),v.Z(),mass,pol,norm,+1);
}

//_____________________________________________________________________
TSpinor::TSpinor(const FourVector<double>& v, double mass, const Polarization& pol,
		 Normalization norm)
  : TObject(), fComponent(new Matrix<4,1>)
{
  // Constructor
  //
  // 'v' is the 4momentum of the particle with mass 'mass'. The particle
  // needn't be on-mass shell.
  // It's polarization state is given by 'pol'
  // 'norm' defines the normalization convention of the spinors
  assert(fComponent); // check allocation
  InitializeSpinor(v[1],v[2],v[3],mass,pol,norm,+1);
}

//_____________________________________________________________________
TSpinor::TSpinor(TRootIOCtor* rio)
  : TObject(), fComponent(0)
{
  // ROOT I/0 Constructor
}

//_____________________________________________________________________
TSpinor::TSpinor(const TSpinor& rhs)
  : fComponent(new Matrix<4,1>(*rhs.fComponent))
{
  // Copy constructor

  assert(fComponent); // check allocation
}

//_____________________________________________________________________
TSpinor::~TSpinor()
{
  // Destructor
  delete fComponent;
}

//_____________________________________________________________________
TSpinor& TSpinor::operator=(const TSpinor& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
    TObject::operator=(rhs);
    *fComponent=*rhs.fComponent;
  }

  return *this;
}

//_____________________________________________________________________
TSpinor* TSpinor::Clone(const char*) const
{
  // Virtual copy constructor
  return new TSpinor(*this);
}

//______________________________________________________________________________
void TSpinor::Streamer(TBuffer &R__b){}
// {
//   // Stream an object of class TSpinor.

//   UInt_t R__s, R__c;
  
//   if (R__b.IsReading()) { // begin reading
//     Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
//     TObject::Streamer(R__b);
//     fComponent = new Matrix<4,1>;
//     double realComponent, imagComponent;
//     for(int i=0; i<4; ++i) {
//       R__b >> realComponent;
//       R__b >> imagComponent;
//       (*fComponent)(i,0) = complex<double>(realComponent,imagComponent);
//     }
//     if( R__v<3 ) {
//       cerr << "ERROR in TSpinor::Streamer(TBuffer&): "
// 	   << "Old version of TSpinor might contain errors.\n";
//       exit(1);
//     }

//     R__b.CheckByteCount(R__s, R__c, TSpinor::IsA());
//   } // end reading 

//   else { // begin writing
//     R__c = R__b.WriteVersion(TSpinor::IsA(), kTRUE);
//     TObject::Streamer(R__b);
//     for(int i=0; i<4; ++i) {
//       R__b << (*fComponent)(i,0).real();
//       R__b << (*fComponent)(i,0).imag();
//     }
//     R__b.SetByteCount(R__c, kTRUE);
//   } // end writing
// }

//______________________________________________________________________________
void TSpinor::InitializeSpinor(double px, double py, double pz, double mass,
			       const Polarization& pol, Normalization norm,
			       int energyState)
{
  // Initialize the spinor following the conventions of Bjorken&Drell and Gross.
  // See eq.(A.15) on p.596 of 
  // Gross,'Relativistic Quantum Mechanics and Field Theory'
  //
  // The sign of energyState determines whether we calculate the positive
  // or negative energy solution.

  Matrix<2,2> pSigma //  = vect(p) * vect(sigma)
    = px * kSigmaPauli[1]
    + py * kSigmaPauli[2] 
    + pz * kSigmaPauli[3];
  double frontFactor // = sqrt(Ep+m)
    = sqrt(sqrt(px*px+py*py+pz*pz+mass*mass)+mass); 
  Matrix<2,1> temp; // an intermediate 2x1 matrix
  
  // the first two components of the spinor (large component)
  temp = frontFactor * pol.GetPauliSpinor(energyState);
  if(energyState>0) {
    (*fComponent)(0,0) = temp(0,0);
    (*fComponent)(1,0) = temp(1,0);
  } else {
    (*fComponent)(2,0) = temp(0,0);
    (*fComponent)(3,0) = temp(1,0);
  }
  
  // the two last components of the spinor (small component)
  temp = (pSigma * pol.GetPauliSpinor(energyState)) *(1./frontFactor);
  if(energyState>0) {
    (*fComponent)(2,0) = temp(0,0);
    (*fComponent)(3,0) = temp(1,0);
  } else {
    (*fComponent)(0,0) = temp(0,0);
    (*fComponent)(1,0) = temp(1,0);
  }

  // Finally make sure the spinor is normalized correctly
  if( norm==kUnity) fComponent->operator*=( 1./sqrt(2.*mass) );
}

//______________________________________________________________________________
Matrix<1,4> TSpinor::Bar(const TSpinor& spinor)
{
  // BEGIN_LATEX
  // Bar(u) = #bar{u} = u^{T}#gamma_{0}
  // END_LATEX
  return (((Matrix<4,1>)spinor).H())*G0;
}

//______________________________________________________________________________
ostream& operator<<(ostream& output, const TSpinor& rhs) 
{
  return output << static_cast<const Matrix<4,1>&>(rhs);
}

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TSpinor::Polarization                                                 //
//                                                                       //
// Polarization state of a spin-1/2 Dirac particle                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TSpinor::Polarization)

//______________________________________________________________________________
TSpinor::Polarization::Polarization(double theta, double phi,
				    State polState)
  : fTheta(theta), fPhi(phi), fPolarizationState(polState)
{
  // Constructor

  SetTheta(theta);
  SetPhi(phi);
}

//______________________________________________________________________________
TSpinor::Polarization::Polarization(TRootIOCtor *rio)
  : TObject(), fTheta(0.), fPhi(0.), fPolarizationState(kUp)
{
  // ROOT I/O Constructor
}

//______________________________________________________________________________
TSpinor::Polarization::Polarization(const Polarization& rhs)
  : fTheta(rhs.fTheta), fPhi(rhs.fPhi), fPolarizationState(rhs.fPolarizationState)
{
  // Copy constructor
}

//______________________________________________________________________________
TSpinor::Polarization::~Polarization()
{
  // Destructor
}

//______________________________________________________________________________
TSpinor::Polarization& TSpinor::Polarization::operator=(const Polarization& rhs)
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
bool TSpinor::Polarization::operator==(const Polarization& rhs)
{
  // Comparison

  return 
    fabs(fTheta-rhs.fTheta)<STRANGEUFLOW &&
    fabs(fPhi-rhs.fPhi)<STRANGEUFLOW &&
    fPolarizationState==rhs.fPolarizationState;
}

//______________________________________________________________________________
bool TSpinor::Polarization::operator==(State polarizationState)
{
  // Compare states
  return fPolarizationState==polarizationState;
}

//______________________________________________________________________________
void TSpinor::Polarization::SetTheta(double theta)
{
  // Set the polar angle of the spin quantization axis (in rad)

  if( theta<0. || theta>TMath::Pi() ) {
    cerr << "ERROR in TSpinor::Polarization::SetTheta(double): "
	 << "polar angle not in range [0,pi].\n";
    exit(1);
  }

  fTheta = theta;
}

//______________________________________________________________________________
void TSpinor::Polarization::SetPhi(double phi)
{
  // Set the azimuthal angle of the spin quantization axis (in rad)
  fPhi = phi;
}

//______________________________________________________________________________
TSpinor::Polarization& TSpinor::Polarization::operator++()
{
  // Move the polarization state one-up like this ++spinorPolarization

  if(fPolarizationState==kDown) fPolarizationState=kUp;

  return *this;
}

//______________________________________________________________________________
TSpinor::Polarization& TSpinor::Polarization::operator--()
{
  // Move the polarization state one-down like this --spinorPolarization

  if(fPolarizationState==kUp) fPolarizationState=kDown;

  return *this;
}

//______________________________________________________________________________
TSpinor::Polarization& TSpinor::Polarization::Invert()
{
  // Inverts the polarization state
  if(fPolarizationState==kUp) fPolarizationState=kDown;
  else fPolarizationState=kUp;

  return *this;
}

//______________________________________________________________________________
TSpinor::Polarization TSpinor::Polarization::Inverse() const
{
  // Returns an object with the inverse polarization state
  Polarization inverse = *this;

  return inverse.Invert();
}

//______________________________________________________________________________
Matrix<2,1> TSpinor::Polarization::GetPauliSpinor(int eSolution) const
{
  // Determine the 2-component spinor which depends upon the 
  // polarization state and the orientation of the quantization 
  // axis defined by theta and phi.
  // if oriented along z-axis this spinor is either (1,0) or (0,1) 
  // depending on the polarization state.
  // For other polarization axes, the spinor is obtained through combinations
  // with Wigner D-matrices. We find
  //BEGIN_LATEX #chi^{1/2} = #left(#splitline{cos#frac{#theta}{2}}{e^{i#phi}sin#frac{#theta}{2}}#right) END_LATEX ,  BEGIN_LATEX #chi^{-1/2} = #left(#splitline{-e^{-i#phi}sin#frac{#theta}{2}}{cos#frac{#theta}{2}}#right) END_LATEX
  //
  // The pauli spinor for particles and anti-particles is different. This is
  // indicated by the sign of the argument. We define
  // BEGIN_LATEX #eta^{s} = #chi^{-s} END_LATEX
 
  Matrix<2,1> pauliSpinor;

  // Positive energy solutions (particles)
  if(eSolution>0) {
    // polarization = -1/2
    if( fPolarizationState==kDown )
      {
	pauliSpinor(0,0) = -sin(fTheta/2.0)
	  *exp(complex<double>(0,-fPhi));
	pauliSpinor(1,0) = cos(fTheta/2.0);
      }
    // polarization = +1/2
    else if( fPolarizationState==kUp )
      {
	pauliSpinor(0,0) = cos(fTheta/2.0);
	pauliSpinor(1,0) = sin(fTheta/2.0)
	  *exp(complex<double>(0,fPhi));
      }
  } // positive energy solutions

  // Negative energy solutions (anti-particles)
  else if( eSolution<0 ) {
    // polarization = 1/2
    if( fPolarizationState==kUp )
      {
	pauliSpinor(0,0) = -sin(fTheta/2.0)
	  *exp(complex<double>(0,-fPhi));
	pauliSpinor(1,0) = cos(fTheta/2.0);
      }
    // polarization = -1/2
    else if( fPolarizationState==kDown )
      {
	pauliSpinor(0,0) = cos(fTheta/2.0);
	pauliSpinor(1,0) = sin(fTheta/2.0)
	  *exp(complex<double>(0,fPhi));
      }
  } // negative energy solutions

  return pauliSpinor;
}

//______________________________________________________________________________
TSpinor::Polarization::operator int() const
{
  // Cast polarization state to int: 
  // returns magnetic quantum number multiplied by two

  if( fPolarizationState==kUp ) return 1;
  else if( fPolarizationState==kDown ) return -1;

  return INT_MAX;
}

//______________________________________________________________________________
TSpinor::Polarization::operator State() const
{
  // Cast to polarization state
  return fPolarizationState;
}

/*//______________________________________________________________________________
TString TSpinor::Polarization::HashName() const
{
  // Returns a unique string that specifies the current state of the
  // spinor's polarization.
  
  // the polarization is fixed by the angles of the quantization axis
  // and the state.
  double key[3] = { fTheta, fPhi, (double)fPolarizationState };
  
  // We create a hash key using the function in the cachetools folder
  char hash[sizeof(double)*3*2 + 1];
  doublestochar(3, key, hash);

  return hash;
}*/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TVSpinor                                                              //
//                                                                       //
// Dirac spinor for negative energy on-mass shell spin-1/2 particles     //
//
//BEGIN_LATEX
//v(#vec{p},s) = C u^{T}(#vec{p},s) = #sqrt{E_{p}+m}#left( #splitline{#frac{#vec{#sigma}#upoint#vec{p}}{E_{p}+m}}{1} #right) #eta^{s},
//END_LATEX
// with BEGIN_LATEX E_{p} = #sqrt{#vec{p}^{2}+m^{2}} END_LATEX and BEGIN_LATEX C=-i#gamma^{0}#gamma^{2} END_LATEX.
// See documentation of TSpinor::Polarization::GetPauliSpinor for the 
// definition of the two-component spinors BEGIN_LATEX #eta^{s} END_LATEX.
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TVSpinor)

//_____________________________________________________________________
TVSpinor::TVSpinor(const TLorentzVector& v, double mass, const Polarization& pol, 
		   Normalization norm)
: TSpinor(v.X(),v.Y(),v.Z(),mass,pol,norm,-1)
{
  // Constructor
}

//_____________________________________________________________________
TVSpinor::TVSpinor(const FourVector<double>& v, double mass, const Polarization& pol,
		   Normalization norm)
  : TSpinor(v[1],v[2],v[3],mass,pol,norm,-1)
{
  // Constructor
}

//_____________________________________________________________________
TVSpinor::TVSpinor(TRootIOCtor* rio)
  : TSpinor(rio)
{
  // ROOT I/0 Constructor
}

//_____________________________________________________________________
TVSpinor::TVSpinor(const TVSpinor& rhs)
  : TSpinor(rhs)
{
  // Copy constructor
}

//_____________________________________________________________________
TVSpinor::~TVSpinor()
{
  // Destructor
}

//_____________________________________________________________________
TVSpinor& TVSpinor::operator=(const TVSpinor& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
    TSpinor::operator=(rhs);
  }

  return *this;
}

//_____________________________________________________________________
TVSpinor* TVSpinor::Clone(const char*) const
{
  // Virtual copy constructor
  return new TVSpinor(*this);
}

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// THelicitySpinor                                                       //
//                                                                       //
// Helicity spinor for positive energy on-mass shell spin-1/2 particles  //
//
//BEGIN_LATEX
//u(#vec{p},#lambda) = #sqrt{E_{p}+m}#left( #splitline{1}{2#lambda#frac{|#vec{p}|}{E_{p}+m}} #right) #chi^{#lambda},
//END_LATEX
// with BEGIN_LATEX E_{p} = #sqrt{#vec{p}^{2}+m^{2}} END_LATEX.
// See documentation of TSpinor::Polarization::GetPauliSpinor for the 
// definition of the two-component spinors BEGIN_LATEX #chi^{#lambda} END_LATEX.
// The two-component spinors are calculated using the polar and azimuthal
// angles of the particle's 3-momentum.
//                                                                       //
///////////////////////////////////////////////////////////////////////////
ClassImp(THelicitySpinor)

//_____________________________________________________________________
THelicitySpinor::THelicitySpinor(const TLorentzVector& v, double mass, 
				 const Polarization::State& polState, 
				 Normalization norm)
: TSpinor(0.,0.,v.P(),mass,
	  (v.P() < STRANGEUFLOW) ? 
	  TSpinor::Polarization(0.,0.,polState) : 
	  TSpinor::Polarization(v.Theta(),v.Phi(),polState),
	  norm,1)
{
  // Constructor
  ToHelicitySpinor(polState,v.Phi());
}

//_____________________________________________________________________
THelicitySpinor::THelicitySpinor(const FourVector<double>& v, double mass, 
				 const Polarization::State& polState, 
				 Normalization norm)
  : TSpinor(0.,0.,sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]), mass,
	    (sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]) < STRANGEUFLOW) ? 
	    TSpinor::Polarization(0.,0.,polState) :
	    TSpinor::Polarization(v.theta(),v.phi(),polState),
	    norm,1)
{
  // Constructor
  ToHelicitySpinor(polState,v.phi());
}

//_____________________________________________________________________
THelicitySpinor::THelicitySpinor(TRootIOCtor* rio)
  : TSpinor(rio)
{
  // ROOT I/0 Constructor
}

//_____________________________________________________________________
THelicitySpinor::THelicitySpinor(const THelicitySpinor& rhs)
  : TSpinor(rhs)
{
  // Copy constructor
}

//_____________________________________________________________________
THelicitySpinor::~THelicitySpinor()
{
  // Destructor
}

//_____________________________________________________________________
THelicitySpinor& THelicitySpinor::operator=(const THelicitySpinor& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
    TSpinor::operator=(rhs);
  }

  return *this;
}

//_____________________________________________________________________
THelicitySpinor* THelicitySpinor::Clone(const char*) const
{
  // Virtual copy constructor
  return new THelicitySpinor(*this);
}

//_____________________________________________________________________
void THelicitySpinor::ToHelicitySpinor(const Polarization::State& polState,
				       const double phi)
{
  // In THelicitySpinor's constructor we obviously make a call
  // to TSpinor's constructor. This means some components of the
  // spinor are incorrect. This function takes care of this.
  
  // We have called TSpinor's constructor. This means the Pauli spinors
  // have not been rotated in the correct way. We need to correct for this
  complex<double> factor;
  if( polState==TSpinor::Polarization::kUp )
    factor = exp(complex<double>(0.,-phi/2.));
  else
    factor = exp(complex<double>(0.,phi/2.));
  
  for(int i=0; i<4; ++i)
    (*fComponent)(i,0) *= factor;

  // In addition some components have the wrong sign.
  if( polState==TSpinor::Polarization::kUp ) {
    (*fComponent)(3,0) *= -1.;
  } else {
    (*fComponent)(2,0) *= -1.;
  }
}

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// THelicityVSpinor                                                      //
//                                                                       //
// Helicity spinor for negative energy on-mass shell spin-1/2 particles  //
//
//BEGIN_LATEX
//v(#vec{p},#lambda) = #sqrt{E_{p}+m}#left( #splitline{-2#lambda#frac{|#vec{p}|}{E_{p}+m}}{1} #right) #eta^{#lambda},
//END_LATEX
// with BEGIN_LATEX E_{p} = #sqrt{#vec{p}^{2}+m^{2}} END_LATEX.
// See documentation of TSpinor::Polarization::GetPauliSpinor for the 
// definition of the two-component spinors BEGIN_LATEX #eta^{#lambda} END_LATEX.
// The two-component spinors are calculated using the polar and azimuthal
// angles of the particle's 3-momentum.
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(THelicityVSpinor)

//_____________________________________________________________________
THelicityVSpinor::THelicityVSpinor(const TLorentzVector& v, double mass, 
				 const Polarization::State& polState, 
				 Normalization norm)
: TSpinor(0.,0.,v.P(),mass,
	  TSpinor::Polarization(v.Theta(),v.Phi(),polState),
	  norm,-1)
{
  // Constructor
  ToHelicityVSpinor(polState,v.Phi());
}

//_____________________________________________________________________
THelicityVSpinor::THelicityVSpinor(const FourVector<double>& v, double mass, 
				 const Polarization::State& polState, 
				 Normalization norm)
  : TSpinor(0.,0.,sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]), mass,
	    TSpinor::Polarization(v.theta(),v.phi(),polState),
	    norm,-1)
{
  // Constructor
  ToHelicityVSpinor(polState,v.phi());
}

//_____________________________________________________________________
THelicityVSpinor::THelicityVSpinor(TRootIOCtor* rio)
  : TSpinor(rio)
{
  // ROOT I/0 Constructor
}

//_____________________________________________________________________
THelicityVSpinor::THelicityVSpinor(const THelicityVSpinor& rhs)
  : TSpinor(rhs)
{
  // Copy constructor
}

//_____________________________________________________________________
THelicityVSpinor::~THelicityVSpinor()
{
  // Destructor
}

//_____________________________________________________________________
THelicityVSpinor& THelicityVSpinor::operator=(const THelicityVSpinor& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
    TSpinor::operator=(rhs);
  }

  return *this;
}

//_____________________________________________________________________
THelicityVSpinor* THelicityVSpinor::Clone(const char*) const
{
  // Virtual copy constructor
  return new THelicityVSpinor(*this);
}

//_____________________________________________________________________
void THelicityVSpinor::ToHelicityVSpinor(const Polarization::State& polState,
					 const double phi)
{
  // In THelicityVSpinor's constructor we obviously make a call
  // to TSpinor's constructor. This means some components of the
  // spinor are incorrect. This function takes care of this.
  
  // We have called TSpinor's constructor. This means the Pauli spinors
  // have not been rotated in the correct way. We need to correct for this
  complex<double> factor;
  if( polState==TSpinor::Polarization::kUp )
    factor = exp(complex<double>(0.,phi/2.));
  else
    factor = exp(complex<double>(0.,-phi/2.));
  
  for(int i=0; i<4; ++i)
    (*fComponent)(i,0) *= factor;

  // In addition some components have the wrong sign.  
  if( polState==TSpinor::Polarization::kUp ) {
    (*fComponent)(0,0) *= -1.;
  } else {
    (*fComponent)(1,0) *= -1.;
  }
}

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TSpinor::LorentzRotation                                              //
//                                                                       //
// A representation of the Lorentz group on the four-dimensional Dirac space.
// Takes care of pure boosts and rotation or any combination. This class is
// closely linked to TLorentzRotation which acts on Lorentz indices.
// 
// A concise treatment of the transformation rules for spinors can be found
// in section 5.9 (p.146-152) in
// "Relativistic Quantum Field Mechanics and Field Theory" by Franz Gross.
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TSpinor::LorentzRotation)

//______________________________________________________________________________
TSpinor::LorentzRotation::LorentzRotation()
: fTransformation(1.)
{
  // Default constructor (unity)
}

//______________________________________________________________________________
TSpinor::LorentzRotation::LorentzRotation(const TRotation& rot)
: fTransformation(1.)
{
  // Constructor
  // Initialize transformation in spinor space corresponding to the given rotation
  double theta;
  TVector3 axis;
  rot.AngleAxis(theta,axis);
  SetRotation( theta*axis.X(), theta*axis.Y(), theta*axis.Z() );
}

//______________________________________________________________________________
TSpinor::LorentzRotation::LorentzRotation(const TVector3& boostVector)
: fTransformation(1.)
{
  // Constructor
  // Initialize transformation in spinor space corresponding to the given boost
  SetBoost(boostVector.X(),boostVector.Y(),boostVector.Z());
}

//______________________________________________________________________________
TSpinor::LorentzRotation::LorentzRotation(double betaX, double betaY, double betaZ)
: fTransformation(1.)
{
  // Constructor
  // Initialize transformation in spinor space corresponding to the given boost
  SetBoost(betaX,betaY,betaZ);
}

//______________________________________________________________________________
TSpinor::LorentzRotation& TSpinor::LorentzRotation::Invert()
{
  // Invert the transformation
  fTransformation = fTransformation.Conjugate();
  return *this;
}

//______________________________________________________________________________
TSpinor::LorentzRotation TSpinor::LorentzRotation::Inverse() const
{
  // Calculate the inverse transformation
  LorentzRotation inverse = *this;
  return inverse.Invert();
}

//______________________________________________________________________________
TSpinor TSpinor::LorentzRotation::operator*(const TSpinor& rhs) const
{
  // Apply the current transformation to the spinor indices
  TSpinor rotated = rhs;
  (*rotated.fComponent) = fTransformation * (*rotated.fComponent);
  return rotated;
}

//______________________________________________________________________________
TSpinor::LorentzRotation TSpinor::LorentzRotation::operator*(const LorentzRotation& m) const
{
  // Compound two transformations
  LorentzRotation product = m;
  product.Transform(*this);
  return product;
}

//______________________________________________________________________________
TSpinor::LorentzRotation& TSpinor::LorentzRotation::operator*=(const LorentzRotation& m)
{
  // Compound two transformations
  LorentzRotation product = m;
  product.Transform(*this);
  return *this = product;
}

//______________________________________________________________________________
GammaStructure TSpinor::LorentzRotation::operator*(const GammaStructure& rhs) const
{
  return static_cast<const GammaStructure&>(*this) * rhs;
}

//______________________________________________________________________________
TSpinor::LorentzRotation& TSpinor::LorentzRotation::Rotate(double angle, const TVector3& axis)
{
  // Rotate over a certain angle around an axis
  SetRotation(angle*axis.X(),angle*axis.Y(),angle*axis.Z());
  return *this;
}

//______________________________________________________________________________
TSpinor::LorentzRotation& TSpinor::LorentzRotation::Transform(const LorentzRotation& m)
{
  // Compound two transformations, e.g. a.Transform(b) --> a=b*a
  fTransformation = m.fTransformation * fTransformation;
  return *this;
}

//______________________________________________________________________________
void TSpinor::LorentzRotation::SetBoost(double betaX, double betaY, double betaZ)
{
  // Construct the transformation matrix for a pure boost defined by BEGIN_LATEX #vec{#beta}=(#beta_{x},#beta_{y},#beta_{z})END_LATEX:
  //BEGIN_LATEX S(L_{#vec{#beta}}) = cosh#frac{#xi}{2} + #sum_{i} sinh#frac{#xi}{2} #beta_{i} #gamma_{0} #gamma_{i} END_LATEX
  // with BEGIN_LATEX #xi END_LATEX the rapidity.

  // Parameters of Pure Lorentz boost
  double beta = sqrt( betaX*betaX + betaY*betaY + betaZ*betaZ );
  double coshksi = sqrt( (1./sqrt(1.-beta*beta) + 1.)/2. );
  double sinhksi = sqrt( (1./sqrt(1.-beta*beta) - 1.)/2. );
  
  // Do nothing in case of no boost
  if( beta<STRANGEUFLOW ) return;

  // Define the transformation matrix
  fTransformation = 
    GammaStructure(coshksi, // 1
  		   0., // G5
  		   0., // G0
  		   0., // G1
  		   0., // G2
  		   0., // G3
  		   0., // G5 G0
  		   0., // G5 G1
  		   0., // G5 G2
  		   0., // G5 G3
  		   sinhksi*betaX/beta, // G0 G1
  		   sinhksi*betaY/beta, // G0 G2
  		   sinhksi*betaZ/beta, // G0 G3
  		   0.,  // G5 G0 G1
  		   0.,  // G5 G0 G2 
  		   0.)  // G5 G0 G3
    * fTransformation;
}

//______________________________________________________________________________
void TSpinor::LorentzRotation::SetRotation(double thetaX, double thetaY, double thetaZ)
{
  // Construct the transformation matrix for a pure rotation defined by BEGIN_LATEX #vec{#theta}=(#theta_{x},#theta_{y},#theta_{z})END_LATEX:
  //BEGIN_LATEX S(L_{#vec{#theta}}) = cosh#frac{#theta}{2} -i #sum_{i} sinh#frac{#theta}{2} #theta_{i} #gamma_{5} #gamma_{0} #gamma_{i} END_LATEX

  // Parameters of Pure rotation
  double theta = sqrt( thetaX*thetaX + thetaY*thetaY + thetaZ*thetaZ );
  double costheta = cos(theta/2.);
  complex<double> sintheta = complex<double>(0., -sin(theta/2.) );
 
  // Do nothing in case of no rotation
  if( theta<STRANGEUFLOW ) return;

  // Define the transformation matrix
  fTransformation = 
    GammaStructure(costheta, // 1
		   0., // G5
		   0., // G0
		   0., // G1
		   0., // G2
		   0., // G3
		   0., // G5 G0
		   0., // G5 G1
		   0., // G5 G2
		   0., // G5 G3
		   0., // G0 G1
		   0., // G0 G2
		   0., // G0 G3
		   sintheta*thetaX/theta,  // G5 G0 G1
		   sintheta*thetaY/theta,  // G5 G0 G2 
		   sintheta*thetaZ/theta)  // G5 G0 G3
    * fTransformation;
}

//______________________________________________________________________________
void TSpinor::LorentzRotation::Streamer(TBuffer &R__b){}
// {
//   // Stream an object of class TSpinor::LorentzRotation.

//   //This works around a msvc bug and should be harmless on other platforms
//   typedef ::TSpinor::LorentzRotation thisClass;
//   UInt_t R__s, R__c;
  
//   if (R__b.IsReading()) { // begin reading
//     Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
//     TObject::Streamer(R__b);

//     // Retrieve all GammaStructure coefficients from the buffer
//     complex<double> compCoeff[16];
//     double real,imag;
//     for(int i=0; i<16; ++i) {
//       R__b >> real;
//       R__b >> imag;
//       compCoeff[i]=complex<double>(real,imag);
//     }
//     fTransformation 
//       = GammaStructure(compCoeff[0],compCoeff[1],compCoeff[2],compCoeff[3],
// 		       compCoeff[4],compCoeff[5],compCoeff[6],compCoeff[7],
// 		       compCoeff[8],compCoeff[9],compCoeff[10],compCoeff[11],
// 		       compCoeff[12],compCoeff[13],compCoeff[14],compCoeff[15]);

//     R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
//   } // end reading 

//   else { // begin writing
//     R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
//     TObject::Streamer(R__b);
    
//     // We wills store the coefficients of the GammaStructure
//     complex<double> compCoeff[16];
//     compCoeff[0] = Trace(fTransformation.value()) /4.;
//     compCoeff[1] = Trace(fTransformation.value()*G5) /4.;
//     compCoeff[2] = Trace(fTransformation.value()*G0) /4.;
//     compCoeff[3] = Trace(fTransformation.value()*G1) /-4.;
//     compCoeff[4] = Trace(fTransformation.value()*G2) /-4.;
//     compCoeff[5] = Trace(fTransformation.value()*G3) /-4.;
//     compCoeff[6] = Trace(fTransformation.value()*G5G0) /-4.;
//     compCoeff[7] = Trace(fTransformation.value()*G5G1) /4.;
//     compCoeff[8] = Trace(fTransformation.value()*G5G2) /4.;
//     compCoeff[9] = Trace(fTransformation.value()*G5G3) /4.;
//     compCoeff[10] = Trace(fTransformation.value()*A1) /4.;
//     compCoeff[11] = Trace(fTransformation.value()*A2) /4.;
//     compCoeff[12] = Trace(fTransformation.value()*A3) /4.;
//     compCoeff[13] = Trace(fTransformation.value()*S1) /4.;
//     compCoeff[14] = Trace(fTransformation.value()*S2) /4.;
//     compCoeff[15] = Trace(fTransformation.value()*S3) /4.;
//     for(int i=0; i<16; ++i) {
//       R__b << compCoeff[i].real();
//       R__b << compCoeff[i].imag();
//     }

//     R__b.SetByteCount(R__c, kTRUE);
//   } // end writing
// }

//______________________________________________________________________________
ostream& operator<<(ostream& output, const TSpinor::LorentzRotation& rhs) 
{
  return output << static_cast<const GammaStructure&>(rhs);
}

//______________________________________________________________________________
Matrix<4,1> operator*(const Matrix<4,4>& lhs, const TSpinor& rhs)
{ 
  return lhs * (Matrix<4,1>)rhs; 
}

//______________________________________________________________________________
complex<double> operator*(const Matrix<1,4>& lhs, const TSpinor& rhs)
{ 
  return lhs * (Matrix<4,1>)rhs; 
}

//______________________________________________________________________________
Matrix<4,1> operator*(const GammaStructure& lhs, const TSpinor& rhs)
{
  return (Matrix<4,4>)lhs * (Matrix<4,1>)rhs;
}

//______________________________________________________________________________
TLorentzRotation THelicitySpinor::HelicityTransformation(const TLorentzVector& p)
{
  // Compute the helicity transformation h(p) as defined in section 1.2.2 of
  // 'Spin in particle physics' by E.Leader.
  //
  //BEGIN_LATEX h(p) = r_{z}(#phi) #times r_{y}(#theta) #times l_{z}(#frac{|#vec{p}|}{p^{0}})END_LATEX
  return 
    TLorentzRotation(0.,0.,p.P()/p.E())
    .RotateY(p.Theta())
    .RotateZ(p.Phi());
}

//______________________________________________________________________________
TSpinor::LorentzRotation THelicitySpinor::HelicitySpinorTransformation(const TLorentzVector& p)
{
  // Compute the helicity transformation h(p) as defined in section 1.2.2 of
  // 'Spin in particle physics' by E.Leader.
  //
  //BEGIN_LATEX h(p) = r_{z}(#phi) #times r_{y}(#theta) #times l_{z}(#frac{|#vec{p}|}{p^{0}})END_LATEX
  return 
    TSpinor::LorentzRotation(0.,0.,p.P()/p.E())
    .RotateY(p.Theta()) 
    .RotateZ(p.Phi());
}

//______________________________________________________________________________
TLorentzQuaternion THelicitySpinor::HelicityQuaternionTransformation(const TLorentzVector& p)
{
  // Compute the helicity transformation h(p) as defined in section 1.2.2 of
  // 'Spin in particle physics' by E.Leader.
  //
  //BEGIN_LATEX h(p) = r_{z}(#phi) #times r_{y}(#theta) #times l_{z}(#frac{|#vec{p}|}{p^{0}})END_LATEX
  return 
    TLorentzQuaternion(0.,0.,p.P()/p.E())
    .RotateY(p.Theta())
    .RotateZ(p.Phi());
}

//______________________________________________________________________________
TLorentzRotation THelicitySpinor::WickRotation(const TLorentzRotation& l,
						 const TLorentzVector& p)
{
  // Compute the Wick helicity rotation as defined in (2.1.7) in section 2.1
  // of 'Spin in particle physics' by E.Leader.
  // The Lorentz transformation 'l', given as argument, acts on 4-vectors,...
  // however not on the reference frame.
  //
  //BEGIN_LATEX r(l,p) = h(l*p)^{-1} #times l #times h(l*p)END_LATEX
  return 
    HelicityTransformation(l*p).Inverse() 
    * l
    * HelicityTransformation(p);
}

//______________________________________________________________________________
TLorentzQuaternion THelicitySpinor::WickRotation(const TLorentzQuaternion& l,
						 const TLorentzVector& p)
{
  // Compute the Wick helicity rotation as defined in (2.1.7) in section 2.1
  // of 'Spin in particle physics' by E.Leader.
  // The Lorentz transformation 'l', given as argument, acts on 4-vectors,...
  // however not on the reference frame.
  //
  //BEGIN_LATEX r(l,p) = h(l*p)^{-1} #times l #times h(l*p)END_LATEX
  return 
    HelicityQuaternionTransformation(l*p).Inverse() 
    * l
    * HelicityQuaternionTransformation(p);
}

//______________________________________________________________________________
void THelicitySpinor::WickRotationEulerAngles(const TLorentzRotation& l,
					      const TLorentzVector& p,
					      double& alpha, double& beta, 
					      double& gamma)
{
  // Compute the Wick helicity rotation as defined in (2.1.7) in section 2.1
  // of 'Spin in particle physics' by E.Leader.
  // The Lorentz transformation 'l', given as argument, acts on 4-vectors,...
  // however not on the reference frame.
  //
  //BEGIN_LATEX r(l,p) = h(l*p)^{-1} #times l #times h(l*p) = r(#alpha,#beta,#gamma)END_LATEX
  // We use the Y-convention Euler angles as explained in the documentation 
  // of TRotation.

  TLorentzRotation wickRotation = WickRotation(l,p);

  // We know the wickRotation is a pure rotation, but ROOT doesn't give
  // us the possibility to turn a TLorentzRotation in a TRotation. We will
  // copy TRotation::GetYTheta(), TRotation::GetYPhi() and TRotation::GetYPsi()
  // verbatim here, but calling wickRotation's members instead.

  // We define an underflow parameter to deal with sines and cosines near
  // the edges.
  const double underflow = 1.e-7;

  Double_t s2 =  1.0 - wickRotation.ZZ()*wickRotation.ZZ();
  if (s2 < 0) {
    if( (fabs(wickRotation.ZZ())-1.) < underflow ) {
      s2 = 0;
    } 
    else {
      cerr << "WARNING in  THelicitySpinor::WickRotationEulerAngles("
	   << "const TLorentzRotation&,const TLorentzVector&,"
	   << "double&,double&,double&): "
	   << "|fzz| > 1" << endl;
      assert(1==0);
    }
  }
  const Double_t sinTheta = TMath::Sqrt(s2);

  if ( fabs(sinTheta) > underflow  ) {
    Double_t cscTheta = 1/sinTheta;
    Double_t cosAbsPhi =  wickRotation.ZY() * cscTheta;
    if( TMath::Abs(cosAbsPhi-1.) < underflow ) {
      cosAbsPhi = 1.;
    } else if( TMath::Abs(cosAbsPhi+1.) < underflow ) {
      cosAbsPhi = -1.;
    } else if ( TMath::Abs(cosAbsPhi) > 1 ) {        // NaN-proofing
      cerr << "Error in  THelicitySpinor::WickRotationEulerAngles("
	   << "const TLorentzRotation&,const TLorentzVector&,"
	   << "double&,double&,double&): "
	   << "finds |cos(phi)| > 1" << endl;
      assert(1==0);
    }
    Double_t absPhi;
    if( fabs(cosAbsPhi-1.) < underflow ) absPhi = 0.;
    else if( fabs(cosAbsPhi+1.) < underflow ) absPhi = TMath::Pi();
    else absPhi= TMath::ACos(cosAbsPhi);
    if (wickRotation.ZX() > 0) {
      alpha = absPhi;
    } else if (wickRotation.ZX() < 0) {
      alpha = -absPhi;
    } else if (wickRotation.ZY() > 0) {
      alpha = 0.0;
    } else {
      alpha = TMath::Pi();
    }
  } else {              // sinTheta == 0 so |Fzz| = 1
    Double_t cosAbsPhi = wickRotation.XX();
    Double_t absPhi;
    if( fabs(cosAbsPhi-1.) < underflow ) absPhi = 0.;
    else if( fabs(cosAbsPhi+1.) < underflow ) absPhi = TMath::PiOver2();
    else absPhi= .5 * TMath::ACos(cosAbsPhi);
    if (wickRotation.XY() > 0) {
      alpha =  -absPhi;
    } else if (wickRotation.XY() < 0) {
      alpha =   absPhi;
    } else if (wickRotation.XX()>0) {
      alpha = 0.0;
    } else {
      alpha = wickRotation.ZZ() * TMath::PiOver2();
    }
  }

  // GetYPhi()
  alpha += TMath::Pi()/2.0;

  // GetYTheta()
  double cosBeta = wickRotation.ZZ();
  if( fabs(cosBeta-1.) < underflow ) beta = 0.;
  else if( fabs(cosBeta+1.) < underflow ) beta = TMath::Pi();
  else beta = TMath::ACos(cosBeta);

  // GetXPsi(): store the finalPsi in gamma
  if ( fabs(sinTheta) > underflow ) {
    Double_t cscTheta = 1/sinTheta;
    Double_t cosAbsPsi =  - wickRotation.YZ() * cscTheta;
    if ( TMath::Abs(cosAbsPsi) > 1 ) {        // NaN-proofing
      cerr << "WARNING in  THelicitySpinor::WickRotationEulerAngles("
	   << "const TLorentzRotation&,const TLorentzVector&,"
	   << "double&,double&,double&): "
	   << "GetPsi(): | cos psi | > 1" << endl;
      cosAbsPsi = 1;
    }
    Double_t absPsi;
    if( fabs(cosAbsPsi-1.) < underflow ) absPsi = 0.;
    else if( fabs(cosAbsPsi+1.) < underflow ) absPsi = TMath::Pi();
    else absPsi= TMath::ACos(cosAbsPsi);
    if (wickRotation.XZ() > 0) {
      gamma = absPsi;
    } else if (wickRotation.XZ() < 0) {
      gamma = -absPsi;
    } else {
      gamma = (wickRotation.YZ() < 0) ? 0 : TMath::Pi();
    }
  } else {              // sinTheta == 0 so |Fzz| = 1
    Double_t cosAbsPsi = wickRotation.XX();
    if ( TMath::Abs(cosAbsPsi) > 1 ) {        // NaN-proofing
      cerr << "WARNING in  THelicitySpinor::WickRotationEulerAngles("
	   << "const TLorentzRotation&,const TLorentzVector&,"
	   << "double&,double&,double&): "
	   << "GetPsi(): | fxx | > 1" << endl;
      cosAbsPsi = 1;
    }
    Double_t absPsi;
    if( fabs(cosAbsPsi-1.) < underflow ) absPsi = 0.;
    else if( fabs(cosAbsPsi+1.) < underflow ) absPsi = TMath::PiOver2();
    else absPsi= .5 * TMath::ACos(cosAbsPsi);
    if (wickRotation.YX() > 0) {
      gamma = absPsi;
    } else if (wickRotation.YX() < 0) {
      gamma = -absPsi;
    } else {
      gamma = (wickRotation.XX() > 0) ? 0 : TMath::PiOver2();
    }
  }

  // GetYPsi()
  gamma -= TMath::Pi()/2.0;
}

//______________________________________________________________________________
Matrix<2,2> THelicitySpinor::RotationMatrix(const TLorentzQuaternion& l,
					    const TLorentzVector& p)
{
  // Compute the Wick helicity rotation as defined in (2.1.7) in section 2.1
  // of 'Spin in particle physics' by E.Leader and determine the corresponding
  // Wigner rotation matrix for spin-1/2 particles.
  //BEGIN_LATEX
  //U(l)|p,#lambda> = #sum_{#lambda'} M_{#lambda,#lambda'} |l*p,#lambda'>
  //END_LATEX
  // with
  //BEGIN_LATEX M_{#lambda,#lambda'} = D_{#lambda,#lambda'}#left(r(l,p)#right)END_LATEX
  //
  // The Lorentz transformation 'l', given as argument, acts on 4-vectors,...
  // however not on the reference frame.
  double omega;				    
  TVector3 axis,boost;
  WickRotation(l,p).Decompose(omega,axis,boost);
  const double theta = axis.Theta();
  const double phi = axis.Phi();
  const complex<double> expIphi = exp(complex<double>(0.,phi+TMath::PiOver2()));
  const complex<double> expMIphi = exp(complex<double>(0.,-phi+TMath::PiOver2()));

  // This formula is taken from Eq.(18) in section (6.2.4) of "Quantum Theory
  // of Angular Momentum" by D.Varshalovich, A.Moskalev, V.Khersonskii.
  return Matrix<2,2>(complex<double>(cos(omega/2.),-sin(omega/2.)*cos(theta)),
		     (-sin(omega/2.)*sin(theta))*expIphi,
		     (-sin(omega/2.)*sin(theta))*expMIphi,
		     complex<double>(cos(omega/2.),sin(omega/2.)*cos(theta)));
}

//______________________________________________________________________________
Matrix<2,2> THelicitySpinor::RotationMatrix(const TLorentzRotation& l,
					    const TLorentzVector& p)
{
  // Compute the Wick helicity rotation as defined in (2.1.7) in section 2.1
  // of 'Spin in particle physics' by E.Leader and determine the corresponding
  // Wigner rotation matrix for spin-1/2 particles.
  //BEGIN_LATEX
  //U(#Lambda)|p,#lambda> = #sum_{#lambda'} M_{#lambda,#lambda'} |#Lambda p,#lambda'>
  //END_LATEX
  // with
  //BEGIN_LATEX M_{#lambda,#lambda'} = D_{#lambda,#lambda'}#left(r(l,p)#right)END_LATEX
  //
  // The Lorentz transformation 'l', given as argument, acts on 4-vectors,...
  // however not on the reference frame.
  //
  // Due to the 2-to-1 homomorphism between SU(2) and O(3) this method is not
  // fail-safe. It is recommended to use the
  // RotationMatrix(const TLorentzQuaternion&,const TLorentzVector&) method.
  double alpha,beta,gamma;
  WickRotationEulerAngles(l,p,alpha,beta,gamma);

  // This formula is taken from Appendix A of 'Spin in particle physics' by E.Leader
  return Matrix<2,2>(exp(complex<double>(0.,-.5*(alpha+gamma)))*cos(beta/2.),
		     exp(complex<double>(0.,.5*(alpha-gamma)))*sin(beta/2.),
		     -exp(complex<double>(0.,-.5*(alpha-gamma)))*sin(beta/2.),
		     exp(complex<double>(0.,.5*(alpha+gamma)))*cos(beta/2.));
}
