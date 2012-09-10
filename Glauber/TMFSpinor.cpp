#include "TMFSpinor.hpp"
#include <constants.hpp>
//#include <numtoa.h>
// #include <TString.h>
#include <FourVector.h>
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <cstdlib>
// #include <climits>
#include <iostream>
#include <complex>
#include <TDeuteron.h>

using TMath::Sign;
using namespace std;

ClassImp(TMFSpinor)

TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double r, double costheta, double phi)
: TObject(), fComponent(new Matrix<4,1>){

  double Gr=nucleus.getWave_G(shellindex,r);
  double Fr=nucleus.getWave_F(shellindex,r);
  complex<double> spinup=exp(I_UNIT*0.5*phi*double(m-1));
  complex<double> spindown=exp(I_UNIT*0.5*phi*double(m+1));
  (*fComponent)(0,0)=I_UNIT*Gr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(1,0)=I_UNIT*Gr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m+1)/2, costheta);
  (*fComponent)(2,0)=-Fr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(3,0)=-Fr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m+1)/2, costheta);
  
}


TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, double costheta, double phi)
: TObject(), fComponent(new Matrix<4,1>){
  complex<double> spinup=exp(I_UNIT*0.5*phi*double(m-1));
  complex<double> spindown=exp(I_UNIT*0.5*phi*double(m+1));
  (*fComponent)(0,0)=I_UNIT*Gr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(1,0)=I_UNIT*Gr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m+1)/2, costheta);
  (*fComponent)(2,0)=-Fr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(3,0)=-Fr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m+1)/2, costheta);
  
}


TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, const double *Spher_harm, double phi)
: TObject(), fComponent(new Matrix<4,1>){
  complex<double> spinup=exp(I_UNIT*0.5*phi*double(m-1));
  complex<double> spindown=exp(I_UNIT*0.5*phi*double(m+1));
  (*fComponent)(0,0)=I_UNIT*Gr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*Spher_harm[0];
  (*fComponent)(1,0)=I_UNIT*Gr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*Spher_harm[1];
  (*fComponent)(2,0)=-Fr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*Spher_harm[2];
  (*fComponent)(3,0)=-Fr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*Spher_harm[3];
  
}

TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double r, double costheta, double phi, Matrix<4,1> CG)
: TObject(), fComponent(new Matrix<4,1>){

  double Gr=nucleus.getWave_G(shellindex,r);
  double Fr=nucleus.getWave_F(shellindex,r);
  complex<double> spinup=exp(I_UNIT*0.5*phi*double(m-1));
  complex<double> spindown=exp(I_UNIT*0.5*phi*double(m+1));
  (*fComponent)(0,0)=I_UNIT*Gr*spinup*CG(0,0)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(1,0)=I_UNIT*Gr*spindown*CG(1,0)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m+1)/2, costheta);
  (*fComponent)(2,0)=-Fr*spinup*CG(2,0)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(3,0)=-Fr*spindown*CG(3,0)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m+1)/2, costheta);
  
}


TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, double costheta, double phi, Matrix<4,1> CG)
: TObject(), fComponent(new Matrix<4,1>){
  complex<double> spinup=exp(I_UNIT*0.5*phi*double(m-1));
  complex<double> spindown=exp(I_UNIT*0.5*phi*double(m+1));
  (*fComponent)(0,0)=I_UNIT*Gr*spinup*CG(0,0)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(1,0)=I_UNIT*Gr*spindown*CG(1,0)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m+1)/2, costheta);
  (*fComponent)(2,0)=-Fr*spinup*CG(2,0)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(3,0)=-Fr*spindown*CG(3,0)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m+1)/2, costheta);
}


TMFSpinor::TMFSpinor(double Gr, double Fr, int m, const double *Spher_harm, double phi, Matrix<4,1> CG)
: TObject(), fComponent(new Matrix<4,1>){
  complex<double> spinup=exp(I_UNIT*0.5*phi*double(m-1));
  complex<double> spindown=exp(I_UNIT*0.5*phi*double(m+1));
  (*fComponent)(0,0)=I_UNIT*Gr*spinup*Spher_harm[0]*CG(0,0);
  (*fComponent)(1,0)=I_UNIT*Gr*spindown*Spher_harm[1]*CG(1,0);
  (*fComponent)(0,0)=-Fr*spinup*Spher_harm[2]*CG(2,0);
  (*fComponent)(1,0)=-Fr*spindown*Spher_harm[3]*CG(3,0);
}


//_____________________________________________________________________
TMFSpinor::TMFSpinor(TRootIOCtor* rio)
  : TObject(), fComponent(0)
{
  // ROOT I/0 Constructor
}


TMFSpinor::TMFSpinor(const TMFSpinor& rhs)
  : fComponent(new Matrix<4,1>(*rhs.fComponent))
{
  // Copy constructor

  assert(fComponent); // check allocation
}

//_____________________________________________________________________
TMFSpinor::~TMFSpinor()
{
  // Destructor
  delete fComponent;
}

//_____________________________________________________________________
TMFSpinor& TMFSpinor::operator=(const TMFSpinor& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
    TObject::operator=(rhs);
    *fComponent=*rhs.fComponent;
  }

  return *this;
}

//_____________________________________________________________________
TMFSpinor* TMFSpinor::Clone(const char*) const
{
  // Virtual copy constructor
  return new TMFSpinor(*this);
}

//______________________________________________________________________________
void TMFSpinor::Streamer(TBuffer &R__b)
{
  // Stream an object of class TSpinor.

  UInt_t R__s, R__c;
  
  if (R__b.IsReading()) { // begin reading
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
    TObject::Streamer(R__b);
    fComponent = new Matrix<4,1>;
    double realComponent, imagComponent;
    for(int i=0; i<4; ++i) {
      R__b >> realComponent;
      R__b >> imagComponent;
      (*fComponent)(i,0) = complex<double>(realComponent,imagComponent);
    }
    if( R__v<3 ) {
      cerr << "ERROR in TSpinor::Streamer(TBuffer&): "
	   << "Old version of TSpinor might contain errors.\n";
      exit(1);
    }

    R__b.CheckByteCount(R__s, R__c, TMFSpinor::IsA());
  } // end reading 

  else { // begin writing
    R__c = R__b.WriteVersion(TMFSpinor::IsA(), kTRUE);
    TObject::Streamer(R__b);
    for(int i=0; i<4; ++i) {
      R__b << (*fComponent)(i,0).real();
      R__b << (*fComponent)(i,0).imag();
    }
    R__b.SetByteCount(R__c, kTRUE);
  } // end writing
}

Matrix<1,4> TMFSpinor::H()
{
  return (((Matrix<4,1>)*this).H());
}

//______________________________________________________________________________
ostream& operator<<(ostream& output, const TMFSpinor& rhs) 
{
  return output << static_cast<const Matrix<4,1>&>(rhs);
}




//______________________________________________________________________________
Matrix<4,1> operator*(const Matrix<4,4>& lhs, const TMFSpinor& rhs)
{ 
  return lhs * (Matrix<4,1>)rhs; 
}

//______________________________________________________________________________
complex<double> operator*(const Matrix<1,4>& lhs, const TMFSpinor& rhs)
{ 
  return lhs * (Matrix<4,1>)rhs; 
}

//______________________________________________________________________________
Matrix<4,1> operator*(const GammaStructure& lhs, const TMFSpinor& rhs)
{
  return (Matrix<4,4>)lhs * (Matrix<4,1>)rhs;
}

