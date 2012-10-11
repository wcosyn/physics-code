/*
 * TSpinor.h
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: August,10 2009
 *
 */

#ifndef TSPINOR_H
#define TSPINOR_H

// Include statements for classes that could be merely forward declared
// but need a full declaration to allow the Dict.cpp file to be compiled. 
#ifdef ROOT_Rtypes
#include "TLorentzQuaternion.h"
#include <TLorentzRotation.h>
#else
class TLorentzQuaternion;
class TLorentzRotation;
#endif

#include "Matrix.h"
#include "FourVector.h"
#include <TObject.h>
#include <complex>
#include "GammaStructure.h"
#include <TVector3.h>
#include <iostream>

// Forward declarations
class TLorentzVector;
class TRotation;
class TString;

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TSpinor                                                               //
//                                                                       //
// Dirac spinor for positive energy on-mass shell spin-1/2 particles     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

class TSpinor : public TObject
{
 public:
  // forward declarations
  class Polarization;

  enum Normalization { kUnity=1,     // Bar(U)U = 1
		       kDoubleMass=2 // Bar(U)U = 2M
  };

 protected:
  TSpinor(double px,double py,double pz,double mass, const Polarization&, Normalization,int energyState);
 public:
  TSpinor(const TLorentzVector&, double mass, const Polarization&, Normalization=kDoubleMass);
  TSpinor(const FourVector<double>&, double mass, const Polarization&, Normalization=kDoubleMass);
  TSpinor(TRootIOCtor*); // ROOT I/0 Constructor
  TSpinor(const TSpinor&);
  virtual ~TSpinor();
  TSpinor& operator=(const TSpinor&);
  virtual TSpinor* Clone(const char* ="") const;

  operator const Matrix<4,1>&() const { return *fComponent; }
  template<int K>
    Matrix<4,K> operator*(const Matrix<1,K>& rhs) { return (Matrix<4,1>)(*this) * rhs; }
  
  friend 
    std::ostream& operator<<(std::ostream&, const TSpinor&);
  
  static Matrix<1,4> Bar(const TSpinor&);
  static const FourVector< Matrix<2,2> > kSigmaPauli; // Pauli matrix 4vector
  
 protected:
  void InitializeSpinor(double px,double py,double pz,double mass,const Polarization&, Normalization,int energyState);
  
 protected:
  Matrix<4,1> *fComponent; //! a spinor is a 4-component vector
 private:

  ClassDef(TSpinor,3); // Dirac spinor for on-mass shell spin-1/2 particles
  
  ///////////////////////////////////////////////////////////////////////////
  //                                                                       //
  // TSpinor::Polarization                                                 //
  //                                                                       //
  // Polarization state of a spin-1/2 Dirac particle                       //
  //                                                                       //
  ///////////////////////////////////////////////////////////////////////////

 public:
  class Polarization : public TObject
  {
  public:
    enum State { kUp=1,
		 kDown=-1 };

  public:   
    Polarization(double,double,State);
    Polarization(TRootIOCtor*); // ROOT I/O Constructor
    Polarization(const Polarization&);
    virtual ~Polarization();
    Polarization& operator=(const Polarization&);

    bool operator==(const Polarization&);
    bool operator!=(const Polarization& pol) { return !((*this)==pol); }
    bool operator==(State polarizationState);
    bool operator!=(State polarizationState) { return !((*this)==polarizationState); }
    operator int() const;
    operator State() const;
    //TString HashName() const;

    void SetTheta(double theta);
    void SetPhi(double phi);
    void SetState(State polState) { fPolarizationState=polState; }

    Polarization& operator++();
    Polarization& operator--();

    Polarization Inverse() const;
    Polarization& Invert();

    Matrix<2,1> GetPauliSpinor(int eSolution=1) const;

  private:
    double            fTheta; // polar angle in rad
    double            fPhi;   // azimuthal angle in rad
    State fPolarizationState; // polarization state

    ClassDef(Polarization,1); // Polarization state of a spin-1/2 Dirac particle

  }; // class TSpinor::Polarization

  ///////////////////////////////////////////////////////////////////////////
  //                                                                       //
  // TSpinor::LorentzRotation                                              //
  //                                                                       //
  // Lorentz transformations (boosts and rotations) in spinor space        //
  //                                                                       //
  ///////////////////////////////////////////////////////////////////////////

 public:
  class LorentzRotation : public TObject
  {
  public:
    LorentzRotation();
    LorentzRotation(const TRotation&); // Pure rotation
    LorentzRotation(const TVector3&); // Pure Lorentz boost
    LorentzRotation(double betaX, double betaY, double betaZ); // Pure Lorentz boost

    // Const Cast to GammaStructure
    operator const GammaStructure&() const { return fTransformation; }

    // Cast to GammaStructure
    operator GammaStructure&() { return fTransformation; }

    friend 
      std::ostream& operator<<(std::ostream&, const LorentzRotation&);

    LorentzRotation& Boost(const TVector3& boostVector) { return Transform(LorentzRotation(boostVector)); }
    LorentzRotation& Boost(double betaX, double betaY, double betaZ) { return Transform(LorentzRotation(betaX,betaY,betaZ)); }
    LorentzRotation  Inverse() const;
    LorentzRotation& Invert();
    TSpinor          operator*(const TSpinor&) const;
    LorentzRotation  operator*(const LorentzRotation&) const;
    LorentzRotation& operator*=(const LorentzRotation&);
    GammaStructure   operator*(const GammaStructure&) const;
    LorentzRotation& Rotate(double angle, const TVector3& axis);
    LorentzRotation& Rotate(double angle, const TVector3* axis) { return Rotate(angle,*axis); }
    LorentzRotation& RotateX(double angle) { return Rotate(angle,TVector3(1.,0.,0.)); }
    LorentzRotation& RotateY(double angle) { return Rotate(angle,TVector3(0.,1.,0.)); }
    LorentzRotation& RotateZ(double angle) { return Rotate(angle,TVector3(0.,0.,1.)); }
    LorentzRotation& Transform(const LorentzRotation&);

  protected:
    void SetBoost(double betaX, double betaY, double betaZ);
    void SetRotation(double thetaX, double thetaY, double thetaZ);
		     
  private:
    GammaStructure fTransformation; // transformation matrix

    ClassDef(LorentzRotation,1); // Lorentz transformations (boosts and rotations) in spinor space

  }; // class TSpinor::LorentzRotation  

}; // class TSpinor

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TVSpinor                                                              //
//                                                                       //
// Dirac spinor for negative energy on-mass shell spin-1/2 particles     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

class TVSpinor : public TSpinor
{
 public:
  TVSpinor(const TLorentzVector&, double mass, const Polarization&, Normalization=kDoubleMass);
  TVSpinor(const FourVector<double>&, double mass, const Polarization&, Normalization=kDoubleMass);
  TVSpinor(TRootIOCtor*); // ROOT I/0 Constructor
  TVSpinor(const TVSpinor&);
  virtual ~TVSpinor();
  TVSpinor& operator=(const TVSpinor&);
  virtual TVSpinor* Clone(const char* ="") const;

  ClassDef(TVSpinor,1); // Dirac spinor for on-mass shell spin-1/2 particles (neg.energy)

}; // class TVSpinor

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// THelicitySpinor                                                       //
//                                                                       //
// Helicity spinor for positive energy on-mass shell spin-1/2 particles  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

class THelicitySpinor : public TSpinor
{
 public:
  THelicitySpinor(const TLorentzVector&, double mass, const Polarization::State&, Normalization=kDoubleMass);
  THelicitySpinor(const FourVector<double>&, double mass, const Polarization::State&, Normalization=kDoubleMass);
  THelicitySpinor(TRootIOCtor*); // ROOT I/0 Constructor
  THelicitySpinor(const THelicitySpinor&);
  virtual ~THelicitySpinor();
  THelicitySpinor& operator=(const THelicitySpinor&);
  virtual THelicitySpinor* Clone(const char* ="") const;

  static TLorentzRotation HelicityTransformation(const TLorentzVector&);
  static TSpinor::LorentzRotation HelicitySpinorTransformation(const TLorentzVector&);
  static TLorentzQuaternion HelicityQuaternionTransformation(const TLorentzVector&);
  static TLorentzRotation WickRotation(const TLorentzRotation&, const TLorentzVector&);
  static TLorentzQuaternion WickRotation(const TLorentzQuaternion&, const TLorentzVector&);
  static void WickRotationEulerAngles(const TLorentzRotation&, const TLorentzVector&,double& alpha, double& beta, double& gamma);
  static Matrix<2,2> RotationMatrix(const TLorentzRotation&, const TLorentzVector&);
  static Matrix<2,2> RotationMatrix(const TLorentzQuaternion&, const TLorentzVector&);

 private:
  void ToHelicitySpinor(const Polarization::State&, const double phi);

  ClassDef(THelicitySpinor,1); // Helicity spinor for on-mass shell spin-1/2 particles (pos.energy)

}; // class THelicitySpinor

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// THelicityVSpinor                                                      //
//                                                                       //
// Helicity spinor for negative energy on-mass shell spin-1/2 particles  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

class THelicityVSpinor : public TSpinor
{
 public:
  THelicityVSpinor(const TLorentzVector&, double mass, const Polarization::State&, Normalization=kDoubleMass);
  THelicityVSpinor(const FourVector<double>&, double mass, const Polarization::State&, Normalization=kDoubleMass);
  THelicityVSpinor(TRootIOCtor*); // ROOT I/0 Constructor
  THelicityVSpinor(const THelicityVSpinor&);
  virtual ~THelicityVSpinor();
  THelicityVSpinor& operator=(const THelicityVSpinor&);
  virtual THelicityVSpinor* Clone(const char* ="") const;

 private:
  void ToHelicityVSpinor(const Polarization::State&, const double phi);

  ClassDef(THelicityVSpinor,1); // Helicity spinor for on-mass shell spin-1/2 particles (neg.energy)

}; // class THelicityVSpinor

//______________________________________________________________________________
Matrix<4,1> operator*(const Matrix<4,4>& lhs, const TSpinor& rhs);
std::complex<double> operator*(const Matrix<1,4>& lhs, const TSpinor& rhs);
Matrix<4,1> operator*(const GammaStructure& lhs, const TSpinor& rhs);

#endif // TSPINOR_H
