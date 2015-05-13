/*
 * TLorentzQuaternion.h
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * Work started on November 13, 2010
 * 
 */

#ifndef LORENTZQUATERNION_H
#define LORENTZQUATERNION_H

// Include statements for classes that could be merely forward declared
// but need a full declaration to allow the Dict.cpp file to be compiled. 
#ifdef ROOT_Rtypes
#include <TLorentzVector.h>
#include <TVector3.h>
#else
class TLorentzVector;
class TVector3;
#endif

#include <TObject.h>
#include <complex>

// Forward declarations
class TRotation;

//////////////////////////////////////////////////////////////////////////
//
// TLorentzQuaternion
//
//////////////////////////////////////////////////////////////////////////

class TLorentzQuaternion : public TObject
{
 public:
  TLorentzQuaternion();
  TLorentzQuaternion(TRootIOCtor*); // ROOT I/O Constructor
  TLorentzQuaternion(const double, const TVector3&); // Pure rotation
  TLorentzQuaternion(const TRotation&); // Pure rotation
  TLorentzQuaternion(const TVector3&); // Pure Lorentz boost
  TLorentzQuaternion(const double betaX, const double betaY, const double betaZ); // Pure Lorentz boost
  
  TLorentzQuaternion& Boost(const TVector3&);
  TLorentzQuaternion& Boost(const double,const double,const double);
  TLorentzQuaternion& Rotate(const double, const TVector3&);
  TLorentzQuaternion& Rotate(const double, const TVector3*);
  TLorentzQuaternion& RotateX(const double);
  TLorentzQuaternion& RotateY(const double);
  TLorentzQuaternion& RotateZ(const double);

  TLorentzQuaternion& Invert();
  TLorentzQuaternion Inverse() const;
  TLorentzQuaternion& Conjugate();
  TLorentzQuaternion& ComplexConjugate();
  TLorentzQuaternion& Hermitian();
  static TLorentzQuaternion Conjugate(const TLorentzQuaternion&);
  static TLorentzQuaternion ComplexConjugate(const TLorentzQuaternion&);
  static TLorentzQuaternion Hermitian(const TLorentzQuaternion&);

  TLorentzQuaternion operator*(const TLorentzQuaternion&) const;
  TLorentzVector operator*(const TLorentzVector&) const;
  TVector3 operator*(const TVector3&) const;
  TLorentzQuaternion& operator*=(const TLorentzQuaternion&);
  TLorentzQuaternion& Transform(const TLorentzQuaternion&);

  bool IsIdentity() const;
  void Decompose(TLorentzQuaternion& rotation, TLorentzQuaternion& boost) const;
void Decompose(double& theta, TVector3& axis, TVector3& boost) const;

  friend std::ostream& operator<<(std::ostream&, const TLorentzQuaternion&);

 protected:
  void SetBoost(const double betaX, const double betaY, const double betaZ);
  void SetRotation(const double, const TVector3&);

 private:
  TLorentzQuaternion(const std::complex<double>& q0,const std::complex<double>& q1,const std::complex<double>& q2,const std::complex<double>& q3);
 
 private:
  std::complex<double> fQ0; // Scalar part biquaternion
  std::complex<double> fQ1; // 1st component vector part biquaternion
  std::complex<double> fQ2; // 2nd component vector part biquaternion
  std::complex<double> fQ3; // 3rd component vector part biquaternion
 
  ClassDef(TLorentzQuaternion,1); // Lorentz transformations (boosts and rotations) represented by a biquaternion

}; // class TLorentzQuaternion

#endif // LORENTZQUATERNION_H
