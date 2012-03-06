/*
 * TLorentzQuaternion.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * Work started on November 13, 2010
 * 
 */

#include "TLorentzQuaternion.h"
 #include "constants.hpp"
#include <complex>
#include <TRotation.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <cmath>
#include <cassert>
#include <iostream>
using std::complex;
using std::sqrt;
using std::cosh;
using std::sinh;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::cosh;
using std::sinh;
using std::atan2;

//////////////////////////////////////////////////////////////////////////
//
// TLorentzQuaternion
//
// Representation of Lorentz transformationS using biquaternions.
//
// |------------|
// | Definition |
// |------------|
// A biquaternion Q is a set of four complex numbers:
//BEGIN_LATEX Q = [Q_{0},Q_{1},Q_{2},Q_{3}] = [Q_{0}, #vec{Q} ]END_LATEX
// for which we define addition and multiplication:
//BEGIN_LATEX
//Q + Q' = [ Q_{0}+Q_{0}, #vec{Q}+#vec{Q}' ]
//QQ' = [ Q_{0}Q'_{0} - #vec{Q}#upoint#vec{Q}', Q_{0}#vec{Q}' + Q'_{0}#vec{Q} + #vec{Q}#times#vec{Q}']
//END_LATEX
// In addition we define conjugation (Conjugate()), complex conjugation (ComplexConjugate())
// and hermitian conjugation (Hermitian()):
//BEGIN_LATEX
//#bar{Q} = [Q_{0},-Q_{1},-Q_{2},-Q_{3}]
//Q* = [Q*_{0},Q*_{1},Q*_{2},Q*_{3}]
//Q^{+} = #bar{Q}*
//END_LATEX
//
// |------------------------|
// | Lorentz transformation |
// |------------------------|
// A general Lorentz transformation can be represented by a biquaterion of
// unit norm, i.e.
//BEGIN_LATEX Q#bar{Q}=1 END_LATEX
// Consecutive transformations can be chained by multiplying from the left.
// Representing a general 4-vector v as a quaternion V
//BEGIN_LATEX v#equiv(t,#vec{x}) #Rightarrow V#equiv[t,i#vec{x}]END_LATEX
// the 4-vector can be Lorentz transformed applying
//BEGIN_LATEX V' =  Q V Q^{+}END_LATEX
// with Q the transformation.
//
// A pure rotation over angle theta around a unit axis n is represented by
//BEGIN_LATEX
//Q_{rotation}(#theta,#vec{n}) = [ cos#frac{#theta}{2}, #vec{n}sin#frac{#theta}{2}]
//END_LATEX
// A pure boost with boost vector beta and rapidity ksi is represented by
//BEGIN_LATEX
//Q_{boost}(#vec{#beta}) = [ cosh#frac{#xi}{2}, #frac{#vec{#beta}}{#beta}sinh#frac{#xi}{2}]
//END_LATEX
//
// |---------------|
// | Decomposition |
// |---------------|
// Any Lorentz transformation can be decomposed in a pure rotation r followed by
// a pure boost b. When the transformation is represented by a biquaternion Q
// this decomposition can easily be performed:
//BEGIN_LATEX
//Q #equiv U + i V  (U,V ordinary quaternions)
//#mu = #sqrt{UU^{+}}
//r := #frac{U}{#mu}
//s := [ #mu, -i #frac{UV^{+}}{#mu} ] 
//END_LATEX
// Decomposition of a TLorentzQuaternion is done with the methods Decompose(..).
//
//////////////////////////////////////////////////////////////////////////

ClassImp(TLorentzQuaternion)

//_____________________________________________________________________
TLorentzQuaternion::TLorentzQuaternion(const std::complex<double>& q0,
				       const std::complex<double>& q1,
				       const std::complex<double>& q2,
				       const std::complex<double>& q3)
: TObject(), fQ0(q0), fQ1(q1), fQ2(q2), fQ3(q3)
{
  // Private constructor
}

//_____________________________________________________________________
TLorentzQuaternion::TLorentzQuaternion()
  : TObject(), fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize identity transformation
}

//_____________________________________________________________________
TLorentzQuaternion::TLorentzQuaternion(TRootIOCtor* rio)
  : TObject(), fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
TLorentzQuaternion::TLorentzQuaternion(const double theta, const TVector3& axis)
  : TObject(), fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize rotation over angle 'theta' around axis defined by 'axis'
  SetRotation(theta,axis);
}

//_____________________________________________________________________
TLorentzQuaternion::TLorentzQuaternion(const TRotation& rotation)
  : TObject(), fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize rotation
  double theta;
  TVector3 axis;
  rotation.AngleAxis(theta,axis);

  SetRotation(theta,axis);
}

//_____________________________________________________________________
TLorentzQuaternion::TLorentzQuaternion(const TVector3& beta)
  : TObject(), fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize pure Lorentz boost defined by boost vector 'beta'
  SetBoost(beta.X(),beta.Y(),beta.Z());
}

//_____________________________________________________________________
TLorentzQuaternion::TLorentzQuaternion(const double betaX, 
				       const double betaY, 
				       const double betaZ)
  : TObject(), fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize pure Lorentz boost defined by components of boost vector
  SetBoost(betaX,betaY,betaZ);
}

//_____________________________________________________________________
void TLorentzQuaternion::SetBoost(const double betaX, 
				  const double betaY, 
				  const double betaZ)
{
  // Create a pure boost (defined by the components of the boost vector).
  const double beta = sqrt(betaX*betaX
			   +betaY*betaY
			   +betaZ*betaZ);
  const double coshksi2 = sqrt( (1./sqrt(1.-beta*beta) + 1.)/2. );
  const double sinhksi2 = sqrt( (1./sqrt(1.-beta*beta) - 1.)/2. );

  fQ0 = complex<double>(coshksi2);
  if( fabs(beta)>STRANGEUFLOW ) {
    fQ1 = complex<double>(0.,betaX/beta*sinhksi2);
    fQ2 = complex<double>(0.,betaY/beta*sinhksi2);
    fQ3 = complex<double>(0.,betaZ/beta*sinhksi2);
  }
}

//_____________________________________________________________________
void TLorentzQuaternion::SetRotation(const double theta, const TVector3& axis)
{
  // Create a rotation over angle 'theta' around axis defined by 'axis'
  const double norm = axis.Mag();
  if( fabs(norm)<STRANGEUFLOW ) {
    cerr << "ERROR in TLorentzQuaternion::SetRotation(const double, "
	 << "const TVector3&): "
	 << "invalid rotation axis." << endl;
    assert(1==0);
  }
  const double sinThetaOver2 = sin(theta/2.);

  fQ0 = complex<double>(cos(theta/2.));
  fQ1 = complex<double>(axis.X()/norm*sinThetaOver2);
  fQ2 = complex<double>(axis.Y()/norm*sinThetaOver2);
  fQ3 = complex<double>(axis.Z()/norm*sinThetaOver2);
}

//_____________________________________________________________________
ostream& operator<<(ostream& stream, const TLorentzQuaternion& q)
{
  stream << "[" << q.fQ0 << "," << q.fQ1 << "," << q.fQ2 << "," << q.fQ3 << "]";
  return stream;
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::Conjugate()
{
  // Conjugate the biquaternion
  //BEGIN_LATEX #bar{Q} = [Q_{0},-Q_{1},-Q_{2},-Q_{3}]END_LATEX
  static const complex<double> minus(-1.,0.);
  
  fQ1 *= minus;
  fQ2 *= minus;
  fQ3 *= minus;

  return *this;
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::ComplexConjugate()
{
  // Complex conjugate the biquaternion
  //BEGIN_LATEX Q* = [Q*_{0},Q*_{1},Q*_{2},Q*_{3}]END_LATEX
  fQ0 = conj(fQ0);
  fQ1 = conj(fQ1);
  fQ2 = conj(fQ2);
  fQ3 = conj(fQ3);

  return *this;
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::Hermitian()
{
  // Hermitian conjugate the biquaternion
  //BEGIN_LATEX Q^{+} = #bar{Q}*END_LATEX
  static const complex<double> minus(-1.,0.);

  fQ0 = conj(fQ0);
  fQ1 = conj(fQ1) * minus;
  fQ2 = conj(fQ2) * minus;
  fQ3 = conj(fQ3) * minus;

  return *this;
}

//_____________________________________________________________________
TLorentzQuaternion TLorentzQuaternion::Conjugate(const TLorentzQuaternion& rhs)
{
  // Conjugate the biquaternion
  //BEGIN_LATEX #bar{Q} = [Q_{0},-Q_{1},-Q_{2},-Q_{3}]END_LATEX
  TLorentzQuaternion lhs = rhs;
  return lhs.Conjugate();
}

//_____________________________________________________________________
TLorentzQuaternion TLorentzQuaternion::ComplexConjugate(const TLorentzQuaternion& rhs)
{
  // Complex conjugate the biquaternion
  //BEGIN_LATEX Q* = [Q*_{0},Q*_{1},Q*_{2},Q*_{3}]END_LATEX
  TLorentzQuaternion lhs = rhs;
  return lhs.ComplexConjugate();
}

//_____________________________________________________________________
TLorentzQuaternion TLorentzQuaternion::Hermitian(const TLorentzQuaternion& rhs)
{
  // Hermitian conjugate the biquaternion
  //BEGIN_LATEX Q^{+} = #bar{Q}*END_LATEX
  TLorentzQuaternion lhs = rhs;
  return lhs.Hermitian();
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::Invert()
{
  // Invert the Lorentz transformation
  return Conjugate();
}

//_____________________________________________________________________
TLorentzQuaternion TLorentzQuaternion::Inverse() const
{
  // Returns the inverse Lorentz transformation
  return Conjugate(*this);
}

//_____________________________________________________________________
TLorentzQuaternion TLorentzQuaternion::operator*(const TLorentzQuaternion& rhs) const
{
  // Multiply with a biquaternion from the right:
  //BEGIN_LATEX QQ' = [ Q_{0}Q'_{0} - #vec{Q}#upoint#vec{Q}', Q_{0}#vec{Q}' + Q'_{0}#vec{Q} + #vec{Q}#times#vec{Q}'] END_LATEX
  TLorentzQuaternion result;
  
  result.fQ0 = fQ0*rhs.fQ0 - fQ1*rhs.fQ1 - fQ2*rhs.fQ2 - fQ3*rhs.fQ3;
  result.fQ1 = fQ0*rhs.fQ1 + rhs.fQ0*fQ1 + fQ2*rhs.fQ3 - fQ3*rhs.fQ2;
  result.fQ2 = fQ0*rhs.fQ2 + rhs.fQ0*fQ2 + fQ3*rhs.fQ1 - fQ1*rhs.fQ3;
  result.fQ3 = fQ0*rhs.fQ3 + rhs.fQ0*fQ3 + fQ1*rhs.fQ2 - fQ2*rhs.fQ1;

  return result;
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::operator*=(const TLorentzQuaternion& rhs)
{
  // Matrix multiplication
  // Note a *= b; <=> a = a * b; while a.transform(b); <=> a = b * a;
  return operator=( (*this) * rhs );
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::Transform(const TLorentzQuaternion& rhs)
{
  // Chain LorentzTransformations.
  // Note a.transform(b); <=> a = b * a; while a *= b; <=> a = a * b;
  return operator=( rhs * (*this) );
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::Boost(const TVector3& boostVector)
{
  // Chain pure Lorentz boost to transformation.
  return Transform( TLorentzQuaternion(boostVector) );
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::Boost(const double betaX, 
					      const double betaY, 
					      const double betaZ)
{
  // Chain pure Lorentz boost to transformation.
  return Transform( TLorentzQuaternion(betaX,betaY,betaZ) );
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::Rotate(const double theta,
					       const TVector3& axis)
{
  // Chain a rotation over 'theta' around 'axis' to transformation
  return Transform( TLorentzQuaternion(theta,axis) );
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::Rotate(const double theta,
					       const TVector3* axis)
{
  // Chain a rotation over 'theta' around 'axis' to transformation
  return Transform( TLorentzQuaternion(theta,*axis) );
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::RotateX(const double theta)
{
  // Chain a rotation over 'theta' around x-axis to transformation
  return Transform( TLorentzQuaternion(theta,TVector3(1.,0.,0.)) );
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::RotateY(const double theta)
{
  // Chain a rotation over 'theta' around y-axis to transformation
  return Transform( TLorentzQuaternion(theta,TVector3(0.,1.,0.)) );
}

//_____________________________________________________________________
TLorentzQuaternion& TLorentzQuaternion::RotateZ(const double theta)
{
  // Chain a rotation over 'theta' around z-axis to transformation
  return Transform( TLorentzQuaternion(theta,TVector3(0.,0.,1.)) );
}

//_____________________________________________________________________
TLorentzVector TLorentzQuaternion::operator*(const TLorentzVector& v) const
{
  // Transform the 4-vector

  // Define quaternion equivalent of 4-vector
  TLorentzQuaternion q(complex<double>(v.E()),
		       complex<double>(0.,v.X()),
		       complex<double>(0.,v.Y()),
		       complex<double>(0.,v.Z()));

  // Transform the quaternion
  TLorentzQuaternion qp = (*this) * q * Hermitian(*this);

  // Go back from quaternion to 4-vector
  return TLorentzVector(qp.fQ1.imag(),
			qp.fQ2.imag(),
			qp.fQ3.imag(),
			qp.fQ0.real());
}

//_____________________________________________________________________
TVector3 TLorentzQuaternion::operator*(const TVector3& v) const
{
  // Transform the 3-vector
  TLorentzVector vp = (*this) * TLorentzVector(v.X(),v.Y(),v.Z(),0.);
  return TVector3(vp.X(),vp.Y(),vp.Z());
}

//_____________________________________________________________________
bool TLorentzQuaternion::IsIdentity() const
{
  return (fabs(fQ1.real())<STRANGEUFLOW) && (fabs(fQ1.imag())<STRANGEUFLOW)
    && (fabs(fQ2.real())<STRANGEUFLOW) && (fabs(fQ2.imag())<STRANGEUFLOW)
    && (fabs(fQ3.real())<STRANGEUFLOW) && (fabs(fQ3.imag())<STRANGEUFLOW)
    && (fabs(fQ0.imag())<STRANGEUFLOW) && (fabs(fQ0.real()-1.)<STRANGEUFLOW);
}

//_____________________________________________________________________
void TLorentzQuaternion::Decompose(TLorentzQuaternion& rotation, 
				   TLorentzQuaternion& boost) const
{
  // Decompose the current Lorentz transformation into a rotation + boost.
  // Details can be found in J.Lambek, The mathematical intelligencer, 
  // Volume 17, Number 4, 7-15 (1995) (http://dx.doi.org/10.1007/BF03024783)

  // Decompose the quaternion in a real and imaginary part
  TLorentzQuaternion u(fQ0.real(),
		       fQ1.real(),
		       fQ2.real(),
		       fQ3.real());
  TLorentzQuaternion v(fQ0.imag(),
		       fQ1.imag(),
		       fQ2.imag(),
		       fQ3.imag());
  
  // Find the norm of the real part
  const double mu = sqrt( (u*Conjugate(u)).fQ0.real() );

  // Define the rotation
  rotation= u;
  rotation.fQ0 /= mu;
  rotation.fQ1 /= mu;
  rotation.fQ2 /= mu;
  rotation.fQ3 /= mu;

  // Define the boost
  boost = u*Conjugate(v);
  const complex<double> factor(0.,-1./mu);
  boost.fQ0 = mu;
  boost.fQ1 *= factor;
  boost.fQ2 *= factor;
  boost.fQ3 *= factor;
}

//_____________________________________________________________________
void TLorentzQuaternion::Decompose(double& theta, TVector3& axis,
				   TVector3& boost) const
{
  // Decompose the current Lorentz transformation into a rotation + boost.
  // The rotation is given as an angle 'theta' with its rotation axis.
  // The boost is given by its boost vector.

  // Decompose the transformation
  TLorentzQuaternion r,b;
  Decompose(r,b);

  // Find the angle+axis of the rotation
  const double sintheta2 = sqrt(r.fQ1.real()*r.fQ1.real()
				+r.fQ2.real()*r.fQ2.real()
				+r.fQ3.real()*r.fQ3.real());
  theta = 2.*atan2(sintheta2,r.fQ0.real());
  if( fabs(sintheta2)>STRANGEUFLOW )
    axis = TVector3(r.fQ1.real()/sintheta2,
		    r.fQ2.real()/sintheta2,
		    r.fQ3.real()/sintheta2);
  else 
    axis = TVector3(0.,0.,1.);
    
  // Find the boost vector
  const double sinhksi2 = sqrt(b.fQ1.imag()*b.fQ1.imag()
			       +b.fQ2.imag()*b.fQ2.imag()
			       +b.fQ3.imag()*b.fQ3.imag());
  const double beta = tanh(2.*log(b.fQ0.real()+sinhksi2));
  if( fabs(sinhksi2)>STRANGEUFLOW)
    boost = TVector3(beta*b.fQ1.imag()/sinhksi2,
		     beta*b.fQ2.imag()/sinhksi2,
		     beta*b.fQ3.imag()/sinhksi2 );
  else 
    boost = TVector3(0.,0.,0.);
}
