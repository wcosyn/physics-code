/*
 * TElectronKinematics.h
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: March,24 2010
 *
 */

#ifndef TELECTRONKINEMATICS_H
#define TELECTRONKINEMATICS_H

#include "TObject.h"
#include "FourVector.h"
// Forward declarations
class TKinematics2to2;
class TKinematics2to3;

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TElectronKinematics                                                   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

class TElectronKinematics : public TObject
{
 public:
  TElectronKinematics(TRootIOCtor*); // ROOT I/O Constructor
  TElectronKinematics(const TElectronKinematics&);
  TElectronKinematics& operator=(const TElectronKinematics&);
  virtual ~TElectronKinematics() {} // Destructor

  virtual TElectronKinematics* Clone(const char * ="") const { return new TElectronKinematics(*this); }

  static TElectronKinematics* CreateWithBeamEnergy(double);
  static TElectronKinematics* CreateWithEpsilon(double);
  static TElectronKinematics* CreateWithCosScatterAngle(double);

  void SetBeamEnergy(double); // in [MeV]
  void SetEpsilon(double);
  void SetCosScatterAngle(double); // in ]-1,1]

  double GetBeamEnergy(const TKinematics2to2&) const;
  double GetEpsilon(const TKinematics2to2&) const;
  double GetCosScatterAngle(const TKinematics2to2&) const;
  double GetTan2HalfAngle(const TKinematics2to2&) const;
  double GetBeamEnergy(const double Q2, const double nu) const;
  double GetEpsilon(const double Q2, const double nu) const;
  double GetCosScatterAngle(const double Q2, const double nu) const;
  double GetTan2HalfAngle(const double Q2, const double nu) const;
  double GetBeamEnergy(const TKinematics2to3&) const;
  double GetEpsilon(const TKinematics2to3&) const;
  double GetCosScatterAngle(const TKinematics2to3&) const;
  double GetTan2HalfAngle(const TKinematics2to3&) const;

  /*! constructs the four vectors for the incoming and outgoing lepton
   * \param kin kinematics of the hadron piece
   * \param[out] k_in [MeV] incoming lepton fourvector
   * \param[out] k_out [MeV] outgoing lepton fourvector
   */
  void GetLeptonVectors(const TKinematics2to2& kin, FourVector<double> &k_in, FourVector<double> &k_out); 
  double GetLeptonMass() const{return 0.511;} /*!< [MeV] returns lepton mass */
  
  
 protected:
  TElectronKinematics();
  bool SolveKinematics(const TKinematics2to2&) const;
  bool SolveKinematics(const double Q2, const double nu) const;
  bool SolveKinematics(const TKinematics2to3&) const;
  
 private:
  enum InputVariable { kBeamEnergy = 0,
		       kEpsilon = 1,
		       kScatterAngle = 2 };

  InputVariable  fInput;        // specifies the input variable
  mutable double fBeamEnergy;   // energy incoming electron beam in [MeV]
  mutable double fEpsilon;      // virtual photon's transverse linear polarization
  mutable double fScatterAngle; // cos of electron scattering angle in ]-1,1]
  mutable double ftan2HalfAngle;

  ClassDef(TElectronKinematics,1); // Abstraction of an eD->e'KYN observable

}; // class TElectronKinematics

#endif // TELECTRONKINEMATICS_H
