/*
 * TKinematics2to3.h
 *
 * TKinematics2to3 solves kinematics for \gamma D -> K Y N
 *
 * The user can pick from different choices of reference frame. 
 * The z axis is always along the photon momentum. One of the final
 * state particles then fixes the x and y axis. The user can 
 * specify this the EAxisConvention type.
 *
 * Internally we take the y axis along P_gamma x P_N
 *
 * All internal variables are in MeV or rad
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 */


#ifndef TKINEMATICS2TO3_H
#define TKINEMATICS2TO3_H

#include "TSpinor.h"
#include "TKinematics2to2.h"
#include "FourVector.h"
#include <TNamed.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TMath.h>

// Forward declarations
class TString;
class TVector3;
class TLorentzQuaternion;

// Global functions for conversion 
// from ROOT 4vectors to strangecalc 4vectors.
// template<typename T>
// FourVector<T> operator*(const TLorentzRotation&, const FourVector<T>&);
// 
// FourVector<double> ToFourVector(const TLorentzVector&);
// TLorentzVector ToLorentzVector(const FourVector<double>&);

class TKinematics2to3 : public TNamed
{
  friend class TKinematics2to3WithLabAngles;

 public:
  enum EAxisConvention { 
    // code to specify the axis convention.
    // z ~ p_gamma and x ~ y x z
    kS = 0, // standard input/output axis convention stored in fAC
    kK = 1, // y ~ p_gamma x p_K
    kY = 2, // y ~ p_gamma x p_Y
    kN = 3  // y ~ p_gamma x p_N (internal convention)
  };

  TKinematics2to3(const char *name, const char *title,
		  int isospin, EAxisConvention ac, 
		  const char *format,
		  double var1, double var2, double var3,
		  double var4, double var5, double var6,
		  double Md=0., double Mn=0., double Mk=0., double My=0.);
  TKinematics2to3(const TKinematics2to3&);
  TKinematics2to3(TRootIOCtor*); // ROOT I/O Constructor
  virtual ~TKinematics2to3();
  TKinematics2to3& operator=(const TKinematics2to3&);

  void           SetVarRange(int varNum,double low,double high,int steps);
  void           FixVariable(int varNum);
  void           FixVariables();
  int            Next();
  void           GoTo(int totalStep);
  void           GoTo(int var1step,int var2step,int var3step,
		      int var4step,int var5step,int var6step);
  int            VariableInfo() const;
  static void    Help();
  bool           operator==(const TKinematics2to3&) const;
//   TString        HashName() const;

  void           SetIsospin(int isospin);
  void           SetFormat(const char *format);
  void           SetVar(int varNum,double value);
  void           SetVar(double,double,double,double,double,double);
  
  const char     *GetFormat() const         { return fFormat; }
  EAxisConvention GetAxisConvention() const { return fAC; }
  bool           IsPhysical() const     { return fIsPhysical; }
  int            GetIsospin() const     { return fIsospin; }
  double         GetMk() const          { return fMk; }
  double         GetMy() const          { return fMy; }
  double         GetMn() const          { return fMn; }
  double         GetMd() const          { return fMd; }
  TVector3       GetK3(EAxisConvention =kS) const;
  TVector3       GetY3(EAxisConvention =kS) const;
  TVector3       GetN3(EAxisConvention =kS) const;
  TVector3       GetG3(EAxisConvention =kS) const;
  TVector3       GetD3(EAxisConvention =kS) const;
  TLorentzVector GetK4(EAxisConvention =kS) const;
  TLorentzVector GetY4(EAxisConvention =kS) const;
  TLorentzVector GetN4(EAxisConvention =kS) const;
  TLorentzVector GetG4(EAxisConvention =kS) const;
  TLorentzVector GetD4(EAxisConvention =kS) const;
  double         GetPk() const          { return fKvector->P(); }
  double         GetPy() const          { return fYvector->P(); }
  double         GetPn() const          { return fNvector->P(); }
  double         GetKlab() const        { return fGvector->P(); }
  double         GetEk() const          { return fKvector->E(); }
  double         GetEy() const          { return fYvector->E(); }
  double         GetEn() const          { return fNvector->E(); }
  double         GetWlab() const        { return fGvector->E(); }
  double         GetThk_deg() const;                     // in degrees [0,180]
  double         GetThy_deg() const;                     // in degrees [0,180]
  double         GetThn_deg() const;                     // in degrees [0,180]
  double         GetThk() const;                         // in rad [0,Pi]
  double         GetThy() const;                         // in rad [0,Pi]
  double         GetThn() const;                         // in rad [0,Pi]
  double         GetCosthk() const      { return fKvector->CosTheta(); }
  double         GetCosthy() const      { return fYvector->CosTheta(); }
  double         GetCosthn() const      { return fNvector->CosTheta(); }
  double         GetPhik_deg(EAxisConvention =kS) const; // in degrees [0,360[
  double         GetPhiy_deg(EAxisConvention =kS) const; // in degrees [0,360[
  double         GetPhin_deg(EAxisConvention =kS) const; // in degrees [0,360[
  double         GetPhik(EAxisConvention =kS) const;     // in rad [0,2Pi[
  double         GetPhiy(EAxisConvention =kS) const;     // in rad [0,2Pi[
  double         GetPhin(EAxisConvention =kS) const;     // in rad [0,2Pi[
  double         GetThky_deg() const;                    // in degrees [0,180]
  double         GetThkn_deg() const;                    // in degrees [0,180]
  double         GetThyn_deg() const;                    // in degrees [0,180]
  double         GetThky() const;                        // in rad [0,Pi]
  double         GetThkn() const;                        // in rad [0,Pi]
  double         GetThyn() const;                        // in rad [0,Pi]
  double         GetCosthky() const     { return TMath::Cos(GetThky()); }
  double         GetCosthkn() const     { return TMath::Cos(GetThkn()); }
  double         GetCosthyn() const     { return TMath::Cos(GetThyn()); }
  double         GetQsquared() const    { return -1*fGvector->Mag2(); }
  double         GetWtot() const        { return TMath::Sqrt(GetStot()); }
  double         GetWky() const         { return fKYsystem->GetW(); }
  double         GetWkn() const         { return fKNsystem->GetW(); }
  double         GetWyn() const         { return fYNsystem->GetW(); }
  double         GetStot() const;
  double         GetSky() const         { return fKYsystem->GetS(); }
  double         GetSkn() const         { return fKNsystem->GetS(); }
  double         GetSyn() const         { return fYNsystem->GetS(); }
  double         GetTkg() const;
  double         GetTyg() const;
  double         GetTng() const;
  double         GetTkd() const;
  double         GetTyd() const;
  double         GetTnd() const;
  double         GetKYsystem_Pk() const { return fKYsystem->GetPk(); }
  double         GetKYsystem_Py() const { return fKYsystem->GetPk(); }
  double         GetKYsystem_Ek() const;
  double         GetKYsystem_Ey() const;
  double         GetKYsystem_S() const  { return fKYsystem->GetS(); }
  double         GetKYsystem_W() const  { return fKYsystem->GetW(); }
  double         GetKYsystem_costhkcm() const;
  double         GetKYsystem_costhycm() const;
  double         GetKYsystem_thkcm() const;
  double         GetKYsystem_thycm() const;
  double         GetKYsystem_phikcm(EAxisConvention =kS) const;
  double         GetKYsystem_phiycm(EAxisConvention =kS) const;
  double         GetKYsystem_phikcm_deg(EAxisConvention =kS) const;
  double         GetKYsystem_phiycm_deg(EAxisConvention =kS) const;
  double         GetKNsystem_Pk() const { return fKNsystem->GetPk(); }
  double         GetKNsystem_Pn() const { return fKNsystem->GetPk(); }
  double         GetKNsystem_Ek() const;
  double         GetKNsystem_En() const;
  double         GetKNsystem_S() const  { return fKNsystem->GetS(); }
  double         GetKNsystem_W() const  { return fKNsystem->GetW(); }
  double         GetKNsystem_costhkcm() const;
  double         GetKNsystem_costhncm() const;
  double         GetKNsystem_thkcm() const;
  double         GetKNsystem_thncm() const;
  double         GetKNsystem_phikcm(EAxisConvention =kS) const;
  double         GetKNsystem_phincm(EAxisConvention =kS) const;
  double         GetKNsystem_phikcm_deg(EAxisConvention =kS) const;
  double         GetKNsystem_phincm_deg(EAxisConvention =kS) const;
  double         GetYNsystem_Py() const { return fYNsystem->GetPk(); }
  double         GetYNsystem_Pn() const { return fYNsystem->GetPk(); }
  double         GetYNsystem_Ey() const;
  double         GetYNsystem_En() const;
  double         GetYNsystem_S() const  { return fYNsystem->GetS(); }
  double         GetYNsystem_W() const  { return fYNsystem->GetW(); }
  double         GetYNsystem_costhycm() const;
  double         GetYNsystem_costhncm() const; 
  double         GetYNsystem_thycm() const;
  double         GetYNsystem_thncm() const;
  double         GetYNsystem_phiycm(EAxisConvention =kS) const;
  double         GetYNsystem_phincm(EAxisConvention =kS) const;
  double         GetYNsystem_phiycm_deg(EAxisConvention =kS) const;
  double         GetYNsystem_phincm_deg(EAxisConvention =kS) const;
  TLorentzVector BoostToKYcm(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzVector BoostToKNcm(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzVector BoostToYNcm(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzVector BoostToGNcm(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzRotation BoostToKYcm(EAxisConvention =kS) const;
  TLorentzRotation BoostToKNcm(EAxisConvention =kS) const;
  TLorentzRotation BoostToYNcm(EAxisConvention =kS) const;
  TLorentzRotation BoostToGNcm(EAxisConvention =kS) const;
  TLorentzVector BoostFromKYcm(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzVector BoostFromKNcm(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzVector BoostFromYNcm(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzVector BoostFromGNcm(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzRotation BoostFromKYcm(EAxisConvention =kS) const;
  TLorentzRotation BoostFromKNcm(EAxisConvention =kS) const;
  TLorentzRotation BoostFromYNcm(EAxisConvention =kS) const;
  TLorentzRotation BoostFromGNcm(EAxisConvention =kS) const;
  TSpinor        BoostToKYcm(const TSpinor&, EAxisConvention =kS) const;
  TSpinor        BoostToKNcm(const TSpinor&, EAxisConvention =kS) const;
  TSpinor        BoostToYNcm(const TSpinor&, EAxisConvention =kS) const;
  TSpinor        BoostToGNcm(const TSpinor&, EAxisConvention =kS) const;
  TSpinor::LorentzRotation BoostSpinorToKYcm(EAxisConvention =kS) const;
  TSpinor::LorentzRotation BoostSpinorToKNcm(EAxisConvention =kS) const;
  TSpinor::LorentzRotation BoostSpinorToYNcm(EAxisConvention =kS) const;
  TSpinor::LorentzRotation BoostSpinorToGNcm(EAxisConvention =kS) const;
  TSpinor        BoostFromKYcm(const TSpinor&, EAxisConvention =kS) const;
  TSpinor        BoostFromKNcm(const TSpinor&, EAxisConvention =kS) const;
  TSpinor        BoostFromYNcm(const TSpinor&, EAxisConvention =kS) const;
  TSpinor        BoostFromGNcm(const TSpinor&, EAxisConvention =kS) const;
  TSpinor::LorentzRotation BoostSpinorFromKYcm(EAxisConvention =kS) const;
  TSpinor::LorentzRotation BoostSpinorFromKNcm(EAxisConvention =kS) const;
  TSpinor::LorentzRotation BoostSpinorFromYNcm(EAxisConvention =kS) const;
  TSpinor::LorentzRotation BoostSpinorFromGNcm(EAxisConvention =kS) const;
  TLorentzQuaternion BoostQuaternionToKYcm(EAxisConvention =kS) const;
  TLorentzQuaternion BoostQuaternionToKNcm(EAxisConvention =kS) const;
  TLorentzQuaternion BoostQuaternionToYNcm(EAxisConvention =kS) const;
  TLorentzQuaternion BoostQuaternionToGNcm(EAxisConvention =kS) const;
  TLorentzQuaternion BoostQuaternionFromKYcm(EAxisConvention =kS) const;
  TLorentzQuaternion BoostQuaternionFromKNcm(EAxisConvention =kS) const;
  TLorentzQuaternion BoostQuaternionFromYNcm(EAxisConvention =kS) const;
  TLorentzQuaternion BoostQuaternionFromGNcm(EAxisConvention =kS) const;
  TLorentzRotation RotateToKac(EAxisConvention =kS) const;
  TLorentzRotation RotateToYac(EAxisConvention =kS) const;
  TLorentzRotation RotateToNac(EAxisConvention =kS) const;
  TLorentzVector RotateToKac(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzVector RotateToYac(const TLorentzVector&, EAxisConvention =kS) const;
  TLorentzVector RotateToNac(const TLorentzVector&, EAxisConvention =kS) const;
  TSpinor::LorentzRotation RotateSpinorToKac(EAxisConvention =kS) const;
  TSpinor::LorentzRotation RotateSpinorToYac(EAxisConvention =kS) const;
  TSpinor::LorentzRotation RotateSpinorToNac(EAxisConvention =kS) const;
  TSpinor RotateToKac(const TSpinor&, EAxisConvention =kS) const;
  TSpinor RotateToYac(const TSpinor&, EAxisConvention =kS) const;
  TSpinor RotateToNac(const TSpinor&, EAxisConvention =kS) const;
  double         GetVar(int) const;
  double         GetVar(const TString& varName) const;
  virtual double *GetVarArray(const TString& varName) const; // user owns the double*
  virtual double *GetPhysicalVarArray(const TString& varName) const; // user owns the double*
  TString        GetVarName(int varNum) const;
  int            GetNumberOfSteps(int varNum) const;
  int            GetNumberOfSteps() const;
  virtual int    GetNumberOfPhysicalSteps() const;
  int            GetNumberOfVariables() const;
  double         GetStepSize(int varNum) const;
  int            GetStep(int varNum) const;
  int            GetStep() const { return fStep; }
  bool           IsFixed() const;
  bool           IsFixed(int varNum) const;

  double GetPhotonThreshold() const;
  void   GetBoundsMomentumParticle1(double th1, double& lower, double& upper) const;
  double GetLowerBoundMomentumParticle1(double th1) const;
  double GetUpperBoundMomentumParticle1(double th1) const;
  double GetMaximumScatterAngleParticle1() const;

 protected:
  virtual void   UpdateKinematics();
  void           SetMd(double md) { fMd = md; }
  void           SetMk(double mk) { fMk = mk; }
  void           SetMy(double my) { fMy = my; }
  void           SetMn(double mn) { fMn = mn; }

 private:
  virtual void   ReadFormatString(const char*);
  void           ChangeIsospinChannel(int);
  void           InitializeSubsystems();
  void           UpdateSubsystems();
  void           UpdateVariables();
  double         GetVarByGetter(const TString& variable) const;
  void           SetAxisConvention(EAxisConvention);
  const char    *GetSolutionType() const { return fKinSolution; }

  // Wrap kinematics getters that take an argument
  double         GetPhik_deg_kS() const            { return GetPhik_deg(); }
  double         GetPhiy_deg_kS() const            { return GetPhiy_deg(); }
  double         GetPhin_deg_kS() const            { return GetPhin_deg(); }
  double         GetPhik_kS() const                { return GetPhik(); }
  double         GetPhiy_kS() const                { return GetPhiy(); }
  double         GetPhin_kS() const                { return GetPhin(); }
  double         GetKYsystem_phikcm_kS() const     { return GetKYsystem_phikcm(); }
  double         GetKYsystem_phiycm_kS() const     { return GetKYsystem_phiycm(); }
  double         GetKYsystem_phikcm_deg_kS() const { return GetKYsystem_phikcm_deg(); }
  double         GetKYsystem_phiycm_deg_kS() const { return GetKYsystem_phiycm_deg(); }
  double         GetKNsystem_phikcm_kS() const     { return GetKNsystem_phikcm(); }
  double         GetKNsystem_phincm_kS() const     { return GetKNsystem_phincm(); }
  double         GetKNsystem_phikcm_deg_kS() const { return GetKNsystem_phikcm_deg(); }
  double         GetKNsystem_phincm_deg_kS() const { return GetKNsystem_phincm_deg(); }
  double         GetYNsystem_phiycm_kS() const     { return GetYNsystem_phiycm(); }
  double         GetYNsystem_phincm_kS() const     { return GetYNsystem_phincm(); }
  double         GetYNsystem_phiycm_deg_kS() const { return GetYNsystem_phiycm_deg(); }
  double         GetYNsystem_phincm_deg_kS() const { return GetYNsystem_phincm_deg(); }

  // Data Members
  // ------------
 private:
  static 
    const double   kUflow;         // underflow parameter (=1e-7)
  char            *fFormat;        //[100] Format string
  char            *fKinSolution;   //[4] How to solve kinematics
  int              fIsospin;       // isospin channel [1-6]
  EAxisConvention  fAC;            // Standard input/output axis convention
  double           fMd;            // deuteron mass
  double           fMk;            // kaon mass
  double           fMy;            // hyperon mass
  double           fMn;            // nucleon mass
  
  // Kinematics variables (in the deuteron rest frame)
  bool             fIsPhysical;    // Does kinematic point lie in physical plane
  double           fQsquared;      // Q^2
  double           fWlab;          // photon energy
  double           fKlab;          // photon momentum
  double           fEk;            // Kaon Energy
  double           fEy;            // Hyperon Energy
  double           fEn;            // Nucleon Energy
  double           fPk;            // Kaon Momentum
  double           fPy;            // Hyperon Momentum
  double           fPn;            // Nucleon Momentum
  double           fThk;           // Kaon Polar angle [rad]
  double           fThy;           // Hyperon Polar angle [rad]
  double           fThn;           // Nucleon Polar angle [rad]
  double           fCosthk;        // Cosine of kaon polar angle
  double           fCosthy;        // Cosine of hyperon polar angle
  double           fCosthn;        // Cosine of nucleon polar angle
  double           fPhik;          // Kaon azimuthal angle [rad]
  double           fPhiy;          // Hyperon azimuthal angle [rad]
  double           fPhin;          // Nucleon azimuthal angle [rad]
  double           fWtot;          // Total invariant mass
  double           fWky;           // Invariant mass KY-system
  double           fWkn;           // Invariant mass KN-system
  double           fWyn;           // Invariant mass YN-system
  double           fStot;          // Mandelstam s
  double           fSky;           // Mandelstam s KY-system
  double           fSkn;           // Mandelstam s KN-system
  double           fSyn;           // Mandelstam s YN-system
  double           fTgk;           // Momentum transfer photon to kaon
  double           fTgy;           // Momentum transfer photon to hyperon
  double           fTgn;           // Momentum transfer photon to nucleon
  double           fTdk;           // Momentum transfer deuteron to kaon
  double           fTdy;           // Momentum transfer deuteron to hyperon
  double           fTdn;           // Momentum transfer deuteron to nucleon
  double           fKYthkcm;       // kaon scatt. angle in KY CM [rad]
  double           fKYthycm;       // hyperon scatt. angle in KY CM [rad]
  double           fKYcosthkcm;    // kaon scatt. angle in KY CM
  double           fKYcosthycm;    // hyperon scatt. angle in KY CM
  double           fKYphikcm;      // kaon azim. angle in KY CM [rad]
  double           fKYphiycm;      // hyperon azim. angle in KY CM [rad]
  double           fKNthkcm;       // kaon scatt. angle in KN CM [rad]
  double           fKNthncm;       // nucleon scatt. angle in KN CM [rad]
  double           fKNcosthkcm;    // kaon scatt. angle in KN CM
  double           fKNcosthncm;    // nucleon scatt. angle in KN CM
  double           fKNphikcm;      // kaon azim. angle in KN CM [rad]
  double           fKNphincm;      // nucleon azim. angle in KN CM [rad]
  double           fYNthycm;       // hyperon scatt. angle in YN CM [rad]
  double           fYNthncm;       // nucleon scatt. angle in YN CM [rad]
  double           fYNcosthycm;    // hyperon scatt. angle in YN CM
  double           fYNcosthncm;    // nucleon scatt. angle in YN CM
  double           fYNphiycm;      // hyperon azim. angle in YN CM [rad]
  double           fYNphincm;      // nucleon azim. angle in YN CM [rad]
     
  // Subsystems
  TKinematics2to2 *fKYsystem;      //-> kaon-hyperon subkinematics
  TKinematics2to2 *fKNsystem;      //-> kaon-nucleon subkinematics
  TKinematics2to2 *fYNsystem;      //-> hyperon-nucleon subkinematics

  // 4vectors for all particles (stored in internal reference frame)
  TLorentzVector  *fKvector;       //-> kaon 4vector
  TLorentzVector  *fYvector;       //-> hyperon 4vector
  TLorentzVector  *fNvector;       //-> nucleon 4vector
  TLorentzVector  *fGvector;       //-> photon 4vector
  TLorentzVector  *fDvector;       //-> deuteron 4vector
  
  // 6 independent variables
  double         **fVar;           //! pointer to independent kinematic variable
  double          *fE1Var;         //! points to independent variable
  double          *fTh1Var;        //! points to independent variable
  double          *fE2Var;         //! points to independent variable
  double          *fTh2Var;        //! points to independent variable
  double          *fPhi2Var;       //! points to independent variable
  double          *fEVar;          //! points to independent variable
  bool            *fIsVar;         //[6] is non-constant kinematic variable?
  double          *fLowerLimit;    //[6] Lower limit non-constant kinematic variable
  double          *fUpperLimit;    //[6] Upper limit non-constant kinematic variable
  double          *fStepSize;      //[6] Step size non-constant kinematic variable
  int             *fNumberOfSteps; //[6] Number of steps non-constant kinematic variable
  int              fStep;          // Step counter (0 or higher)
  mutable int      fNrOfPhysicals; // Number of physical steps (-1 if not determined)
  
  
  ClassDef(TKinematics2to3,1)     // Kinematics for gamma D -> K Y N
    
};

// //_____________________________________________________________________
// template<typename T>
// FourVector<T> operator*(const TLorentzRotation& rot, const FourVector<T>& v)
// {
//   // Perform a Lorentz transformation on a strangecalc FourVector object.
//   FourVector<T> rotated;
// 
//   // In principle we need a simple matrix multiplication:
//   // rotated(i) = sum(j) rot(i,j) * v(j)
//   //
//   // Unfortunately, FourVector and ROOT label vectors differently
//   //    | FourVector  |  ROOT  |  (ROOT+1)%4
//   //    | ------------|--------|-------------
//   //  T |     0       |    3   |      0
//   //  X |     1       |    0   |      1
//   //  Y |     2       |    1   |      2
//   //  Z |     3       |    2   |      3
// 
//   for(int rooti=0; rooti<4; ++rooti) {
//     for(int rootj=0; rootj<4; ++rootj) {
//       rotated[(rooti+1)%4] += rot(rooti,rootj) * v[(rootj+1)%4];
//     }
//   }
//   
//   return rotated;
// }

#endif
