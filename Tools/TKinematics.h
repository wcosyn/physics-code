/*
 * TKinematics.h
 *
 * a TKinematics object represents a series of kinematic points for
 * photon or electron induced production of 2 particles off a nucleon.
 *
 * We currently implement:
 *   (1)   ===>   g + p   ->  K+ + L0
 *   (2)   ===>   g + p   ->  K+ + S0
 *   (3)   ===>   g + p   ->  K0 + S+
 *   (4)   ===>   g + n   ->  K0 + L0
 *   (5)   ===>   g + n   ->  K0 + S0
 *   (6)   ===>   g + n   ->  K+ + S-
 *   (7)   ===>   g + p   ->  P0 + p
 *   (8)   ===>   g + p   ->  P+ + n
 *   (9)   ===>   g + n   ->  P0 + n
 *   (10)  ===>   g + n   ->  P- + p
 *   (11)  ===>   K- + p  ->  g  + L0
 *   (12)  ===>   K- + p  ->  g  + S0
 *
 * These can be set with TKinematics::SetIsospin(int)  
 *
 * Author:         Pieter Vancraeyveld (pieter.vancraeyveld@UGent.be)
 * 
 */

#ifndef TKINEMATICS_H
#define TKINEMATICS_H

#include <TNamed.h>
#include <TString.h>
#include <TBuffer.h>
// #include <Structures.h>
#include "constants.hpp"


class TKinematics : public TNamed
{
 public:
  TKinematics(TRootIOCtor*); // ROOT I/0 constructor
  TKinematics(const char *name,const char *title,
	      int isospin,const char *formatString,
	      double var1,double var2,double var3);
  TKinematics(const TKinematics&);
  virtual ~TKinematics();
  TKinematics& operator=(const TKinematics&);
  bool          operator==(const TKinematics&) const;

  void          SetVarRange(int varNum,double lowerLimit,double upperLimit,int steps);
  void          FixVariable(int varNumber);
  void          FixVariables();
  int           Next();
  void          GoTo(int totalStep);
  void          GoTo(int stepVar1,int stepVar2,int stepVar3);
  virtual void  Print(Option_t* option="") const;
  int           VariableInfo() const;
  static void   Help();

  void          SetIsospin(int isospin);
  void          SetFormat(const char *formatString);
  void          SetVar(int varNum,double value);

  int           GetIsospin() const     { return fIsospin; }
  const char   *GetFormat() const      { return fFormat; }
  double        GetNucleonMass() const { return fMp; }
  double        GetMesonMass() const   { return fMk; }
  double        GetHyperonMass() const { return fMy; }
  double        GetWlab() const        { return fWlab; }
  double        GetWcm() const         { return fWcm; }
  double        GetKlab() const        { return fKlab; }
  double        GetKcm() const         { return fKcm; }
  double        GetPk() const          { return fPk; }
  double        GetPkcm() const        { return fPk; }
  double        GetPklab() const       { return fPklab; }
  double        GetCosthkcm() const    { return fCosthkcm; }
  double        GetCosthklab() const    { return fCosthklab; }
  double        GetPYlab() const       { return fPYlab; }
  double        GetCosthYlab() const    { return fCosthYlab; }
  double        GetQsquared() const    { return fQsquared; }
  double        GetW() const           { return fW; }
  double        GetS() const           { return fS; }
  double        GetT() const           { return fT; }
  double        GetU() const           { return fU; }
  double        GetXb() const          { return fXb; }
  bool          IsPhysical() const     { return fIsPhysical; }
  double        GetVar(int) const;
  double       *GetVarArray(const TString& varName) const; // user owns the double*
  double       *GetPhysicalVarArray(const TString& varName) const; // user owns the double*
  TString       GetVarName(int varNumber) const;
  int           GetNumberOfSteps(int varNumber) const;
  int           GetNumberOfSteps() const;
  int           GetNumberOfPhysicalSteps() const;
  int           GetNumberOfVariables() const;
  double        GetStepSize(int varNumber) const;
  int           GetStep(int varNumber) const;
  int           GetStep() const        { return fStep; }
  bool          IsFixed() const;
  bool          IsFixed(int varNumber) const;

 protected:
  void          SetMn(double mn) { fMp = mn; }
  void          SetMk(double mk) { fMk = mk; }
  void          SetMy(double my) { fMy = my; }
  void          UpdateKinematics();
  // Underflow parameter for TKinematics
  virtual 
    double      Underflow() const { return STRANGEUFLOW; }

 private:
  void          ReadFormatString(const char*);
  void          ChangeIsospinChannel(int);
  double        EnergyConservation(double);
  
  // Data Members
  // ------------
  char    *fFormat;        //[100] Format string
  int      fIsospin;       // isospin channel [1-6]
  double   fMp;            // nucleon mass
  double   fMk;            // meson mass
  double   fMy;            // hyperon mass

  // Kinematic variables
  double   fWlab;          // photon energy (LAB)
  double   fWcm;           // photon energy (CM)
  double   fKlab;          // photon momentum (LAB)
  double   fKcm;           // photon momentum (CM)
  double   fPk;            // kaon momentum (CM)
  double   fPklab;         // kaon momentum (LAB)
  double   fCosthkcm;      // cos of kaon scattering angle in CM
  double fCosthklab;
  double   fPYlab;	//hyperon mom in lab
  double fCosthYlab;
  double   fQsquared;      // Q^2
  double   fXb;            // Bjorken x
  double   fW;             // invariant mass
  double   fS;             // s (Mandelstam variable)
  double   fT;             // t (Mandelstam variable)
  double   fU;             // u (Mandelstam variable)
  bool     fIsPhysical;    // Does kinematic point lie in physical plane
  mutable 
    int    fNrOfPhysicals; // Number of physical points (-1 when not determined)

  // 3 independent variables
  double **fVar;           //! pointer to independent kinematic variable
  double  *fEVar;          //! points to independent energy variable
  double  *fAVar;          //! points to independent angular variable
  bool    *fIsVar;         //[3] is non-constant kinematic variable?
  double  *fLowerLimit;    //[3] Lower limit non-constant kinematic variable
  double  *fUpperLimit;    //[3] Upper limit non-constant kinematic variable
  double  *fStepSize;      //[3] Step size non-constant kinematic variable
  int     *fNumberOfSteps; //[3] Number of steps non-constant kinematic variable
  int      fStep;          // Step counter (0 or higher)

  ClassDef(TKinematics,4)  // Kinematics for g N -> K Y
};


#endif
