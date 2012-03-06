/*
 * TCstPWF.h
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on July,23 2010
 *
 */

#ifndef TCSTPWF_H
#define TCSTPWF_H

#include "TWavefunctionImplementation.h"

class TCstPWF : public TWavefunctionImplementation
{
 private:
  enum ESpace { // wavefunction in r-space or p-space?
    kR = 1,     // configuration space
    kP = 2      // momentum space
  };
  
  enum EWave { // type of wave
    kU = 0,    // S-wave
    kW = 1,    // D-wave
    kVT = 2,   // triplet P-wave
    kVS = 3    // singlet P-wave
  };

 public:
  TCstPWF(int *order, double *leadingMass, double *stepMass, double *tailMass,
	  int *exponent, double *parU, double *parW, double *parVT, double *parVS);
  TCstPWF(TRootIOCtor*); // ROOT I/O Constructor
  TCstPWF(const TCstPWF&);
  TCstPWF& operator=(const TCstPWF&);
  virtual ~TCstPWF();

  virtual TCstPWF *Clone(const char* ="") const {return new TCstPWF(*this); }

  virtual double GetUr(double r) const; // r in [fm]
  virtual double GetWr(double r) const; // r in [fm]
  virtual double GetVTr(double r) const; // r in [fm]
  virtual double GetVSr(double r) const; // r in [fm]

  virtual double GetUp(double p) const; // p in [MeV]
  virtual double GetWp(double p) const; // p in [MeV]
  virtual double GetVTp(double p) const; // p in [MeV]
  virtual double GetVSp(double p) const; // p in [MeV]
  virtual double GetUpoff(const TVector3& p) const; // p in [MeV]
  virtual double GetWpoff1(const TVector3& p) const; // p in [MeV]
  virtual double GetWpoff2(double pperp2) const; // p in [MeV]

 protected:
  static int AngularMomentum(EWave wave);
  double Mass(EWave wave, int i) const;
  double GetG(ESpace space, EWave wave, int term, double q) const;
 
 private:
  int     *fOrder; //[4] Number of terms in parametrization
  double **fCoeff; //[4][fOrder[i]] expansion coefficients [MeV^-3/2]
  double  *fAlpha; //[4] leading mass [MeV]
  double  *fMx;    //[4] step mass [MeV]
  double  *fMn;    //[4] tail mass [MeV]
  int     *fExp;   //[4] exponent in expansion functions
  
  ClassDef(TCstPWF,1); //  Covariant spectator theory parametrized wavefunction

}; // class TCstPWF

#endif // TCSTPWF_H
