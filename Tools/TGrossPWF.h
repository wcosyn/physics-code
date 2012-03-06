/*
 * TGrossPWF.h
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on February,12 2009
 *
 */

#ifndef TGROSSPWF_H
#define TGROSSPWF_H

#include "TWavefunctionImplementation.h"

class TGrossPWF : public TWavefunctionImplementation
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
  TGrossPWF(int order, double *alpha, double m0, 
	    double *parU, double *parW, double *parVT, double *parVS);
  TGrossPWF(TRootIOCtor*); // ROOT I/O Constructor
  TGrossPWF(const TGrossPWF&);
  TGrossPWF& operator=(const TGrossPWF&);
  virtual ~TGrossPWF();

  virtual TGrossPWF *Clone(const char* ="") const {return new TGrossPWF(*this); }

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
  double Mass(int L, int i) const;
  double GetK(int L, int term, int index) const;
  double GetG(ESpace space, int L, int term, double q) const;
  double GetF(ESpace space, int L, int index, double q) const;

 private:
  int     fOrder; // Number of terms in parametrization
  double *fAlpha; //[3] first mass parameter for each L-wave [MeV]
  double  fM0;    // second mass parameter [MeV]
  double *fPar;   //[4*fOrder] parameters [MeV^1/2]


  ClassDef(TGrossPWF,1); // Gross parametrized wavefunction
};

#endif
