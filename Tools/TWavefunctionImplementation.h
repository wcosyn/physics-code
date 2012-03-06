/*
 * TWavefunctionImplementation.h
 *
 * Author: Pieter Vancraeyveld
 *         pieter.vancraeyveld@ugent.be
 *
 * Worked started February,5 2009
 */

#ifndef TWAVEFUNCTIONIMPLEMENTATION_H
#define TWAVEFUNCTIONIMPLEMENTATION_H

#include <TROOT.h>
#include <TObject.h>
#include <TVector3.h>

class TWavefunctionImplementation : public TObject
{
 public:
  static const double fgHbarc; // hbar*c [fm MeV]

 public:
  TWavefunctionImplementation(TRootIOCtor*); // ROOT I/O Constructor
  virtual ~TWavefunctionImplementation();

  // pure virtual copy constructor
  virtual TWavefunctionImplementation *Clone(const char* ="") const =0;

  // L=0 wave in configuration space (pure virtual)
  virtual double GetUr(double r) const =0; // r in [fm]

  // L=2 wave in configuration space (pure virtual)
  virtual double GetWr(double r) const =0; // r in [fm]

  // L=1 triplet wave in configuration space (pure virtual)
  virtual double GetVTr(double r) const =0; // r in [fm]

  // L=1 singlet wave in configuration space (pure virtual)
  virtual double GetVSr(double r) const =0; // r in [fm]

  // L=0 wave in momentum space (pure virtual)
  virtual double GetUp(double p) const =0; // p in [MeV]

  // L=2 wave in momentum space (pure virtual)
  virtual double GetWp(double p) const =0; // p in [MeV]

  // L=1 triplet wave in momentum space (pure virtual)
  virtual double GetVTp(double p) const =0; // p in [MeV]

  // L=1 singlet wave in momentum space (pure virtual)
  virtual double GetVSp(double p) const =0; // p in [MeV]

  // L=0 off-shell wave in momentum space (pure virtual)
  virtual double GetUpoff(const TVector3& p) const =0; // p in [MeV]

  // L=2 off-shell wave in momentum space (pure virtual)
  virtual double GetWpoff1(const TVector3& p) const =0; // p in [MeV]

  // L=2 off-shell wave in momentum space (pure virtual)
  virtual double GetWpoff2(const double pperp2) const =0; // p in [MeV]


  virtual double Radial_p(int l, double p) const; // p in [MeV]
  virtual double Radial_r(int l, double r) const; // r in [fm]


  ClassDef(TWavefunctionImplementation,1); // Interface for wavefunction implementation
};

#endif



