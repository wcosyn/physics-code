/*
 * TYukawaPWF.h
 *
 * Author: Pieter Vancraeyveld
 *         pieter.vancraeyveld@ugent.be
 *
 * Worked started February,5 2009
 *
 */

#ifndef TYUKAWAPWF_H
#define TYUKAWAPWF_H

#include <TROOT.h>
#include "TWavefunctionImplementation.h"

class TYukawaPWF : public TWavefunctionImplementation
{ 
 public:
  TYukawaPWF(int order, double alpha, double m0, double *c, double *d);
  TYukawaPWF(int order, double *c, double *d, double *m);
  TYukawaPWF(TRootIOCtor*); // ROOT I/O Constructor
  TYukawaPWF(const TYukawaPWF&);
  TYukawaPWF& operator=(const TYukawaPWF&);
  virtual ~TYukawaPWF();

  virtual TYukawaPWF *Clone(const char* ="") const { return new TYukawaPWF(*this); }

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
  double getResidu() const; //get pole residu for on-shell extrapolation
  
 protected:
  double Mass(int i) const;

 private:
  int     fOrder;  // Number of terms in parametrization
  double  fAlpha;  // mass parameter in parametrization [fm^-1]
  double  fM0;     // second mass parameter in parametrization [fm^-1]
  double *fC;      //[fOrder] parameters for s-wave [fm^1/2]
  double *fD;      //[fOrder] parameters for d-wave [fm^1/2]
  double *fM;	//[fOrder] parameters for masses [fm^-1]	

  ClassDef(TYukawaPWF,1); // Yukawa parametrized wavefunction
  
};




#endif
