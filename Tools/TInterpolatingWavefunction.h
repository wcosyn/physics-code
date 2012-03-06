/*
 * TInterpolatingWavefunction.h
 *
 * Author: Pieter Vancraeyveld
 *         pieter.vancraeyveld@ugent.be
 *
 * Worked started February,16 2009
 *
 */

#ifndef TINTERPOLATINGWAVEFUNCTION_H
#define TINTERPOLATINGWAVEFUNCTION_H

#include <map>
#include "TWavefunctionImplementation.h"

class TInterpolatingWavefunction : public TWavefunctionImplementation
{
 public:
  TInterpolatingWavefunction();
  TInterpolatingWavefunction(const TInterpolatingWavefunction&);
  TInterpolatingWavefunction& operator=(const TInterpolatingWavefunction&);
  virtual ~TInterpolatingWavefunction();

  virtual TInterpolatingWavefunction *Clone(const char* ="") const
  { return new TInterpolatingWavefunction(*this); }

  void           SetUr(int n,double *r, double *list); // r in [fm], wf in [fm^-1/2]
  void           SetWr(int n,double *r, double *list); // r in [fm], wf in [fm^-1/2]
  void           SetVTr(int n,double *r, double *list); // r in [fm], wf in [fm^-1/2]
  void           SetVSr(int n,double *r, double *list); // r in [fm], wf in [fm^-1/2]
  
  void           SetUp(int n,double *p, double *list); // p in [MeV], wf in [MeV^-3/2]
  void           SetWp(int n,double *p, double *list); // p in [MeV], wf in [MeV^-3/2]
  void           SetVTp(int n,double *p, double *list); // p in [MeV], wf in [MeV^-3/2]
  void           SetVSp(int n,double *p, double *list); // p in [MeV], wf in [MeV^-3/2]

  void           AddUr(double r, double value); // r in [fm], wf in [fm^-1/2]
  void           AddWr(double r, double value); // r in [fm], wf in [fm^-1/2]
  void           AddVTr(double r, double value); // r in [fm], wf in [fm^-1/2]
  void           AddVSr(double r, double value); // r in [fm], wf in [fm^-1/2]
  
  void           AddUp(double p, double value); // p in [MeV], wf in [MeV^-3/2]
  void           AddWp(double p, double value); // p in [MeV], wf in [MeV^-3/2]
  void           AddVTp(double p, double value); // p in [MeV], wf in [MeV^-3/2]
  void           AddVSp(double p, double value); // p in [MeV], wf in [MeV^-3/2]

  void           RemoveUr();
  void           RemoveWr();
  void           RemoveVTr();
  void           RemoveVSr();

  void           RemoveUp();
  void           RemoveWp();
  void           RemoveVTp();
  void           RemoveVSp();

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
  double Interpolate(std::map<double,double>* list, double q) const;

 private:
  std::map<double,double> *fUr; // list of L=0 wave in r-space
  std::map<double,double> *fWr; // list of L=2 wave in r-space
  std::map<double,double> *fVTr; // list of L=1 triplet wave in r-space
  std::map<double,double> *fVSr; // list of L=1 singlet wave in r-space

  std::map<double,double> *fUp; // list of L=0 wave in p-space
  std::map<double,double> *fWp; // list of L=2 wave in p-space
  std::map<double,double> *fVTp; // list of L=1 triplet wave in p-space
  std::map<double,double> *fVSp; // list of L=1 singlet wave in p-space

  ClassDef(TInterpolatingWavefunction,1); // Interpolated wavefunction
};

#endif
