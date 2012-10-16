#ifndef NUCLEONEMOPERATOR_HPP
#define NUCLEONEMOPERATOR_HPP

#include <iostream>
// #include <complex>
#include "constants.hpp"
#include "FourVector.h"
#include "GammaStructure.h"

class MeanFieldNucleusThick;

class NucleonEMOperator{
public:
  NucleonEMOperator();
  NucleonEMOperator(const double Q2, const bool proton, const int para);
  ~NucleonEMOperator();
  double getGE() const{ return GE;}
  double getGM() const{return GM;}
  double getF1() const{return F1;}
  double getF2() const{return F2;}
  double getGE(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{ return GE*getEmod(r,medium,nucleus);}
  double getGM(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{ return GM*getMmod(r,medium,nucleus);}
  double getF1(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  double getF2(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  FourVector<GammaStructure> getCC1(const FourVector<double> &pi, const FourVector<double> &pf, 
				    const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  FourVector<GammaStructure> getCC2(const FourVector<double> &q, const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  FourVector<GammaStructure> getCC3(const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				    const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  FourVector<GammaStructure> getCC(const int current, const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf, 
				   const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  static const FourVector<GammaStructure> gamma_mu;
  static const GammaStructure Id;
  
private:
  bool proton;
  int parametrization;
  double Q2;
  int Q2index;
  double Q2_interp;
  double tau;
  double GE_null;
  double GM_null;
  double GE;
  double GM;
  double F1;
  double F2;  //kappa*F2!!!!  with kappa anomalous magneticmoment
  
  void setGE();
  void setGM();
  double Get_Gdipole(const double Q2);
  void setF1();
  void setF2();
  
  double getEmod(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  double getMmod(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const;
  
  
  static const double QMCGE[61][15];
  static const double QMCGM[61][15];
  static const double CQSMGE[61][16];
  static const double CQSMGM[61][16];
  
};

#endif
