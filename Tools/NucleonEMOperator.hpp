#ifndef NUCLEONEMOPERATOR_HPP
#define NUCLEONEMOPERATOR_HPP

#include <iostream>
// #include <complex>
#include "constants.hpp"
#include "FourVector.h"
#include "GammaStructure.h"


class NucleonEMOperator{
public:
  NucleonEMOperator();
  NucleonEMOperator(double Q2, bool proton, int para);
  ~NucleonEMOperator();
  double getGE();
  double getGM();
  double getF1();
  double getF2();
  FourVector<GammaStructure> getCC1(FourVector<double> pi, FourVector<double> pf);
  FourVector<GammaStructure> getCC2(FourVector<double> q);
  FourVector<GammaStructure> getCC3(FourVector<double> q, FourVector<double> pi, FourVector<double> pf);
  FourVector<GammaStructure> getCC(int current, FourVector<double> q, FourVector<double> pi, FourVector<double> pf);
  static const FourVector<GammaStructure> gamma_mu;
  static const GammaStructure Id;
  
private:
  bool proton;
  int parametrization;
  double Q2;
  double tau;
  double GE_null;
  double GM_null;
  double GE;
  double GM;
  double F1;
  double F2;  //kappa*F2!!!!  with kappa anomalous magneticmoment
  
  void setGE();
  void setGM();
  double Get_Gdipole(double Q2);
  void setF1();
  void setF2();
};

#endif
