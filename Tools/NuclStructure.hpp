#ifndef NUCLSTRUCTURE_HPP
#define NUCLSTRUCTURE_HPP

#include <string>

using namespace std;

class NuclStructure{
  
public:
  NuclStructure(bool proton, double var1, double var2, int switchvar, string name);
  NuclStructure(bool proton, double Q2in, double xin, double Wsqin, string name);
  void getF_Alekhin(double &F1, double &F2);
  double getF1_Alekhin();
  double getF2_Alekhin();
  void getF_SLAC(double &F1, double &F2) const;
  double getF1_SLAC() const;
  double getF2_SLAC() const;
  void getF_CB(double &F1, double &F2) const;
  double getF1_CB() const;
  double getF2_CB() const;
  void getF(double &F1, double &F2);
  double getF1();
  double getF2();
  const string getName() const{return name;}
  
  
private:
  string name;
  bool proton;
  double mass;
  string dir;
  double x;
  double Q2;
  double Wsq;
  
  double f2p_b(double massi, double xp, double q2, double fm) const;
  double f2n_b(double massi, double xp, double q2, double fm) const;
  
  double bodek(double wm, double qsq) const;
};

#endif