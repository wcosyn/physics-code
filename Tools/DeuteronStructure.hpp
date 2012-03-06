#ifndef DEUTERONSTRUCTURE_HPP
#define DEUTERONSTRUCTURE_HPP

#include "TKinematics2to2.h"
#include "TElectronKinematics.h"
#include <string>

using namespace std;

class DeuteronStructure{

public:
  DeuteronStructure(TElectronKinematics &el, int proton, string name);
  void getStructureFunctions(TKinematics2to2 &kin, double &FL, double &FT, double &FTT, double &FTL, double Einoff) const;
  double getStructure(TKinematics2to2 &kin, double phi, double Einoff) const;
  double getavgStructure(TKinematics2to2 &kin, double Einoff) const;
  double getInclStructure(TKinematics2to2 &kin, double Einoff) const;
  
private:
  TElectronKinematics electron;
  int proton;
  string name;
  double massi; 
};

#endif