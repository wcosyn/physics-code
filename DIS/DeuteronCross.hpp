#ifndef DEUTERONCROSS_HPP
#define DEUTERONCROSS_HPP

#include "DeuteronMomDistr.hpp"
#include <TKinematics2to2.h>
#include <DeuteronStructure.hpp>
#include <TElectronKinematics.h>

#include <string>

class DeuteronCross{
  
public:
  DeuteronCross(TElectronKinematics elec, std::string name, bool proton, std::string strucname,
    double sigmain, double betain, double epsilonin, double betaoffin, double lambdain, int offshellset);
  ~DeuteronCross();
  double getavgCross(TKinematics2to2 &kin, int pw, double Einoff);
  
private:
  double massi;
  TElectronKinematics electron;
  DeuteronMomDistr momdistr;
  DeuteronStructure structure;  
};

#endif