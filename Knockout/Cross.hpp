#ifndef CROSS_HPP
#define CROSS_HPP

#include <TElectronKinematics.h>
#include <MeanFieldNucleusThick.hpp>
#include <TKinematics2to2.h>
#include "Model.hpp"

class Cross{
public:
  Cross(TElectronKinematics &elec, MeanFieldNucleusThick *pnucl, string dir);
  ~Cross();
  double getDiffCross(TKinematics2to2 &kin, int current, int thick, int SRC, int CT, int pw, int shellindex, double phi);
  double getElCross(TKinematics2to2 &kin, int current, double phi);
private:
  string homedir;
  TElectronKinematics electron;
  MeanFieldNucleusThick *pnucl;
  double kinfactors[6];
  double response[6];
  double frontfactor;
  double mott;
  Model *reacmodel;
  
};

#endif