#ifndef DOUBLENCROSS_HPP
#define DOUBLENCROSS_HPP

#include <TElectronKinematics.h>
#include <MeanFieldNucleusThick.hpp>
#include <TKinematics2to3WithLabAngles.h>
#include "DoubleNModel.hpp"

class DoubleNCross{
public:
  DoubleNCross(TElectronKinematics &elec, MeanFieldNucleusThick *pnucl, string dir);
  ~DoubleNCross();
  double getDiffCross(const TKinematics2to3 &kin, bool SRC, bool CT, bool pw, bool corr, int shellindex1, int shellindex2,
		      double phi);
private:
  string homedir;
  TElectronKinematics electron;
  MeanFieldNucleusThick *pnucl;
  double kinfactors[6];
  double response[9];
  double frontfactor;
  double mott;
  DoubleNModel *reacmodel;
  
};



#endif