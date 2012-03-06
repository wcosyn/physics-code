/*
 * TKinematics2to3WithLabAngles.h
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: August, 20 2009
 *
 */

#ifndef TKINEMATICS2TO3WITHLABANGLES_H
#define TKINEMATICS2TO3WITHLABANGLES_H

#include "TKinematics2to3.h"
#include <fstream>

class TKinematics2to3WithLabAngles : public TKinematics2to3
{
 public:
  TKinematics2to3WithLabAngles(const char *name, const char *title,
			       int isospin, EAxisConvention ac, 
			       const char *format,
			       double var1, double var2, double var3,
			       double var4, double var5, double var6,
			       double Md =0.,double Mn=0., double Mk=0.,double My=0.);
  TKinematics2to3WithLabAngles(const TKinematics2to3WithLabAngles&);
  TKinematics2to3WithLabAngles(TRootIOCtor*); // ROOT I/O Constructor
  TKinematics2to3WithLabAngles& operator=(const TKinematics2to3WithLabAngles&);
  virtual ~TKinematics2to3WithLabAngles();

  int GetNumberOfSolutions() const { return fNrOfSolutions; }
  const TKinematics2to3& GetSolution(int) const;
  virtual int     GetNumberOfPhysicalSteps() const;
  virtual double *GetVarArray(const TString& varName) const; // user owns the double*
  virtual double *GetPhysicalVarArray(const TString& varName) const; // user owns the double*
 
 protected:
  virtual void ReadFormatString(const char*);
  virtual void UpdateKinematics();
  static void TurnOutputOff(bool=true);

 private:
  int fNrOfSolutions;           //! 0, 1 or 2 solutions for current Omega_2
  TKinematics2to3 **fSolutions; //! the real solutions
  static const std::ofstream fgDevNull; // stream to /dev/null

  ClassDef(TKinematics2to3WithLabAngles,1);
  
}; // class TKinematics2to3WithLabAngles

#endif // TKINEMATICS2TO3WITHLABANGLES_H
