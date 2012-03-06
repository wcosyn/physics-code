/*
 * TKinematics2to2.h
 *
 * Author:        Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * 
 */

#ifndef TKINEMATICS2TO2_H
#define TKINEMATICS2TO2_H

#include "TKinematics.h"

class TKinematics2to2 : public TKinematics
{
 public:
  TKinematics2to2(const char *name,const char *title,
		  double mN, double mK, double mY,
		  const char *formatString,
		  double var1,double var2,double var3);
  TKinematics2to2(TRootIOCtor*); // ROOT I/O Constructor
  
  virtual void SetMasses(double mN, double mK, double mY); // in [MeV]

 protected:
  virtual double Underflow() const { return 1e-6; }
  
  ClassDef(TKinematics2to2,1)    // Kinematics for gamma N -> K Y
    
};


#endif
