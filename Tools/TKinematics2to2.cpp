/*
 * TKinematics2to2.cpp
 *
 * Author:        Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * 
 */

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TKinematics2to2                                                       //
//                                                                       //
// TKinematics2to2 allows to calculate kinematics for 2 -> 2 scattering. //
//                                                                       //
// All quantities are in MeV or rad.                                     //
//                                                                       //
// TKinematics2to2 is an extension of the TKinematics class defined      //
// in the strangecalc-wrapper project.                                   //
// The user can go beyond the 6 standard 'isospin' channels by setting   //
// the nucleon, kaon and hyperon mass and choosing Q^2 intelligently.    //
//                                                                       // 
// Because TKinematics was not intended for general 2->2 scattering, the //
// getters and setters can have awkward names. This should not be an     //
// insurmountable problem however.                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "TKinematics2to2.h"

ClassImp(TKinematics2to2)

//_____________________________________________________________________
TKinematics2to2::TKinematics2to2(const char *name,const char *title,
				 double mN, double mK, double mY,
				 const char *formatString,
				 double var1,double var2,double var3)
: TKinematics(name,title,1,formatString,var1,var2,var3)
{
  // Constructor
  //
  // Q^2,mN = masses of initial state particles
  // mK,mY  = masses of final state particles
  //
  // Remember all quantities should be in MeV.
  //
  // See TKinematics::TKinematics for more info

  SetMasses(mN,mK,mY);
}

//_____________________________________________________________________
TKinematics2to2::TKinematics2to2(TRootIOCtor *rio)
  : TKinematics(rio)
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
void TKinematics2to2::SetMasses(double mN, double mK, double mY)
{
  // mN = mass of initial state particle
  // mK,mY  = masses of final state particles 
  //
  // Masses should be in MeV.

  SetMn(mN);
  SetMk(mK);
  SetMy(mY);
  UpdateKinematics();
}
