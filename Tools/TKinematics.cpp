/*
 * TKinematics.cpp
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 
 * \copyright
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
// TKinematics                                                                      //
//                                                                                  //
// TKinematics allows to calculate kinematics for g N -> K Y.                       //
// All quantities are in MeV or rad.   
//
// The kinematics are solved in the reference frame where
// * the z-axis lies along the photon's momentum
// * the x-axis defined by the component of the kaon momentum perpendicular to the
//   z-axis
// * y ~ z times x
//                                                                                  //
// TKinematics is designed to work with a fixed number of final states, confusingly //
// labeled isospin channels:                                                        //
//   (1)   ===>   g + p   ->  K+ + L0                                               //
//   (2)   ===>   g + p   ->  K+ + S0                                               //
//   (3)   ===>   g + p   ->  K0 + S+                                               //
//   (4)   ===>   g + n   ->  K0 + L0                                               //
//   (5)   ===>   g + n   ->  K0 + S0                                               //
//   (6)   ===>   g + n   ->  K+ + S-                                               //
//   (7)   ===>   g + p   ->  P0 + p                                                //
//   (8)   ===>   g + p   ->  P+ + n                                                //
//   (9)   ===>   g + n   ->  P0 + n                                                //
//   (10)  ===>   g + n   ->  P- + p                                                //
//   (11)  ===>   K- + p  ->  g  + L0                                               //
//   (12)  ===>   K- + p  ->  g  + S0                                               //
//
// The user needs to specify 2, 3 or 4 variables to unambiguously fix the kinematics.
// Obviously, he/she needs to make sure the variabels of choice are independent.
// More info can be found in SetFormat(const char*).
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include "TKinematics.h"
// #include <numtoa.h>
#include <vector>
#include <cassert>

using std::sqrt;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

/*!
 * \class TKinematics
 *
 * A TKinematics object represents a series of kinematic points for
 * photon or electron induced production of 2 particles off a nucleon.
 *
 * See http://rprmodel.ugent.be/api/TKinematics.html.
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 * 
 */

ClassImp(TKinematics)

//_____________________________________________________________________
TKinematics::TKinematics(TRootIOCtor* rio)
: TNamed::TNamed()
  ,fFormat(0)
  ,fIsospin(1)
  ,fNVars(0)
  ,fMp(0.0)
  ,fMk(0.0)
  ,fMy(0.0)
  ,fWlab(0.0)
  ,fWcm(0.0)
  ,fKlab(0.0)
  ,fKcm(0.0)
  ,fPk(0.0)
  ,fCosthkcm(0.0), fCosthklab(0.),fPYlab(0.),fCosthYlab(0.)
  ,fQsquared(0.0)
  ,fXb(0.0)
  ,fW(0.0)
  ,fS(0.0)
  ,fT(0.0)
  ,fU(0.0)
  ,fPhi(0)
  ,fPhiMin(0)
  ,fIsPhysical(false)
  ,fNrOfPhysicals(-1)
  ,fVar()
  ,fEVar(0)
  ,fAVar(0)
  ,fIsVar()
  ,fLowerLimit()
  ,fUpperLimit()
  ,fStepSize()
  ,fNumberOfSteps()
  ,fStepVar()
{
  // ROOT I/O Constructor

  // This constructor leaves the object in an inconsistent state. 
  // After calling the default constructor you should make sure to at 
  // least call SetFormat() and SetIsospin(int)!
  // Failing to do this will result in a serious ERROR when attempting to
  // TKinematics::Write() the object.
}

//_____________________________________________________________________
TKinematics::TKinematics(const char* name, const char* title,
			 int isospin, const char *format, 
			 double var1, double var2)
  : TNamed::TNamed(name,title)
  ,fFormat(0)
  ,fIsospin(isospin)
  ,fNVars(2)
  ,fMp(0.0)
  ,fMk(0.0)
  ,fMy(0.0)
  ,fWlab(0.0)
  ,fWcm(0.0)
  ,fKlab(0.0)
  ,fKcm(0.0)
  ,fPk(0.0)
  ,fPklab(0.0)
  ,fCosthkcm(0.0), fCosthklab(0.),fPYlab(0.),fCosthYlab(0.)
  ,fQsquared(0.0)
  ,fXb(0.0)
  ,fW(0.0)
  ,fS(0.0)
  ,fT(0.0)
  ,fU(0.0)
  ,fPhi(0)
  ,fPhiMin(0)
  ,fIsPhysical(false)
  ,fNrOfPhysicals(-1)
  ,fVar()
  ,fEVar(0)
  ,fAVar(0)
  ,fIsVar()
  ,fLowerLimit()
  ,fUpperLimit()
  ,fStepSize()
  ,fNumberOfSteps()
  ,fStepVar()
{
  // Constructor
  //
  // The user needs to specify 
  // * the isospin channel (see SetIsospin(int) for more info).
  // * 2 independent variables with the 'format' string. See SetFormat(const char*)
  //   for info about what constitutes a valid format string.
  // * the initial values (var1-2) of the independent variables in the same order as listed
  //  in the format string.
  //
  // The variables Q^2 and phi will be set to zero.

  ChangeIsospinChannel(isospin);

  ReadFormatString(format);
  if( fNVars != 2 ) {
    cerr << "ERROR in TKinematics::TKinematics(const char*,const char*,int,const char*,"
	 << "double,double): "
	 << "Constructor requires 2 independent variables in format string." << endl;
    assert( fNVars==2 );
  }

  *fVar[0] = var1;
  *fVar[1] = var2;
  UpdateKinematics();

  InitializeArrays();
}

//_____________________________________________________________________
TKinematics::TKinematics(const char* name, const char* title,
			 int isospin, const char *format, 
			 double var1, double var2, double var3)
  : TNamed::TNamed(name,title)
  ,fFormat(0)
  ,fIsospin(isospin)
  ,fNVars(3)
  ,fMp(0.0)
  ,fMk(0.0)
  ,fMy(0.0)
  ,fWlab(0.0)
  ,fWcm(0.0)
  ,fKlab(0.0)
  ,fKcm(0.0)
  ,fPk(0.0)
  ,fPklab(0.0)
  ,fCosthkcm(0.0), fCosthklab(0.),fPYlab(0.),fCosthYlab(0.)
  ,fQsquared(0.0)
  ,fXb(0.0)
  ,fW(0.0)
  ,fS(0.0)
  ,fT(0.0)
  ,fU(0.0)
  ,fPhi(0)
  ,fPhiMin(0)
  ,fIsPhysical(false)
  ,fNrOfPhysicals(-1)
  ,fVar()
  ,fEVar(0)
  ,fAVar(0)
  ,fIsVar()
  ,fLowerLimit()
  ,fUpperLimit()
  ,fStepSize()
  ,fNumberOfSteps()
  ,fStepVar()
{
  // Constructor
  //
  // The user needs to specify 
  // * the isospin channel (see SetIsospin(int) for more info).
  // * the 3 independent variables with the 'format' string. See SetFormat(const char*)
  //   for info about what constitutes a valid format string.
  // * the initial values (var1-3) of the independent variables in the same order as listed
  //  in the format string.
  //
  // The unspecified variable Q^2 and/or phi will be set to zero.

  ChangeIsospinChannel(isospin);

  ReadFormatString(format);
  if( fNVars != 3 ) {
    cerr << "ERROR in TKinematics::TKinematics(const char*,const char*,int,const char*,"
	 << "double,double,double): "
	 << "Constructor requires 3 independent variables in format string." << endl;
    assert( fNVars==3 );
  }

  *fVar[0] = var1;
  *fVar[1] = var2;
  *fVar[2] = var3;
  UpdateKinematics();

  InitializeArrays();
}

//_____________________________________________________________________
TKinematics::TKinematics(const char* name, const char* title,
			 int isospin, const char *format, 
			 double var1, double var2, double var3, double var4)
  : TNamed::TNamed(name,title)
  ,fFormat(0)
  ,fIsospin(isospin)
  ,fNVars(4)
  ,fMp(0.0)
  ,fMk(0.0)
  ,fMy(0.0)
  ,fWlab(0.0)
  ,fWcm(0.0)
  ,fKlab(0.0)
  ,fKcm(0.0)
  ,fPk(0.0)
  ,fPklab(0.0)
  ,fCosthkcm(0.0)
  ,fQsquared(0.0)
  ,fXb(0.0)
  ,fW(0.0)
  ,fS(0.0)
  ,fT(0.0)
  ,fU(0.0)
  ,fPhi(0)
  ,fPhiMin(0)
  ,fIsPhysical(false)
  ,fNrOfPhysicals(-1)
  ,fVar()
  ,fEVar(0)
  ,fAVar(0)
  ,fIsVar()
  ,fLowerLimit()
  ,fUpperLimit()
  ,fStepSize()
  ,fNumberOfSteps()
  ,fStepVar()
{
  // Constructor
  //
  // The user needs to specify 
  // * the isospin channel (see SetIsospin(int) for more info).
  // * the 4 independent variables with the 'format' string. See SetFormat(const char*)
  //   for info about what constitutes a valid format string.
  // * the initial values (var1-4) of the independent variables in the same order as listed
  //  in the format string.

  ChangeIsospinChannel(isospin);

  ReadFormatString(format);
  if( fNVars != 4 ) {
    cerr << "ERROR in TKinematics::TKinematics(const char*,const char*,int,const char*,"
	 << "double,double,double,double): "
	 << "Constructor requires 4 independent variables in format string." << endl;
    assert( fNVars==4 );
  }

  *fVar[0] = var1;
  *fVar[1] = var2;
  *fVar[2] = var3;
  *fVar[3] = var4;
  UpdateKinematics();

  InitializeArrays();
}

//_____________________________________________________________________
TKinematics::TKinematics(const TKinematics& toCopy)
  : TNamed(toCopy)
  ,fFormat(0)
  ,fIsospin(toCopy.fIsospin)
  ,fNVars(toCopy.fNVars)
  ,fMp(toCopy.fMp)
  ,fMk(toCopy.fMk)
  ,fMy(toCopy.fMy)
  ,fWlab(toCopy.fWlab)
  ,fWcm(toCopy.fWcm)
  ,fKlab(toCopy.fKlab)
  ,fKcm(toCopy.fKcm)
  ,fPk(toCopy.fPk)
  ,fPklab(toCopy.fPklab)
  ,fCosthkcm(toCopy.fCosthkcm),fCosthklab(toCopy.fCosthklab),fPYlab(toCopy.fPYlab),fCosthYlab(toCopy.fCosthYlab)
  ,fQsquared(toCopy.fQsquared)
  ,fXb(toCopy.fXb)
  ,fW(toCopy.fW)
  ,fS(toCopy.fS)
  ,fT(toCopy.fT)
  ,fU(toCopy.fU)
  ,fPhi(toCopy.fPhi)
  ,fPhiMin(toCopy.fPhiMin)
  ,fIsPhysical(toCopy.fIsPhysical)
  ,fNrOfPhysicals(toCopy.fNrOfPhysicals)
  ,fVar()
  ,fEVar(0)
  ,fAVar(0)
  ,fIsVar(toCopy.fIsVar)
  ,fLowerLimit(toCopy.fLowerLimit)
  ,fUpperLimit(toCopy.fUpperLimit)
  ,fStepSize(toCopy.fStepSize)
  ,fNumberOfSteps(toCopy.fNumberOfSteps)
  ,fStepVar(toCopy.fStepVar)
{
  // Copy constructor
  if(toCopy.fFormat) // when toCopy.fFormat is non-NULL
    ReadFormatString(toCopy.fFormat);
}

//_____________________________________________________________________
void TKinematics::InitializeArrays()
{
  // Resize all vector members to fNVars and 
  // assume the independent variables to be constant.

  fIsVar.resize(fNVars);
  fLowerLimit.resize(fNVars);
  fUpperLimit.resize(fNVars);
  fStepSize.resize(fNVars);
  fNumberOfSteps.resize(fNVars);
  fStepVar.resize(fNVars);

  for(int i=0; i<fNVars; ++i)
    {
      fIsVar[i] = false;
      fLowerLimit[i] = fUpperLimit[i] = *fVar[i];
      fStepSize[i] = 0.0;
      fNumberOfSteps[i] = 1;
      fStepVar[i] = 0;
    }
}

//_____________________________________________________________________
TKinematics* TKinematics::Clone(const char *newname) const
{
  // Virtual copy constructor
  // If newname is specified, this will be the name of the new object.
  TKinematics *clone = new TKinematics(*this);
  if (newname && std::strlen(newname)) clone->SetName(newname);
  return clone;
}

//_____________________________________________________________________
TKinematics::~TKinematics()
{
  // Destructor
  delete[] fFormat;
}

//_____________________________________________________________________
TKinematics& TKinematics::operator=(const TKinematics& toCopy)
{
  // Assignment
  if( this!=&toCopy ) { // avoid self-assignment

    TNamed::operator=(toCopy);
    
    fIsospin = toCopy.fIsospin;
    fNVars = toCopy.fNVars;
    fMp = toCopy.fMp;
    fMk = toCopy.fMk;
    fMy = toCopy.fMy;
    ReadFormatString(toCopy.fFormat);

    fWlab = toCopy.fWlab;
    fWcm = toCopy.fWcm;
    fKlab = toCopy.fKlab;
    fKcm = toCopy.fKcm;
    fPk = toCopy.fPk;
    fPklab = toCopy.fPklab;
    fCosthkcm = toCopy.fCosthkcm;
    fQsquared = toCopy.fQsquared;
    fXb = toCopy.fXb;
    fW = toCopy.fW;
    fS = toCopy.fS;
    fT = toCopy.fT;
    fU = toCopy.fU;
    fPhi = toCopy.fPhi;
    fPhiMin = toCopy.fPhiMin;
    fIsPhysical = toCopy.fIsPhysical;
    fNrOfPhysicals = toCopy.fNrOfPhysicals;

    fIsVar = toCopy.fIsVar;
    fLowerLimit = toCopy.fLowerLimit;
    fUpperLimit = toCopy.fUpperLimit;
    fStepSize = toCopy.fStepSize;
    fNumberOfSteps = toCopy.fNumberOfSteps;
    fStepVar = toCopy.fStepVar;
  
  }
  
  return *this;
}

//_____________________________________________________________________
void TKinematics::SetIsospin(int isospin)
{
  // Change isospin channel. Kinematics will immediately be updated.
  //
  // List of known isospin channels:
  //   (1)   ===>   g + p   ->  K+ + L0
  //   (2)   ===>   g + p   ->  K+ + S0
  //   (3)   ===>   g + p   ->  K0 + S+
  //   (4)   ===>   g + n   ->  K0 + L0
  //   (5)   ===>   g + n   ->  K0 + S0
  //   (6)   ===>   g + n   ->  K+ + S-
  //   (7)   ===>   g + p   ->  P0 + p 
  //   (8)   ===>   g + p   ->  P+ + n 
  //   (9)   ===>   g + n   ->  P0 + n 
  //   (10)  ===>   g + n   ->  P- + p 
  //   (11)  ===>   K- + p  ->  g  + L0
  //   (12)  ===>   K- + p  ->  g  + S0

  if( ((isospin==11 || isospin==12) && fIsospin<11) ||
      (isospin<11 && (fIsospin==11 || fIsospin==12)) ) {
    cerr << "WARNING in TKinematics::SetIsospin(int): "
	 << "Format string for current isospin channel "
	 << fIsospin << " can't work for requested isospin channel "
	 << isospin << ". Ignoring request!\n";
    return;
  }
  
  ChangeIsospinChannel(isospin);
  UpdateKinematics();
}

//_____________________________________________________________________
void TKinematics::ChangeIsospinChannel(int isospin)
{
  // Set the particle masses according to the isospin channel

  /*   (1)   ===>   g + p   ->  K+ + L0
   *   (2)   ===>   g + p   ->  K+ + S0
   *   (3)   ===>   g + p   ->  K0 + S+
   *   (4)   ===>   g + n   ->  K0 + L0
   *   (5)   ===>   g + n   ->  K0 + S0
   *   (6)   ===>   g + n   ->  K+ + S-
   *   (7)   ===>   g + p   ->  P0 + p
   *   (8)   ===>   g + p   ->  P+ + n
   *   (9)   ===>   g + n   ->  P0 + n
   *   (10)  ===>   g + n   ->  P- + p
   *   (11)  ===>   K- + p  ->  g  + L0
   *   (12)  ===>   K- + p  ->  g  + S0
   *
   * This implies changing the masses of the target
   * and outgoing meson and hyperon.
   */

  fIsospin = isospin;

  switch( isospin ) 
    {
    case 1:
      fMp = M_P;
      fMk = M_KP;
      fMy = M_L;
      break;
    case 2:
      fMp = M_P;
      fMk = M_KP;
      fMy = M_S0;
      break;
    case 3:
      fMp = M_P;
      fMk = M_K0;
      fMy = M_SP;
      break;
    case 4:
      fMp = M_N;
      fMk = M_K0;
      fMy = M_L;
      break;
    case 5:
      fMp = M_N;
      fMk = M_K0;
      fMy = M_S0;
      break;
    case 6:
      fMp = M_N;
      fMk = M_KP;
      fMy = M_SM;
      break;
    case 7:
      fMp = M_P;
      fMk = M_PI0;
      fMy = M_P;
      break;
    case 8:
      fMp = M_P;
      fMk = M_PIP;
      fMy = M_N;
      break;
    case 9:
      fMp = M_N;
      fMk = M_PI0;
      fMy = M_N;
      break;
    case 10:
      fMp = M_N;
      fMk = M_PIM;
      fMy = M_P;
      break;
    case 11:
      fMp = M_P;
      fMk = M_KP;
      fMy = M_L;
      break;
    case 12:
      fMp = M_P;
      fMk = M_KP;
      fMy = M_S0;
      break;
    default:
      cerr << "ERROR in TKinematics::ChangeIsospinChannel(int):\n"
	   << "Isospin channel " << isospin << " is not implemented!\n";
      break;
    }  
  
  // Number of physical points may change
  fNrOfPhysicals = -1;
}

//_____________________________________________________________________
void TKinematics::SetFormat(const char* format)
{
  /* Begin_Html
     The user needs to specify 2, 3 or 4 variables to unambiguously fix the kinematics.
     Obviously, he/she needs to make sure the variables of choice are independent.
     TKinematics needs:
     <ul>
     <li> 1 energy variable
     <li> 1 angular variable
     <li> photon's virtuality Q<sup>2</sup> <i>(optional)</i>
     <li> Angle &Phi; between lepton and reaction plane <i>(optional)</i>
     </ul>
     <p>
     The 2, 3 or 4 independent variables are set by providing a ':' separated list of variables.
     In case Q<sup>2</sup> and/or &Phi; are not listed as independent variables, they are set to
     zero.
     <p>
     Below is a list of the 'code names' of variables:
     <table>
     <tr><td>Q<sup>2</sup></td><td><tt>qsquared</tt></td></tr>
     <tr><td>&Phi;</td><td><tt>phi</tt></td></tr>
     <tr><td>photon energy (LAB)</sub></td><td><tt>wlab</tt></td></tr>
     <tr><td>photon momentum (LAB)</td><td><tt>klab</tt></td></tr>
     <tr><td>photon energy (CM)</td><td><tt>wcm</tt></td></tr>
     <tr><td>photon momentum (CM)</td><td><tt>kcm</tt></td></tr>
     <tr><td>kaon momentum (LAB)</td><td><tt>pklab</tt></td></tr>
     <tr><td>kaon momentum (CM)</td><td><tt>pk</tt> or <tt>pkcm</tt></td></tr>
     <tr><td>kaon scattering angle (CM)</td><td><tt>costhkcm</tt></td></tr>
     <tr><td>invariant mass</td><td><tt>w</tt></td></tr>
     <tr><td>Mandelstam s</td><td><tt>s</tt></td></tr>
     <tr><td>Mandelstam t</td><td><tt>t</tt></td></tr>
     <tr><td>Mandelstam u</td><td><tt>u</tt></td></tr>
     <tr><td>Bjorken x</td><td><tt>xb</tt></td></tr>
     </table>
     <p>
     Attention! Calling this function will fix the kinematics in its current
     state. All variables will be fixed.
     End_Html */

  if(!IsFixed()) {
    cout << "INFO in TKinematics::SetFormat(const char*): "
	 << "Fixing all variables.\n";
    for(int i=0; i<fNVars; ++i)
      FixVariable(i+1);
  }

  ReadFormatString(format);
  
  // Finally, we reset all the arrays
  InitializeArrays();
}

//_____________________________________________________________________
void TKinematics::ReadFormatString(const TString format)
{
  // Process the format string
  //
  // The number of independent variables can change. It is advised
  // to re-initialize the member arrays.
  //
  // In case Q^2 and/or phi are not listed as independent variables, 
  // they are set to zero.

  // Store format in TKinematics::fFormat
  // Read TKinematics::fFormat and let fVar[] point
  // to the 2, 3 or 4 independent kinematical points as specified
  // by the user.
  // Check that the choosen variables make sense.

  // Tokenize the format string
  TObjArray *listOfTokens = format.Tokenize(":");
  fNVars = listOfTokens->GetEntries();

  // and check the number of variables
  if( (fNVars<2) && (fNVars>4) ) {
    cerr << "ERROR in TKinematics::ReadFormatString(const TString):\n"
	 << "Format string should provide 2, 3 or 4 independent variables: "
	 << "\"var1:var2[:var3[:var4]]\""  << "\"\n";
    assert( (fNVars>=2) && (fNVars<=4) );
  }

  // allocate memory for fFormat and fVar
  if(!fFormat) fFormat = new char[100];
  fVar.resize( fNVars );

  // Put format string in member and in temporary *char.
  std::strcpy(fFormat,format);
  
  // Loop over all tokens
  // and determine independent kinematic variables
  int qsqCount=0, energyCount=0, angleCount=0, phiCount=0; // count type of variables

  for(int itoken=0; itoken<fNVars; ++itoken)
    {
      TString token = dynamic_cast<TObjString*>(listOfTokens->At(itoken))->GetString();

      if( !token.CompareTo("wlab",TString::kIgnoreCase) ) {
	fVar[itoken] = fEVar = &fWlab;
	++energyCount;
      }
      else if( !token.CompareTo("costhkcm",TString::kIgnoreCase) ) {
	fVar[itoken] = fAVar = &fCosthkcm;
	++angleCount;
      }
      else if( !token.CompareTo("qsquared",TString::kIgnoreCase) ) {
	fVar[itoken] = &fQsquared;
	++qsqCount;
      }
      else if( !token.CompareTo("w",TString::kIgnoreCase) ) {
	fVar[itoken] = fEVar = &fW;
	++energyCount;
      }
      else if( !token.CompareTo("s",TString::kIgnoreCase) ) {
	fVar[itoken] = fEVar = &fS;
	++energyCount;
      }
      else if( !token.CompareTo("t",TString::kIgnoreCase) ) {
	fVar[itoken] = fAVar = &fT;
	++angleCount;
      }
      else if( !token.CompareTo("wcm",TString::kIgnoreCase) ) {
	fVar[itoken] = fEVar = &fWcm;
	++energyCount;
      }
      else if( !token.CompareTo("klab",TString::kIgnoreCase) ) {
	fVar[itoken] = fEVar = &fKlab;
	++energyCount;
      }
      else if( !token.CompareTo("kcm",TString::kIgnoreCase) ) {
	fVar[itoken] = fEVar = &fKcm;
	++energyCount;
      }
      else if( !token.CompareTo("pk",TString::kIgnoreCase) ) {
	fVar[itoken] = fEVar = &fPk;
	++energyCount;
      }
      else if( !token.CompareTo("pkcm",TString::kIgnoreCase) ) {
	fVar[itoken] = fEVar = &fPk;
	++energyCount;
      }
      else if( !token.CompareTo("pklab",TString::kIgnoreCase) ) {
	//hacked by me, coz I use pklab as an angle variable!!!
	fVar[itoken] = fAVar = &fPklab;
	++angleCount;
      }
      else if( !token.CompareTo("u",TString::kIgnoreCase) ) {
	fVar[itoken] = fAVar = &fU;
	++angleCount;
      }
      else if( !token.CompareTo("xb",TString::kIgnoreCase) ) {
	fVar[itoken] = fEVar = &fXb;
	++energyCount;
      }
      else if( !token.CompareTo("phi",TString::kIgnoreCase) ) {
	fVar[itoken] = &fPhi;
	++phiCount;
      }
      else {
	fVar[itoken] = NULL;
	cerr << "ERROR in  TKinematics::ReadFormatString(const TString):\n"
	     << "Kinematic variable \"" << token << "\" unknown!\n";
	exit(1);
      }
    }  // end loop over tokens

  // Get rid of the list of tokens
  delete listOfTokens;

  // Check whether the user provided a meaningful format string
  // First the optional arguments: Q^2 and phi (set them to zero when necessary)
  if( fNVars == 3 ) {
    if( qsqCount != 1 ) {
      if( phiCount != 1 ) {
	cerr << "ERROR in TKinematics::ReadFormatString(const TString):\n"
	     << "The angle phi or Q^2 need to be one of the 3 independent kinematic variables!\n";
	exit(1);
      }
      else {
	fQsquared = 0.;
      }
    }
    else {
      if( phiCount == 1 ) {
	cerr << "ERROR in TKinematics::ReadFormatString(const TString):\n"
	     << "The angle phi and Q^2 cannot simultaneously be one of the 3 independent kinematic variables!\n";
	exit(1);
      }
      else {
	fPhi = 0.;
      }
    }
  }
  else if( fNVars == 4 && qsqCount != 1 && phiCount != 1 ) {
    cerr << "ERROR in TKinematics::ReadFormatString(const TString):\n"
	 << "The angle phi and Q^2 need to be one of the 4 independent kinematic variables!\n";
    exit(1);
  }
  else {
    fQsquared = fPhi = 0.;
  }

  // Check for 1 and only 1 energy variable
  if( energyCount == 0 ) {
    cerr << "ERROR in TKinematics::ReadFormatString(const TString):\n"
	 << "wlab/wcm/klab/kcm/pk/pklab/w/s/xb needs to be one of " << fNVars
	 << " independent kinematic variables.\n";
    exit(1); }
  else if( energyCount != 1 ) {
    cerr << "ERROR in TKinematics::ReadFormatString(const TString):\n"
	 << "Only one of following kinematic variables can be independent: "
	 << "wlab/wcm/klab/kcm/pk/pklab/w/s/xb.\n";
    exit(1);
  }

  // Check for 1 and only 1 angular variable
  if( angleCount == 0 ) {
    cerr << "ERROR in TKinematics::ReadFormatString(const TString):\n"
	 << "costhkcm/t/u needs to be one of " << fNVars
	 << " independent kinematic variables.\n";
    exit(1); }
  else if( angleCount != 1 ) {
    cerr << "ERROR in TKinematics::ReadFormatString(const TString):\n"
	 << "Only one of following kinematic variables can be independent: "
	 << "costhkcm/t/u.\n";
    exit(1);
  }

  // Some variables are not possible for certain isospin channels
  if( fIsospin == 11 || fIsospin == 12 ) {
    if( fEVar==&fKlab || fEVar==&fWlab || fEVar==&fXb ) {
      cerr << "ERROR in  TKinematics::ReadFormatString(const TString):\n"
	   << "Kinematic variable wlab/klab/xb is not allowed "
	   << "as independent variable in radiative capture kinematics!\n";
      exit(1);
    }
  }
  else {
    if( fEVar==&fPk || fEVar==&fPklab ) {
      cerr << "ERROR in  TKinematics::ReadFormatString(const TString):\n"
	   << "Kinematic variable pk/pkcm/pklab is not allowed "
	   << "as independent variable in EM meson production!\n";
      exit(1);
    }
  }
}

//_____________________________________________________________________
TString TKinematics::GetVarName(int varNumber) const
{
  // Returns the name of the 'varNumber'th independent variable 
  // as specified in the format string.

  // First we do some checks
  if( varNumber<1 || varNumber>fNVars ) {
    cerr << "ERROR in TKinematics::GetVarName(int):\n"
	 << "Index out-of-bounds!\n";
    return TString();
  }

  if( fVar.size()==0 ) {
    cerr << "ERROR in TKinematics::GetVarName(int):\n"
	 << "Independent kinematic variables have not been initialized!\n";
    return TString();
  }

  TString varName;
  
  // Lookup name
  if( fVar[varNumber-1] == &fWlab )
    varName = "wlab";
  else if( fVar[varNumber-1] == &fWcm )
    varName = "wcm";
  else if( fVar[varNumber-1] == &fKlab )
    varName = "klab";
  else if( fVar[varNumber-1] == &fKcm )
    varName = "kcm";
  else if( fVar[varNumber-1] == &fPk )
    varName = "pk";
  else if( fVar[varNumber-1] == &fPklab )
    varName = "pklab";
  else if( fVar[varNumber-1] == &fCosthkcm )
    varName = "costhkcm";
  else if( fVar[varNumber-1] == &fQsquared )
    varName = "qsquared";
  else if( fVar[varNumber-1] == &fW )
    varName = "w";
  else if( fVar[varNumber-1] == &fS )
    varName = "s";
  else if( fVar[varNumber-1] == &fT )
    varName = "t";
  else if( fVar[varNumber-1] == &fU )
    varName = "u";
  else if( fVar[varNumber-1] == &fXb )
    varName = "xb";
  else if( fVar[varNumber-1] == &fPhi )
    varName = "phi";
  else {
    cerr << "Error in TKinematics::GetVarName(int):\n"
	 << "TKinematics::fVar[" << varNumber-1 << "] "
	 << "points to unknown variable.\n";
  }

  return varName;
}

//_____________________________________________________________________
void TKinematics::SetVar(int varNumber, double value)
{
  // Set the value of the 'varNumber'th independent variable given
  // in the format string.
  // This function can't be called when a range has been given to
  // this variable.

  if( varNumber<1 || varNumber>fNVars )
    cerr << "ERROR in TKinematics::SetVar(int,double): "
	 << "Index out-of-bounds!\n";

  else if( fVar.size()==0 )
    cerr << "ERROR in TKinematics::SetVar(int,double): "
	 << "Independent kinematic variables have not been initialized!\n";
    
  else if( fIsVar[varNumber-1] )
    cerr << "ERROR in TKinematics::SetVar(int,double): "
	 << "This method can not be used because the variable is not fixed.\n";

  else {
    fLowerLimit[varNumber-1] = value;
    fUpperLimit[varNumber-1] = value;
    GoTo(0);
    fNrOfPhysicals = -1; // Number of physical points may change
  }
}

//_____________________________________________________________________
void TKinematics::SetPhiLimits(double phiMin, double phiMax)
{
  // Set fPhi = phiMax and fPhiMin = fMin, and fix the variable fPhi

  if(phiMin == phiMax)
    cerr << "ERROR in TKinematics::SetPhiLimits(double,double): "
	 << "Upper and lower limit for phi are the same.\n";

  else if(phiMin > phiMax)
    cerr << "ERROR in TKinematics::SetPhiLimits(double,double): "
	 << "phiMax has to be larger than phiMin.\n";

  else
    {
      int i = 1;
      while(i <= fNVars && GetVarName(i) != "phi")
	i++;

      if(i != fNVars + 1)
	{
	  SetVar(i,phiMax);
	  FixVariable(i);
	}
      else
	fPhi = phiMax;
      
      fPhiMin = phiMin;
    }
}

//_____________________________________________________________________
double TKinematics::GetVar(int varNumber) const
{
  // Returns the current value of the 'varNumber'th independent
  // variable as specified in the format string.

  // First we do some checks
  if( varNumber<1 || varNumber>fNVars ) {
    cerr << "ERROR in TKinematics::GetVar(int):\n"
	 << "Index out-of-bounds!\n";
    return -1e16;
  }

  if( fVar.size()==0 ) {
    cerr << "ERROR in TKinematics::GetVar(int):\n"
	 << "Independent kinematic variables have not been initialized!\n";
    return -1e16;
  }

  return *fVar[varNumber-1];
  
}

//_____________________________________________________________________
double TKinematics::GetTMin() const
{
  // Returns the minimum value of -t.

  return fMk*fMk - fQsquared -2.0*fWcm*sqrt(fPk*fPk+fMk*fMk) + 2.0*fPk*fKcm;
}

//_____________________________________________________________________
template <typename T>
void TKinematics::Streamer(TBuffer &R__b, vector<T>& vector, Version_t R__v)
{
  // Stream an object of type vector<>.

  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
    // No need to check different versions.
    // They are always written to file in the same way.
    vector.clear();
    vector.resize(fNVars);
    T *tmp = new T[fNVars];
    R__b.ReadFastArray(tmp,fNVars);
    vector.assign( tmp , tmp+fNVars );
    delete[] tmp;
  }
  else { // writing
    // We exploit the fact that vectors are continuous in memory
    R__b.WriteFastArray(&vector[0],fNVars);
  }
}

//_____________________________________________________________________
template <>
void TKinematics::Streamer<bool>(TBuffer &R__b, vector<bool>& vector, Version_t R__v)
{
  // Stream an object of type vector<bool>.

  // A template specialization is needed for bools, because the implementation
  // of vector<bool> is not what you would expect.
  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
    vector.clear();
    vector.resize(fNVars);
    bool *tmp = new bool[fNVars];
    if (R__v < 5) {
      R__b.ReadFastArray(tmp,fNVars);
    } else {
      for(int i=0; i<fNVars; ++i) R__b >> tmp[i];
    }
    vector.assign( tmp , tmp+fNVars );
    delete[] tmp;
  }
  else { // writing
    for(int i=0; i<fNVars; ++i) R__b << vector[i];
  }
}

//______________________________________________________________________________
void TKinematics::Streamer(TBuffer &R__b)
{
  // Stream an object of class TKinematics.

  // This function was generated by rootcint with #pragma link c++ class TKinematics+
  // We implement it ourselfs to be able to initialize the persistent member **fVar.
  // This comes down to calling ReadFormatString(fFormat) after reading the buffer.
  // See ROOT user guide chpt.11 -> Streamers for more info

  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
    TNamed::Streamer(R__b);
    delete [] fFormat;
    fFormat = new char[100];
    char *tempFormat = new char[100]; // Added by PVC
    R__b.ReadFastArray(tempFormat,100); // Added by PVC
    R__b >> fIsospin;
    if (R__v > 4 ) R__b >> fNVars; // added in version 5
    else fNVars = 3;
    R__b >> fMp;
    R__b >> fMk;
    R__b >> fMy;
    R__b >> fWlab;
    R__b >> fWcm;
    R__b >> fKlab;
    R__b >> fKcm;
    R__b >> fPk;
    if (R__v > 1) R__b >> fPklab; // added version 2
    R__b >> fCosthkcm;
    R__b >> fQsquared;
    if (R__v > 3) R__b >> fXb; // added in version 4
    else fXb = fQsquared / 2. / fMp / fWlab;
    R__b >> fW;
    R__b >> fS;
    R__b >> fT;
    R__b >> fU;
    if( R__v > 4) { // added in version 5
      R__b >> fPhi;
      R__b >> fPhiMin;
    } else 
      fPhi = fPhiMin = 0.;
    R__b >> fIsPhysical;
    if(R__v<3) fNrOfPhysicals = -1; // added version 3
    else R__b >> fNrOfPhysicals;
    Streamer(R__b, fIsVar, R__v);
    Streamer(R__b, fLowerLimit, R__v);
    Streamer(R__b, fUpperLimit, R__v);
    Streamer(R__b, fStepSize, R__v);
    Streamer(R__b, fNumberOfSteps, R__v);
    if( R__v > 4 ) Streamer(R__b, fStepVar, R__v);
    else {
      // Prior to version 5, we store fStep (a single integer).
      // This variable represented to total step
      int step;
      R__b >> step;
      
      // Prepare the vector
      fStepVar.clear();
      fStepVar.resize(3);
      
      // Calculate the individual steps
      fStepVar[2] = step / (fNumberOfSteps[1]*fNumberOfSteps[0]);
      fStepVar[1] = (step % (fNumberOfSteps[1]*fNumberOfSteps[0]))/ fNumberOfSteps[0];
      fStepVar[0] = ((step % (fNumberOfSteps[1]*fNumberOfSteps[0])) % fNumberOfSteps[0]);
    }
    R__b.CheckByteCount(R__s, R__c, TKinematics::IsA());
    ReadFormatString(tempFormat); // Added by PVC
    delete[] tempFormat; // Added by PVC
    if (R__v < 2) UpdateKinematics();  // added version 2
  } else {
    R__c = R__b.WriteVersion(TKinematics::IsA(), kTRUE);
    TNamed::Streamer(R__b);
    R__b.WriteFastArray(fFormat,100);
    R__b << fIsospin;
    R__b << fNVars;
    R__b << fMp;
    R__b << fMk;
    R__b << fMy;
    R__b << fWlab;
    R__b << fWcm;
    R__b << fKlab;
    R__b << fKcm;
    R__b << fPk;
    R__b << fPklab;
    R__b << fCosthkcm;
    R__b << fQsquared;
    R__b << fXb;
    R__b << fW;
    R__b << fS;
    R__b << fT;
    R__b << fU;
    R__b << fPhi;
    R__b << fPhiMin;
    R__b << fIsPhysical;
    R__b << fNrOfPhysicals;
    Streamer(R__b, fIsVar, R__c);
    Streamer(R__b, fLowerLimit, R__c);
    Streamer(R__b, fUpperLimit, R__c);
    Streamer(R__b, fStepSize, R__c);
    Streamer(R__b, fNumberOfSteps, R__c);
    Streamer(R__b, fStepVar, R__c);
    R__b.SetByteCount(R__c, kTRUE);
  }
}

//_____________________________________________________________________
void TKinematics::SetVarRange(int varNumber, double lowerLimit,
			      double upperLimit, int numberOfSteps)
{
  // Give independent variable number 'varNumber' a range in the 
  // interval [lowerLimit,upperLimit].

  // First we do some checks
  if( varNumber<1 || varNumber>fNVars )
    cerr << "ERROR in TKinematics::SetVarRange(int,double,double,int):\n"
	 << "Index out-of-bounds!\n";

  else if( fVar.size()==0 ) 
    cerr << "ERROR in TKinematics::SetVarRange(int,double,double,int)\n"
	 << "Independent kinematic variables have not been initialized!\n";

  else if( numberOfSteps<2 )
    cerr << "ERROR in TKinematics::SetVarRange(int,double,double,int)\n"
	 << "Number of steps need to be >=2.\n";

  else {
    fIsVar[varNumber-1] = true;

    // Set limits
    fLowerLimit[varNumber-1] = lowerLimit;
    fUpperLimit[varNumber-1] = upperLimit;

    // Set number of steps
    fNumberOfSteps[varNumber-1] = numberOfSteps;

    // Calculate stepsize
    fStepSize[varNumber-1] = (upperLimit-lowerLimit) / ((numberOfSteps-1)*1.0);

    // Go to first step
    GoTo(0);
    
    // Number of physical points may change
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
void TKinematics::FixVariable(int varNumber)
{
  // Unset the range of variable number 'varNumber'. Set it to its current value.

  // First we do some checks
  if( varNumber<1 || varNumber>fNVars )
    cerr << "ERROR in TKinematics::FixVariable(int):\n"
	 << "Index out-of-bounds!\n";

  else if( fVar.size()==0 ) 
    cerr << "ERROR in TKinematics::FixVariable(int)\n"
	 << "Independent kinematic variables have not been initialized!\n";

  else if( fIsVar[varNumber-1] ) {
    fIsVar[varNumber-1] = false;
    
    // Reset all variables
    fLowerLimit[varNumber-1] = *fVar[varNumber-1];
    fUpperLimit[varNumber-1] = *fVar[varNumber-1];
    fStepSize[varNumber-1] = 0.0;
    fNumberOfSteps[varNumber-1] = 1;
    
    // Go to first step
    fStepVar[varNumber-1] = 0;
    
    // Number of physical points may change
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
void TKinematics::FixVariables()
{
  // Calls FixVariable(int) on all variables
  for(int var=1; var<=fNVars; ++var)
    FixVariable(var);
}

//_____________________________________________________________________
void TKinematics::GoTo(int totalStep)
{
  //  Proceed to the 'totalStep'th point in the grid of independent variables.

  // Total number of steps
  int totNrOfSteps = GetNumberOfSteps();

  // First we do some checks
  if( totalStep<0 || totalStep>=totNrOfSteps )
    cerr << "ERROR in TKinematics::GoTo(int):\n"
	 << "Invalid step. Should be in range [0,"
	 << totNrOfSteps-1 << "].\n";

  else if( fVar.size()==0 ) 
    cerr << "ERROR in TKinematics::GoTo(int)\n"
	 << "Independent kinematic variables have not been initialized!\n";

  else {  
    int prodN;
    for(int i = fNVars - 1; i >= 0; i--)
      {
	prodN = GetNumberOfSteps()/fNumberOfSteps[fNVars - 1];
	fStepVar[i] = totalStep;
	for(int j = fNVars - 1; j > i; j--)
	  {
	    fStepVar[i] -= fStepVar[j]*prodN;
	    prodN = prodN/fNumberOfSteps[j - 1];
	  }
	fStepVar[i] = fStepVar[i]/prodN;
      }

    // Set variables
    for(int i=0; i<fNVars; ++i)
      *fVar[i] = fLowerLimit[i] + fStepVar[i] * fStepSize[i];

    UpdateKinematics();
  }
}

//_____________________________________________________________________
void TKinematics::GoTo(int* steps)
{
  // Proceed to point in the grid of independent variables specified by the
  // arguments. The i'th argument fixes the position in the grid of the 'i'th variable.
  //
  // CAUTION: It is assumed that the size of the steps array equals fNVars.


  bool updateSteps = true;
  int i = 0;
  
  if( fVar.size()==0 )
    {
      cerr << "ERROR in TKinematics::GoTo(int[]):\n"
	   << "Independent kinematic variables have not been initialized!\n";
      updateSteps = false;
    }
 
  while(updateSteps && i < fNVars)
    {
      if(steps[i] < 0 || steps[i] >= fNumberOfSteps[i])
	{
	  cerr << "ERROR in TKinematics::GoTo(int[]):\n"
	       << "Invalid \"" << GetVarName(i + 1) << "\" step. "
	       << "Should be in range [0," << fNumberOfSteps[i] - 1 << "].\n";
	  updateSteps = false;
	}
      i++;
    }

  if(updateSteps)
    {
      for(int i = 0; i < fNVars; i++) {
	fStepVar[i] = steps[i];
	*fVar[i] = fLowerLimit[i] + fStepVar[i] * fStepSize[i];
      }
      UpdateKinematics();
    }
}

//_____________________________________________________________________
void TKinematics::GoTo(int stepVar1, int stepVar2)
{
  // Proceed to point in the grid of independent variables specified by the
  // arguments. The i'th argument fixes the position in the grid of the 'i'th variable.
  if( fNVars != 2 ) {
    cerr << "ERROR in TKinematics::GoTo(int,int): "
	 << "Method only valid in case of 2 independent variables." << endl;
    assert( fNVars == 2 );
  }

  int list[] = {stepVar1,stepVar2};
  GoTo( list );
}

//_____________________________________________________________________
void TKinematics::GoTo(int stepVar1, int stepVar2, int stepVar3)
{
  // Proceed to point in the grid of independent variables specified by the
  // arguments. The i'th argument fixes the position in the grid of the 'i'th variable.
  if( fNVars != 3 ) {
    cerr << "ERROR in TKinematics::GoTo(int,int,int): "
	 << "Method only valid in case of 3 independent variables." << endl;
    assert( fNVars == 3 );
  }

  int list[] = {stepVar1,stepVar2,stepVar3};
  GoTo( list );
}

//_____________________________________________________________________
void TKinematics::GoTo(int stepVar1, int stepVar2, int stepVar3, int stepVar4)
{
  // Proceed to point in the grid of independent variables specified by the
  // arguments. The i'th argument fixes the position in the grid of the 'i'th variable.
  if( fNVars != 4 ) {
    cerr << "ERROR in TKinematics::GoTo(int,int,int,int): "
	 << "Method only valid in case of 4 independent variables." << endl;
    assert( fNVars == 4 );
  }

  int list[] = {stepVar1,stepVar2,stepVar3,stepVar4};
  GoTo( list );
}

//_____________________________________________________________________
bool TKinematics::IsFixed() const
{
  // Check whether any of the independent variables has been given a range.
  for(int i = 0; i < fNVars; i++)
    if(fIsVar[i])
      return false;

  return true;
}

//_____________________________________________________________________
bool TKinematics::IsFixed(int varNumber) const
{
  //  Check whether the 'varNumber'th variable has been given a range.

  // First we do some checks
  if( varNumber<1 || varNumber>fNVars ) {
    cerr << "ERROR in TKinematics::IsFixed(int):\n"
	 << "Index out-of-bounds!\n";
    return true;
  }

  return !fIsVar[varNumber-1];
}

//_____________________________________________________________________
int TKinematics::Next()
{
  // Proceed to the next step in the grid of independent variables.
  // Return value is zero when the end of the grid is reached.

  int step = GetStep();

  // Check whether we are at the end of the grid.
  if( step == (GetNumberOfSteps()-1) ) return 0;

  GoTo(++step);
  return 1;
}

//_____________________________________________________________________
int TKinematics::GetNumberOfSteps() const
{
  // Number of points in the grid of independent variables
  int NSteps = 1;

  for(int i = 0; i < fNVars; i++)
    NSteps *= fNumberOfSteps[i];

  return NSteps;
}

//_____________________________________________________________________
int TKinematics::GetNumberOfPhysicalSteps() const
{
  // Number of physical points in the grid of independent variables
  if( fNrOfPhysicals==-1 ) {
  
    TKinematics tk(*this); // temporary copy
    fNrOfPhysicals = 0;    // final result
    
    for(int step=0; step<tk.GetNumberOfSteps(); ++step) {
      tk.GoTo(step);
      if(tk.IsPhysical()) ++fNrOfPhysicals;
    }
  }

  return fNrOfPhysicals;
}    

//_____________________________________________________________________
int TKinematics::GetNumberOfSteps(int varNumber) const
{
  // Number of points in the grid over the 'varNumber'th independent
  // variable as specified in the format string.

  // First we do some checks
  if( varNumber<1 || varNumber>fNVars ) {
    cerr << "ERROR in TKinematics::GetNumberOfSteps(int):\n"
	 << "Index out-of-bounds!\n";
    return 0;
  }

  return fNumberOfSteps[varNumber-1];
}

//_____________________________________________________________________
double TKinematics::GetStepSize(int varNumber) const
{
  // Returns the size of the steps in between points in the grid of
  // the 'varNumber'th independent variable as specified in the format string.

  // First we do some checks
  if( varNumber<1 || varNumber>fNVars ) {
    cerr << "ERROR in TKinematics::GetStepSize(int):\n"
	 << "Index out-of-bounds!\n";
    return 0;
  }

  return fStepSize[varNumber-1];
}

//_____________________________________________________________________
int TKinematics::GetStep(int varNumber) const
{
  // Return the current position in the grid of the 'varNumber'th 
  // independent variable as specified in the format string.

  if( varNumber<1 || varNumber>fNVars)
    cerr << "ERROR in TKinematics::GetStep(int):\n"
	 << "Index out-of-bounds!\n";
  else 
    return fStepVar[varNumber-1];
  
  return -1;
}

//_____________________________________________________________________
int TKinematics::GetStep() const
{
  // Return the current position in the grid of the 'varNumber'th 
  // independent variable as specified in the format string.

  int 
    step = 0,
    prodN = 1;

  for(int i = 0; i < fNVars; i++)
    {
      step += fStepVar[i]*prodN;
      prodN *= fNumberOfSteps[i];
    }

  return step;
}


//_____________________________________________________________________
double TKinematics::EnergyConservation(double pk)
{
  // Returns the deviation from conservation of energy in the CM for the
  // given kaon momentum with the photon's energy and momentum fixed.

  // See (74) on p.18 of Stijn's notes.
  return sqrt(fMk*fMk+pk*pk) + sqrt(fMy*fMy+pk*pk)
    - sqrt(fMp*fMp+fKcm*fKcm) - fWcm;
}

//_____________________________________________________________________
void TKinematics::UpdateKinematics()
{
  // Solve the kinematics!

  // The kinematics for radiative capture (iso 11 & 12) and
  // EM meson production ( iso 1 - 10 ) need to be solved separately.

  if( fIsospin==11 || fIsospin==12) { // radiative capture

    // photon should be real!
    if( fQsquared!=0.0 ) {
      cerr << "ERROR in TKinematics::UpdateKinematics(): "
	   << "In radiative capture, the emitted photon is real!\n";
      exit(1);
    }

    // Determine all energy variables 
    if( fEVar == &fW ) {
      fS = fW*fW;
      fPklab = sqrt(pow((fS - fMk*fMk - fMp*fMp)/(2.*fMp),2) - fMk*fMk);
      fPk = fMp*fPklab/fW;
      fWcm = fKcm = fW/2. - fMy*fMy/2./fW;
      fWlab = fKlab = fWcm/fW*(fMp+fPklab+sqrt(fPklab*fPklab+fMk*fMk));
    }
    if( fEVar == &fS ) {
      fW = sqrt(fS);
      fPklab = sqrt(pow((fS - fMk*fMk - fMp*fMp)/(2.*fMp),2) - fMk*fMk);
      fPk = fMp*fPklab/fW;
      fWcm = fKcm = fW/2. - fMy*fMy/2./fW;
      fWlab = fKlab = fWcm/fW*(fMp+fPklab+sqrt(fPklab*fPklab+fMk*fMk));
    }
    else if( fEVar == &fPk ) {
      fW = sqrt(fPk*fPk+fMk*fMk) + sqrt(fPk*fPk+fMp*fMp);
      fS = fW*fW;
      fPklab = fW*fPk/fMp;
      fWcm = fKcm = fW/2. - fMy*fMy/2./fW;
      fWlab = fKlab = fWcm/fW*(fMp+fPklab+sqrt(fPklab*fPklab+fMk*fMk));      
    }
    else if( fEVar == &fPklab ) {
      fS = fMk*fMk + fMp*fMp + 2.*fMp*sqrt(fPklab*fPklab+fMk*fMk);
      fW = sqrt(fS);
      fPk = fMp*fPklab/fW;
      fWcm = fKcm = fW/2. - fMy*fMy/2./fW;
      fWlab = fKlab = fWcm/fW*(fMp+fPklab+sqrt(fPklab*fPklab+fMk*fMk));
    }
    else if( fEVar == &fWcm ) {
      fKcm = fWcm;
      fW = fWcm + sqrt(fWcm*fWcm+fMy*fMy);
      fS = fW*fW;
      fPklab = sqrt(pow((fS - fMk*fMk - fMp*fMp)/(2.*fMp),2) - fMk*fMk);
      fPk = fMp*fPklab/fW;
      fWlab = fKlab = fWcm/fW*(fMp+fPklab+sqrt(fPklab*fPklab+fMk*fMk));
    }
    else if( fEVar == &fKcm ) {
      fWcm = fKcm;
      fW = fWcm + sqrt(fWcm*fWcm+fMy*fMy);
      fS = fW*fW;
      fPklab = sqrt(pow((fS - fMk*fMk - fMp*fMp)/(2.*fMp),2) - fMk*fMk);
      fPk = fMp*fPklab/fW;
      fWlab = fKlab = fWcm/fW*(fMp+fPklab+sqrt(fPklab*fPklab+fMk*fMk));
    }
    fXb = 0.;

    if( fPk<0. || fWcm<0. || fWlab<0. )
      fIsPhysical = false;
    else
      fIsPhysical = true;

    // Determine the angular variables
    if( fAVar == &fCosthkcm ) {
      fT = fMk*fMk -2.0*fWcm*sqrt(fPk*fPk+fMk*fMk) + 2.0*fPk*fKcm*fCosthkcm;
      fU = fMy*fMy +fMk*fMk -2.0*sqrt(fKcm*fKcm+fMy*fMy)*sqrt(fPk*fPk+fMk*fMk) 
	- 2.0*fPk*fKcm*fCosthkcm;
    }
    else if( fAVar == &fT ) {
      fCosthkcm = ( fT + 2.0*fWcm*sqrt(fPk*fPk+fMk*fMk) - fMk*fMk )
	/ 2.0 / fPk / fKcm;
      fU = fMy*fMy +fMk*fMk -2.0*sqrt(fKcm*fKcm+fMy*fMy)*sqrt(fPk*fPk+fMk*fMk) 
	- 2.0*fPk*fKcm*fCosthkcm;
    }
    else {
      fCosthkcm = ( fU + 2.0*sqrt(fKcm*fKcm+fMy*fMy)*sqrt(fPk*fPk+fMk*fMk)
		    - fMy*fMy - fMk*fMk ) / 2. / fPk / fKcm;
      fT = fMk*fMk -2.0*fWcm*sqrt(fPk*fPk+fMk*fMk) + 2.0*fPk*fKcm*fCosthkcm;
    }

    // Check whether kinematics lie on the physical plane
    if( std::fabs(fCosthkcm) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( std::fabs(fCosthkcm)-1.0 < Underflow() ) {
      fIsPhysical = true && fIsPhysical;
      fCosthkcm = ( fCosthkcm>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;

  } // end radiative capture

  else { // meson production

    // Determine all energy variables (s and wlab first)
    if( fEVar == &fWlab ) {
      fS = fMp*fMp + 2.0 * fMp * fWlab - fQsquared;
      fW = sqrt(fS);
      fKlab = sqrt(fWlab*fWlab + fQsquared);
      fWcm = (fWlab*fMp-fQsquared)/sqrt(fMp*fMp+2.0*fWlab*fMp-fQsquared);
      fKcm = sqrt( fWcm*fWcm + fQsquared );
      fXb = fQsquared / 2. / fMp / fWlab;
    }
    else if( fEVar == &fXb ) {
      fWlab = fQsquared / 2. / fMp / fXb;
      fS = fMp*fMp + 2.0 * fMp * fWlab - fQsquared;
      fW = sqrt(fS);
      fKlab = sqrt(fWlab*fWlab + fQsquared);
      fWcm = (fWlab*fMp-fQsquared)/sqrt(fMp*fMp+2.0*fWlab*fMp-fQsquared);
      fKcm = sqrt( fWcm*fWcm + fQsquared );
    }
    else if( fEVar == &fW ) {
      fS = fW*fW;
      fWlab = ( fS + fQsquared - fMp*fMp ) / 2.0 / fMp;
      fKlab = sqrt(fWlab*fWlab + fQsquared);
      fWcm = (fWlab*fMp-fQsquared)/sqrt(fMp*fMp+2.0*fWlab*fMp-fQsquared);
      fKcm = sqrt(fWcm*fWcm + fQsquared);
      fXb = fQsquared / 2. / fMp / fWlab;
    }
    else if( fEVar == &fS ) {
      fWlab = ( fS + fQsquared - fMp*fMp ) / 2.0 / fMp;
      fW = sqrt(fS);
      fWcm = (fWlab*fMp-fQsquared)/sqrt(fMp*fMp+2.0*fWlab*fMp-fQsquared);
      fKcm = sqrt( fWcm*fWcm + fQsquared );
      fKlab = sqrt(fWlab*fWlab + fQsquared);
      fXb = fQsquared / 2. / fMp / fWlab;
    }
    else if( fEVar == &fWcm ) {
      fS = std::pow( fWcm + sqrt(fWcm*fWcm+fQsquared+fMp*fMp) ,2.0);
      fWlab = ( fS + fQsquared - fMp*fMp ) / 2.0 / fMp;
      fW = sqrt(fS);
      fKcm = sqrt( fWcm*fWcm + fQsquared );
      fKlab = sqrt(fWlab*fWlab + fQsquared);
      fXb = fQsquared / 2. / fMp / fWlab;
    }
    else if( fEVar == &fKlab ) {
      fWlab = sqrt( fKlab*fKlab - fQsquared );
      fS = fMp*fMp + 2.0 * fMp * fWlab - fQsquared;
      fW = sqrt(fS);
      fWcm = (fWlab*fMp-fQsquared)/sqrt(fMp*fMp+2.0*fWlab*fMp-fQsquared);
      fKcm = sqrt( fWcm*fWcm + fQsquared );
      fXb = fQsquared / 2. / fMp / fWlab;
    }
    else if( fEVar == &fKcm ) {
      fS = std::pow( sqrt(fKcm*fKcm-fQsquared) + sqrt(fKcm*fKcm+fMp*fMp) ,2.0);
      fWlab = ( fS + fQsquared - fMp*fMp ) / 2.0 / fMp;
      fW = sqrt(fS);
      fWcm = (fWlab*fMp-fQsquared)/sqrt(fMp*fMp+2.0*fWlab*fMp-fQsquared);
      fKlab = sqrt(fWlab*fWlab + fQsquared);
      fXb = fQsquared / 2. / fMp / fWlab;
    }

    // Determine pk
    double ek = (fW + (fMk*fMk-fMy*fMy)/fW)/2.;

    if( ek>fMk ) {
      fIsPhysical = true;
      fPk = sqrt(ek*ek-fMk*fMk);
    }
    else {
      fIsPhysical = false;
      fPk = -1.;
    }

    if(fAVar == &fPklab){
      fPYlab = sqrt(pow(fWlab+fMp-sqrt(fMk*fMk+fPklab*fPklab),2.)-fMy*fMy);
      fCosthYlab = (fPklab*fPklab-fPYlab*fPYlab-fKlab*fKlab)/(-2.*fPYlab*fKlab);
      fCosthklab = (fPYlab*fPYlab-fPklab*fPklab-fKlab*fKlab)/(-2.*fPklab*fKlab);
      fT = -fQsquared+fMk*fMk-2.*fWlab*sqrt(fPklab*fPklab+fMk*fMk) + 2.*fPklab*fKlab*fCosthklab;
      fU = -fQsquared+fMy*fMy-2.*fWlab*sqrt(fPYlab*fPYlab+fMy*fMy) + 2.*fPYlab*fKlab*fCosthYlab;
      fCosthkcm = (fT-(fMk*fMk - fQsquared -2.0*fWcm*sqrt(fPk*fPk+fMk*fMk)))/ 2.0 / fPk / fKcm;
      
    }
    else{
      // Determine the angular variables
      if( fAVar == &fCosthkcm ) {
	fT = fMk*fMk - fQsquared -2.0*fWcm*sqrt(fPk*fPk+fMk*fMk) + 2.0*fPk*fKcm*fCosthkcm;
	fU = fMy*fMy - fQsquared -2.0*fWcm*sqrt(fPk*fPk+fMy*fMy) - 2.0*fPk*fKcm*fCosthkcm;
      }
      else if( fAVar == &fT ) {
	fCosthkcm = ( fT + fQsquared + 2.0*fWcm*sqrt(fPk*fPk+fMk*fMk) - fMk*fMk )
	  / 2.0 / fPk / fKcm;
	fU = fMp*fMp + fMk*fMk + fMy*fMy - fQsquared - fS -fT;
      }
      else {
	fT = fMp*fMp + fMk*fMk + fMy*fMy - fQsquared - fS -fU;
	fCosthkcm = ( fT + fQsquared + 2.0*fWcm*sqrt(fPk*fPk+fMk*fMk) - fMk*fMk )
	  / 2.0 / fPk / fKcm;
      }
    
      // Determine kaon momentum in LAB frame by boosting know kaon cm energy to lab
      fPklab = sqrt(pow(((fWlab+fMp)*sqrt(fPk*fPk+fMk*fMk)+fKlab*fPk*fCosthkcm)/fW,2)-fMk*fMk);
      fPYlab = sqrt(pow(fWlab+fMp-sqrt(fMk*fMk+fPklab*fPklab),2.)-fMy*fMy);
      fCosthYlab = (fPklab*fPklab-fPYlab*fPYlab-fKlab*fKlab)/(-2.*fPYlab*fKlab);
      fCosthklab = (fPYlab*fPYlab-fPklab*fPklab-fKlab*fKlab)/(-2.*fPklab*fKlab);
    }
    if( std::fabs(fCosthkcm) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( std::fabs(fCosthkcm)-1.0 < Underflow() ) {
      fIsPhysical = true && fIsPhysical;
      fCosthkcm = ( fCosthkcm>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;

    
  } // end meson production

}

//_____________________________________________________________________
double* TKinematics::GetVarArray(const TString& variable) const
{
  // Return an array of the variable of type 'variable'. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. pklab.
  //
  // The returned pointer has a size GetNumberOfSteps() and is owned by the user.
  TKinematics tk(*this); // temporary copy
  double(TKinematics::*ptrToVarGetter)() const; // pointer to TKinematics getters

  // Assign the correct Getter
  if( !std::strcmp(variable,"wlab") ) {
    ptrToVarGetter = &TKinematics::GetWlab;
  }
  else if( !std::strcmp(variable,"costhkcm") ) {
    ptrToVarGetter = &TKinematics::GetCosthkcm;
  }
  else if( !std::strcmp(variable,"qsquared") ) {
    ptrToVarGetter = &TKinematics::GetQsquared;
  }
  else if( !std::strcmp(variable,"w") ) {
    ptrToVarGetter = &TKinematics::GetW;
  }
  else if( !std::strcmp(variable,"s") ) {
    ptrToVarGetter = &TKinematics::GetS;
  }
  else if( !std::strcmp(variable,"t") || 
	   !std::strcmp(variable,"-t") ) {
    ptrToVarGetter = &TKinematics::GetT;
  }
  else if( !std::strcmp(variable,"wcm") ) {
    ptrToVarGetter = &TKinematics::GetWcm;
  }
  else if( !std::strcmp(variable,"klab") ) {
    ptrToVarGetter = &TKinematics::GetKlab;
  }
  else if( !std::strcmp(variable,"kcm") ) {
    ptrToVarGetter = &TKinematics::GetKcm;
  }
  else if( !std::strcmp(variable,"pk") ) {
    ptrToVarGetter = &TKinematics::GetPk;
  }
  else if( !std::strcmp(variable,"pkcm") ) {
    ptrToVarGetter = &TKinematics::GetPk;
  }  
  else if( !std::strcmp(variable,"pklab") ) {
    ptrToVarGetter = &TKinematics::GetPklab;
  }
  else if( !std::strcmp(variable,"u") ) {
    ptrToVarGetter = &TKinematics::GetU;
  }
  else if( !std::strcmp(variable,"xb") ) {
    ptrToVarGetter = &TKinematics::GetXb;
  }
  else if( !std::strcmp(variable,"phi") ) {
    ptrToVarGetter = &TKinematics::GetPhi;
  }
  else {
    cerr << "ERROR in  TKinematics::GetVarArray(TString): "
	 << "Kinematic variable \"" << variable << "\" unknown!\n";
    exit(1);
  }

  double *array = new double[tk.GetNumberOfSteps()];
  
  // Fill the array
  for(int step=0; step<tk.GetNumberOfSteps(); ++step) {
    tk.GoTo(step);
    array[step] = (tk.*ptrToVarGetter)();
    if(!std::strcmp(variable,"-t")) array[step] *= -1.;
  }    

  return array;
}

//_____________________________________________________________________
double*  TKinematics::GetPhysicalVarArray(const TString& variable) const
{
  // Return an array of the variable of type 'variable' in all physical points. 
  // The variable can effectively by any one for which their exists a getter. 
  // It need not be listed in the info for SetFormat(const char*), eg. pklab.
  //
  // The returned pointer has a size GetNumberOfPhysicalSteps() and is owned by the user.

  TKinematics tk(*this); // temporary copy
  double(TKinematics::*ptrToVarGetter)() const; // pointer to TKinematics getters

  // Assign the correct Getter
  if( !std::strcmp(variable,"wlab") ) {
    ptrToVarGetter = &TKinematics::GetWlab;
  }
  else if( !std::strcmp(variable,"costhkcm") ) {
    ptrToVarGetter = &TKinematics::GetCosthkcm;
  }
  else if( !std::strcmp(variable,"qsquared") ) {
    ptrToVarGetter = &TKinematics::GetQsquared;
  }
  else if( !std::strcmp(variable,"w") ) {
    ptrToVarGetter = &TKinematics::GetW;
  }
  else if( !std::strcmp(variable,"s") ) {
    ptrToVarGetter = &TKinematics::GetS;
  }
  else if( !std::strcmp(variable,"t")  || 
	   !std::strcmp(variable,"-t") ) {
    ptrToVarGetter = &TKinematics::GetT;
  }
  else if( !std::strcmp(variable,"wcm") ) {
    ptrToVarGetter = &TKinematics::GetWcm;
  }
  else if( !std::strcmp(variable,"klab") ) {
    ptrToVarGetter = &TKinematics::GetKlab;
  }
  else if( !std::strcmp(variable,"kcm") ) {
    ptrToVarGetter = &TKinematics::GetKcm;
  }
  else if( !std::strcmp(variable,"pk") ) {
    ptrToVarGetter = &TKinematics::GetPk;
  }
  else if( !std::strcmp(variable,"pkcm") ) {
    ptrToVarGetter = &TKinematics::GetPk;
  }  
  else if( !std::strcmp(variable,"pklab") ) {
    ptrToVarGetter = &TKinematics::GetPklab;
  }
  else if( !std::strcmp(variable,"u") ) {
    ptrToVarGetter = &TKinematics::GetU;
  }
  else if( !std::strcmp(variable,"xb") ) {
    ptrToVarGetter = &TKinematics::GetXb;
  }
  else if( !std::strcmp(variable,"phi") ) {
    ptrToVarGetter = &TKinematics::GetPhi;
  }
  else {
    cerr << "ERROR in  TKinematics::GetVarArray(TString): "
	 << "Kinematic variable \"" << variable << "\" unknown!\n";
    exit(1);
  }

  double *array = new double[tk.GetNumberOfPhysicalSteps()];
  
  // Fill the array
  int physicalStep = 0;
  for(int step=0; step<tk.GetNumberOfSteps(); ++step) {
    tk.GoTo(step);
    if(tk.IsPhysical()) {
      array[physicalStep++] = (tk.*ptrToVarGetter)();
      if(!std::strcmp(variable,"-t")) array[physicalStep-1] *= -1.;
    }
  }    

  return array;     
}

//_____________________________________________________________________
int TKinematics::GetNumberOfVariables() const
{
  // Determine how many variables are not fixed
  int nrOfVar = 0;
  for(int var=1; var<=fNVars; ++var)
    if( !IsFixed(var) ) ++nrOfVar;

  return nrOfVar;
}

//_____________________________________________________________________
void TKinematics::Print(Option_t* option) const
{
  // Print info about the object and the current kinematics.
  cout << GetName() << ": " << GetTitle()
       << "\n****************************************************\n"
       << "Point " << (fIsPhysical ? "on " : "out of ") << "physical plane\n\n"
       << "Q^2   =\t" << std::setw(12) << fQsquared            
       << "\t\tw_lab =\t" << std::setw(12) << fWlab
       << "\nW     =\t" << std::setw(12) << fW
       << "\t\tw_cm  =\t" << std::setw(12) << fWcm
       << "\ncos   =\t" << std::setw(12) << fCosthkcm
       << "\t\tk_lab =\t" << std::setw(12) << fKlab
       << "\ns     =\t" << std::setw(12) << fS
       << "\t\tk_cm  =\t" << std::setw(12) << fKcm
       << "\nt     =\t" << std::setw(12) << fT
       << "\t\tp_K   =\t" << std::setw(12) << fPk
       << "\nu     =\t" << std::setw(12) << fU
       << "\t\tphi   =\t" << std::setw(12) << fPhi
       << "\n****************************************************\n";
}

//_____________________________________________________________________
int TKinematics::VariableInfo() const
{
  // Print info about the range that has been given to the variables.
  cout << GetName() << ": " << GetTitle()
       << "\n****************************************************\n";

  for(int var=0; var<fNVars; ++var) {
    cout << var+1 << ": " << GetVarName(var+1) << " is ";

    if( fIsVar[var] )
      cout << "variable.\n"
	   << "   Range:\t\t" << fLowerLimit[var] << " <-> "
	   << fUpperLimit[var]
	   << "\n   Number of steps:\t" << fNumberOfSteps[var]
	   << "\n   Current step:\t" << GetStep(var+1) << "\n";

    else
      cout << "fixed.\n";
  }
  
  cout << "****************************************************\n";

  return GetNumberOfVariables();
}

//_____________________________________________________________________
bool TKinematics::operator==(const TKinematics& rhs) const
{
  // Evaluate whether both objects currently refer to the same point in the
  // kinematical plane.

  // when the objects are the same, the answer is trivial
  if( this==&rhs ) return true;

  // Each kinematic point is fixed by the isospin channel, 
  // the Mandelstam variables s and t and Q^2
  double key1[5] = {GetIsospin()*1., GetS(), GetT(), GetQsquared(), GetPhi()};
  double key2[5] = {rhs.GetIsospin()*1., rhs.GetS(), rhs.GetT(), rhs.GetQsquared(), rhs.GetPhi()};

  // We create a hash key using the function in the cachetools folder
  char hash1[sizeof(double)*5*2 + 1];
  char hash2[sizeof(double)*5*2 + 1];

  //doublestochar(5, key1, hash1);
  //doublestochar(5, key2, hash2);

  return !std::strcmp(hash1,hash2);
}

//_____________________________________________________________________
void TKinematics::Help()
{
  cout << "You seek help, young Padawan?\n"
       << "1) Open emacs\n"
       << "2) type Alt+x\n"
       << "3) type doctor\n"
       << "4) press return\n";
}
