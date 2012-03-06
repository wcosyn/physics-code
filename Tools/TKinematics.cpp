/*
 * TKinematics.cpp
 *
 * Author:         Pieter Vancraeyveld (pieter.vancraeyveld@UGent.be)
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
// The user needs to specify 3 variables to unambiguously fix the kinematics.
// Obviously, he/she needs to make sure the variabels of choice are independent.
// More info can be found in SetFormat(const char*).
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "TKinematics.h"
// #include <numtoa.h>

using std::sqrt;

ClassImp(TKinematics)

//_____________________________________________________________________
TKinematics::TKinematics(TRootIOCtor* rio)
: TNamed::TNamed(),fFormat(0),fIsospin(1),fMp(0.0),fMk(0.0),
  fMy(0.0),fWlab(0.0),fWcm(0.0),fKlab(0.0),fKcm(0.0),fPk(0.0),fCosthkcm(0.0),
  fCosthklab(0.),fPYlab(0.),fCosthYlab(0.),
  fQsquared(0.0),fXb(0.0), fW(0.0),fS(0.0),fT(0.0),fU(0.0),fIsPhysical(false), 
  fNrOfPhysicals(-1),
  fVar(0),fIsVar(0),fLowerLimit(0),fUpperLimit(0),fStepSize(0),
  fNumberOfSteps(0),fStep(0)
{
  // ROOT I/O Constructor

  // This constructor leaves the object in an inconsistent state. 
  // After calling the default constructor you should make sure to at 
  // least call ReadFormatString() and SetIsospin(int)!
  // Failing to do this will result in a serious ERROR when attempting to
  // TKinematics::Write() the object.

  // Initialize kinematic variables
  fIsVar = new bool[3];
  fLowerLimit = new double[3];
  fUpperLimit = new double[3];
  fStepSize = new double[3];
  fNumberOfSteps = new int[3];

  // Assume the independent variables to be constant
  for(int i=0; i<3; ++i) {
    fIsVar[i] = false;
    fLowerLimit[i] = 0.0;
    fUpperLimit[i] = 0.0;
    fStepSize[i] = 0.0;
    fNumberOfSteps[i] = 1;
  }
}

//_____________________________________________________________________
TKinematics::TKinematics(const char* name, const char* title,
			 int isospin, const char *format, 
			 double var1, double var2, double var3)
  : TNamed::TNamed(name,title),fFormat(0),fIsospin(isospin),fMp(0.0),fMk(0.0),
    fMy(0.0),fWlab(0.0),fWcm(0.0),fKlab(0.0),fKcm(0.0),fPk(0.0),fPklab(0.0),fCosthkcm(0.0),
    fCosthklab(0.),fPYlab(0.),fCosthYlab(0.),
    fQsquared(0.0),fXb(0.0), fW(0.0),fS(0.0),fT(0.0),fU(0.0),fIsPhysical(false), 
    fNrOfPhysicals(-1),
    fVar(0),fIsVar(0),fLowerLimit(0),fUpperLimit(0),fStepSize(0),
    fNumberOfSteps(0),fStep(0)
{
  // Constructor
  //
  // The user needs to specify 
  // * the isospin channel (see SetIsospin(int) for more info).
  // * the 3 independent variables with the 'format' string. See SetFormat(const char*)
  //   for info about what constitutes a valid format string.
  // * the initial values (var1-3) of the independent variables in the same order as listed
  //  in the format string.

  ChangeIsospinChannel(isospin);

  ReadFormatString(format);
  *fVar[0] = var1;
  *fVar[1] = var2;
  *fVar[2] = var3;
  UpdateKinematics();

  // Initialize kinematic variables
  fIsVar = new bool[3];
  fLowerLimit = new double[3];
  fUpperLimit = new double[3];
  fStepSize = new double[3];
  fNumberOfSteps = new int[3];
  
  // Assume the independent variables to be constant
  for(int i=0; i<3; ++i) {
    fIsVar[i] = false;
    fLowerLimit[i] = *fVar[i];
    fUpperLimit[i] = *fVar[i];
    fStepSize[i] = 0.0;
    fNumberOfSteps[i] = 1;
  }
}

//_____________________________________________________________________
TKinematics::TKinematics(const TKinematics& toCopy)
  : TNamed(toCopy),fFormat(0),fIsospin(toCopy.fIsospin),
    fMp(toCopy.fMp),fMk(toCopy.fMk),fMy(toCopy.fMy),
    fWlab(toCopy.fWlab),fWcm(toCopy.fWcm),fKlab(toCopy.fKlab),
    fKcm(toCopy.fKcm),fPk(toCopy.fPk),fPklab(toCopy.fPklab),fCosthkcm(toCopy.fCosthkcm),
    fCosthklab(toCopy.fCosthklab),fPYlab(toCopy.fPYlab),fCosthYlab(toCopy.fCosthYlab),
    fQsquared(toCopy.fQsquared), fXb(toCopy.fXb), fW(toCopy.fW),fS(toCopy.fS),
    fT(toCopy.fT),fU(toCopy.fU),fIsPhysical(toCopy.fIsPhysical),
    fNrOfPhysicals(toCopy.fNrOfPhysicals),
    fVar(0),fIsVar(0),fLowerLimit(0),fUpperLimit(0),fStepSize(0),
    fNumberOfSteps(0),fStep(toCopy.fStep)
{
  // Copy constructor
  if(toCopy.fFormat) // when toCopy.fFormat is non-NULL
    ReadFormatString(toCopy.fFormat);

  if(toCopy.fIsVar) { // when kinematic variables have been initialized
    fIsVar = new bool[3];
    fLowerLimit = new double[3];
    fUpperLimit = new double[3];
    fStepSize = new double[3];
    fNumberOfSteps = new int[3];

    for(int i=0; i<3; ++i) {
      fIsVar[i] = toCopy.fIsVar[i];
      fLowerLimit[i] = toCopy.fLowerLimit[i];
      fUpperLimit[i] = toCopy.fUpperLimit[i];
      fStepSize[i] = toCopy.fStepSize[i];
      fNumberOfSteps[i] = toCopy.fNumberOfSteps[i];
    }
  }
}

//_____________________________________________________________________
TKinematics::~TKinematics()
{
  // Destructor
  delete[] fFormat;
  delete[] fVar;
  delete[] fIsVar;
  delete[] fLowerLimit;
  delete[] fUpperLimit;
  delete[] fStepSize;
  delete[] fNumberOfSteps;
}

//_____________________________________________________________________
TKinematics& TKinematics::operator=(const TKinematics& toCopy)
{
  // Assignment
  if( this!=&toCopy ) { // avoid self-assignment

    TNamed::operator=(toCopy);
    
    fIsospin = toCopy.fIsospin;
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
    fIsPhysical = toCopy.fIsPhysical;
    fNrOfPhysicals = toCopy.fNrOfPhysicals;
    fStep = toCopy.fStep;

    if(toCopy.fIsVar) {
      if( !fIsVar) fIsVar = new bool[3];
      if( !fLowerLimit) fLowerLimit = new double[3];
      if( !fUpperLimit) fUpperLimit = new double[3];
      if( !fStepSize) fStepSize = new double[3];
      if( !fNumberOfSteps) fNumberOfSteps = new int[3];
      
      for(int i=0; i<3; ++i) {
	fIsVar[i] = toCopy.fIsVar[i];
	fLowerLimit[i] = toCopy.fLowerLimit[i];
	fUpperLimit[i] = toCopy.fUpperLimit[i];
	fStepSize[i] = toCopy.fStepSize[i];
	fNumberOfSteps[i] = toCopy.fNumberOfSteps[i];
      }
    }
    else {
      if( fIsVar) { delete[] fIsVar; fIsVar = 0; }
      if( fLowerLimit) { delete[] fLowerLimit; fLowerLimit = 0; }
      if( fUpperLimit) { delete[] fUpperLimit; fUpperLimit = 0; }
      if( fStepSize) { delete[] fStepSize; fStepSize = 0; }
      if( fNumberOfSteps) { delete[] fNumberOfSteps; fNumberOfSteps = 0; }
    }
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
    std::cerr << "WARNING in TKinematics::SetIsospin(int): "
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
      std::cerr << "ERROR in TKinematics::ChangeIsospinChannel(int):\n"
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
  The user needs to specify 3 variables to unambiguously fix the kinematics. Obviously,
  he/she needs to make sure the variables of choice are independent.
  TKinematics needs:
  <ul>
  <li> 1 energy variable
  <li> 1 angular variable
  <li> Q<sup>2</sup>
  </ul>
  <p>
  The 3 independent variables are set by providing a ':' separated list of variables.
  <p>
  Below is a list of the 'code names' of variables:
  <table>
  <tr><td>Q<sup>2</sup></td><td><tt>qsquared</tt></td></tr>
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
  End_Html */

  if(!IsFixed()) {
    std::cout << "INFO in TKinematics::SetFormat(const char*): "
	      << "Fixing all variables.\n";
    for(int i=0; i<3; ++i)
      FixVariable(i+1);
  }

  ReadFormatString(format);
  
  // Finally, set the correct upper and lower limits
  for(int i=0; i<3; ++i) {
    fLowerLimit[i] = *fVar[i];
    fUpperLimit[i] = *fVar[i];
  }
}

//_____________________________________________________________________
void TKinematics::ReadFormatString(const char* format)
{
  // Process the format string

  // Store format in TKinematics::fFormat
  // Read TKinematics::fFormat and let fVar[] point
  // to the 3 independent kinematical points as specified
  // by the user.
  // Check that the choosen variables make sense.


  // allocate memory for fFormat and fVar
  if(!fFormat) fFormat = new char[100];
  if(!fVar)    fVar    = new double*[3];

  // Put format string in member and in temporary *char.
  std::strcpy(fFormat,format);
  char tmpchar[100]; std::strcpy(tmpchar,format);

  // Tokenize format string
  // and determine independent kinematic variables
  char* token; int nrOfTokens =0;
  int qsqCount=0, energyCount=0, angleCount=0; // count type of variables

  token = std::strtok(tmpchar,":"); // read first token
  while( token != NULL  &&  nrOfTokens < 3 ) 
    {
      if( !std::strcmp(token,"wlab") ) {
	fVar[nrOfTokens] = fEVar = &fWlab;
	++energyCount;
      }
      else if( !std::strcmp(token,"costhkcm") ) {
	fVar[nrOfTokens] = fAVar = &fCosthkcm;
	++angleCount;
      }
      else if( !std::strcmp(token,"qsquared") ) {
	fVar[nrOfTokens] = &fQsquared;
	++qsqCount;
      }
      else if( !std::strcmp(token,"w") ) {
	fVar[nrOfTokens] = fEVar = &fW;
	++energyCount;
      }
      else if( !std::strcmp(token,"s") ) {
	fVar[nrOfTokens] = fEVar = &fS;
	++energyCount;
      }
      else if( !std::strcmp(token,"t") ) {
	fVar[nrOfTokens] = fAVar = &fT;
	++angleCount;
      }
      else if( !std::strcmp(token,"wcm") ) {
	fVar[nrOfTokens] = fEVar = &fWcm;
	++energyCount;
      }
      else if( !std::strcmp(token,"klab") ) {
	fVar[nrOfTokens] = fEVar = &fKlab;
	++energyCount;
      }
      else if( !std::strcmp(token,"kcm") ) {
	fVar[nrOfTokens] = fEVar = &fKcm;
	++energyCount;
      }
      else if( !std::strcmp(token,"pk") ) {
	fVar[nrOfTokens] = fEVar = &fPk;
	++energyCount;
      }
      else if( !std::strcmp(token,"pkcm") ) {
	fVar[nrOfTokens] = fEVar = &fPk;
	++energyCount;
      }
      else if( !std::strcmp(token,"pklab") ) {
	fVar[nrOfTokens] = fAVar = &fPklab;
	++angleCount;
      }
      else if( !std::strcmp(token,"u") ) {
	fVar[nrOfTokens] = fAVar = &fU;
	++angleCount;
      }
      else if( !std::strcmp(token,"xb") ) {
	fVar[nrOfTokens] = fEVar = &fXb;
	++energyCount;
      }
      else {
	fVar[nrOfTokens] = NULL;
	std::cerr << "ERROR in  TKinematics::ReadFormatString(const char*):\n"
		  << "Kinematic variable \"" << token << "\" unknown!\n";
	exit(1);
      }
      nrOfTokens++;
      token = strtok(NULL,":"); // read next token
    }  // end tokenization

  if( nrOfTokens < 3 ) {
    std::cerr << "ERROR in TKinematics::ReadFormatString(const char*):\n"
	      << "Format string should provide 3 independent variables: "
	      << "\"var1:var2:var3\"\n";
    exit(1);
  }

  // Check whether the user provided a meaningful format string
  if( qsqCount != 1 ) {
    std::cerr << "ERROR in TKinematics::ReadFormatString(const char*):\n"
	      << "Q^2 needs to be one of 3 independent kinematic variables.\n";
    exit(1);
  }
  if( energyCount == 0 ) {
    std::cerr << "ERROR in TKinematics::ReadFormatString(const char*):\n"
	      << "wlab/wcm/klab/kcm/pk/pklab/w/s/xb needs to be one of 3 "
	      << "independent kinematic variables.\n";
    exit(1); }
  else if( energyCount != 1 ) {
    std::cerr << "ERROR in TKinematics::ReadFormatString(const char*):\n"
	      << "Only one of following kinematic variables can be independent: "
	      << "wlab/wcm/klab/kcm/pk/pklab/w/s/xb.\n";
    exit(1);
  }
  if( angleCount == 0 ) {
    std::cerr << "ERROR in TKinematics::ReadFormatString(const char*):\n"
	      << "costhkcm/t/u needs to be one of 3 "
	      << "independent kinematic variables.\n";
    exit(1); }
  else if( angleCount != 1 ) {
    std::cerr << "ERROR in TKinematics::ReadFormatString(const char*):\n"
	      << "Only one of following kinematic variables can be independent: "
	      << "costhkcm/t/u.\n";
    exit(1);
  }
  if( fIsospin == 11 || fIsospin == 12 ) {
    if( fEVar==&fKlab || fEVar==&fWlab || fEVar==&fXb ) {
      std::cerr << "ERROR in  TKinematics::ReadFormatString(const char*):\n"
		<< "Kinematic variable wlab/klab/xb is not allowed "
		<< "as independent variable in radiative capture kinematics!\n";
      exit(1);
    }
  }
  else {
    if( fEVar==&fPk || fEVar==&fPklab ) {
      std::cerr << "ERROR in  TKinematics::ReadFormatString(const char*):\n"
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
  if( varNumber<1 || varNumber>3 ) {
    std::cerr << "ERROR in TKinematics::GetVarName(int):\n"
	      << "Index out-of-bounds!\n";
    return TString();
  }

  if( !fVar ) {
    std::cerr << "ERROR in TKinematics::GetVarName(int):\n"
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
  else {
    std::cerr << "Error in TKinematics::GetVarName(int):\n"
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

  if( varNumber<1 || varNumber>3 )
    std::cerr << "ERROR in TKinematics::SetVar(int,double): "
	      << "Index out-of-bounds!\n";

  else if( !fVar )
    std::cerr << "ERROR in TKinematics::SetVar(int,double): "
	      << "Independent kinematic variables have not been initialized!\n";
    
  else if( fIsVar[varNumber-1] )
    std::cerr << "ERROR in TKinematics::SetVar(int,double): "
	      << "This method can not be used because the variable is not fixed.\n";

  else {
    fLowerLimit[varNumber-1] = value;
    fUpperLimit[varNumber-1] = value;
    GoTo(0);
    fNrOfPhysicals = -1; // Number of physical points may change
  }
}

//_____________________________________________________________________
double TKinematics::GetVar(int varNumber) const
{
  // Returns the current value of the 'varNumber'th independent
  // variable as specified in the format string.

  // First we do some checks
  if( varNumber<1 || varNumber>3 ) {
    std::cerr << "ERROR in TKinematics::GetVar(int):\n"
	      << "Index out-of-bounds!\n";
    return -1e16;
  }

  if( !fVar ) {
    std::cerr << "ERROR in TKinematics::GetVar(int):\n"
	      << "Independent kinematic variables have not been initialized!\n";
    return -1e16;
  }

  return *fVar[varNumber-1];
  
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
    R__b >> fIsPhysical;
    if(R__v<3) fNrOfPhysicals = -1; // added version 3
    else R__b >> fNrOfPhysicals;
    delete [] fIsVar;
    fIsVar = new bool[3];
    R__b.ReadFastArray(fIsVar,3);
    delete [] fLowerLimit;
    fLowerLimit = new double[3];
    R__b.ReadFastArray(fLowerLimit,3);
    delete [] fUpperLimit;
    fUpperLimit = new double[3];
    R__b.ReadFastArray(fUpperLimit,3);
    delete [] fStepSize;
    fStepSize = new double[3];
    R__b.ReadFastArray(fStepSize,3);
    delete [] fNumberOfSteps;
    fNumberOfSteps = new int[3];
    R__b.ReadFastArray(fNumberOfSteps,3);
    R__b >> fStep;
    R__b.CheckByteCount(R__s, R__c, TKinematics::IsA());
    ReadFormatString(tempFormat); // Added by PVC
    delete[] tempFormat; // Added by PVC
    if (R__v < 2) UpdateKinematics();  // added version 2
  } else {
    R__c = R__b.WriteVersion(TKinematics::IsA(), kTRUE);
    TNamed::Streamer(R__b);
    R__b.WriteFastArray(fFormat,100);
    R__b << fIsospin;
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
    R__b << fIsPhysical;
    R__b << fNrOfPhysicals;
    R__b.WriteFastArray(fIsVar,3);
    R__b.WriteFastArray(fLowerLimit,3);
    R__b.WriteFastArray(fUpperLimit,3);
    R__b.WriteFastArray(fStepSize,3);
    R__b.WriteFastArray(fNumberOfSteps,3);
    R__b << fStep;
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
  if( varNumber<1 || varNumber>3 )
    std::cerr << "ERROR in TKinematics::SetVarRange(int,double,double,int):\n"
	      << "Index out-of-bounds!\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics::SetVarRange(int,double,double,int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else if( numberOfSteps<2 )
    std::cerr << "ERROR in TKinematics::SetVarRange(int,double,double,int)\n"
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
  if( varNumber<1 || varNumber>3 )
    std::cerr << "ERROR in TKinematics::FixVariable(int):\n"
	      << "Index out-of-bounds!\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics::FixVariable(int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else if( fIsVar[varNumber-1] ) {
    // Determine current step
    int step1 = GetStep(1);
    int step2 = GetStep(2);
    int step3 = GetStep(3);

    fIsVar[varNumber-1] = false;
    
    // Reset all variables
    fLowerLimit[varNumber-1] = *fVar[varNumber-1];
    fUpperLimit[varNumber-1] = *fVar[varNumber-1];
    fStepSize[varNumber-1] = 0.0;
    fNumberOfSteps[varNumber-1] = 1;

    // Go to first step
    switch(varNumber) {
    case 1: GoTo(0,step2,step3); break;
    case 2: GoTo(step1,0,step3); break;
    case 3: GoTo(step1,step2,0); break;
    }

    // Number of physical points may change
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
void TKinematics::FixVariables()
{
  // Calls FixVariable(int) on all variables
  for(int var=1; var<=3; ++var)
    FixVariable(var);
}

//_____________________________________________________________________
void TKinematics::GoTo(int totalStep)
{
  //  Proceed to the 'totalStep'th point in the grid of independent variables.

  // Total number of steps
  int totNrOfSteps = 1;
  for(int i=0; i<3; ++i)
    totNrOfSteps *= fNumberOfSteps[i];

  // First we do some checks
  if( totalStep<0 || totalStep>=totNrOfSteps )
    std::cerr << "ERROR in TKinematics::GoTo(int):\n"
	      << "Invalid step. Should be in range [0,"
	      << totNrOfSteps-1 << "].\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics::GoTo(int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else {
    fStep = totalStep;
    
    // Determine step per variable
    int stepVar[3];
    stepVar[2] = totalStep / (fNumberOfSteps[1]*fNumberOfSteps[0]);
    stepVar[1] = (totalStep % (fNumberOfSteps[1]*fNumberOfSteps[0]))/ fNumberOfSteps[0];
    stepVar[0] = ((totalStep % (fNumberOfSteps[1]*fNumberOfSteps[0])) % fNumberOfSteps[0]);

    // Set variables
    for(int i=0; i<3; ++i)
      *fVar[i] = fLowerLimit[i] + stepVar[i] * fStepSize[i];

    UpdateKinematics();
  }
}

//_____________________________________________________________________
void TKinematics::GoTo(int stepVar1, int stepVar2, int stepVar3)
{
  // Proceed to point in the grid of independent variables specified by the
  // arguments.
  // Each stepVari fixes the position in the grid of the 'i'th variable.

  // First we do some checks
  if( stepVar1<0 || stepVar1>=fNumberOfSteps[0] )
    std::cerr << "ERROR in TKinematics::GoTo(int,int,int):\n"
	      << "Invalid \"" << GetVarName(1) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[0]-1 << "].\n";

  else if( stepVar2<0 || stepVar2>=fNumberOfSteps[1] )
    std::cerr << "ERROR in TKinematics::GoTo(int,int,int):\n"
	      << "Invalid \"" << GetVarName(2) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[1]-1 << "].\n";

  else if( stepVar3<0 || stepVar3>=fNumberOfSteps[2] )
    std::cerr << "ERROR in TKinematics::GoTo(int,int,int):\n"
	      << "Invalid \"" << GetVarName(3) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[2]-1 << "].\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics::GoTo(int,int,int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else {
    fStep = stepVar1 
      + stepVar2 * fNumberOfSteps[0]
      + stepVar3 * fNumberOfSteps[0] * fNumberOfSteps[1];
    
    // Set variables
    *fVar[0] = fLowerLimit[0] + stepVar1 * fStepSize[0];
    *fVar[1] = fLowerLimit[1] + stepVar2 * fStepSize[1];
    *fVar[2] = fLowerLimit[2] + stepVar3 * fStepSize[2];

    UpdateKinematics();
  }
}

//_____________________________________________________________________
bool TKinematics::IsFixed() const
{
  // Check whether any of the independent variables has been given a range.
  return !fIsVar[0] && !fIsVar[1] && !fIsVar[2];
}

//_____________________________________________________________________
bool TKinematics::IsFixed(int varNumber) const
{
  //  Check whether the 'varNumber'th variable has been given a range.

  // First we do some checks
  if( varNumber<1 || varNumber>3 ) {
    std::cerr << "ERROR in TKinematics::IsFixed(int):\n"
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

  // Total number of steps
  int totNrOfSteps = 1;
  for(int i=0; i<3; ++i) totNrOfSteps *= fNumberOfSteps[i];

  // Check whether we are at the end of the grid.
  if( fStep == (totNrOfSteps-1) ) return 0;

  GoTo(++fStep);
  return 1;
}

//_____________________________________________________________________
int TKinematics::GetNumberOfSteps() const
{
  // Number of points in the grid of independent variables
  return fNumberOfSteps[0] * fNumberOfSteps[1] * fNumberOfSteps[2];
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
  if( varNumber<1 || varNumber>3 ) {
    std::cerr << "ERROR in TKinematics::GetNumberOfSteps(int):\n"
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
  if( varNumber<1 || varNumber>3 ) {
    std::cerr << "ERROR in TKinematics::GetStepSize(int):\n"
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

  switch( varNumber ) 
    {
    case 1:
      return ((fStep % (fNumberOfSteps[1]*fNumberOfSteps[0])) % fNumberOfSteps[0]);
      break;

    case 2:
      return (fStep % (fNumberOfSteps[1]*fNumberOfSteps[0]))/ fNumberOfSteps[0];
      break;

    case 3:
      return fStep / (fNumberOfSteps[1]*fNumberOfSteps[0]);
      break;
      
    default:
      std::cerr << "ERROR in TKinematics::GetStep(int):\n"
		<< "Index out-of-bounds!\n";
    }

  return -1;
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
      std::cerr << "ERROR in TKinematics::UpdateKinematics(): "
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
    // pk can be found by solving the implicit TKinematics::EnergyConservation()
    // equation. We use the bisection method, which looks for a solution in the
    // interval [lowerBound=0.0,upperBound=fW].
    double lowerBound=0.0, upperBound=fW;
  
    // First we need to make sure that the function changes sign in this region.
    if( EnergyConservation(lowerBound) > 0.0 || EnergyConservation(upperBound) < 0.0 ) {
      fPk = -1.0;
      fIsPhysical = false;
    }
    else {
      double delta = upperBound - lowerBound;
      while( delta > 0.00001 ) {
	if( EnergyConservation(lowerBound+delta/2.0) > 0.0 )
	  upperBound = lowerBound + delta/2.0;
	else
	  lowerBound += delta/2.0;

	delta = upperBound - lowerBound;
      }

      fPk = lowerBound + delta/2.0;
      if( fPk < 0.0 )
	fIsPhysical = false;
      else
	fIsPhysical = true;
    }

    if(fAVar == &fPklab){
      fPYlab = sqrt(pow(fWlab+fMp-sqrt(fMk*fMk+fPklab*fPklab),2.)-fMy*fMy);
      fCosthYlab = (fPklab*fPklab-fPYlab*fPYlab-fKlab*fKlab)/(-2.*fPYlab*fKlab);
      fCosthklab = (fPYlab*fPYlab-fPklab*fPklab-fKlab*fKlab)/(-2.*fPklab*fKlab);
      fT = -fQsquared+fMk*fMk-2.*fWlab*sqrt(fPklab*fPklab+fMk*fMk) + 2.*fPklab*fKlab*fCosthklab;
      fU = -fQsquared+fMy*fMy-2.*fWlab*sqrt(fPYlab*fPYlab+fMy*fMy) + 2.*fPYlab*fKlab*fCosthYlab;
      fCosthkcm = fT-(fMk*fMk - fQsquared -2.0*fWcm*sqrt(fPk*fPk+fMk*fMk));
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
  else {
    std::cerr << "ERROR in  TKinematics::GetVarArray(TString): "
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
  else {
    std::cerr << "ERROR in  TKinematics::GetVarArray(TString): "
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
  for(int var=1; var<=3; ++var)
    if( !IsFixed(var) ) ++nrOfVar;

  return nrOfVar;
}

//_____________________________________________________________________
void TKinematics::Print(Option_t* option) const
{
  // Print info about the object and the current kinematics.
  std::cout << GetName() << ": " << GetTitle()
	    << "\n****************************************************\n"
	    << "Point " << (fIsPhysical ? "on " : "out of ") << "physical plane\n\n"
	    << "Q^2   =\t" << std::setw(12) << fQsquared
	    << "\nw_lab =\t" << std::setw(12) << fWlab
	    << "\t\tk_lab =\t" << std::setw(12) << fKlab
	    << "\nw_cm  =\t" << std::setw(12) << fWcm
	    << "\t\tk_cm  =\t" << std::setw(12) << fKcm
	    << "\nW     =\t" << std::setw(12) << fW
	    << "\t\tp_K   =\t" << std::setw(12) << fPk
	    << "\ns     =\t" << std::setw(12) << fS
	    << "\t\tu     =\t" << std::setw(12) << fU
	    << "\ncos   =\t" << std::setw(12) << fCosthkcm
	    << "\t\tt     =\t" << std::setw(12) << fT
	    << "\n****************************************************\n";
}

//_____________________________________________________________________
int TKinematics::VariableInfo() const
{
  // Print info about the range that has been given to the variables.
  std::cout << GetName() << ": " << GetTitle()
	    << "\n****************************************************\n";

  for(int var=0; var<3; ++var) {
    std::cout << var+1 << ": " << GetVarName(var+1) << " is ";

    if( fIsVar[var] )
      std::cout << "variable.\n"
		<< "   Range:\t\t" << fLowerLimit[var] << " <-> "
		<< fUpperLimit[var]
		<< "\n   Number of steps:\t" << fNumberOfSteps[var]
		<< "\n   Current step:\t" << GetStep(var+1) << "\n";

    else
      std::cout << "fixed.\n";
  }
  
  std::cout << "****************************************************\n";

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
  float key1[4] = {GetIsospin(), GetS(), GetT(), GetQsquared()};
  float key2[4] = {rhs.GetIsospin(), rhs.GetS(), rhs.GetT(), rhs.GetQsquared()};

  // We create a hash key using the function in the cachetools folder
  char hash1[sizeof(float)*4*2 + 1];
  char hash2[sizeof(float)*4*2 + 1];

  //floatstochar(4, key1, hash1);
  //floatstochar(4, key2, hash2);

  return !std::strcmp(hash1,hash2);
}

//_____________________________________________________________________
void TKinematics::Help()
{
  std::cout << "You seek help, young Padawan?\n"
	    << "1) Open emacs\n"
	    << "2) type Alt+x\n"
	    << "3) type doctor\n"
	    << "4) press return\n";
}
