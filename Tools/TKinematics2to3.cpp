/*
 * TKinematics2to3.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * 
 */

///////////////////////////////////////////////////////////////////////////
/* Begin_Html
   <center><h2>TKinematics2to3</h2></center>
   TKinematics2to3 allows to calculate kinematics for g D -> K Y N.
   <p>
   All quantities are in MeV or rad.
   <p>


   <h3>Reference frame</h3>
   This class solves the kinematics in the LAB frame, where the deuteron
   is at rest. The <b>axis convention</b>, however, can be chosen in a number of
   different ways. The user can specify this with the <b>EAxisConvention</b> type.
   <p>
   The z axis is always along the photon momentum and we define x ~ y x z.
   This leaves us the choice of y axis:
   <ul>
   <li><b>kN</b>:&nbsp  y ~ p<sub>g</sub> x p<sub>N</sub>
   <li><b>kK</b>:&nbsp  y ~ p<sub>g</sub> x p<sub>K</sub>
   <li><b>kY</b>:&nbsp  y ~ p<sub>g</sub> x p<sub>Y</sub>
   </ul>
   Internally, TKinematics2to3 uses the kN convention. Getters that depend
   on the type of axis convention, take an EAxisConvention as argument.  
   The constant <b>kS</b> is defined as the axis convention chosen by the user  
   in the constructor. This allows the user to omit the EAxisConvention  
   argument in getters, which defaults to kS.
   <p>
   Some getters return quantities in the centre-of-mass (CM) of a 2-particle
   subsystem. These getters don't use the same axis as the LAB frame getters.
   The y axis is still fixed with the axis convention, but the x and z axis are
   rotated such that the z axis lies along the total momentum of the 2 particles.


   <h3>Isospin channels</h3>
   TKinematics2to3 is designed to work with a fixed number of final states, confusingly labeled <b>isospin channels</b>.
   <ul>
   <li> (1)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>+</sup> + L<sup>0</sup> + n
   <li> (2)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>+</sup> + S<sup>0</sup> + n
   <li> (3)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>0</sup> + S<sup>+</sup> + n
   <li> (4)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>0</sup> + L<sup>0</sup> + p
   <li> (5)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>0</sup> + S<sup>0</sup> + p
   <li> (6)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>+</sup> + S<sup>-</sup> + p
   </ul>


   <a name="indepvar"></a><h3>Independent variables</h3>
   The user needs to specify 6 variables to unambiguously fix the kinematics.
   Obviously, he/she needs to made sure the variables of choice are independent. 
   TKinematics2to3 is capable of dealing with 2 situations.
   <a name="indepvarCaseA"></a><h3><h4>Case A: Mother with children</h4>
   We consider the final state as the sum of one particle (mother) and the subsystem
   of the two other particles (children). The mother should be the same particle
   as the one that fixes the x and y axis (axis convention).
   <ul>
   <li>The user provides info about the energy and the theta angle (either cos or angle) of the mother.
   <li>A 3rd variable should either specify
   <ul><li>the photon's energy
   <li>the invariant mass of the total system
   <li>the invariant mass of the children subsystem</ul>
   <li>Finally, the user specifies the solid angle (theta+phi) of one of the children in the CM of the subsystem.
   <li>The virtualiy Q<sup>2</sup> of the photon is always one of the independent variables.
   </ul> 
   Note: the user can provide the solid angle of one of the children in the LAB frame 
   as well. Such a set of independent variables, however, can have two solutions. This
   class is not capable of dealing with this and will print a warning before arbitrarily
   choosing one of both solutions. See TKinematics2to3WithLabAngles for a kinematics
   class that is designed to deal with this situation.
   <h4>Case B: Papa Mandelstam</h4>
   All independent variables are Mandelstam variables, with the exception of Q<sup>2</sup> of course.
   Again, we define particle 1 as the one that fixes the x and y axis.
   <ul><li>s<sub>tot</sub> should always be given
   <li>For the other variables the user has 7 options:
   <ul><li>s<sub>12</sub>, s<sub>23</sub>, t<sub>g1</sub>, t<sub>d3</sub>
   <li>s<sub>12</sub>, s<sub>23</sub>, t<sub>g1</sub>, t<sub>d2</sub>
   <li>s<sub>12</sub>, s<sub>23</sub>, t<sub>d2</sub>, t<sub>d3</sub>
   <li>s<sub>12</sub>, t<sub>d2</sub>, t<sub>g1</sub>, t<sub>d3</sub>
   <li>s<sub>12</sub>, s<sub>23</sub>, t<sub>g1</sub>, t<sub>g2</sub>
   <li>s<sub>12</sub>, s<sub>23</sub>, t<sub>g2</sub>, t<sub>d3</sub>
   <li>t<sub>g2</sub>, s<sub>23</sub>, t<sub>g1</sub>, t<sub>d3</sub>
   </ul></ul>
   <p>
   End_Html */
///////////////////////////////////////////////////////////////////////////

#include "TKinematics2to3.h"
#include "TMPI.h"
#include "TLorentzQuaternion.h"
#include "constants.hpp"
#include <iostream>
#include <cstring>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <Math/Polynomial.h>
#include <vector>
#include <TLorentzRotation.h>
#include <TVector2.h>
#include <TVector3.h>
// #include <numtoa.h>
using std::cout; using std::cerr;
using ROOT::Math::Polynomial;
using std::vector;
using std::strcmp;
using TMath::Sqrt; using TMath::Cos; using TMath::Sin;
using TMath::ACos; using TMath::Abs; using TMath::Power;

ClassImp(TKinematics2to3)

//_____________________________________________________________________
const double TKinematics2to3::kUflow = 1.e-7;

//_____________________________________________________________________
TKinematics2to3::TKinematics2to3(const char *name, const char *title,
				 int isospin, EAxisConvention ac, 
				 const char *format,
				 double var1, double var2, double var3,
				 double var4, double var5, double var6,
				 double Md, double Mn, double Mk, double My)
  : TNamed(name,title), fKYsystem(0), fKNsystem(0), fYNsystem(0),
    fKvector(0), fYvector(0), fNvector(0), fGvector(0), fDvector(0),
    fIsospin(0), fMd(Md), fMk(Mk), fMy(My), fMn(Mn), fIsPhysical(false),
    fQsquared(0), fWlab(0), fKlab(0), fEk(0), fEy(0), fEn(0),
    fPk(0), fPy(0), fPn(0), fThk(0), fThy(0), fThn(0), fCosthk(0),
    fCosthy(0), fCosthn(0), fPhik(0), fPhiy(0), fPhin(0), fWtot(0),
    fWky(0), fWkn(0), fWyn(0), fStot(0), fSky(0), fSkn(0), fSyn(0), fTgk(0),
    fTgy(0), fTgn(0), fKYcosthkcm(0), fKYcosthycm(0), fKNcosthkcm(0),
    fKNcosthncm(0), fYNcosthycm(0), fYNcosthncm(0),
    fVar(0), fE1Var(0), fTh1Var(0), fE2Var(0), fTh2Var(0), fPhi2Var(0),
    fTdk(0), fTdy(0), fTdn(0), fKYphikcm(0), fKYphiycm(0), fKNphikcm(0),
    fKNphincm(0), fYNphiycm(0), fYNphincm(0), fKYthkcm(0), fKYthycm(0),
    fKNthkcm(0), fKNthncm(0), fYNthycm(0), fYNthncm(0),
    fEVar(0), fIsVar(0), fLowerLimit(0), fUpperLimit(0),
    fStepSize(0), fNumberOfSteps(0), fStep(0), fFormat(0),
    fAC(ac), fKinSolution(0), fNrOfPhysicals(-1)
{
  // Constructor
  //
  // See SetFormat(const char*) for info about what constitutes a valid
  // format string. 
  //
  // var1-6 will set the initial value of the independent variables in the 
  // same order as listed in the format string.
  //
  // Make sure the independent variables you choose in the format string are
  // compatible with your axis convention.

  if( ac==kS ) {
    cerr << "ERROR in TKinematics2to3::TKinematics2to3(..): "
	 << "TKinematics2to3::kS is NOT a valid axis convention.\n";
    exit(1);
  }

  ChangeIsospinChannel(isospin);

  // allocate memory for fFormat, fVar and f*vector
  fFormat = new char[100];
  fKinSolution = new char[4];
  fVar    = new double*[6];
  fKvector = new TLorentzVector;
  fYvector = new TLorentzVector;
  fNvector = new TLorentzVector;
  fGvector = new TLorentzVector;
  fDvector = new TLorentzVector;

  ReadFormatString(format);
  *fVar[0] = var1;
  *fVar[1] = var2;
  *fVar[2] = var3;
  *fVar[3] = var4;
  *fVar[4] = var5;
  *fVar[5] = var6;

  InitializeSubsystems();

  UpdateKinematics();
  
  // Initialize kinematic variables
  fIsVar = new bool[6];
  fLowerLimit = new double[6];
  fUpperLimit = new double[6];
  fStepSize = new double[6];
  fNumberOfSteps = new int[6];
  
  // Assume the independent variables to be constant
  for(int i=0; i<6; ++i) {
    fIsVar[i] = false;
    fLowerLimit[i] = *fVar[i];
    fUpperLimit[i] = *fVar[i];
    fStepSize[i] = 0.0;
    fNumberOfSteps[i] = 1;
  }
}


//_____________________________________________________________________
TKinematics2to3::TKinematics2to3(const TKinematics2to3& toCopy)
  : TNamed(toCopy), fKYsystem(0), fKNsystem(0), fYNsystem(0),
    fKvector(0), fYvector(0), fNvector(0), fGvector(0), fDvector(0),
    fIsospin(toCopy.fIsospin), fMd(toCopy.fMd), fMk(toCopy.fMk), 
    fMy(toCopy.fMy), fMn(toCopy.fMn), fIsPhysical(toCopy.fIsPhysical),
    fQsquared(toCopy.fQsquared), fWlab(toCopy.fWlab), fKlab(toCopy.fKlab), 
    fEk(toCopy.fEk), fEy(toCopy.fEy), fEn(toCopy.fEn), fPk(toCopy.fPk), 
    fPy(toCopy.fPy), fPn(toCopy.fPn), fThk(toCopy.fThk), fThy(toCopy.fThy), 
    fThn(toCopy.fThn), fCosthk(toCopy.fCosthk), fCosthy(toCopy.fCosthy), 
    fCosthn(toCopy.fCosthn), fPhik(toCopy.fPhik), fPhiy(toCopy.fPhiy), 
    fPhin(toCopy.fPhin), fWtot(toCopy.fWtot), fWky(toCopy.fWky), fWkn(toCopy.fWkn), 
    fWyn(toCopy.fWyn), fStot(toCopy.fStot), fSky(toCopy.fSky), fSkn(toCopy.fSkn), 
    fSyn(toCopy.fSyn), fTgk(toCopy.fTgk), fTgy(toCopy.fTgy), fTgn(toCopy.fTgn), 
    fKYcosthkcm(toCopy.fKYcosthkcm), fKYcosthycm(toCopy.fKYcosthycm), 
    fKNcosthkcm(toCopy.fKNcosthkcm), fKNcosthncm(toCopy.fKNcosthncm), 
    fYNcosthycm(toCopy.fYNcosthycm), fYNcosthncm(toCopy.fYNcosthncm),
    fVar(0), fE1Var(0), fTh1Var(0), fE2Var(0), fTh2Var(0), fPhi2Var(0),
    fTdk(toCopy.fTdk), fTdy(toCopy.fTdy), fTdn(toCopy.fTdy), 
    fKYphikcm(toCopy.fKYphikcm), fKYphiycm(toCopy.fKYphiycm), 
    fKNphikcm(toCopy.fKNphikcm), fKNphincm(toCopy.fKNphincm), 
    fYNphiycm(toCopy.fYNphiycm), fYNphincm(toCopy.fYNphincm), 
    fKYthkcm(toCopy.fKYthkcm), fKYthycm(toCopy.fKYthycm), 
    fKNthkcm(toCopy.fKNthkcm), fKNthncm(toCopy.fKNthncm), 
    fYNthycm(toCopy.fYNthycm), fYNthncm(toCopy.fYNthncm),
    fEVar(0), fIsVar(0), fLowerLimit(0), fUpperLimit(0),
    fStepSize(0), fNumberOfSteps(0), fStep(toCopy.fStep), fFormat(0),
    fAC(toCopy.fAC), fKinSolution(0), fNrOfPhysicals(toCopy.fNrOfPhysicals)
{
  // Copy constructor

  if(toCopy.fVar) fVar =  new double*[6];
  if(toCopy.fFormat) {
    fFormat = new char[100];
    fKinSolution = new char[4];
    ReadFormatString(toCopy.fFormat);
  }

  fKvector = new TLorentzVector(*toCopy.fKvector);
  fYvector = new TLorentzVector(*toCopy.fYvector);
  fNvector = new TLorentzVector(*toCopy.fNvector);
  fGvector = new TLorentzVector(*toCopy.fGvector);
  fDvector = new TLorentzVector(*toCopy.fDvector);

  if(toCopy.fKYsystem) fKYsystem = new TKinematics2to2(*toCopy.fKYsystem);
  if(toCopy.fKNsystem) fKNsystem = new TKinematics2to2(*toCopy.fKNsystem);
  if(toCopy.fYNsystem) fYNsystem = new TKinematics2to2(*toCopy.fYNsystem);

  if(toCopy.fIsVar) {  
    fIsVar = new bool[6];
    fLowerLimit = new double[6];
    fUpperLimit = new double[6];
    fStepSize = new double[6];
    fNumberOfSteps = new int[6];
  
    for(int i=0; i<6; ++i) {
      fIsVar[i] = toCopy.fIsVar[i];
      fLowerLimit[i] = toCopy.fLowerLimit[i];
      fUpperLimit[i] = toCopy.fUpperLimit[i];
      fStepSize[i] = toCopy.fStepSize[i];
      fNumberOfSteps[i] = toCopy.fNumberOfSteps[i];
    }
  }
}

//_____________________________________________________________________
TKinematics2to3::TKinematics2to3(TRootIOCtor *rio)
  : TNamed(), fKYsystem(0), fKNsystem(0), fYNsystem(0),
    fKvector(0), fYvector(0), fNvector(0), fGvector(0), fDvector(0),
    fIsospin(0), fMd(0), fMk(0), fMy(0), fMn(0), fIsPhysical(false),
    fQsquared(0), fWlab(0), fKlab(0), fEk(0), fEy(0), fEn(0),
    fPk(0), fPy(0), fPn(0), fThk(0), fThy(0), fThn(0), fCosthk(0),
    fCosthy(0), fCosthn(0), fPhik(0), fPhiy(0), fPhin(0), fWtot(0),
    fWky(0), fWkn(0), fWyn(0), fStot(0), fSky(0), fSkn(0), fSyn(0), fTgk(0),
    fTgy(0), fTgn(0), fKYcosthkcm(0), fKYcosthycm(0), fKNcosthkcm(0),
    fKNcosthncm(0), fYNcosthycm(0), fYNcosthncm(0),
    fVar(0), fE1Var(0), fTh1Var(0), fE2Var(0), fTh2Var(0), fPhi2Var(0),
    fTdk(0), fTdy(0), fTdn(0), fKYphikcm(0), fKYphiycm(0), fKNphikcm(0),
    fKNphincm(0), fYNphiycm(0), fYNphincm(0), fKYthkcm(0), fKYthycm(0),
    fKNthkcm(0), fKNthncm(0), fYNthycm(0), fYNthncm(0),
    fEVar(0), fIsVar(0), fLowerLimit(0), fUpperLimit(0),
    fStepSize(0), fNumberOfSteps(0), fStep(0), fFormat(0),
    fAC(kN), fKinSolution(0), fNrOfPhysicals(-1)
{
  // ROOT I/O constructor
}

//_____________________________________________________________________
TKinematics2to3::~TKinematics2to3()
{
  // Destructor

  delete[] fFormat;
  delete[] fKinSolution;
  delete[] fVar;
  delete[] fIsVar;
  delete[] fLowerLimit;
  delete[] fUpperLimit;
  delete[] fStepSize;
  delete[] fNumberOfSteps;

  delete fKYsystem;
  delete fKNsystem;
  delete fYNsystem;
  delete fKvector;
  delete fYvector;
  delete fNvector;
  delete fGvector;
  delete fDvector;
}

//_____________________________________________________________________
TKinematics2to3& TKinematics2to3::operator=(const TKinematics2to3& toCopy)
{
  // Assignment

  if( this!=&toCopy ) { // avoid self-assignment
    
    TNamed::operator=(toCopy);
    
    ChangeIsospinChannel(toCopy.fIsospin);
    fAC = toCopy.fAC;
    ReadFormatString(toCopy.fFormat);

    fIsPhysical = toCopy.fIsPhysical;
    fQsquared = toCopy.fQsquared;
    fWlab = toCopy.fWlab;
    fKlab = toCopy.fKlab;
    fEk = toCopy.fEk;
    fEy = toCopy.fEy;
    fEn = toCopy.fEn;
    fPk = toCopy.fPk;
    fPy = toCopy.fPy;
    fPn = toCopy.fPn;
    fThk = toCopy.fThk;
    fThy = toCopy.fThy;
    fThn = toCopy.fThn;
    fCosthk = toCopy.fCosthk;
    fCosthy = toCopy.fCosthy;
    fCosthn = toCopy.fCosthn;
    fPhik = toCopy.fPhik;
    fPhiy = toCopy.fPhiy;
    fPhin = toCopy.fPhin;
    fWtot = toCopy.fWtot;
    fWky = toCopy.fWky;
    fWkn = toCopy.fWkn;
    fWyn = toCopy.fWyn;
    fStot = toCopy.fStot;
    fSky = toCopy.fSky;
    fSkn = toCopy.fSkn;
    fSyn = toCopy.fSyn;
    fTgk = toCopy.fTgk;
    fTgy = toCopy.fTgy;
    fTgn = toCopy.fTgn;
    fTdk = toCopy.fTdk;
    fTdy = toCopy.fTdy;
    fTdn = toCopy.fTdn;
    fKYthkcm = toCopy.fKYthkcm;
    fKYthycm = toCopy.fKYthycm;
    fKNthkcm = toCopy.fKNthkcm;
    fKNthncm = toCopy.fKNthncm;
    fYNthycm = toCopy.fYNthycm;
    fYNthncm = toCopy.fYNthncm;
    fKYcosthkcm = toCopy.fKYcosthkcm;
    fKYcosthycm = toCopy.fKYcosthycm;
    fKNcosthkcm = toCopy.fKNcosthkcm;
    fKNcosthncm = toCopy.fKNcosthncm;
    fYNcosthycm = toCopy.fYNcosthycm;
    fYNcosthncm = toCopy.fYNcosthncm;
    fKYphikcm = toCopy.fKYphikcm;
    fKYphiycm = toCopy.fKYphiycm;
    fKNphikcm = toCopy.fKNphikcm;
    fKNphincm = toCopy.fKNphincm;
    fYNphiycm = toCopy.fYNphiycm;
    fYNphincm = toCopy.fYNphincm;

    *fKYsystem = *toCopy.fKYsystem;
    *fKNsystem = *toCopy.fKNsystem;
    *fYNsystem = *toCopy.fYNsystem;

    *fKvector = *toCopy.fKvector;
    *fYvector = *toCopy.fYvector;
    *fNvector = *toCopy.fNvector;
    *fGvector = *toCopy.fGvector;
    *fDvector = *toCopy.fDvector;

    fStep = toCopy.fStep;
    fNrOfPhysicals = toCopy.fNrOfPhysicals;

    if(toCopy.fIsVar) {
      if( !fIsVar) fIsVar = new bool[6];
      if( !fLowerLimit) fLowerLimit = new double[6];
      if( !fUpperLimit) fUpperLimit = new double[6];
      if( !fStepSize) fStepSize = new double[6];
      if( !fNumberOfSteps) fNumberOfSteps = new int[6];
      
      for(int i=0; i<6; ++i) {
	fIsVar[i] = toCopy.fIsVar[i];
	fLowerLimit[i] = toCopy.fLowerLimit[i];
	fUpperLimit[i] = toCopy.fUpperLimit[i];
	fStepSize[i] = toCopy.fStepSize[i];
	fNumberOfSteps[i] = toCopy.fNumberOfSteps[i];
      }
    } else {
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
void TKinematics2to3::SetIsospin(int isospin)
{
  // Change isospin channel. Kinematics will immediately be updated.

  ChangeIsospinChannel(isospin);
  UpdateKinematics();
}

//_____________________________________________________________________
void TKinematics2to3::ChangeIsospinChannel(int isospin)
{
  // Set the particle masses according to the isospin channel

  /*   (1)   ===>   g + d   ->  K+ + L0 + n
   *   (2)   ===>   g + d   ->  K+ + S0 + n
   *   (3)   ===>   g + d   ->  K0 + S+ + n
   *   (4)   ===>   g + d   ->  K0 + L0 + p
   *   (5)   ===>   g + d   ->  K0 + S0 + p
   *   (6)   ===>   g + d   ->  K+ + S- + p
   */

  fIsospin = isospin;

  //fMd = M_D;  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  switch( isospin ) 
    {
    case 1:
      fMd = MASSD;
      fMn = M_N;
      fMk = M_KP;
      fMy = M_L;
      break;
    case 2:
      fMd = MASSD;
      fMn = M_N;
      fMk = M_KP;
      fMy = M_S0;
      break;
    case 3:
      fMd = MASSD;
      fMn = M_N;
      fMk = M_K0;
      fMy = M_SP;
      break;
    case 4:
      fMd = MASSD;
      fMn = M_P;
      fMk = M_K0;
      fMy = M_L;
      break;
    case 5:
      fMd = MASSD;
      fMn = M_P;
      fMk = M_K0;
      fMy = M_S0;
      break;
    case 6:
      fMd = MASSD;
      fMn = M_P;
      fMk = M_KP;
      fMy = M_SM;
      break;
    case 7:
      break;
    default:
      std::cerr << "ERROR in TKinematics2to3::ChangeIsospinChannel(int):\n"
		<< "Isospin channel " << isospin << " is not implemented!\n";
      break;
    }  

  // if the subsystems have been created, update their masses
  if(fKYsystem) {   
    fKYsystem->SetMasses(fMd,fMk,fMy);
    fKNsystem->SetMasses(fMd,fMk,fMn);
    fYNsystem->SetMasses(fMd,fMn,fMy);
  }
  
  // Reset the number of physical steps
  fNrOfPhysicals = -1;
}

//_____________________________________________________________________
void TKinematics2to3::SetFormat(const char* format)
{
  /* Begin_Html
     Set the 6 independent variables by providing a ':' separated list of variables.
  
     All variables are in the LAB frame, except when their name explicitly
     suggests differently.
  
     Look <a href="#indepvar">here</a> to find out what your format string should
     look like.
     <p>
     Below is a list of the 'code names' of variables:
     <table>
     <tr><td>Q<sup>2</sup></td><td><tt>qsquared</tt></td></tr>
     <tr><td>photon energy</sub></td><td><tt>wlab</tt></td></tr>
     <tr><td>photon momentum</td><td><tt>klab</tt></td></tr>
     <tr><td>kaon energy</td><td><tt>ek</tt></td></tr>
     <tr><td>hyperon energy</td><td><tt>ey</tt></td></tr>
     <tr><td>nucleon energy</td><td><tt>en</tt></td></tr>
     <tr><td>kaon momentum</td><td><tt>pk</tt></td></tr>
     <tr><td>hyperon momentum</td><td><tt>py</tt></td></tr>
     <tr><td>nucleon momentum</td><td><tt>pn</tt></td></tr>
     <tr><td>kaon theta</td><td><tt>thk</tt></td></tr>
     <tr><td>hyperon theta</td><td><tt>tky</tt></td></tr>
     <tr><td>nucleon theta</td><td><tt>thn</tt></td></tr>
     <tr><td>cos(theta<sub>K</sub>)</td><td><tt>costhk</tt></td></tr>
     <tr><td>cos(theta<sub>Y</sub>)</td><td><tt>costhy</tt></td></tr>
     <tr><td>cos(theta<sub>N</sub>)</td><td><tt>costhn</tt></td></tr>
     <tr><td>kaon phi</td><td><tt>phik</tt></td></tr>
     <tr><td>hyperon phi</td><td><tt>phiy</tt></td></tr>
     <tr><td>nucleon phi</td><td><tt>phin</tt></td></tr>
     <tr><td>total invariant mass</td><td><tt>wtot</tt></td></tr>
     <tr><td>invariant mass KY-system</td><td><tt>wky</tt></td></tr>
     <tr><td>invariant mass KN-system</td><td><tt>wkn</tt></td></tr>
     <tr><td>invariant mass YN-system</td><td><tt>wyn</tt></td></tr>
     <tr><td>total s</td><td><tt>stot</tt></td></tr>
     <tr><td>s of KY-system</td><td><tt>sky</tt></td></tr>
     <tr><td>s of KN-system</td><td><tt>skn</tt></td></tr>
     <tr><td>s of YN-system</td><td><tt>syn</tt></td></tr>
     <tr><td>Mom.transf. photon-kaon</td><td><tt>tgk</tt></td></tr>
     <tr><td>Mom.transf. photon-hyperon</td><td><tt>tgy</tt></td></tr>
     <tr><td>Mom.transf. photon-nucleon</td><td><tt>tgn</tt></td></tr>
     <tr><td>Mom.transf. deut.-kaon</td><td><tt>tdk</tt></td></tr>
     <tr><td>Mom.transf. deut.-hyperon</td><td><tt>tdy</tt></td></tr>
     <tr><td>Mom.transf. deut.-nucleon</td><td><tt>tdn</tt></td></tr>
     <tr><td>kaon theta in KY-system</td><td><tt>kythkcm</tt></td></tr>
     <tr><td>hyperon theta in KY-system</td><td><tt>kythycm</tt></td></tr>
     <tr><td>kaon theta in KN-system</td><td><tt>knthkcm</tt></td></tr>
     <tr><td>nucleon theta in KN-system</td><td><tt>knthncm</tt></td></tr>
     <tr><td>hyperon theta in YN-system</td><td><tt>ynthycm</tt></td></tr>
     <tr><td>nucleon theta in YN-system</td><td><tt>ynthncm</tt></td></tr>
     <tr><td>cos(<tt>kythkcm</tt>)</td><td><tt>kycosthkcm</tt></td></tr>
     <tr><td>cos(<tt>kythycm</tt>)</td><td><tt>kycosthycm</tt></td></tr>
     <tr><td>cos(<tt>knthkcm</tt>)</td><td><tt>kncosthkcm</tt></td></tr>
     <tr><td>cos(<tt>knthncm</tt>)</td><td><tt>kncosthncm</tt></td></tr>
     <tr><td>cos(<tt>ynthycm</tt>)</td><td><tt>yncosthycm</tt></td></tr>
     <tr><td>cos(<tt>ynthncm</tt>)</td><td><tt>yncosthncm</tt></td></tr>
     <tr><td>kaon phi in KY-system</td><td><tt>kyphikcm</tt></td></tr>
     <tr><td>hyperon phi in KY-system</td><td><tt>kyphiycm</tt></td></tr>
     <tr><td>kaon phi in KN-system</td><td><tt>knphikcm</tt></td></tr>
     <tr><td>nucleon phi in KN-system</td><td><tt>knphincm</tt></td></tr>
     <tr><td>hyperon phi in YN-system</td><td><tt>ynphiycm</tt></td></tr>
     <tr><td>nucleon phi in YN-system</td><td><tt>ynphincm</tt></td></tr>
     </table>
     End_Html */ 

  if(!IsFixed()) {
    TMPI::Cout() << "INFO in TKinematics2to3::SetFormat(const char*): "
		 << "Fixing all variables.\n";
    for(int i=0; i<6; ++i)
      FixVariable(i+1);
  }

  UpdateVariables();

  ReadFormatString(format);
  
  // Finally, set the upper and lower limit
  for(int i=0; i<6; ++i) {
    fLowerLimit[i] = *fVar[i];
    fUpperLimit[i] = *fVar[i];
  }
}

//_____________________________________________________________________
void TKinematics2to3::SetAxisConvention(EAxisConvention ac)
{
  // Change the axis convention. Be sure you know what you're doing before
  // you call this function. This function doesn't check whether what 
  // you are doing makes sense.
  fAC = ac;
}

//_____________________________________________________________________
void TKinematics2to3::ReadFormatString(const char *format)
{
  // Process the format string.

  // Store format in TKinematics2to3::fFormat
  // Read TKinematics2to3::fFormat and let fVar[] point
  // to the 6 independent kinematical points as specified
  // by the user.
  
  // Check that the chosen variables make sense.
  // The user should always provide Q^2 of the photon

  // We identify 3 distinct scenarios:
  // 1) 1+(23)
  // 2) 1+2
  // 3) mandelstam
  // Given the format string this function will determine what
  // scenario is adequate to solve the kinematics. Some specifics
  // will be stored in fKinSolution to assist UpdateKinematics().
  std::strcpy(fKinSolution,"---");

  // Some more info on the different scenarios:
  // 1) The user provides info about the energy (fE1Var)
  //    and the theta-angle (fTh1Var) of the particle that
  //    fixes the x-axis (fAC).
  //    One variable should either specify the photon energy
  //    or the invariant mass of the total system or the
  //    (2+3) subsystem (fEVar).
  //    In addition angular info needs to be given about a
  //    second final particle. Both theta (fTh2Var) and phi (fPhi2Var)
  //    should be either in the cm of the (2+3) subsystem or in the
  //    lab frame. This info will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "1YL": scenario 1, 2nd particle = hyperon with angles in LAB
  //    "1K*": scenario 1, 2nd particle = kaon with angles in CM
  //
  // 2) The user provides info about the energy (fE1Var)
  //    and the theta-angle (fTh1Var) of the particle that
  //    fixes the x-axis (fAC).
  //    In addition a second final particle is fully specified,
  //    energy, theta and phi (fE2Var,fTh2Var,fPhi2Var) in the lab frame.
  //    The name of the 2nd particle will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "2N-": scenario 2, 2nd particle = nucleon
  //
  // 3) All variables are mandelstam variables, except for Q^2
  //    obviously. 
  //    S_tot = (P_gamma + P_deut)^2 should always be given.
  //    fAC fixes the 1st particle. This leaves 7 combinations:
  //    (3a) s_12 + s_23 + t_g1 + t_d3
  //    (3b) s_12 + s_23 + t_g1 + t_d2
  //    (3c) s_12 + s_23 + t_d2 + t_d3
  //    (3d) s_12 + t_d2 + t_g1 + t_d3
  //    (3e) s_12 + s_23 + t_g1 + t_g2
  //    (3f) s_12 + s_23 + t_g2 + t_d3
  //    (3g) t_g2 + s_23 + t_g1 + t_d3
  //    The name of the 2nd particle will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "3Kb": scenario 3b with 2nd particle = kaon

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!                                                      !!
  // !! ATTENTION                                            !!
  // !! some cases are not valid anymore, until I find a way !!
  // !! to solve them without ambiguities.                   !!
  // !!                                                      !!
  // !! Not valid: Case 1, with angles in the LAB            !!
  // !!            Case 2                                    !!
  // !!                                                      !!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // Set all pointers to independent variables to NULL                 
  fE1Var = fE2Var = fEVar = fTh1Var = fTh2Var = fPhi2Var = 0;

  // Put format string in member and in temporary TString.
  std::strcpy(fFormat,format);
  TString formatString = format;

  // Tokenize format string
  TObjArray *tokens = formatString.Tokenize(":");
  char token[11];

  // Check the number of tokens
  int nrOfTokens = tokens->GetEntries();
  if( nrOfTokens != 6 ) {
    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	 << "Format string should provide 6 independent variables: "
	 << "\"var1:var2:var3:var4:var5:var6\"\n";
    exit(1);
  }

  // Transform tokens to lower case
  for(int i=0; i<nrOfTokens; ++i)
    ((TObjString*)tokens->At(i))->String().ToLower();

  // Check that Q^2 is present
  if( !formatString.Contains("qsquared",TString::kIgnoreCase) ) {
    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	 << "Q^2 needs to be one of 6 independent kinematic variables.\n";
    exit(1);
  }

  // Check whether all tokens are mandelstam variables. We do this by
  // counting the s's and t's
  TString *mandel = 0;
  if( formatString.CountChar('s')+formatString.CountChar('S')-1
      +formatString.CountChar('t')+formatString.CountChar('T')-2 ==5 
      && formatString.Contains("stot",TString::kIgnoreCase) ) {

    fKinSolution[0] = '3';

    mandel = new TString("abcdefg"); // list all subscenarios
    
    // find out name of '2nd particle'
    if( fAC == kK ) {
      if( formatString.Contains("sky",TString::kIgnoreCase) ||
	  formatString.Contains("syk",TString::kIgnoreCase) ||
	  formatString.Contains("tgy",TString::kIgnoreCase) ||
	  formatString.Contains("tyg",TString::kIgnoreCase) )
	fKinSolution[1] = 'Y';
      else
	fKinSolution[1] = 'N';
    } // fAC == kK
    else if( fAC == kY ) {
      if( formatString.Contains("sky",TString::kIgnoreCase) ||
	  formatString.Contains("syk",TString::kIgnoreCase) ||
	  formatString.Contains("tgk",TString::kIgnoreCase) ||
	  formatString.Contains("tkg",TString::kIgnoreCase) )
	fKinSolution[1] = 'K';
      else
	fKinSolution[1] = 'N';
    } // fAC == kY
    else if( fAC == kN ) {
      if( formatString.Contains("sny",TString::kIgnoreCase) ||
	  formatString.Contains("syn",TString::kIgnoreCase) ||
	  formatString.Contains("tgy",TString::kIgnoreCase) ||
	  formatString.Contains("tyg",TString::kIgnoreCase) )
	fKinSolution[1] = 'Y';
      else
	fKinSolution[1] = 'K';
    } // fAC == kN
  } // end mandelstam case

  // Determine independent kinematic variables
  for(int i=0; i<nrOfTokens; ++i) 
    {
      std::strcpy(token,((TObjString*)tokens->At(i))->String());

      //__________________________________
      if( !strcmp(token,"qsquared") ) {
	fVar[i] = &fQsquared;
      }
      //__________________________________
      else if( !strcmp(token,"wlab") ) {
	if(fEVar || fKinSolution[0] == '2') {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	fVar[i] = fEVar = &fWlab;
	fKinSolution[0] = '1';
      }
      //__________________________________
      else if( !strcmp(token,"klab") ) {
	if(fEVar || fKinSolution[0] == '2') {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	fVar[i] = fEVar = &fKlab;
	fKinSolution[0] = '1';
      }
      //__________________________________
      else if( !strcmp(token,"ek") ) {
	if( fAC == kK ) {
	  if(fE1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE1Var = &fEk;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE2Var = &fEk;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'K';
	}
      }
      //__________________________________
      else if( !strcmp(token,"ey") ) {
	if( fAC == kY ) {
	  if(fE1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE1Var = &fEy;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE2Var = &fEy;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'Y';
	}	
      }
      //__________________________________
      else if( !strcmp(token,"en") ) {
	if( fAC == kN ) {
	  if(fE1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE1Var = &fEn;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE2Var = &fEn;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'N';
	}	
      }
      //__________________________________
      else if( !strcmp(token,"pk") ) {
	if( fAC == kK ) {
	  if(fE1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE1Var = &fPk;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE2Var = &fPk;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'K';
	}
      }
      //__________________________________
      else if( !strcmp(token,"py") ) {
	if( fAC == kY ) {
	  if(fE1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE1Var = &fPy;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE2Var = &fPy;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'Y';
	}	
      }
      //__________________________________
      else if( !strcmp(token,"pn") ) {
	if( fAC == kN ) {
	  if(fE1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE1Var = &fPn;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fE2Var = &fEn;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'P';
	}	
      }
      //__________________________________
      else if( !strcmp(token,"thk") ) {
	if( fAC == kK ) {
	  if(fTh1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh1Var = &fThk;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fThk;
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( !strcmp(token,"thy") ) {
	if( fAC == kY ) {
	  if(fTh1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh1Var = &fThy;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fThy;
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( !strcmp(token,"thn") ) {
	if( fAC == kN ) {
	  if(fTh1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh1Var = &fThn;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fThn;
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( !strcmp(token,"costhk") ) {
	if( fAC == kK ) {
	  if(fTh1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh1Var = &fCosthk;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fCosthk;
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( !strcmp(token,"costhy") ) {
	if( fAC == kY ) {
	  if(fTh1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh1Var = &fCosthy;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fCosthy;
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( !strcmp(token,"costhn") ) {
	if( fAC == kN ) {
	  if(fTh1Var) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh1Var = &fCosthn;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fCosthn;
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( !strcmp(token,"phik") ) {
	if( fAC == kK ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane! Because of the axis convention "
	       << "phik is zero by definition.\n";
	  exit(1);
	}
	else {
	  if(fPhi2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fPhi2Var = &fPhik;
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( !strcmp(token,"phiy") ) {
	if( fAC == kY ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane! Because of the axis convention "
	       << "phiy is zero by definition.\n";
	  exit(1);
	}
	else {
	  if(fPhi2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fPhi2Var = &fPhiy;
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( !strcmp(token,"phin") ) {
	if( fAC == kN ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane! Because of the axis convention "
	       << "phin is zero by definition.\n";
	  exit(1);
	}
	else {
	  if(fPhi2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fPhi2Var = &fPhin;
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( !strcmp(token,"wtot") ) {
	if(fEVar || fKinSolution[0] == '2') {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	fVar[i] = fEVar = &fWtot;
	fKinSolution[0] = '1';
      }
      //__________________________________
      else if( !strcmp(token,"wky") || !strcmp(token,"wyk") ) {
	if( fAC != kN ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fEVar || fKinSolution[0] == '2') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fEVar = &fWky;
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( !strcmp(token,"wkn") || !strcmp(token,"wnk") ) {
	if( fAC != kY ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fEVar || fKinSolution[0] == '2') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fEVar = &fWkn;
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( !strcmp(token,"wyn") || !strcmp(token,"wny") ) {
	if( fAC != kK ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fEVar || fKinSolution[0] == '2') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fEVar = &fWyn;
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( !strcmp(token,"stot") ) {
	if(fEVar || fKinSolution[0] == '2') {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else if(fKinSolution[0] == '3') {
	  fVar[i] = &fStot;
	}
	else {
	  fVar[i] = fEVar = &fStot;
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( !strcmp(token,"sky") || !strcmp(token,"syk") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kK && fKinSolution[1]=='Y') ||
	      (fAC==kY && fKinSolution[1]=='K') ) {
	    mandel->ReplaceAll("g","-");
	  }
	  else if( fAC==kN ) {
	    mandel->ReplaceAll("d","-");
	  }
	  else {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = &fSky;
	} // end mandelstam case
	else {
	  if( fAC != kN ) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  else {
	    if(fEVar || fKinSolution[0] == '2') {
	      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!\n";
	      exit(1);
	    }
	    fVar[i] = fEVar = &fSky;
	    fKinSolution[0] = '1';
	  }
	}
      }
      //__________________________________
      else if( !strcmp(token,"skn") || !strcmp(token,"snk") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kK && fKinSolution[1]=='N') ||
	      (fAC==kN && fKinSolution[1]=='K') ) {
	    mandel->ReplaceAll("g","-");
	  }
	  else if( fAC==kY ) {
	    mandel->ReplaceAll("d","-");
	  }
	  else {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = &fSkn;
	} // end mandelstam case
	else {
	  if( fAC != kY ) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  else {
	    if(fEVar || fKinSolution[0] == '2') {
	      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!\n";
	      exit(1);
	    }
	    fVar[i] = fEVar = &fSkn;
	    fKinSolution[0] = '1';
	  }
	}
      }
      //__________________________________
      else if( !strcmp(token,"syn") || !strcmp(token,"sny") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kY && fKinSolution[1]=='N') ||
	      (fAC==kN && fKinSolution[1]=='Y') ) {
	    mandel->ReplaceAll("g","-");
	  }
	  else if( fAC==kK ) {
	    mandel->ReplaceAll("d","-");
	  }
	  else {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = &fSyn;
	} // end mandelstam case
	else {
	  if( fAC != kK ) {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  else {
	    if(fEVar || fKinSolution[0] == '2') {
	      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!\n";
	      exit(1);
	    }
	    fVar[i] = fEVar = &fSyn;
	    fKinSolution[0] = '1';
	  }
	}
      }
      //__________________________________
      else if( !strcmp(token,"tgk") || !strcmp(token,"tkg") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kY || fAC==kN ) && fKinSolution[1]=='K' ) {
	    mandel->ReplaceAll("a","-");
	    mandel->ReplaceAll("b","-");
	    mandel->ReplaceAll("c","-");
	    mandel->ReplaceAll("d","-");
	  }
	  else if( fAC==kK ) {
	    mandel->ReplaceAll("c","-");
	    mandel->ReplaceAll("f","-");
	  }
	  else {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = &fTgk;
	} // end mandelstam case
	else if (fKinSolution[0] == '2') {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if( fAC == kK ) {
	    if(fTh1Var) {
	      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!\n";
	      exit(1);
	    }
	    fVar[i] = fTh1Var = &fTgk;
	  }
	  else {
	    if(fTh2Var || fKinSolution[2] == 'L' 
	       || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!\n";
	      exit(1);
	    }
	    fVar[i] = fTh2Var = &fTgk;
	    fKinSolution[1] = 'K';
	    fKinSolution[2] = '*';
	  }
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( !strcmp(token,"tgy") || !strcmp(token,"tyg") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kK || fAC==kN ) && fKinSolution[1]=='Y' ) {
	    mandel->ReplaceAll("a","-");
	    mandel->ReplaceAll("b","-");
	    mandel->ReplaceAll("c","-");
	    mandel->ReplaceAll("d","-");
	  }
	  else if( fAC==kY ) {
	    mandel->ReplaceAll("c","-");
	    mandel->ReplaceAll("f","-");
	  }
	  else {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = &fTgy;
	} // end mandelstam case
	else if (fKinSolution[0] == '2') {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if( fAC == kY ) {
	    if(fTh1Var) {
	      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!\n";
	      exit(1);
	    }
	    fVar[i] = fTh1Var = &fTgy;
	  }
	  else {
	    if(fTh2Var || fKinSolution[2] == 'L' 
	       || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!\n";
	      exit(1);
	    }
	    fVar[i] = fTh2Var = &fTgy;
	    fKinSolution[1] = 'Y';
	    fKinSolution[2] = '*';
	  }
	  fKinSolution[0] = '1'; 
	}
      }
      //__________________________________
      else if( !strcmp(token,"tgn") || !strcmp(token,"tng") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kK || fAC==kY ) && fKinSolution[1]=='N' ) {
	    mandel->ReplaceAll("a","-");
	    mandel->ReplaceAll("b","-");
	    mandel->ReplaceAll("c","-");
	    mandel->ReplaceAll("d","-");
	  }
	  else if( fAC==kN ) {
	    mandel->ReplaceAll("c","-");
	    mandel->ReplaceAll("f","-");
	  }
	  else {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = &fTgn;
	} // end mandelstam case
	else if (fKinSolution[0] == '2') {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if( fAC == kN ) {
	    if(fTh1Var) {
	      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!\n";
	      exit(1);
	    }
	    fVar[i] = fTh1Var = &fTgn;
	  }
	  else {
	    if(fTh2Var || fKinSolution[2] == 'L' 
	       || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!\n";
	      exit(1);
	    }
	    fVar[i] = fTh2Var = &fTgn;
	    fKinSolution[1] = 'N';
	    fKinSolution[2] = '*';
	  }
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( !strcmp(token,"tdk") || !strcmp(token,"tkd") ) {
	if( fKinSolution[0] != '3' || fAC==kK ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else if( fKinSolution[1] == 'K' ) {
	  mandel->ReplaceAll("a","-");
	  mandel->ReplaceAll("e","-");
	  mandel->ReplaceAll("f","-");
	  mandel->ReplaceAll("g","-");
	}
	else {
	  mandel->ReplaceAll("b","-");
	  mandel->ReplaceAll("e","-");
	}
	fVar[i] = &fTdk;
      }
      //__________________________________
      else if( !strcmp(token,"tdy") || !strcmp(token,"tyd") ) {
	if( fKinSolution[0] != '3' || fAC==kY ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else if( fKinSolution[1] == 'Y' ) {
	  mandel->ReplaceAll("a","-");
	  mandel->ReplaceAll("e","-");
	  mandel->ReplaceAll("f","-");
	  mandel->ReplaceAll("g","-");
	}
	else {
	  mandel->ReplaceAll("b","-");
	  mandel->ReplaceAll("e","-");
	}
	fVar[i] = &fTdy;
      }
      //__________________________________
      else if( !strcmp(token,"tdn") || !strcmp(token,"tnd") ) {
	if( fKinSolution[0] != '3' || fAC==kN ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else if( fKinSolution[1] == 'N' ) {
	  mandel->ReplaceAll("a","-");
	  mandel->ReplaceAll("e","-");
	  mandel->ReplaceAll("f","-");
	  mandel->ReplaceAll("g","-");
	}
	else {
	  mandel->ReplaceAll("b","-");
	  mandel->ReplaceAll("e","-");
	}
	fVar[i] = &fTdn;
      }
      //__________________________________
      else if( !strcmp(token,"kythkcm") || !strcmp(token,"ykthkcm") ) {
	if( fAC != kN ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fKYthkcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"kythycm") || !strcmp(token,"ykthycm") ) {
	if( fAC != kN ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fKYthycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"knthkcm") || !strcmp(token,"nkthkcm") ) {
	if( fAC != kY ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fKNthkcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"knthncm") || !strcmp(token,"nkthncm") ) {
	if( fAC != kY ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fKNthncm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"ynthycm") || !strcmp(token,"nythycm") ) {
	if( fAC != kK ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fYNthycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"ynthncm") || !strcmp(token,"nythncm") ) {
	if( fAC != kK ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fYNthncm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"kycosthkcm") || !strcmp(token,"ykcosthkcm") ) {
	if( fAC != kN ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fKYcosthkcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"kycosthycm") || !strcmp(token,"ykcosthycm") ) {
	if( fAC != kN ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fKYcosthycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"kncosthkcm") || !strcmp(token,"nkcosthkcm") ) {
	if( fAC != kY ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fKNcosthkcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"kncosthncm") || !strcmp(token,"nkcosthncm") ) {
	if( fAC != kY ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fKNcosthncm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"yncosthycm") || !strcmp(token,"nycosthycm") ) {
	if( fAC != kK ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fYNcosthycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"yncosthncm") || !strcmp(token,"nycosthncm") ) {
	if( fAC != kK ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fTh2Var = &fYNcosthncm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"kyphikcm") || !strcmp(token,"ykphikcm") ) {
	if( fAC != kN ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fPhi2Var = &fKYphikcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"kyphiycm") || !strcmp(token,"ykphiycm") ) {
	if( fAC != kN ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fPhi2Var = &fKYphiycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"knphikcm") || !strcmp(token,"nkphikcm") ) {
	if( fAC != kY ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fPhi2Var = &fKNphikcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"knphincm") || !strcmp(token,"nkphincm") ) {
	if( fAC != kY ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fPhi2Var = &fKNphincm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"ynphiycm") || !strcmp(token,"nyphiycm") ) {
	if( fAC != kK ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fPhi2Var = &fYNphiycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( !strcmp(token,"ynphincm") || !strcmp(token,"nyphincm") ) {
	if( fAC != kK ) {
	  cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!\n";
	  exit(1);
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!\n";
	    exit(1);
	  }
	  fVar[i] = fPhi2Var = &fYNphincm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }

      //__________________________________
      else {
	fVar[i] = NULL;
	cerr << "ERROR in  TKinematics2to3::ReadFormatString(const char*): "
	     << "Kinematic variable \"" << token << "\" unknown!\n";
	exit(1);
      }
    } // end loop over tokens
  
  // Check whether the user provided a meaningful format string
  // in the mandelstam case
  if(mandel) {
    if(mandel->CountChar('-') != 6) {
      cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	   << "Format string is not sane!\n";
      exit(1);
    }
    fKinSolution[2] = mandel->Strip(TString::kBoth,'-')[0];
  }

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!                                                      !!
  // !! ATTENTION                                            !!
  // !! some cases are not valid anymore, until I find a way !!
  // !! to solve them without ambiguities.                   !!
  // !!                                                      !!
  // !! Not valid: Case 1, with angles in the LAB            !!
  // !!            Case 2                                    !!
  // !!                                                      !!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if( fKinSolution[0] == '2' ) {
    cerr << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	 << "I cannot solve kinematics with these variables.\n";
    exit(1);
  }
  else if( fKinSolution[0] == '1' && fKinSolution[2] == 'L' )
    cerr << "WARNING in TKinematics2to3::ReadFormatString(const char*): "
	 << "CAVEAT: risk of running into multiple solutions.\n"
	 << "Perhaps TKinematics2to3WithLabAngles is something up your alley.\n";

  // delete dynamics memory
  delete tokens; delete mandel;
}

//_____________________________________________________________________
void TKinematics2to3::InitializeSubsystems()
{
  // Create and set the KY-, KN- and YN-subsystems.

  if(!fKYsystem) fKYsystem = new TKinematics2to2("fKYsystem","KY CM subsystem",
						 fMd,fMk,fMy,"s:t:qsquared",
						 0.0,0.0,0.0);
  if(!fKNsystem) fKNsystem = new TKinematics2to2("fKNsystem","KN CM subsystem",
						 fMd,fMk,fMn,"s:t:qsquared",
						 0.0,0.0,0.0);
  if(!fYNsystem) fYNsystem = new TKinematics2to2("fYNsystem","YN CM subsystem",
						 fMd,fMn,fMy,"s:t:qsquared",
						 0.0,0.0,0.0);
}

//_____________________________________________________________________
void TKinematics2to3::UpdateKinematics()
{
  // Solve the kinematics!

  fIsPhysical = true;

  // declare 4vectors in the fAC axis convention.
  TLorentzVector vector1,vector2,vector3;
  TLorentzVector vector23;
  
  // pointer to subsystem
  TKinematics2to2 *subSystem;
  if( fAC==kK ) subSystem = fYNsystem;
  else if( fAC==kY ) subSystem = fKNsystem;
  else subSystem = fKYsystem;

  // declare final state particle masses
  double m1,m2,m3;
  if( fAC==kK ) { m1 = fMk;
    if( fKinSolution[1]=='Y' ) {
      m2 = fMy;
      m3 = fMn;
    } else {
      m2 = fMn;
      m3 = fMy;
    }
  } else if( fAC==kY ) { m1 = fMy;
    if( fKinSolution[1]=='K' ) {
      m2 = fMk;
      m3 = fMn;
    } else {
      m2 = fMn;
      m3 = fMk;
    }
  } else { m1 = fMn;
    if( fKinSolution[1]=='K' ) {
      m2 = fMk;
      m3 = fMy;
    } else {
      m2 = fMy;
      m3 = fMk;
    }
  }

  fDvector->SetXYZM(0.,0.,0.,fMd); // trivial

  // 1) The user provides info about the energy (fE1Var)
  //    and the theta-angle (fTh1Var) of the particle that
  //    fixes the x-axis (fAC).
  //    One variable should either specify the photon energy
  //    or the invariant mass of the total system or the
  //    (2+3) subsystem (fEVar).
  //    In addition angular info needs to be given about a
  //    second final particle. Both theta (fTh2Var) and phi (fPhi2Var)
  //    should be either in the cm of the (2+3) subsystem or in the
  //    lab frame. This info will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "1YL": scenario 1, 2nd particle = hyperon with angles in LAB
  //    "1K*": scenario 1, 2nd particle = kaon with angles in CM
  //
  //    The first option with the angles in the LAB system can give
  //    multiple solutions. We advise the user not to use this option
  //    in ReadFormatString(const char*).
  //    TKinematics2to3WitLabAngles can be used for these cases.
  
  if( fKinSolution[0] == '1' ) {

    // Determine 4vector of photon & 1st particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Energy 1st particle
    double e1,p1;
    if( fE1Var==&fPk || fE1Var==&fPy || fE1Var==&fPn ) {
      p1 = *fE1Var;
      e1 = Sqrt( p1*p1 + m1*m1 );
    } else {
      e1 = *fE1Var;
      p1 = Sqrt( e1*e1 - m1*m1 );
    }
    if(e1<0. || p1 <0.) fIsPhysical = false;

    // Determine photon energy
    //
    // Special attention is needed when energy is given as s23 or w23
    double eg=-1, s23=-1;
    bool egCalculated = true;
    if( fEVar==&fWlab )
      eg = fWlab;
    else if( fEVar==&fKlab )
      eg = Sqrt( fKlab*fKlab - fQsquared );
    else if( fEVar==&fWtot )
      eg = ( fWtot*fWtot + fQsquared - fMd*fMd ) / 2. / fMd;
    else if( fEVar==&fStot )
      eg = ( fStot + fQsquared - fMd*fMd ) / 2. / fMd;
    else if( (fEVar==&fWky || fEVar==&fWkn || fEVar==&fWyn) &&
	     (fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn) )
      eg = ( *fEVar**fEVar - *fTh1Var -fMd*fMd + 2.*fMd*e1 ) / 2. / fMd;
    else if( (fEVar==&fSky || fEVar==&fSkn || fEVar==&fSyn) &&
	     (fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn) )
      eg = ( *fEVar - *fTh1Var -fMd*fMd + 2.*fMd*e1 ) / 2. / fMd;
    else {
      // if energy is given by s_{23} or W_{23}, it is less straightforward
      // to find the photon's energy.
      // First, we set s23
      if( fEVar==&fWky || fEVar==&fWkn || fEVar==&fWyn )
	s23 = *fEVar * *fEVar;
      else
	s23 = *fEVar;

      if( s23<(pow(m2+m3,2.)-STRANGEUFLOW) ) fIsPhysical = false;

      // we can calculate eg now, if t_{g1} is the angular variable
      if( fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn )
	eg = (s23 - *fTh1Var - fMd*fMd + 2.*e1*fMd) /2. /fMd;
      else
	egCalculated = false;
    }

    // Polar angle 1st particle
    double costh1;
    if( fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn )
      costh1 = ( *fTh1Var + fQsquared - m1*m1 + 2.*eg*e1 ) 
	/ 2. / p1 / Sqrt(eg*eg+fQsquared);
    else if( fTh1Var==&fThk || fTh1Var==&fThy || fTh1Var==&fThn )
      costh1 = Cos( *fTh1Var );
    else
      costh1 = *fTh1Var;
    
    if( Abs(costh1) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( Abs(costh1)-1.0 < kUflow ) {
      fIsPhysical = true && fIsPhysical;
      costh1 = ( costh1>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;

    // sin polar angle always positive
    // cos -> sin calculation is very unstable
    // when cos close to +-1.
    double sinth1=-2.;
    if( 1.-Abs(costh1) <= kUflow) sinth1=0.;
    else sinth1 = Sqrt(1. - costh1*costh1); 

    // Photon energy for the remaining case,
    // that is when energy is given by s_{23} or W_{23}
    // and t_{g1} isn't the angular variable.
    //
    // We do this by solving s23 = (p_g + p_D - p_1)^2
    if( !egCalculated ) {
      if( fabs(costh1)<kUflow ) {
	eg = ( s23 + fQsquared - fMd*fMd - m1*m1 + 2.*fMd*e1 ) / 2. / ( fMd - e1 );
      } 
      else {
	Polynomial eq(pow((fMd-e1)/p1/costh1,2)-1.,
		      2.*(fMd-e1)/p1/costh1*(fMd*fMd+m1*m1-2.*fMd*e1-s23-fQsquared)/2./p1/costh1,
		      pow((fMd*fMd+m1*m1-2.*fMd*e1-s23-fQsquared)/2./p1/costh1,2)-fQsquared);
	vector<double> roots = eq.FindRealRoots();

	if(roots.size() == 0 ) {
	  fIsPhysical = false;
	} 
	else { 
	  // test the solutions:
	  // * should be positive
	  // * sign of sqrt(eg^2+Q^2) should be positive
	  bool sol1=false, sol2=false;
	  
	  if( roots[0]>0. && (s23 + fQsquared -fMd*fMd -m1*m1 
			      +2.*e1*fMd +2.*roots[0]*(e1-fMd))/costh1>0. ) sol1=true;
	  
	  if( roots.size()==2 ) {
	    if( roots[1]>0. && (s23 + fQsquared -fMd*fMd -m1*m1 
				+2.*e1*fMd +2.*roots[1]*(e1-fMd))/costh1>0. ) sol2=true;
	  }
	  
	  if( roots.size()==1 ) {
	    eg = sol1 ? roots[0] : -1;
	  }
	  else {
	    if( sol1 && sol2 ) {
	      if( fabs(roots[0]-roots[1]) > STRANGEUFLOW*fabs(roots[0]))
		TMPI::Cout() << "WARNING in TKinematics2to3::UpdateKinematics(): "
			     << "found multiple solutions. Choosing one arbitrarily.\n";
	      eg = roots[0];
	    }
	    else {
	      eg = sol1 ? roots[0] : roots[1];
	    }
	  }
	} // solution found
      } // non-zero costh1
    } // if( !egCalculated )

    double threshold = ((m1+m2+m3)*(m1+m2+m3)-fMd*fMd+fQsquared) / 2. / fMd;
    if( eg<threshold-STRANGEUFLOW ) fIsPhysical = false;

    // Set 4vectors
    fGvector->SetXYZT(0.,0.,Sqrt(eg*eg+fQsquared),eg);
    vector1.SetXYZT(p1*sinth1,0.,p1*costh1,e1);
    
    // Determine (23) subsystem
    //~~~~~~~~~~~~~~~~~~~~~~~~~

    vector23 = *fGvector + *fDvector - vector1;
    subSystem->SetVar(1,vector23.Mag2()); // set S_23
    subSystem->SetVar(3,-1.*(*fGvector-vector1).Mag2()); // set -t_g1
       

    // Determine lorentz rotations
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define rotation matrix that changes reference frame to
    // z' ~ p23,  y' ~ y and x' ~ y' x z'
    TLorentzRotation rotateTo23; 
    double rotationAngle=0.; 
    if( 1.-fabs(vector23.CosTheta())>kUflow )
      rotationAngle = -1.*TMath::Sign(1.,vector23.X())
	*vector23.Vect().Angle(fGvector->Vect());
    rotateTo23.RotateY(rotationAngle);
    // Define Lorentzboost that boosts in reference frame defined above
    // to 23 CM system.
    TVector3 boostVector = -1.0 * (rotateTo23*vector23).BoostVector();
    TLorentzRotation boostCM23 = boostVector;

 
    // Determine 4vector of 2nd particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CM energy particle 2
    double p2cm = subSystem->GetPk();
    double e2cm = Sqrt( p2cm*p2cm + m2*m2 );

    if( fKinSolution[2] == '*' ) {
      // CM angle particle 2
      double costh2cm;
      if( fTh2Var==&fKYthkcm || fTh2Var==&fKYthycm || fTh2Var==&fKNthncm ||
	  fTh2Var==&fKNthkcm || fTh2Var==&fYNthycm || fTh2Var==&fYNthncm )
	costh2cm = Cos(*fTh2Var);
      else if( fTh2Var==&fKYcosthkcm || fTh2Var==&fKYcosthycm || fTh2Var==&fKNcosthncm ||
	       fTh2Var==&fKNcosthkcm || fTh2Var==&fYNcosthycm || fTh2Var==&fYNcosthncm )
	costh2cm = *fTh2Var;
      else {
	double wcm = ( boostCM23*rotateTo23 * *fGvector ).E();
	double kcm = ( boostCM23*rotateTo23 * *fGvector ).P();
	costh2cm = ( *fTh2Var + fQsquared - m2*m2 + 2.*wcm*e2cm ) / 2. / kcm / p2cm;
      }
      
      if( Abs(costh2cm) <= 1.0 )
	fIsPhysical = true && fIsPhysical;
      else if( Abs(costh2cm)-1.0 < kUflow ) {
	fIsPhysical = true && fIsPhysical;
	costh2cm = ( costh2cm>0.0 ? 1.0 : -1.0 );
      }
      else
	fIsPhysical = false;

      // sin polar angle always positive
      // cos -> sin calculation is very unstable
      // when cos close to +-1.
      double sinth2cm=-2.;
      if( 1.-Abs(costh2cm) <= kUflow ) sinth2cm=0.;
      else sinth2cm = Sqrt(1. - costh2cm*costh2cm); 

      if( Abs(sinth2cm) <= 1.0 )
	fIsPhysical = true && fIsPhysical;
      else if( Abs(sinth2cm)-1.0 < kUflow ) {
	fIsPhysical = true && fIsPhysical;
	sinth2cm = ( sinth2cm>0.0 ? 1.0 : -1.0 );
      }
      else
	fIsPhysical = false;

      // set 4vector for particle 2 in (23) CM system
      vector2.SetXYZM(p2cm*sinth2cm*Cos(*fPhi2Var),
		      p2cm*sinth2cm*Sin(*fPhi2Var),
		      p2cm*costh2cm, m2);

      // calculate 4vector particle 3
      // Remark:
      // one would logically now boost/rotate vector2 back to
      // the LAB frame and then calculate vector3 by taking the
      // difference with vector23. This however (the god of the
      // floating point knows why), introduces a small error (order 1e-9)
      // that keeps propagating. The procedure used below is more
      // accurate. If the reader of this comment understands why,
      // please let me (=Pieter) know.
      vector3.SetXYZM(-vector2.X(),-vector2.Y(),-vector2.Z(),m3);
      vector3 = rotateTo23.Inverse() * (boostCM23.Inverse() * vector3);

      // calculate 4vector particle 2 in the LAB system
      vector2 = vector23 - vector3;

    } // angles 2nd particle in CM system
    else {
      // LAB angle particle 2
      double costh2;
      if( fTh2Var==&fThk || fTh2Var==&fThy || fTh2Var==&fThn )
	costh2 = Cos( *fTh2Var );
      else
	costh2 = *fTh2Var;
      
      if( Abs(costh2) <= 1.0 )
	fIsPhysical = true && fIsPhysical;
      else if( Abs(costh2)-1.0 < kUflow ) {
	fIsPhysical = true && fIsPhysical;
	costh2 = ( costh2>0.0 ? 1.0 : -1.0 );
      }
      else
	fIsPhysical = false;
      
      // sin polar angle always positive
      // cos -> sin calculation is very unstable
      // when cos close to +-1.
      double sinth2=-2.;
      if( 1.-Abs(costh2) <= kUflow ) sinth2=0.;
      else sinth2 = Sqrt(1. - costh2*costh2); 
      
      if( Abs(sinth2) <= 1.0 )
	fIsPhysical = true && fIsPhysical;
      else if( Abs(sinth2)-1.0 < kUflow ) {
	fIsPhysical = true && fIsPhysical;
	sinth2 = ( sinth2>0.0 ? 1.0 : -1.0 );
      }
      else
	fIsPhysical = false;

      // We wish to determine the 3-momentum of the 2nd particle in the LAB frame
      // Because we have the LAB angles we know the LAB 4-vector except for the
      // 3-momentum. This 4-vector can be boosted to the CM frame.
      // Equating the energy component of this boosted 4-vector with the energy of
      // the 2nd particle in the CM frame (which we calculated supra), leads to a
      // quadratic equation. 
      // The equation can have 2 physical solutions!
      double e23 = vector23.E();
      double w23 = vector23.M();
      double p23 = (*fGvector-vector1).P();
      double costh23 = vector23.CosTheta();
      if( 1.-fabs(costh23)<kUflow ) costh23 = TMath::Sign(1.,costh23);
      double sinth23 = -2.;
      if( 1.-Abs(costh23) <= kUflow ) sinth23=0.;
      else sinth23 = Sqrt(1. - costh23*costh23);

      Polynomial eq(// coefficient p_2^2
		    p23*p23*sinth23*sinth23*sinth2*sinth2*Cos(*fPhi2Var)*Cos(*fPhi2Var)
		    +p23*p23*costh23*costh23*costh2*costh2
		    -2*p23*p23*sinth23*costh23*costh2*sinth2*Cos(*fPhi2Var)
		    -e23*e23,
		    // coefficient p_2
		    2*e2cm*w23*p23*costh23*costh2
		    -2*e2cm*w23*p23*sinth23*sinth2*Cos(*fPhi2Var),
		    // constant term
		    e2cm*e2cm*w23*w23
		    -e23*e23*m2*m2);
      vector<double> roots = eq.FindRealRoots();
      
      double p2 = -1.;
      if(roots.size() == 0 ) {
	fIsPhysical = false;
      } 
      else { 
	// test the solutions:
	// * should be positive
	// * sign of sqrt(p2^2+m2^2) should be positive
	bool sol1=false, sol2=false;
	  
	if( roots[0]>0. && (e2cm*w23-roots[0]
			    *(p23*sinth23*sinth2*Cos(*fPhi2Var)
			      -p23*costh23*costh2))>0. ) sol1=true;
	  
	if( roots.size()==2 ) {
	  if( roots[1]>0. && (e2cm*w23-roots[1]
			      *(p23*sinth23*sinth2*Cos(*fPhi2Var)
				-p23*costh23*costh2))>0. ) sol2=true;
	}
	  
	if( roots.size()==1 ) {
	  p2 = sol1 ? roots[0] : -1;
	}
	else {
	  if( sol1 && sol2 ) {
	    if( fabs(roots[0]-roots[1]) > STRANGEUFLOW*fabs(roots[0]))
	      TMPI::Cout() << "WARNING in TKinematics2to3::UpdateKinematics(): "
			   << "found multiple solutions. Choosing one arbitrarily.\n";
	    p2 = roots[0];
	    TMPI::Cout() << "sol1= " << roots[0] << std::endl
			 << "sol2= " << roots[1] << std::endl;
	  }
	  else {
	    p2 = sol1 ? roots[0] : roots[1];
	  }
	}
      } // solution found
      
      if(p2<0.) fIsPhysical = false;

      // set 4vector for particle 2
      vector2.SetXYZM(sinth2*Cos(*fPhi2Var)*p2,
		      sinth2*Sin(*fPhi2Var)*p2,
		      costh2*p2, m2);

      vector3 = vector23 - vector2;
   
    } // angles 2nd particle in LAB system
   
  } // end 1 + (2+3) case

  // 2) The user provides info about the energy (fE1Var)
  //    and the theta-angle (fTh1Var) of the particle that
  //    fixes the x-axis (fAC).
  //    In addition a second final particle is fully specified,
  //    energy, theta and phi (fE2Var,fTh2Var,fPhi2Var) in the lab frame.
  //    The name of the 2nd particle will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "2N-": scenario 2, 2nd particle = nucleon

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!                                                            !!
  // !! No longer considered a valid way of solving the kinematics !!
  // !!                                                            !!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  else if( fKinSolution[0] == '2' ) {
   
    //     // Determine 4vector of 1st particle
    //     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //     // Energy 1st particle
    //     double e1,p1;
    //     if( fE1Var==&fPk || fE1Var==&fPy || fE1Var==&fPn ) {
    //       p1 = *fE1Var;
    //       e1 = Sqrt( p1*p1 + m1*m1 );
    //     } else {
    //       e1 = *fE1Var;
    //       p1 = Sqrt( e1*e1 - m1*m1 );
    //     }
    //     if(e1<0. || p1 <0.) fIsPhysical = false;

    //     // Polar angle 1st particle
    //     double costh1;
    //     if( fTh1Var==&fThk || fTh1Var==&fThy || fTh1Var==&fThn )
    //       costh1 = Cos( *fTh1Var );
    //     else
    //       costh1 = *fTh1Var;
    
    //     if( Abs(costh1) <= 1.0 )
    //       fIsPhysical = true && fIsPhysical;
    //     else if( Abs(costh1)-1.0 < kUflow ) {
    //       fIsPhysical = true && fIsPhysical;
    //       costh1 = ( costh1>0.0 ? 1.0 : -1.0 );
    //     }
    //     else
    //       fIsPhysical = false;

    //     // set 4vector
    //     vector1.SetXYZT(p1*Sqrt(1.-costh1*costh1),0.,p1*costh1,e1);

 
    //     // Determine 4vector of 2nd particle
    //     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //     // Energy 2nd particle
    //     double e2,p2;
    //     if( fE2Var==&fPk || fE2Var==&fPy || fE2Var==&fPn ) {
    //       p2 = *fE2Var;
    //       e2 = Sqrt( p2*p2 + m2*m2 );
    //     } else {
    //       e2 = *fE2Var;
    //       p2 = Sqrt( e2*e2 - m2*m2 );
    //     }
    //     if(e2<0. || p2 <0.) fIsPhysical = false;

    //     // Polar angle 2nd particle
    //     double costh2;
    //     if( fTh2Var==&fThk || fTh2Var==&fThy || fTh2Var==&fThn )
    //       costh2 = Cos( *fTh2Var );
    //     else
    //       costh2 = *fTh2Var;
    //     double sinth2 = Sqrt(1. -costh2*costh2);
    
    //     if( Abs(costh2) <= 1.0 )
    //       fIsPhysical = true && fIsPhysical;
    //     else if( Abs(costh2)-1.0 < kUflow ) {
    //       fIsPhysical = true && fIsPhysical;
    //       costh2 = ( costh2>0.0 ? 1.0 : -1.0 );
    //     }
    //     else
    //       fIsPhysical = false;

    //     // set 4vector
    //     vector2.SetXYZT(sinth2*Cos(*fPhi2Var)*p2,
    // 		    sinth2*Sin(*fPhi2Var)*p2,
    // 		    costh2*p2, e2);
    
    //     // Determine 4vector of 3rd particle and photon
    //     // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //     // We do this by writing conservation of 4momentum: p_g+p_D=p_1+p_2+p_3
    //     // in components:
    //     //        wlab + fMd  = e1 + e2 + e3                     (1)
    //     //                 0  = p1x + p2x + p3*sinth3*cos(phi3)  (2)
    //     //                 0  = p2y + p3*sinth3*sin(phi3)        (3)
    //     // sqrt(wlab^2 + Q^2) = p1z + p2z + p3*costh3            (4)
    //     //
    //     // costh3 = +-sqrt(1 - ((p1x+p2x)^2+p2y^2)/p3^2 )        (5) = (2)^2+(3)^2

    //     // Azimuthal angle: tan(phi3) = (3)/(2) = -p2y / (-p1x-p2x)
    //     double phi3 = 
    //       TVector2::Phi_0_2pi(atan2(-1.*vector2.Y(),-1.*(vector1.X()+vector2.X())));

    //     // momentum 3rd particle. Iteratively solve (4) after inserting (1) and (5)
    //     double p3low2=pow(vector1.X()+vector2.X(),2)+vector2.Y()*vector2.Y(),
    //       p3high2=p3low2;

    //     // first we look for a lower and upper limit
    //     while( fabs( -1.*(vector1.Z()+vector2.Z()) 
    // 		 + Sqrt(pow(e1+e2+Sqrt(p3high2+m3*m3)-fMd,2)+fQsquared) )
    // 	 - Sqrt(p3high2-(pow(vector1.X()+vector2.X(),2)+vector2.Y()*vector2.Y()))
    // 	 > 0. ) {
    //       p3low2 = p3high2;
    //       p3high2 *= 2.;
    //     }
    //     if(p3high2 <= p3low2) fIsPhysical = false;

    //     // solve for p3^2
    //     double delta = p3high2 - p3low2;
    
    //     while( delta/p3low2>1e-9 && p3low2<p3high2 ) {
     
    //       if( fabs( -1.*(vector1.Z()+vector2.Z()) 
    // 		+ Sqrt(pow(e1+e2+Sqrt(p3low2+delta/2.+m3*m3)-fMd,2)
    // 		       +fQsquared) )
    // 	  - Sqrt(p3low2+delta/2.
    // 		 -(pow(vector1.X()+vector2.X(),2)+vector2.Y()*vector2.Y()))
    // 	  < 0. )
    // 	p3high2 -= delta/2.;
    //       else
    // 	p3low2 += delta/2.;
      
    //       delta = p3high2 - p3low2;
    //     }
    
    //     double p3=0.;
    //     if( p3low2>=0. ) p3 = Sqrt(p3low2);
    //     else fIsPhysical = false;


    //     // energy 3rd particle
    //     double e3 = Sqrt( p3*p3 + m3*m3 );

    //     // Determine photon 4vector
    //     //~~~~~~~~~~~~~~~~~~~~~~~~~
    //     // photon energy from (1)
    //     double eg = e1 + e2 + e3 - fMd;

    //     // set photon 4vector
    //     fGvector->SetXYZT(0.,0.,Sqrt(eg*eg+fQsquared),eg);

    //     // Polar angle 3rd particle from (4)
    //     double costh3 = ( Sqrt(eg*eg+fQsquared)-vector1.Z()-vector2.Z() ) / p3;
    //     if( Abs(costh3) <= 1.0 )
    //       fIsPhysical = true && fIsPhysical;
    //     else if( Abs(costh3)-1.0 < 1e-14 ) {
    //       fIsPhysical = true && fIsPhysical;
    //       costh3 = ( costh3>0.0 ? 1.0 : -1.0 );
    //     }
    //     else
    //       fIsPhysical = false;
    //     double sinth3 = Sqrt(1. - costh3*costh3);

    //     // set 4vector 3rd particle
    //     vector3.SetXYZT(p3*sinth3*Cos(phi3),
    // 		    p3*sinth3*Sin(phi3),
    // 		    p3*costh3, e3);

  } // end 1 + 2 (+3) case

  // 3) All variables are mandelstam variables, except for Q^2
  //    obviously. 
  //    S_tot = (P_gamma + P_deut)^2 should always be given.
  //    fAC fixes the 1st particle. This leaves 7 combinations:
  //    (3a) s_12 + s_23 + t_g1 + t_d3
  //    (3b) s_12 + s_23 + t_g1 + t_d2
  //    (3c) s_12 + s_23 + t_d2 + t_d3
  //    (3d) s_12 + t_d2 + t_g1 + t_d3
  //    (3e) s_12 + s_23 + t_g1 + t_g2
  //    (3f) s_12 + s_23 + t_g2 + t_d3
  //    (3g) t_g2 + s_23 + t_g1 + t_d3
  //    The name of the 2nd particle will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "3Kb": scenario 3b with 2nd particle = kaon

  else if( fKinSolution[0] == '3' ) {

    double *s12,*s23,*tg1,*tg2,*td2,*td3;

    if( fAC==kK ) {
      tg1 = &fTgk;
      s23 = &fSyn;
      if( fKinSolution[1] == 'Y' ) {
	s12 = &fSky;
	tg2 = &fTgy;
	td2 = &fTdy;
	td3 = &fTdn;
      } else {
	s12 = &fSkn;
	tg2 = &fTgn;
	td2 = &fTdn;
	td3 = &fTdy;
      }
    } else if( fAC==kY ) {
      tg1 = &fTgy;
      s23 = &fSkn;
      if( fKinSolution[1] == 'K' ) {
	s12 = &fSky;
	tg2 = &fTgk;
	td2 = &fTdk;
	td3 = &fTdn;
      } else {
	s12 = &fSyn;
	tg2 = &fTgn;
	td2 = &fTdn;
	td3 = &fTdk;
      }
    } else {
      tg1 = &fTgn;
      s23 = &fSky;
      if( fKinSolution[1] == 'K' ) {
	s12 = &fSkn;
	tg2 = &fTgk;
	td2 = &fTdk;
	td3 = &fTdy;
      } else {
	s12 = &fSyn;
	tg2 = &fTgy;
	td2 = &fTdy;
	td3 = &fTdk;
      }
    }

    // determine the remaining mandelstam variables
    switch(fKinSolution[2]) {
    case 'a':
      *tg2 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg1 - *s12;
      *td2 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td3 - *s23;
      break;
    case 'b':
      *td3 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td2 - *s23;
      *tg2 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg1 - *s12;
      break;
    case 'c':
      *tg1 = *td2 + *td3 + *s23 - fMd*fMd - m2*m2 - m3*m3;
      *tg2 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg1 - *s12;
      break;
    case 'd':
      *s23 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td2 - *td3;
      *tg2 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg1 - *s12;
      break;
    case 'e':
      *td3 = *tg2 + *tg1 + *s12 + fQsquared - m1*m1 - m2*m2;
      *td2 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td3 - *s23;
      break;
    case 'f':
      *tg1 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg2 - *s12;
      *td2 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td3 - *s23;
      break;
    case 'g':
      *s12 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg2 - *tg1;
      *td2 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td3 - *s23;
      break;
    default:
      cerr << "ERROR in TKinematics2to3::UpdateKinematics(): "
	   << "Invalid option.\n";
      exit(1);
    }

    double threshold = (m1+m2+m3)*(m1+m2+m3);
    if( fStot<threshold-STRANGEUFLOW ) fIsPhysical = false;

    // Determine photon 4vector
    //~~~~~~~~~~~~~~~~~~~~~~~~~
    // Photon energy
    double eg = ( fStot + fQsquared - fMd*fMd ) / 2. / fMd;
    if(eg<=0.) fIsPhysical = false;
    double klab = Sqrt(eg*eg+fQsquared);

    // set photon 4vector
    fGvector->SetXYZT(0.,0.,klab,eg);

    
    // Determine 4vector 1st particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Energy 1st particle
    double e1 = (2.*fMd*eg - *s23 + *tg1 + fMd*fMd) / 2. / fMd;
    if(e1<m1) fIsPhysical = false;
    double p1 = Sqrt(e1*e1 - m1*m1);

    // Polar angle 1st particle
    double costh1 = ( *tg1 + fQsquared - m1*m1 + 2.*eg*e1 ) / 2. / klab / p1;
    
    if( Abs(costh1) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( Abs(costh1)-1.0 < kUflow ) {
      fIsPhysical = true && fIsPhysical;
      costh1 = ( costh1>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;
    
    // sin polar angle always positive
    // cos -> sin calculation is very unstable
    // when cos close to +-1.
    double sinth1=-2.;
    if( 1.-Abs(costh1) <= kUflow ) sinth1=0.;
    else sinth1 = Sqrt(1. - costh1*costh1); 

    // set 1st particle 4vector
    vector1.SetXYZT(p1*sinth1,0.,p1*costh1,e1);


    // Determine (23) subsystem energy wise
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vector23 = *fGvector + *fDvector - vector1;
    subSystem->SetVar(1,*s23);
    subSystem->SetVar(3,-1.* *tg1);


    // Determine lorentz rotations
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define rotation matrix that changes reference frame to
    // z' ~ p23,  y' ~ y and x' ~ y' x z'
    TLorentzRotation rotateTo23; 
    double rotationAngle=0.; 
    if( 1.-fabs(vector23.CosTheta())>kUflow )
      rotationAngle = -1.*TMath::Sign(1.,vector23.X())
	*vector23.Vect().Angle(fGvector->Vect());
    rotateTo23.RotateY(rotationAngle);
    // Define Lorentzboost that boosts in reference frame defined above
    // to 23 CM system.
    TVector3 boostVector = -1.0 * (rotateTo23*vector23).BoostVector();
    TLorentzRotation boostCM23 = boostVector;


    // Determine 4vector 2nd and 3rd particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Energy 2nd particle
    double e2 = ( fMd*fMd + m2*m2 - *td2 ) / 2. / fMd;
    if(e2<m2) fIsPhysical = false;
    double p2 = Sqrt(e2*e2 - m2*m2);

    // Energy 3rd particle
    double e3 = eg + fMd - e1 - e2;
    if(e3<m3) fIsPhysical = false;
    double p3 = Sqrt(e3*e3 - m3*m3);

    // Polar angle 2nd particle
    double costh2 = ( *tg2 + fQsquared - m2*m2 + 2.*eg*e2 ) / 2. / klab / p2;
    
    if( Abs(costh2) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( Abs(costh2)-1.0 < kUflow ) {
      fIsPhysical = true && fIsPhysical;
      costh2 = ( costh2>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;

    // sin polar angle always positive
    // cos -> sin calculation is very unstable
    // when cos close to +-1.
    double sinth2=-2.;
    if( 1.-Abs(costh2) <= kUflow ) sinth2=0.;
    else sinth2 = Sqrt(1. - costh2*costh2);
    
    //express conservation of 4momentum: p_g+p_D=p_1+p_2+p_3 in components:
    //
    //    p3*sinth3*cos(phi3) = -p1x - p2*sinth2*cos(phi2)  (1)
    //    p3*sinth3*sin(phi3) = -p2*sinth2*cos(phi2)        (2)
    
    // Polar angle 3rd particle
    double costh3 = ( vector23.Z() - p2*costh2 ) / p3;
    
    if( Abs(costh3) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( Abs(costh3)-1.0 < kUflow ) {
      fIsPhysical = true && fIsPhysical;
      costh3 = ( costh3>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;

    // sin polar angle always positive
    // cos -> sin calculation is very unstable
    // when cos close to +-1.
    double sinth3=-2.;
    if( 1.-Abs(costh3) <= kUflow ) sinth3=0.;
    else sinth3 = Sqrt(1. - costh3*costh3);

    // Find azimuthal angle particle 2
    double phi2;
    if( vector1.X() != 0. ) {
      // Azimuthal angle particle 2 from (1)^2 + (2)^2
      phi2 =ACos( ( p3*p3*sinth3*sinth3 - vector1.X()*vector1.X()
		    - p2*p2*sinth2*sinth2 )
		  / 2. / vector1.X() / p2 / sinth2 ) ;
      
      // with the formula above we can only determine Phi_2 upto its sign.
      // This might seem wrong but it isn't. There's a point symmetry.
    }

    else {
      // if particle 1 has no x-component, the x- and y-axis are not properly defined
      phi2=0;      
    }

    // set 4vector for particle 2 again
    vector2.SetXYZT(sinth2*Cos(phi2)*p2,
		    sinth2*Sin(phi2)*p2,
		    costh2*p2, e2);
    

    // Determine 4vector of 3rd particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vector3 = vector23 - vector2;

  } // end of mandelstam case


  // Put the final state 4vectors in the correct members
  if( fAC==kK ) {
    *fKvector = vector1;
    if( fKinSolution[1] == 'Y' ) {
      *fYvector = vector2;
      *fNvector = vector3;
    } else {
      *fNvector = vector2;
      *fYvector = vector3;
    }
  } else if( fAC==kY ) {
    *fYvector = vector1;
    if( fKinSolution[1] == 'K' ) {
      *fKvector = vector2;
      *fNvector = vector3;
    } else {
      *fNvector = vector2;
      *fKvector = vector3;
    }
  } else {
    *fNvector = vector1;
    if( fKinSolution[1] == 'K' ) {
      *fKvector = vector2;
      *fYvector = vector3;
    } else {
      *fYvector = vector2;
      *fKvector = vector3;
    }
  }

  // The final state 4vectors have been calculated in the user's frame of choice.
  // When necessary we need to rotate them to our choice of frame (=kN).
  if( fAC!=kN ) {
    double angle = -1. * fNvector->Vect().XYvector().Phi();
    fKvector->RotateZ(angle);
    fYvector->RotateZ(angle);
    fNvector->RotateZ(angle);
  }
  
  UpdateSubsystems();

  // Do a final check whether we are in a physical point
  // by checking conservation of energy and momentum
  TLorentzVector sum = *fGvector + *fDvector - *fKvector - *fYvector - *fNvector;
  if( fabs(sum.E()) > STRANGEUFLOW || fabs(sum.X()) > STRANGEUFLOW ||
      fabs(sum.Y()) > STRANGEUFLOW || fabs(sum.Z()) > STRANGEUFLOW )
    fIsPhysical = false;

  fIsPhysical *= fKYsystem->IsPhysical() 
    * fKNsystem->IsPhysical() 
    * fYNsystem->IsPhysical();
}

//_____________________________________________________________________
void TKinematics2to3::UpdateSubsystems()
{
  // Recalculate the KY-, KN- and YN-subsystems.

  fKYsystem->SetVar(1,(*fKvector+*fYvector).Mag2()); // set s_ky
  fKYsystem->SetVar(2,(*fGvector-*fNvector-*fKvector).Mag2()); // set (p_g-p_n-p_k)^2
  fKYsystem->SetVar(3,-1.*(*fGvector-*fNvector).Mag2()); // set -(p_g-p_n)^2

  fKNsystem->SetVar(1,(*fKvector+*fNvector).Mag2()); // set s_kn
  fKNsystem->SetVar(2,(*fGvector-*fYvector-*fKvector).Mag2()); // set (p_g-p_y-p_k)^2
  fKNsystem->SetVar(3,-1.*(*fGvector-*fYvector).Mag2()); // set -(p_g-p_y)^2

  fYNsystem->SetVar(1,(*fYvector+*fNvector).Mag2()); // set s_yn
  fYNsystem->SetVar(2,(*fGvector-*fKvector-*fNvector).Mag2()); // set (p_g-p_k-p_N)^2
  fYNsystem->SetVar(3,-1.*(*fGvector-*fKvector).Mag2()); // set -(p_g-p_k)^2
}

//_____________________________________________________________________
void TKinematics2to3::UpdateVariables()
{
  // After UpdateKinematics() not all data member are up-to-date.
  // This function makes sure they all are.

  fWlab = GetWlab();
  fKlab = GetKlab();
  fEk = GetEk();
  fEy = GetEy();
  fEn = GetEn();
  fPk = GetPk();
  fPy = GetPy();
  fPn = GetPn();
  fThk = GetThk();
  fThy = GetThy();
  fThn = GetThn();
  fCosthk = GetCosthk();
  fCosthy = GetCosthy();
  fCosthn = GetCosthn();
  fPhik = GetPhik(kN);
  fPhiy = GetPhiy(kN);
  fPhin = GetPhin(kN);
  fWtot = GetWtot();
  fWky = GetWky();
  fWkn = GetWkn();
  fWyn = GetWyn();
  fStot = GetStot();
  fSky = GetSky();
  fSkn = GetSkn();
  fSyn = GetSyn();
  fTgk = GetTkg();
  fTgy = GetTyg();
  fTgn = GetTng();
  fTdk = (*fDvector-*fKvector).Mag2();
  fTdy = (*fDvector-*fYvector).Mag2();
  fTdn = (*fDvector-*fNvector).Mag2();
  fKYthkcm = ACos(GetKYsystem_costhkcm());
  fKYthycm = ACos(GetKYsystem_costhycm());
  fKYcosthkcm = GetKYsystem_costhkcm();
  fKYcosthycm = GetKYsystem_costhycm();
  fKYphikcm = GetKYsystem_phikcm();
  fKYphiycm = GetKYsystem_phiycm();
  fKNthkcm = ACos(GetKNsystem_costhkcm());
  fKNthncm = ACos(GetKNsystem_costhncm());
  fKNcosthkcm = GetKNsystem_costhkcm();
  fKNcosthncm = GetKNsystem_costhncm();
  fKNphikcm = GetKNsystem_phikcm();
  fKNphincm = GetKNsystem_phincm();
  fYNthycm = ACos(GetYNsystem_costhycm());
  fYNthncm = ACos(GetYNsystem_costhncm());
  fYNcosthycm = GetYNsystem_costhycm();
  fYNcosthncm = GetYNsystem_costhncm();
  fYNphiycm = GetYNsystem_phiycm();
  fYNphincm = GetYNsystem_phincm();
}

//_____________________________________________________________________
void TKinematics2to3::SetVarRange(int varNumber, double lowerLimit,
				  double upperLimit, int numberOfSteps)
{
  // Give independent variable number 'varNumber' a range in the interval [lowerLimit,upperLimit].

  // First we do some checks
  if( varNumber<1 || varNumber>6 )
    std::cerr << "ERROR in TKinematics2to3::SetVarRange(int,double,double,int):\n"
	      << "Index out-of-bounds!\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics2to3::SetVarRange(int,double,double,int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else if( numberOfSteps<2 )
    std::cerr << "ERROR in TKinematics2to3::SetVarRange(int,double,double,int)\n"
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

    // Reset the number of physical steps
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
void TKinematics2to3::FixVariable(int varNumber)
{
  // Unset the range of variable number 'varNumber'. Set it to its current value.

  // First we do some checks
  if( varNumber<1 || varNumber>6 )
    std::cerr << "ERROR in TKinematics2to3::FixVariable(int):\n"
	      << "Index out-of-bounds!\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics2to3::FixVariable(int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else if( fIsVar[varNumber-1] ){
    // Determine the current step
    int step1 = GetStep(1);
    int step2 = GetStep(2);
    int step3 = GetStep(3);
    int step4 = GetStep(4);
    int step5 = GetStep(5);
    int step6 = GetStep(6);
    fIsVar[varNumber-1] = false;
    
    // Reset all variables
    fLowerLimit[varNumber-1] = *fVar[varNumber-1];
    fUpperLimit[varNumber-1] = *fVar[varNumber-1];
    fStepSize[varNumber-1] = 0.0;
    fNumberOfSteps[varNumber-1] = 1;

    // Go to first step
    switch(varNumber) {
    case 1: GoTo(0,step2,step3,step4,step5,step6); break;
    case 2: GoTo(step1,0,step3,step4,step5,step6); break;
    case 3: GoTo(step1,step2,0,step4,step5,step6); break;
    case 4: GoTo(step1,step2,step3,0,step5,step6); break;
    case 5: GoTo(step1,step2,step3,step4,0,step6); break;
    case 6: GoTo(step1,step2,step3,step4,step5,0); break;
    }

    // Reset the number of physical steps
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
void TKinematics2to3::FixVariables()
{
  // Calls FixVariable(int) on all variables.
  for(int i=1; i<=6; ++i)
    FixVariable(i);
}

//_____________________________________________________________________
int TKinematics2to3::Next()
{
  // Proceed to the next step in the grid of independent variables.
  // Return value is zero when the end of the grid is reached.

  // Total number of steps
  int totNrOfSteps = 1;
  for(int i=0; i<6; ++i) totNrOfSteps *= fNumberOfSteps[i];

  // Check whether we are at the end of the grid.
  if( fStep == (totNrOfSteps-1) ) return 0;

  GoTo(++fStep);
  return 1;
}

//_____________________________________________________________________
void TKinematics2to3::GoTo(int step)
{
  // Proceed to the 'step'th point in the grid of independent variables.

  // Total number of steps
  int totNrOfSteps = 1;
  for(int i=0; i<6; ++i)
    totNrOfSteps *= fNumberOfSteps[i];

  // First we do some checks
  if( step<0 || step>=totNrOfSteps )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int):\n"
	      << "Invalid step. Should be in range [0,"
	      << totNrOfSteps-1 << "].\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics2to3::GoTo(int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else {
    fStep = step;
    
    // Determine step per variable
    int stepVar[6];
    int remainder=step;
    int denum=1;
  
    for(int i=6; i>0; --i)
      {
	denum = 1;
	for(int j=0; j<i-1; ++j)
	  denum *= fNumberOfSteps[j];
	
	stepVar[i-1] = remainder / denum;
	remainder = remainder % denum;
      }
    
    // Set variables
    for(int i=0; i<6; ++i)
      *fVar[i] = fLowerLimit[i] + stepVar[i] * fStepSize[i];

    UpdateKinematics();
  }
}

//_____________________________________________________________________
void TKinematics2to3::GoTo(int stepVar1, int stepVar2, int stepVar3,
			   int stepVar4, int stepVar5, int stepVar6)
{
  // Proceed to point in the grid of independent variables specified by the
  // arguments. 
  // Each stepVari fixes the position in the grid of the 'i'th variable.

  // First we do some checks
  if( stepVar1<0 || stepVar1>=fNumberOfSteps[0] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(1) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[0]-1 << "].\n";

  else if( stepVar2<0 || stepVar2>=fNumberOfSteps[1] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(2) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[1]-1 << "].\n";

  else if( stepVar3<0 || stepVar3>=fNumberOfSteps[2] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(3) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[2]-1 << "].\n";

  else if( stepVar4<0 || stepVar4>=fNumberOfSteps[3] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(4) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[3]-1 << "].\n";

  else if( stepVar5<0 || stepVar5>=fNumberOfSteps[4] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(5) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[4]-1 << "].\n";

  else if( stepVar6<0 || stepVar6>=fNumberOfSteps[5] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(6) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[5]-1 << "].\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6))\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else {
    fStep = stepVar1 
      + stepVar2 * fNumberOfSteps[0]
      + stepVar3 * fNumberOfSteps[0] * fNumberOfSteps[1]
      + stepVar4 * fNumberOfSteps[0] * fNumberOfSteps[1] * fNumberOfSteps[2]
      + stepVar5 * (fNumberOfSteps[0] * fNumberOfSteps[1] 
		    * fNumberOfSteps[2]*fNumberOfSteps[3])
      + stepVar5 * (fNumberOfSteps[0] * fNumberOfSteps[1] * fNumberOfSteps[2] 
		    * fNumberOfSteps[3]*fNumberOfSteps[4]);
    
    // Set variables
    *fVar[0] = fLowerLimit[0] + stepVar1 * fStepSize[0];
    *fVar[1] = fLowerLimit[1] + stepVar2 * fStepSize[1];
    *fVar[2] = fLowerLimit[2] + stepVar3 * fStepSize[2];
    *fVar[3] = fLowerLimit[3] + stepVar4 * fStepSize[3];
    *fVar[4] = fLowerLimit[4] + stepVar5 * fStepSize[4];
    *fVar[5] = fLowerLimit[5] + stepVar6 * fStepSize[5];

    UpdateKinematics();
  }
}

//_____________________________________________________________________
int TKinematics2to3::VariableInfo() const
{
  // Print information on the type of independent variables and their status.

  TMPI::Cout() << GetName() << ": " << GetTitle()
	       << "\n****************************************************\n";

  for(int var=0; var<6; ++var) {
    TMPI::Cout() << var+1 << ": " << GetVarName(var+1) << " is ";
    
    if( fIsVar[var] )
      TMPI::Cout() << "variable.\n"
		   << "   Range:\t\t" << fLowerLimit[var] << " <-> "
		   << fUpperLimit[var]
		   << "\n   Number of steps:\t" << fNumberOfSteps[var]
		   << "\n   Current step:\t" << GetStep(var+1) << "\n";
    
    else
      TMPI::Cout() << "fixed.\n";
  }
  
  TMPI::Cout() << "****************************************************\n";

  return GetNumberOfVariables();
}

//_____________________________________________________________________
void TKinematics2to3::SetVar(int varNumber, double value)
{
  if( varNumber<1 || varNumber>6 )
    std::cerr << "ERROR in TKinematics2to3::SetVar(int,double): "
	      << "Index out-of-bounds!\n";

  else if( !fVar )
    std::cerr << "ERROR in TKinematics2to3::SetVar(int,double): "
	      << "Independent kinematic variables have not been initialized!\n";
    
  else if( fIsVar[varNumber-1] )
    std::cerr << "ERROR in TKinematics2to3::SetVar(int,double): "
	      << "This method can not be used because the variable is not fixed.\n";

  else {
    fLowerLimit[varNumber-1] = value;
    fUpperLimit[varNumber-1] = value;
    GoTo(0);
    
    // Reset the number of physical steps
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
void TKinematics2to3::SetVar(double var1, double var2, double var3,
			     double var4, double var5, double var6)
{
  if( !fVar )
    std::cerr << "ERROR in TKinematics2to3::SetVar(double x6): "
	      << "Independent kinematic variables have not been initialized!\n";
    
  else if( !IsFixed() )
    std::cerr << "ERROR in TKinematics2to3::SetVar(double x6): "
	      << "This method can not be used because some variables are not fixed.\n";

  else {
    fLowerLimit[0] = var1;
    fUpperLimit[0] = var1;
    fLowerLimit[1] = var2;
    fUpperLimit[1] = var2;
    fLowerLimit[2] = var3;
    fUpperLimit[2] = var3;
    fLowerLimit[3] = var4;
    fUpperLimit[3] = var4;
    fLowerLimit[4] = var5;
    fUpperLimit[4] = var5;
    fLowerLimit[5] = var6;
    fUpperLimit[5] = var6;
    GoTo(0);

    // Reset the number of physical steps
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
double TKinematics2to3::GetVar(int varNumber) const
{
  // Returns the value of the independent variable number 'varNumber'
  // as specified in the format string.

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetVar(int):\n"
	      << "Index out-of-bounds!\n";
    return -1e16;
  }

  if( !fVar ) {
    std::cerr << "ERROR in TKinematics2to3::GetVar(int):\n"
	      << "Independent kinematic variables have not been initialized!\n";
    return -1e16;
  }

  return *fVar[varNumber-1];
  
}

//_____________________________________________________________________
double TKinematics2to3::GetVar(const TString& variable) const
{
  // Returns the value of the variable given by name. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.

  double(TKinematics2to3::*ptrToVarGetter)() const; // pointer to TKinematics getters
  
  TString varName = variable;
  varName.ToLower();
  
  // Assign the correct Getter
  if( !strcmp(varName,"pk") ) {
    ptrToVarGetter = &TKinematics2to3::GetPk;
  }
  else if( !strcmp(varName,"py") ) {
    ptrToVarGetter = &TKinematics2to3::GetPy;
  }
  else if( !strcmp(varName,"pn") ) {
    ptrToVarGetter = &TKinematics2to3::GetPn;
  }
  else if( !strcmp(varName,"klab") ) {
    ptrToVarGetter = &TKinematics2to3::GetKlab;
  }
  else if( !strcmp(varName,"ek") ) {
    ptrToVarGetter = &TKinematics2to3::GetEk;
  }
  else if( !strcmp(varName,"ey") ) {
    ptrToVarGetter = &TKinematics2to3::GetEy;
  }
  else if( !strcmp(varName,"en") ) {
    ptrToVarGetter = &TKinematics2to3::GetEn;
  }
  else if( !strcmp(varName,"wlab") ) {
    ptrToVarGetter = &TKinematics2to3::GetWlab;
  }
  else if( !strcmp(varName,"thk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk_deg;
  }
  else if( !strcmp(varName,"thy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy_deg;
  }
  else if( !strcmp(varName,"thn_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn_deg;
  }
  else if( !strcmp(varName,"thk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk;
  }
  else if( !strcmp(varName,"thy") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy;
  }
  else if( !strcmp(varName,"thn") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn;
  }
  else if( !strcmp(varName,"costhk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthk;
  }
  else if( !strcmp(varName,"costhy") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthy;
  }
  else if( !strcmp(varName,"costhn") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthn;
  }
  else if( !strcmp(varName,"phik_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_deg_kS;
  }
  else if( !strcmp(varName,"phiy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_deg_kS;
  }
  else if( !strcmp(varName,"phin_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_deg_kS;
  }
  else if( !strcmp(varName,"phik") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_kS;
  }
  else if( !strcmp(varName,"phiy") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_kS;
  }
  else if( !strcmp(varName,"phin") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_kS;
  }
  else if( !strcmp(varName,"thky_deg") ||
	   !strcmp(varName,"thyk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThky_deg;
  }
  else if( !strcmp(varName,"thkn_deg") ||
	   !strcmp(varName,"thnk_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn_deg;
  }
  else if( !strcmp(varName,"thyn_deg") ||
	   !strcmp(varName,"thny_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn_deg;
  }
  else if( !strcmp(varName,"thky") ||
	   !strcmp(varName,"thyk" )) {
    ptrToVarGetter = &TKinematics2to3::GetThky;
  }
  else if( !strcmp(varName,"thkn") ||
	   !strcmp(varName,"thnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn;
  }
  else if( !strcmp(varName,"thyn") ||
	   !strcmp(varName,"thny") ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn;
  }
  else if( !strcmp(varName,"costhky") ||
	   !strcmp(varName,"costhyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthky;
  }
  else if( !strcmp(varName,"costhkn") ||
	   !strcmp(varName,"costhnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthkn;
  }
  else if( !strcmp(varName,"costhyn") ||
	   !strcmp(varName,"costhny") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthyn;
  }
  else if( !strcmp(varName,"qsquared") ) {
    ptrToVarGetter = &TKinematics2to3::GetQsquared;
  }
  else if( !strcmp(varName,"wtot") ) {
    ptrToVarGetter = &TKinematics2to3::GetWtot;
  }
  else if( !strcmp(varName,"wky") ||
	   !strcmp(varName,"wyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWky;
  }
  else if( !strcmp(varName,"wkn") ||
	   !strcmp(varName,"wnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWkn;
  }
  else if( !strcmp(varName,"wyn") ||
	   !strcmp(varName,"wny") ) {
    ptrToVarGetter = &TKinematics2to3::GetWyn;
  }
  else if( !strcmp(varName,"stot") ) {
    ptrToVarGetter = &TKinematics2to3::GetStot;
  }
  else if( !strcmp(varName,"sky") ||
	   !strcmp(varName,"syk") ) {
    ptrToVarGetter = &TKinematics2to3::GetSky;
  }
  else if( !strcmp(varName,"skn") ||
	   !strcmp(varName,"snk" )) {
    ptrToVarGetter = &TKinematics2to3::GetSkn;
  }
  else if( !strcmp(varName,"syn") ||
	   !strcmp(varName,"sny" )) {
    ptrToVarGetter = &TKinematics2to3::GetSyn;
  }
  else if( !strcmp(varName,"tkg") ||
	   !strcmp(varName,"tgk" )) {
    ptrToVarGetter = &TKinematics2to3::GetTkg;
  }
  else if( !strcmp(varName,"tyg") ||
	   !strcmp(varName,"tgy" )) {
    ptrToVarGetter = &TKinematics2to3::GetTyg;
  }
  else if( !strcmp(varName,"tng") ||
	   !strcmp(varName,"tgn" )) {
    ptrToVarGetter = &TKinematics2to3::GetTng;
  }
  else if( !strcmp(varName,"tkd") ||
	   !strcmp(varName,"tdk" )) {
    ptrToVarGetter = &TKinematics2to3::GetTkd;
  }
  else if( !strcmp(varName,"tyd") ||
	   !strcmp(varName,"tdy" )) {
    ptrToVarGetter = &TKinematics2to3::GetTyd;
  }
  else if( !strcmp(varName,"tnd") ||
	   !strcmp(varName,"tdn" )) {
    ptrToVarGetter = &TKinematics2to3::GetTnd;
  }
  else if( !strcmp(varName,"kypk") ||
	   !strcmp(varName,"ykpk" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Pk;
  }
  else if( !strcmp(varName,"kypy") ||
	   !strcmp(varName,"ykpy" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Py;
  }
  else if( !strcmp(varName,"kyek") ||
	   !strcmp(varName,"ykek" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ek;
  }
  else if( !strcmp(varName,"kyey") ||
	   !strcmp(varName,"ykey" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ey;
  }
  else if( !strcmp(varName,"kys") ||
	   !strcmp(varName,"yks" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_S;
  }
  else if( !strcmp(varName,"kyw") ||
	   !strcmp(varName,"ykw" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_W;
  }
  else if( !strcmp(varName,"kycosthkcm") ||
	   !strcmp(varName,"ykcosthkcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhkcm;
  }
  else if( !strcmp(varName,"kycosthycm") ||
	   !strcmp(varName,"ykcosthycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhycm;
  }
  else if( !strcmp(varName,"kythkcm") ||
	   !strcmp(varName,"ykthkcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thkcm;
  }
  else if( !strcmp(varName,"kythycm") ||
	   !strcmp(varName,"ykthycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thycm;
  }
  else if( !strcmp(varName,"kyphikcm") ||
	   !strcmp(varName,"ykphikcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_kS;
  }
  else if( !strcmp(varName,"kyphiycm") ||
	   !strcmp(varName,"ykphiycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_kS;
  }
  else if( !strcmp(varName,"kyhpikcm_deg") ||
	   !strcmp(varName,"ykphikcm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_deg_kS;
  }
  else if( !strcmp(varName,"kyphiycm_deg") ||
	   !strcmp(varName,"ykphiycm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_deg_kS;
  }
  else if( !strcmp(varName,"knpk") ||
	   !strcmp(varName,"nkpk" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pk;
  }
  else if( !strcmp(varName,"knpn") ||
	   !strcmp(varName,"nkpn" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pn;
  }
  else if( !strcmp(varName,"knek") ||
	   !strcmp(varName,"nkek" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Ek;
  }
  else if( !strcmp(varName,"knen") ||
	   !strcmp(varName,"nken" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_En;
  }
  else if( !strcmp(varName,"kns") ||
	   !strcmp(varName,"nks" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_S;
  }
  else if( !strcmp(varName,"knw") ||
	   !strcmp(varName,"nkw" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_W;
  }
  else if( !strcmp(varName,"kncosthkcm") ||
	   !strcmp(varName,"nkcosthkcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhkcm;
  }
  else if( !strcmp(varName,"kncosthncm") ||
	   !strcmp(varName,"nkcosthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhncm;
  }
  else if( !strcmp(varName,"knthkcm") ||
	   !strcmp(varName,"nkthkcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thkcm;
  }
  else if( !strcmp(varName,"knthncm") ||
	   !strcmp(varName,"nkthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thncm;
  }
  else if( !strcmp(varName,"knphikcm") ||
	   !strcmp(varName,"nkphikcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_kS;
  }
  else if( !strcmp(varName,"knphincm") ||
	   !strcmp(varName,"nkphincm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_kS;
  }
  else if( !strcmp(varName,"knphikcm_deg") ||
	   !strcmp(varName,"nkphikcm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_deg_kS;
  }
  else if( !strcmp(varName,"knphincm_deg") ||
	   !strcmp(varName,"nkphincm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_deg_kS;
  }
  else if( !strcmp(varName,"ynpy") ||
	   !strcmp(varName,"nypy" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Py;
  }
  else if( !strcmp(varName,"ynpn") ||
	   !strcmp(varName,"nypn" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Pn;
  }
  else if( !strcmp(varName,"yney") ||
	   !strcmp(varName,"nyey" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Ey;
  }
  else if( !strcmp(varName,"ynen") ||
	   !strcmp(varName,"nyen" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_En;
  }
  else if( !strcmp(varName,"yns") ||
	   !strcmp(varName,"nys" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_S;
  }
  else if( !strcmp(varName,"ynw") ||
	   !strcmp(varName,"nyw" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_W;
  }
  else if( !strcmp(varName,"yncosthycm") ||
	   !strcmp(varName,"nycosthycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhycm;
  }
  else if( !strcmp(varName,"yncosthncm") ||
	   !strcmp(varName,"nycosthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhncm;
  }
  else if( !strcmp(varName,"ynthycm") ||
	   !strcmp(varName,"nythycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thycm;
  }
  else if( !strcmp(varName,"ynthncm") ||
	   !strcmp(varName,"nythncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thncm;
  }
  else if( !strcmp(varName,"ynphiycm") ||
	   !strcmp(varName,"nyphiycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_kS;
  }
  else if( !strcmp(varName,"ynphincm") ||
	   !strcmp(varName,"nyphincm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_kS;
  }
  else if( !strcmp(varName,"ynphiycm_deg") ||
	   !strcmp(varName,"nyphiycm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_deg_kS;
  }
  else if( !strcmp(varName,"ynphincm_deg") ||
	   !strcmp(varName,"nyphincm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_deg_kS;
  }
  else {
    std::cerr << "ERROR in  TKinematics2to3::GetVar(const TString&): "
	      << "Kinematic variable \"" << variable << "\" unknown!\n";
    exit(1);
  }

  return (this->*ptrToVarGetter)();
}

//_____________________________________________________________________
double* TKinematics2to3::GetVarArray(const TString& variable) const
{
  // Return an array of the variable of type 'variable'. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.
  //
  // The returned pointer has a size GetNumberOfSteps() and is owned by the user.

  TKinematics2to3 tk(*this); // temporary copy
  double(TKinematics2to3::*ptrToVarGetter)() const; // pointer to TKinematics getters
  
  TString varName = variable;
  varName.ToLower();
  
  // Assign the correct Getter
  if( !strcmp(varName,"pk") ) {
    ptrToVarGetter = &TKinematics2to3::GetPk;
  }
  else if( !strcmp(varName,"py") ) {
    ptrToVarGetter = &TKinematics2to3::GetPy;
  }
  else if( !strcmp(varName,"pn") ) {
    ptrToVarGetter = &TKinematics2to3::GetPn;
  }
  else if( !strcmp(varName,"klab") ) {
    ptrToVarGetter = &TKinematics2to3::GetKlab;
  }
  else if( !strcmp(varName,"ek") ) {
    ptrToVarGetter = &TKinematics2to3::GetEk;
  }
  else if( !strcmp(varName,"ey") ) {
    ptrToVarGetter = &TKinematics2to3::GetEy;
  }
  else if( !strcmp(varName,"en") ) {
    ptrToVarGetter = &TKinematics2to3::GetEn;
  }
  else if( !strcmp(varName,"wlab") ) {
    ptrToVarGetter = &TKinematics2to3::GetWlab;
  }
  else if( !strcmp(varName,"thk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk_deg;
  }
  else if( !strcmp(varName,"thy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy_deg;
  }
  else if( !strcmp(varName,"thn_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn_deg;
  }
  else if( !strcmp(varName,"thk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk;
  }
  else if( !strcmp(varName,"thy") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy;
  }
  else if( !strcmp(varName,"thn") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn;
  }
  else if( !strcmp(varName,"costhk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthk;
  }
  else if( !strcmp(varName,"costhy") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthy;
  }
  else if( !strcmp(varName,"costhn") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthn;
  }
  else if( !strcmp(varName,"phik_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_deg_kS;
  }
  else if( !strcmp(varName,"phiy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_deg_kS;
  }
  else if( !strcmp(varName,"phin_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_deg_kS;
  }
  else if( !strcmp(varName,"phik") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_kS;
  }
  else if( !strcmp(varName,"phiy") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_kS;
  }
  else if( !strcmp(varName,"phin") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_kS;
  }
  else if( !strcmp(varName,"thky_deg") ||
	   !strcmp(varName,"thyk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThky_deg;
  }
  else if( !strcmp(varName,"thkn_deg") ||
	   !strcmp(varName,"thnk_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn_deg;
  }
  else if( !strcmp(varName,"thyn_deg") ||
	   !strcmp(varName,"thny_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn_deg;
  }
  else if( !strcmp(varName,"thky") ||
	   !strcmp(varName,"thyk" )) {
    ptrToVarGetter = &TKinematics2to3::GetThky;
  }
  else if( !strcmp(varName,"thkn") ||
	   !strcmp(varName,"thnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn;
  }
  else if( !strcmp(varName,"thyn") ||
	   !strcmp(varName,"thny") ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn;
  }
  else if( !strcmp(varName,"costhky") ||
	   !strcmp(varName,"costhyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthky;
  }
  else if( !strcmp(varName,"costhkn") ||
	   !strcmp(varName,"costhnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthkn;
  }
  else if( !strcmp(varName,"costhyn") ||
	   !strcmp(varName,"costhny") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthyn;
  }
  else if( !strcmp(varName,"qsquared") ) {
    ptrToVarGetter = &TKinematics2to3::GetQsquared;
  }
  else if( !strcmp(varName,"wtot") ) {
    ptrToVarGetter = &TKinematics2to3::GetWtot;
  }
  else if( !strcmp(varName,"wky") ||
	   !strcmp(varName,"wyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWky;
  }
  else if( !strcmp(varName,"wkn") ||
	   !strcmp(varName,"wnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWkn;
  }
  else if( !strcmp(varName,"wyn") ||
	   !strcmp(varName,"wny") ) {
    ptrToVarGetter = &TKinematics2to3::GetWyn;
  }
  else if( !strcmp(varName,"stot") ) {
    ptrToVarGetter = &TKinematics2to3::GetStot;
  }
  else if( !strcmp(varName,"sky") ||
	   !strcmp(varName,"syk") ) {
    ptrToVarGetter = &TKinematics2to3::GetSky;
  }
  else if( !strcmp(varName,"skn") ||
	   !strcmp(varName,"snk" )) {
    ptrToVarGetter = &TKinematics2to3::GetSkn;
  }
  else if( !strcmp(varName,"syn") ||
	   !strcmp(varName,"sny" )) {
    ptrToVarGetter = &TKinematics2to3::GetSyn;
  }
  else if( !strcmp(varName,"tkg") ||
	   !strcmp(varName,"tgk" )) {
    ptrToVarGetter = &TKinematics2to3::GetTkg;
  }
  else if( !strcmp(varName,"tyg") ||
	   !strcmp(varName,"tgy" )) {
    ptrToVarGetter = &TKinematics2to3::GetTyg;
  }
  else if( !strcmp(varName,"tng") ||
	   !strcmp(varName,"tgn" )) {
    ptrToVarGetter = &TKinematics2to3::GetTng;
  }
  else if( !strcmp(varName,"tkd") ||
	   !strcmp(varName,"tdk" )) {
    ptrToVarGetter = &TKinematics2to3::GetTkd;
  }
  else if( !strcmp(varName,"tyd") ||
	   !strcmp(varName,"tdy" )) {
    ptrToVarGetter = &TKinematics2to3::GetTyd;
  }
  else if( !strcmp(varName,"tnd") ||
	   !strcmp(varName,"tdn" )) {
    ptrToVarGetter = &TKinematics2to3::GetTnd;
  }
  else if( !strcmp(varName,"kypk") ||
	   !strcmp(varName,"ykpk" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Pk;
  }
  else if( !strcmp(varName,"kypy") ||
	   !strcmp(varName,"ykpy" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Py;
  }
  else if( !strcmp(varName,"kyek") ||
	   !strcmp(varName,"ykek" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ek;
  }
  else if( !strcmp(varName,"kyey") ||
	   !strcmp(varName,"ykey" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ey;
  }
  else if( !strcmp(varName,"kys") ||
	   !strcmp(varName,"yks" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_S;
  }
  else if( !strcmp(varName,"kyw") ||
	   !strcmp(varName,"ykw" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_W;
  }
  else if( !strcmp(varName,"kycosthkcm") ||
	   !strcmp(varName,"ykcosthkcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhkcm;
  }
  else if( !strcmp(varName,"kycosthycm") ||
	   !strcmp(varName,"ykcosthycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhycm;
  }
  else if( !strcmp(varName,"kythkcm") ||
	   !strcmp(varName,"ykthkcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thkcm;
  }
  else if( !strcmp(varName,"kythycm") ||
	   !strcmp(varName,"ykthycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thycm;
  }
  else if( !strcmp(varName,"kyphikcm") ||
	   !strcmp(varName,"ykphikcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_kS;
  }
  else if( !strcmp(varName,"kyphiycm") ||
	   !strcmp(varName,"ykphiycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_kS;
  }
  else if( !strcmp(varName,"kyhpikcm_deg") ||
	   !strcmp(varName,"ykphikcm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_deg_kS;
  }
  else if( !strcmp(varName,"kyphiycm_deg") ||
	   !strcmp(varName,"ykphiycm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_deg_kS;
  }
  else if( !strcmp(varName,"knpk") ||
	   !strcmp(varName,"nkpk" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pk;
  }
  else if( !strcmp(varName,"knpn") ||
	   !strcmp(varName,"nkpn" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pn;
  }
  else if( !strcmp(varName,"knek") ||
	   !strcmp(varName,"nkek" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Ek;
  }
  else if( !strcmp(varName,"knen") ||
	   !strcmp(varName,"nken" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_En;
  }
  else if( !strcmp(varName,"kns") ||
	   !strcmp(varName,"nks" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_S;
  }
  else if( !strcmp(varName,"knw") ||
	   !strcmp(varName,"nkw" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_W;
  }
  else if( !strcmp(varName,"kncosthkcm") ||
	   !strcmp(varName,"nkcosthkcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhkcm;
  }
  else if( !strcmp(varName,"kncosthncm") ||
	   !strcmp(varName,"nkcosthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhncm;
  }
  else if( !strcmp(varName,"knthkcm") ||
	   !strcmp(varName,"nkthkcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thkcm;
  }
  else if( !strcmp(varName,"knthncm") ||
	   !strcmp(varName,"nkthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thncm;
  }
  else if( !strcmp(varName,"knphikcm") ||
	   !strcmp(varName,"nkphikcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_kS;
  }
  else if( !strcmp(varName,"knphincm") ||
	   !strcmp(varName,"nkphincm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_kS;
  }
  else if( !strcmp(varName,"knphikcm_deg") ||
	   !strcmp(varName,"nkphikcm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_deg_kS;
  }
  else if( !strcmp(varName,"knphincm_deg") ||
	   !strcmp(varName,"nkphincm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_deg_kS;
  }
  else if( !strcmp(varName,"ynpy") ||
	   !strcmp(varName,"nypy" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Py;
  }
  else if( !strcmp(varName,"ynpn") ||
	   !strcmp(varName,"nypn" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Pn;
  }
  else if( !strcmp(varName,"yney") ||
	   !strcmp(varName,"nyey" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Ey;
  }
  else if( !strcmp(varName,"ynen") ||
	   !strcmp(varName,"nyen" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_En;
  }
  else if( !strcmp(varName,"yns") ||
	   !strcmp(varName,"nys" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_S;
  }
  else if( !strcmp(varName,"ynw") ||
	   !strcmp(varName,"nyw" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_W;
  }
  else if( !strcmp(varName,"yncosthycm") ||
	   !strcmp(varName,"nycosthycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhycm;
  }
  else if( !strcmp(varName,"yncosthncm") ||
	   !strcmp(varName,"nycosthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhncm;
  }
  else if( !strcmp(varName,"ynthycm") ||
	   !strcmp(varName,"nythycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thycm;
  }
  else if( !strcmp(varName,"ynthncm") ||
	   !strcmp(varName,"nythncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thncm;
  }
  else if( !strcmp(varName,"ynphiycm") ||
	   !strcmp(varName,"nyphiycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_kS;
  }
  else if( !strcmp(varName,"ynphincm") ||
	   !strcmp(varName,"nyphincm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_kS;
  }
  else if( !strcmp(varName,"ynphiycm_deg") ||
	   !strcmp(varName,"nyphiycm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_deg_kS;
  }
  else if( !strcmp(varName,"ynphincm_deg") ||
	   !strcmp(varName,"nyphincm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_deg_kS;
  }
  else {
    std::cerr << "ERROR in  TKinematics2to3::GetVarArray(TString): "
	      << "Kinematic variable \"" << variable << "\" unknown!\n";
    exit(1);
  }

  double *array = new double[tk.GetNumberOfSteps()];
  
  // Fill the array
  for(int step=0; step<tk.GetNumberOfSteps(); ++step) {
    tk.GoTo(step);
    array[step] = (tk.*ptrToVarGetter)();
  }    

  return array;
}

//_____________________________________________________________________
double* TKinematics2to3::GetPhysicalVarArray(const TString& variable) const
{
  // Return an array of the variable of type 'variable'. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.
  //
  // The returned pointer has a size GetNumberOfPhysicalSteps() and is owned by the user.

  TKinematics2to3 tk(*this); // temporary copy
  double(TKinematics2to3::*ptrToVarGetter)() const; // pointer to TKinematics getters

  TString varName = variable;
  varName.ToLower();

  // Assign the correct Getter
  if( !strcmp(varName,"pk") ) {
    ptrToVarGetter = &TKinematics2to3::GetPk;
  }
  else if( !strcmp(varName,"py") ) {
    ptrToVarGetter = &TKinematics2to3::GetPy;
  }
  else if( !strcmp(varName,"pn") ) {
    ptrToVarGetter = &TKinematics2to3::GetPn;
  }
  else if( !strcmp(varName,"klab") ) {
    ptrToVarGetter = &TKinematics2to3::GetKlab;
  }
  else if( !strcmp(varName,"ek") ) {
    ptrToVarGetter = &TKinematics2to3::GetEk;
  }
  else if( !strcmp(varName,"ey") ) {
    ptrToVarGetter = &TKinematics2to3::GetEy;
  }
  else if( !strcmp(varName,"en") ) {
    ptrToVarGetter = &TKinematics2to3::GetEn;
  }
  else if( !strcmp(varName,"wlab") ) {
    ptrToVarGetter = &TKinematics2to3::GetWlab;
  }
  else if( !strcmp(varName,"thk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk_deg;
  }
  else if( !strcmp(varName,"thy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy_deg;
  }
  else if( !strcmp(varName,"thn_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn_deg;
  }
  else if( !strcmp(varName,"thk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk;
  }
  else if( !strcmp(varName,"thy") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy;
  }
  else if( !strcmp(varName,"thn") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn;
  }
  else if( !strcmp(varName,"costhk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthk;
  }
  else if( !strcmp(varName,"costhy") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthy;
  }
  else if( !strcmp(varName,"costhn") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthn;
  }
  else if( !strcmp(varName,"phik_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_deg_kS;
  }
  else if( !strcmp(varName,"phiy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_deg_kS;
  }
  else if( !strcmp(varName,"phin_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_deg_kS;
  }
  else if( !strcmp(varName,"phik") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_kS;
  }
  else if( !strcmp(varName,"phiy") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_kS;
  }
  else if( !strcmp(varName,"phin") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_kS;
  }
  else if( !strcmp(varName,"thky_deg") ||
	   !strcmp(varName,"thyk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThky_deg;
  }
  else if( !strcmp(varName,"thkn_deg") ||
	   !strcmp(varName,"thnk_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn_deg;
  }
  else if( !strcmp(varName,"thyn_deg") ||
	   !strcmp(varName,"thny_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn_deg;
  }
  else if( !strcmp(varName,"thky") ||
	   !strcmp(varName,"thyk" )) {
    ptrToVarGetter = &TKinematics2to3::GetThky;
  }
  else if( !strcmp(varName,"thkn") ||
	   !strcmp(varName,"thnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn;
  }
  else if( !strcmp(varName,"thyn") ||
	   !strcmp(varName,"thny") ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn;
  }
  else if( !strcmp(varName,"costhky") ||
	   !strcmp(varName,"costhyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthky;
  }
  else if( !strcmp(varName,"costhkn") ||
	   !strcmp(varName,"costhnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthkn;
  }
  else if( !strcmp(varName,"costhyn") ||
	   !strcmp(varName,"costhny") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthyn;
  }
  else if( !strcmp(varName,"qsquared") ) {
    ptrToVarGetter = &TKinematics2to3::GetQsquared;
  }
  else if( !strcmp(varName,"wtot") ) {
    ptrToVarGetter = &TKinematics2to3::GetWtot;
  }
  else if( !strcmp(varName,"wky") ||
	   !strcmp(varName,"wyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWky;
  }
  else if( !strcmp(varName,"wkn") ||
	   !strcmp(varName,"wnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWkn;
  }
  else if( !strcmp(varName,"wyn") ||
	   !strcmp(varName,"wny") ) {
    ptrToVarGetter = &TKinematics2to3::GetWyn;
  }
  else if( !strcmp(varName,"stot") ) {
    ptrToVarGetter = &TKinematics2to3::GetStot;
  }
  else if( !strcmp(varName,"sky") ||
	   !strcmp(varName,"syk") ) {
    ptrToVarGetter = &TKinematics2to3::GetSky;
  }
  else if( !strcmp(varName,"skn") ||
	   !strcmp(varName,"snk" )) {
    ptrToVarGetter = &TKinematics2to3::GetSkn;
  }
  else if( !strcmp(varName,"syn") ||
	   !strcmp(varName,"sny" )) {
    ptrToVarGetter = &TKinematics2to3::GetSyn;
  }
  else if( !strcmp(varName,"tkg") ||
	   !strcmp(varName,"tgk" )) {
    ptrToVarGetter = &TKinematics2to3::GetTkg;
  }
  else if( !strcmp(varName,"tyg") ||
	   !strcmp(varName,"tgy" )) {
    ptrToVarGetter = &TKinematics2to3::GetTyg;
  }
  else if( !strcmp(varName,"tng") ||
	   !strcmp(varName,"tgn" )) {
    ptrToVarGetter = &TKinematics2to3::GetTng;
  }
  else if( !strcmp(varName,"tkd") ||
	   !strcmp(varName,"tdk" )) {
    ptrToVarGetter = &TKinematics2to3::GetTkd;
  }
  else if( !strcmp(varName,"tyd") ||
	   !strcmp(varName,"tdy" )) {
    ptrToVarGetter = &TKinematics2to3::GetTyd;
  }
  else if( !strcmp(varName,"tnd") ||
	   !strcmp(varName,"tdn" )) {
    ptrToVarGetter = &TKinematics2to3::GetTnd;
  }
  else if( !strcmp(varName,"kypk") ||
	   !strcmp(varName,"ykpk" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Pk;
  }
  else if( !strcmp(varName,"kypy") ||
	   !strcmp(varName,"ykpy" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Py;
  }
  else if( !strcmp(varName,"kyek") ||
	   !strcmp(varName,"ykek" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ek;
  }
  else if( !strcmp(varName,"kyey") ||
	   !strcmp(varName,"ykey" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ey;
  }
  else if( !strcmp(varName,"kys") ||
	   !strcmp(varName,"yks" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_S;
  }
  else if( !strcmp(varName,"kyw") ||
	   !strcmp(varName,"ykw" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_W;
  }
  else if( !strcmp(varName,"kycosthkcm") ||
	   !strcmp(varName,"ykcosthkcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhkcm;
  }
  else if( !strcmp(varName,"kycosthycm") ||
	   !strcmp(varName,"ykcosthycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhycm;
  }
  else if( !strcmp(varName,"kythkcm") ||
	   !strcmp(varName,"ykthkcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thkcm;
  }
  else if( !strcmp(varName,"kythycm") ||
	   !strcmp(varName,"ykthycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thycm;
  }
  else if( !strcmp(varName,"kyphikcm") ||
	   !strcmp(varName,"ykphikcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_kS;
  }
  else if( !strcmp(varName,"kyphiycm") ||
	   !strcmp(varName,"ykphiycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_kS;
  }
  else if( !strcmp(varName,"kyhpikcm_deg") ||
	   !strcmp(varName,"ykphikcm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_deg_kS;
  }
  else if( !strcmp(varName,"kyphiycm_deg") ||
	   !strcmp(varName,"ykphiycm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_deg_kS;
  }
  else if( !strcmp(varName,"knpk") ||
	   !strcmp(varName,"nkpk" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pk;
  }
  else if( !strcmp(varName,"knpn") ||
	   !strcmp(varName,"nkpn" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pn;
  }
  else if( !strcmp(varName,"knek") ||
	   !strcmp(varName,"nkek" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Ek;
  }
  else if( !strcmp(varName,"knen") ||
	   !strcmp(varName,"nken" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_En;
  }
  else if( !strcmp(varName,"kns") ||
	   !strcmp(varName,"nks" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_S;
  }
  else if( !strcmp(varName,"knw") ||
	   !strcmp(varName,"nkw" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_W;
  }
  else if( !strcmp(varName,"kncosthkcm") ||
	   !strcmp(varName,"nkcosthkcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhkcm;
  }
  else if( !strcmp(varName,"kncosthncm") ||
	   !strcmp(varName,"nkcosthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhncm;
  }
  else if( !strcmp(varName,"knthkcm") ||
	   !strcmp(varName,"nkthkcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thkcm;
  }
  else if( !strcmp(varName,"knthncm") ||
	   !strcmp(varName,"nkthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thncm;
  }
  else if( !strcmp(varName,"knphikcm") ||
	   !strcmp(varName,"nkphikcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_kS;
  }
  else if( !strcmp(varName,"knphincm") ||
	   !strcmp(varName,"nkphincm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_kS;
  }
  else if( !strcmp(varName,"knphikcm_deg") ||
	   !strcmp(varName,"nkphikcm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_deg_kS;
  }
  else if( !strcmp(varName,"knphincm_deg") ||
	   !strcmp(varName,"nkphincm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_deg_kS;
  }
  else if( !strcmp(varName,"ynpy") ||
	   !strcmp(varName,"nypy" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Py;
  }
  else if( !strcmp(varName,"ynpn") ||
	   !strcmp(varName,"nypn" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Pn;
  }
  else if( !strcmp(varName,"yney") ||
	   !strcmp(varName,"nyey" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Ey;
  }
  else if( !strcmp(varName,"ynen") ||
	   !strcmp(varName,"nyen" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_En;
  }
  else if( !strcmp(varName,"yns") ||
	   !strcmp(varName,"nys" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_S;
  }
  else if( !strcmp(varName,"ynw") ||
	   !strcmp(varName,"nyw" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_W;
  }
  else if( !strcmp(varName,"yncosthycm") ||
	   !strcmp(varName,"nycosthycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhycm;
  }
  else if( !strcmp(varName,"yncosthncm") ||
	   !strcmp(varName,"nycosthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhncm;
  }
  else if( !strcmp(varName,"ynthycm") ||
	   !strcmp(varName,"nythycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thycm;
  }
  else if( !strcmp(varName,"ynthncm") ||
	   !strcmp(varName,"nythncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thncm;
  }
  else if( !strcmp(varName,"ynphiycm") ||
	   !strcmp(varName,"nyphiycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_kS;
  }
  else if( !strcmp(varName,"ynphincm") ||
	   !strcmp(varName,"nyphincm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_kS;
  }
  else if( !strcmp(varName,"ynphiycm_deg") ||
	   !strcmp(varName,"nyphiycm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_deg_kS;
  }
  else if( !strcmp(varName,"ynphincm_deg") ||
	   !strcmp(varName,"nyphincm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_deg_kS;
  }
  else {
    std::cerr << "ERROR in  TKinematics2to3::GetVarArray(TString): "
	      << "Kinematic variable \"" << variable << "\" unknown!\n";
    exit(1);
  }

  double *array = new double[tk.GetNumberOfPhysicalSteps()];
  
  // Fill the array
  int physicalStep = 0;
  for(int step=0; step<tk.GetNumberOfSteps(); ++step) {
    tk.GoTo(step);
    if(tk.IsPhysical()) array[physicalStep++] = (tk.*ptrToVarGetter)();
  }    

  return array;
}

//_____________________________________________________________________
TString TKinematics2to3::GetVarName(int varNumber) const
{

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetVarName(int):\n"
	      << "Index out-of-bounds!\n";
    return TString();
  }

  if( !fVar ) {
    std::cerr << "ERROR in TKinematics2to3::GetVarName(int):\n"
	      << "Independent kinematic variables have not been initialized!\n";
    return TString();
  }

  TString varName;
  
  // Lookup name
  if( fVar[varNumber-1] == &fQsquared )
    varName = "qsquared";
  else if( fVar[varNumber-1] == &fWlab )
    varName = "wlab";
  else if( fVar[varNumber-1] == &fKlab )
    varName = "klab";
  else if( fVar[varNumber-1] == &fEk )
    varName = "ek";
  else if( fVar[varNumber-1] == &fEy )
    varName = "ey";
  else if( fVar[varNumber-1] == &fEn )
    varName = "en";
  else if( fVar[varNumber-1] == &fPk )
    varName = "pk"; 
  else if( fVar[varNumber-1] == &fPy )
    varName = "py"; 
  else if( fVar[varNumber-1] == &fPn )
    varName = "pn"; 
  else if( fVar[varNumber-1] == &fThk )
    varName = "thk";   
  else if( fVar[varNumber-1] == &fThy )
    varName = "thy";   
  else if( fVar[varNumber-1] == &fThn )
    varName = "thn";   
  else if( fVar[varNumber-1] == &fCosthk )
    varName = "costhk";
  else if( fVar[varNumber-1] == &fCosthy )
    varName = "costhy";
  else if( fVar[varNumber-1] == &fCosthn )
    varName = "costhn";
  else if( fVar[varNumber-1] == &fPhik )
    varName = "phik";
  else if( fVar[varNumber-1] == &fPhiy )
    varName = "phiy";  
  else if( fVar[varNumber-1] == &fPhin )
    varName = "phin";  
  else if( fVar[varNumber-1] == &fWtot )
    varName = "wtot";  
  else if( fVar[varNumber-1] == &fWky )
    varName = "wky";   
  else if( fVar[varNumber-1] == &fWkn )
    varName = "wkn";   
  else if( fVar[varNumber-1] == &fWyn )
    varName = "wyn";   
  else if( fVar[varNumber-1] == &fStot )
    varName = "stot";  
  else if( fVar[varNumber-1] == &fSky )
    varName = "sky";   
  else if( fVar[varNumber-1] == &fSkn )
    varName = "skn";   
  else if( fVar[varNumber-1] == &fSyn )
    varName = "syn";   
  else if( fVar[varNumber-1] == &fTgk )
    varName = "tgk";   
  else if( fVar[varNumber-1] == &fTgy )
    varName = "tgy";   
  else if( fVar[varNumber-1] == &fTgn )
    varName = "tgn";   
  else if( fVar[varNumber-1] == &fTdk )
    varName = "tdk";   
  else if( fVar[varNumber-1] == &fTdy )
    varName = "tdy";   
  else if( fVar[varNumber-1] == &fTdn )
    varName = "tdn";   
  else if( fVar[varNumber-1] == &fKYthkcm )
    varName = "kythkcm";
  else if( fVar[varNumber-1] == &fKYthycm )
    varName = "kythycm";
  else if( fVar[varNumber-1] == &fKYcosthkcm )
    varName = "kycosthkcm";
  else if( fVar[varNumber-1] == &fKYcosthycm )
    varName = "kycosthycm";
  else if( fVar[varNumber-1] == &fKYphikcm )
    varName = "kyphikcm";
  else if( fVar[varNumber-1] == &fKYphiycm )
    varName = "kyphiycm";  
  else if( fVar[varNumber-1] == &fKNthkcm )
    varName = "knthkcm";   
  else if( fVar[varNumber-1] == &fKNthncm )
    varName = "knthncm";   
  else if( fVar[varNumber-1] == &fKNcosthkcm )
    varName = "kncosthkcm";
  else if( fVar[varNumber-1] == &fKNcosthncm )
    varName = "kncosthncm";
  else if( fVar[varNumber-1] == &fKNphikcm )
    varName = "knphikcm";  
  else if( fVar[varNumber-1] == &fKNphincm )
    varName = "knphincm";  
  else if( fVar[varNumber-1] == &fYNthycm )
    varName = "ynthycm";   
  else if( fVar[varNumber-1] == &fYNthncm )
    varName = "ynthncm";   
  else if( fVar[varNumber-1] == &fYNcosthycm )
    varName = "yncosthycm";
  else if( fVar[varNumber-1] == &fYNcosthncm )
    varName = "yncosthncm";
  else if( fVar[varNumber-1] == &fYNphiycm )
    varName = "ynphiycm";  
  else if( fVar[varNumber-1] == &fYNphincm )
    varName = "ynphincm";  
  else {
    std::cerr << "Error in TKinematics2to3::GetVarName(int):\n"
	      << "TKinematics2to3::fVar[" << varNumber-1 << "] "
	      << "points to unknown variable.\n";
  }

  return varName;
}

//_____________________________________________________________________
int TKinematics2to3::GetNumberOfSteps() const
{
  return fNumberOfSteps[0] * fNumberOfSteps[1] * fNumberOfSteps[2]
    * fNumberOfSteps[3] * fNumberOfSteps[4] * fNumberOfSteps[5];
}

//_____________________________________________________________________
int TKinematics2to3::GetNumberOfPhysicalSteps() const
{
  if( fNrOfPhysicals==-1 ) {
    TKinematics2to3 tk(*this); // temporary copy
    fNrOfPhysicals = 0;        // final result
    
    for(int step=0; step<tk.GetNumberOfSteps(); ++step) {
      tk.GoTo(step);
      if(tk.IsPhysical()) ++fNrOfPhysicals;
    }
  }

  return fNrOfPhysicals;
}    

//_____________________________________________________________________
int TKinematics2to3::GetNumberOfSteps(int varNumber) const
{
  
  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetNumberOfSteps(int):\n"
	      << "Index out-of-bounds!\n";
    return 0;
  }

  return fNumberOfSteps[varNumber-1];
}

//_____________________________________________________________________
int TKinematics2to3::GetNumberOfVariables() const
{
  // Determine how many independent variables are not fixed

  int nrOfVar = 0;
  for(int var=1; var<=6; ++var)
    if( !IsFixed(var) ) ++nrOfVar;

  return nrOfVar;
}

//_____________________________________________________________________
double TKinematics2to3::GetStepSize(int varNumber) const
{

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetStepSize(int):\n"
	      << "Index out-of-bounds!\n";
    return 0;
  }

  return fStepSize[varNumber-1];
}

//_____________________________________________________________________
int TKinematics2to3::GetStep(int varNumber) const
{

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetStep(int):\n"
	      << "Index out-of-bounds!\n";
    return -1;
  }

  // Determine step per variable
  int stepVar[6];
  int remainder=fStep;
  int denum=1;
  
  for(int i=6; i>0; --i)
    {
      denum = 1;
      for(int j=0; j<i-1; ++j)
	denum *= fNumberOfSteps[j];
      
      stepVar[i-1] = remainder / denum;
      remainder = remainder % denum;
    }
    
  return stepVar[varNumber-1];  
}

//_____________________________________________________________________
bool TKinematics2to3::IsFixed() const
{
  return !fIsVar[0] && !fIsVar[1] && !fIsVar[2]
    && !fIsVar[3] && !fIsVar[4] && !fIsVar[5];
}

//_____________________________________________________________________
bool TKinematics2to3::IsFixed(int varNumber) const
{

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::IsFixed(int):\n"
	      << "Index out-of-bounds!\n";
    return true;
  }

  return !fIsVar[varNumber-1];
}

//______________________________________________________________________________
void TKinematics2to3::Streamer(TBuffer &R__b)
{
  // Stream an object of class TKinematics2to3.

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
    delete [] fKinSolution;
    fKinSolution = new char[4];
    R__b.ReadFastArray(fKinSolution,4);
    R__b >> fIsospin;
    R__b >> (Int_t&)fAC;
    R__b >> fMd;
    R__b >> fMk;
    R__b >> fMy;
    R__b >> fMn;
    R__b >> fIsPhysical;
    R__b >> fQsquared;
    R__b >> fWlab;
    R__b >> fKlab;
    R__b >> fEk;
    R__b >> fEy;
    R__b >> fEn;
    R__b >> fPk;
    R__b >> fPy;
    R__b >> fPn;
    R__b >> fThk;
    R__b >> fThy;
    R__b >> fThn;
    R__b >> fCosthk;
    R__b >> fCosthy;
    R__b >> fCosthn;
    R__b >> fPhik;
    R__b >> fPhiy;
    R__b >> fPhin;
    R__b >> fWtot;
    R__b >> fWky;
    R__b >> fWkn;
    R__b >> fWyn;
    R__b >> fStot;
    R__b >> fSky;
    R__b >> fSkn;
    R__b >> fSyn;
    R__b >> fTgk;
    R__b >> fTgy;
    R__b >> fTgn;
    R__b >> fTdk;
    R__b >> fTdy;
    R__b >> fTdn;
    R__b >> fKYthkcm;
    R__b >> fKYthycm;
    R__b >> fKYcosthkcm;
    R__b >> fKYcosthycm;
    R__b >> fKYphikcm;
    R__b >> fKYphiycm;
    R__b >> fKNthkcm;
    R__b >> fKNthncm;
    R__b >> fKNcosthkcm;
    R__b >> fKNcosthncm;
    R__b >> fKNphikcm;
    R__b >> fKNphincm;
    R__b >> fYNthycm;
    R__b >> fYNthncm;
    R__b >> fYNcosthycm;
    R__b >> fYNcosthncm;
    R__b >> fYNphiycm;
    R__b >> fYNphincm;
    R__b >> fKYsystem;
    R__b >> fKNsystem;
    R__b >> fYNsystem;
    R__b >> fKvector;
    R__b >> fYvector;
    R__b >> fNvector;
    R__b >> fGvector;
    R__b >> fDvector;
    delete [] fIsVar;
    fIsVar = new bool[6];
    R__b.ReadFastArray(fIsVar,6);
    delete [] fLowerLimit;
    fLowerLimit = new double[6];
    R__b.ReadFastArray(fLowerLimit,6);
    delete [] fUpperLimit;
    fUpperLimit = new double[6];
    R__b.ReadFastArray(fUpperLimit,6);
    delete [] fStepSize;
    fStepSize = new double[6];
    R__b.ReadFastArray(fStepSize,6);
    delete [] fNumberOfSteps;
    fNumberOfSteps = new int[6];
    R__b.ReadFastArray(fNumberOfSteps,6);
    R__b >> fStep;
    R__b >> fNrOfPhysicals;
    R__b.CheckByteCount(R__s, R__c, TKinematics2to3::IsA());
    fVar    = new double*[6];  // Added by PVC
    ReadFormatString(tempFormat); // Added by PVC
    delete[] tempFormat;
  } else {
    R__c = R__b.WriteVersion(TKinematics2to3::IsA(), kTRUE);
    TNamed::Streamer(R__b);
    R__b.WriteFastArray(fFormat,100);
    R__b.WriteFastArray(fKinSolution,4);
    R__b << fIsospin;
    R__b << (Int_t)fAC;
    R__b << fMd;
    R__b << fMk;
    R__b << fMy;
    R__b << fMn;
    R__b << fIsPhysical;
    R__b << fQsquared;
    R__b << fWlab;
    R__b << fKlab;
    R__b << fEk;
    R__b << fEy;
    R__b << fEn;
    R__b << fPk;
    R__b << fPy;
    R__b << fPn;
    R__b << fThk;
    R__b << fThy;
    R__b << fThn;
    R__b << fCosthk;
    R__b << fCosthy;
    R__b << fCosthn;
    R__b << fPhik;
    R__b << fPhiy;
    R__b << fPhin;
    R__b << fWtot;
    R__b << fWky;
    R__b << fWkn;
    R__b << fWyn;
    R__b << fStot;
    R__b << fSky;
    R__b << fSkn;
    R__b << fSyn;
    R__b << fTgk;
    R__b << fTgy;
    R__b << fTgn;
    R__b << fTdk;
    R__b << fTdy;
    R__b << fTdn;
    R__b << fKYthkcm;
    R__b << fKYthycm;
    R__b << fKYcosthkcm;
    R__b << fKYcosthycm;
    R__b << fKYphikcm;
    R__b << fKYphiycm;
    R__b << fKNthkcm;
    R__b << fKNthncm;
    R__b << fKNcosthkcm;
    R__b << fKNcosthncm;
    R__b << fKNphikcm;
    R__b << fKNphincm;
    R__b << fYNthycm;
    R__b << fYNthncm;
    R__b << fYNcosthycm;
    R__b << fYNcosthncm;
    R__b << fYNphiycm;
    R__b << fYNphincm;
    R__b << fKYsystem;
    R__b << fKNsystem;
    R__b << fYNsystem;
    R__b << fKvector;
    R__b << fYvector;
    R__b << fNvector;
    R__b << fGvector;
    R__b << fDvector;
    R__b.WriteFastArray(fIsVar,6);
    R__b.WriteFastArray(fLowerLimit,6);
    R__b.WriteFastArray(fUpperLimit,6);
    R__b.WriteFastArray(fStepSize,6);
    R__b.WriteFastArray(fNumberOfSteps,6);
    R__b << fStep;
    R__b << fNrOfPhysicals;
    R__b.SetByteCount(R__c, kTRUE);
  }
}

//_____________________________________________________________________
TVector3 TKinematics2to3::GetK3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the kaon in the LAB.

  if( axisConvention == kS )
    return GetK3(fAC);

  else if( axisConvention == kY ) {
    TVector3 rotatedKvector = fKvector->Vect();
    rotatedKvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedKvector;
  }
  
  else if( axisConvention == kK ) {
    TVector3 rotatedKvector = fKvector->Vect();
    rotatedKvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedKvector;
  }

  return fKvector->Vect();
}

//_____________________________________________________________________
TVector3 TKinematics2to3::GetY3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the hyperon in the LAB.

  if( axisConvention == kS )
    return GetY3(fAC);

  else if( axisConvention == kY ) {
    TVector3 rotatedYvector = fYvector->Vect();
    rotatedYvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedYvector;
  }
  
  else if( axisConvention == kK ) {
    TVector3 rotatedYvector = fYvector->Vect();
    rotatedYvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedYvector;
  }

  return fYvector->Vect();
}

//_____________________________________________________________________
TVector3 TKinematics2to3::GetN3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the nucleon in the LAB.

  if( axisConvention == kS )
    return GetN3(fAC);

  else if( axisConvention == kY ) {
    TVector3 rotatedNvector = fNvector->Vect();
    rotatedNvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedNvector;
  }
  
  else if( axisConvention == kK ) {
    TVector3 rotatedNvector = fNvector->Vect();
    rotatedNvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedNvector;
  }

  return fNvector->Vect();
}

//_____________________________________________________________________
TVector3 TKinematics2to3::GetG3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the photon in the LAB.
  return fGvector->Vect();
}

//_____________________________________________________________________
TVector3 TKinematics2to3::GetD3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the deuteron in the LAB.
  return fDvector->Vect();
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::GetK4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the kaon in the LAB.

  if( axisConvention == kS )
    return GetK4(fAC);

  else if( axisConvention == kY ) {
    TLorentzVector rotatedKvector = *fKvector;
    rotatedKvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedKvector;
  }
  
  else if( axisConvention == kK ) {
    TLorentzVector rotatedKvector = *fKvector;
    rotatedKvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedKvector;
  }

  return *fKvector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::GetY4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the hyperon in the LAB.

  if( axisConvention == kS )
    return GetY4(fAC);

  else if( axisConvention == kY ) {
    TLorentzVector rotatedYvector = *fYvector;
    rotatedYvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedYvector;
  }
  
  else if( axisConvention == kK ) {
    TLorentzVector rotatedYvector = *fYvector;
    rotatedYvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedYvector;
  }

  return *fYvector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::GetN4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the nucleon in the LAB.

  if( axisConvention == kS )
    return GetN4(fAC);

  else if( axisConvention == kY ) {
    TLorentzVector rotatedNvector = *fNvector;
    rotatedNvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedNvector;
  }
  
  else if( axisConvention == kK ) {
    TLorentzVector rotatedNvector = *fNvector;
    rotatedNvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedNvector;
  }

  return *fNvector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::GetG4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the photon in the LAB.
  return *fGvector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::GetD4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the deuteron in the LAB
  return *fDvector;
}

//_____________________________________________________________________
double TKinematics2to3::GetStot() const
{
  return fMd*fMd + 2*GetWlab()*fMd - GetQsquared();
}

//_____________________________________________________________________
double TKinematics2to3::GetThk() const 
{ 
  return fKvector->Angle( fGvector->Vect() ); 
}
  
//_____________________________________________________________________
double TKinematics2to3::GetThy() const 
{ 
  return fYvector->Angle( fGvector->Vect() ); 
}
 
//_____________________________________________________________________
double TKinematics2to3::GetThn() const 
{ 
  return fNvector->Angle( fGvector->Vect() ); 
}

//_____________________________________________________________________
double TKinematics2to3::GetThk_deg() const 
{ 
  return GetThk() * TMath::RadToDeg(); 
}
  
//_____________________________________________________________________
double TKinematics2to3::GetThy_deg() const 
{ 
  return GetThy() * TMath::RadToDeg(); 
}
 
//_____________________________________________________________________
double TKinematics2to3::GetThn_deg() const 
{ 
  return GetThn() * TMath::RadToDeg(); 
}

//_____________________________________________________________________
double TKinematics2to3::GetPhik(EAxisConvention axisConvention) const
{
  if( axisConvention == kS )
    GetPhik(fAC);

  else if( axisConvention == kN )
    return fKvector->Vect().XYvector().Phi();

  else if( axisConvention == kK )
    return 0.0;

  // axisConvention == kY
  return TVector2::Phi_0_2pi( GetPhik(kN) - GetPhiy(kN) );
}

//_____________________________________________________________________
double TKinematics2to3::GetPhiy(EAxisConvention axisConvention) const
{
  if( axisConvention == kS )
    return GetPhiy(fAC);

  else if( axisConvention == kN )
    return fYvector->Vect().XYvector().Phi();

  else if( axisConvention == kY )
    return 0.0;

  // axisConvention == kK
  return TVector2::Phi_0_2pi( GetPhiy(kN) - GetPhik(kN) );
}

//_____________________________________________________________________
double TKinematics2to3::GetPhin(EAxisConvention axisConvention) const
{
  if( axisConvention == kS )
    return GetPhin(fAC);

  else if( axisConvention == kN )
    return 0.0;

  else if( axisConvention == kK )
    return TVector2::Phi_0_2pi( -1.0*GetPhik(kN) );

  // axisConvention == kY
  return TVector2::Phi_0_2pi( -1.0*GetPhiy(kN) );
}

//_____________________________________________________________________
double TKinematics2to3::GetPhik_deg(EAxisConvention axisConvention) const
{
  return GetPhik(axisConvention) * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetPhiy_deg(EAxisConvention axisConvention) const
{
  return GetPhiy(axisConvention) * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetPhin_deg(EAxisConvention axisConvention) const
{
  return GetPhin(axisConvention) * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetThky() const
{
  return fGvector->Angle( fKvector->Vect() + fYvector->Vect() );
}

//_____________________________________________________________________
double TKinematics2to3::GetThkn() const
{
  return fGvector->Angle( fKvector->Vect() + fNvector->Vect() );
}

//_____________________________________________________________________
double TKinematics2to3::GetThyn() const
{
  return fGvector->Angle( fYvector->Vect() + fNvector->Vect() );
}

//_____________________________________________________________________
double TKinematics2to3::GetThky_deg() const
{
  return GetThky() * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetThkn_deg() const
{
  return GetThkn() * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetThyn_deg() const
{
  return GetThyn() * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetTkg() const
{
  return ( *fKvector - *fGvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTyg() const
{
  return ( *fYvector - *fGvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTng() const
{
  return ( *fNvector - *fGvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTkd() const
{
  return ( *fKvector - *fDvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTyd() const
{
  return ( *fYvector - *fDvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTnd() const
{
  return ( *fNvector - *fDvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_Ek() const
{
  return Sqrt( fMk*fMk + fKYsystem->GetPk()*fKYsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_Ey() const
{
  return Sqrt( fMy*fMy + fKYsystem->GetPk()*fKYsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_Ek() const
{
  return Sqrt( fMk*fMk + fKNsystem->GetPk()*fKNsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_En() const
{
  return Sqrt( fMn*fMn + fKNsystem->GetPk()*fKNsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_Ey() const
{
  return Sqrt( fMy*fMy + fYNsystem->GetPk()*fYNsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_En() const
{
  return Sqrt( fMn*fMn + fYNsystem->GetPk()*fYNsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_costhkcm() const
{
  return fKYsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_thkcm() const
{
  return ACos(fKYsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_costhycm() const
{
  return -1.0 * fKYsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_thycm() const
{
  return ACos(-1.0 * fKYsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_costhkcm() const
{
  return fKNsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_thkcm() const
{
  return ACos(fKNsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_costhncm() const
{
  return -1.0 * fKNsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_thncm() const
{
  return ACos(-1.0 * fKNsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_costhycm() const
{
  return -1.0 * fYNsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_thycm() const
{
  return ACos(-1.0 * fYNsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_costhncm() const
{
  return fYNsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_thncm() const
{
  return ACos(fYNsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_phikcm(EAxisConvention ac) const
{
  // Transform Kvector to the KY cm 
  TLorentzVector vectorKcm = BoostToKYcm(GetK4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorKcm.X())<STRANGEUFLOW && fabs(vectorKcm.Y())<STRANGEUFLOW)
    return 0.;
  
  return (vectorKcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_phiycm(EAxisConvention ac) const
{
  // Transform Yvector to the KY cm 
  TLorentzVector vectorYcm = BoostToKYcm(GetY4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorYcm.X())<STRANGEUFLOW && fabs(vectorYcm.Y())<STRANGEUFLOW)
    return 0.;
  
  return (vectorYcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_phikcm(EAxisConvention ac) const
{
  // Transform Kvector to the KN cm 
  TLorentzVector vectorKcm = BoostToKNcm(GetK4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorKcm.X())<STRANGEUFLOW && fabs(vectorKcm.Y())<STRANGEUFLOW)
    return 0.;
  
  return (vectorKcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_phincm(EAxisConvention ac) const
{
  // Transform Nvector to the KN cm 
  TLorentzVector vectorNcm = BoostToKNcm(GetN4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorNcm.X())<STRANGEUFLOW && fabs(vectorNcm.Y())<STRANGEUFLOW)
    return 0.;
  
  return (vectorNcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_phiycm(EAxisConvention ac) const
{
  // Transform Yvector to the YN cm 
  TLorentzVector vectorYcm = BoostToYNcm(GetY4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorYcm.X())<STRANGEUFLOW && fabs(vectorYcm.Y())<STRANGEUFLOW)
    return 0.;
  
  return (vectorYcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_phincm(EAxisConvention ac) const
{
  // Transform Nvector to the YN cm 
  TLorentzVector vectorNcm = BoostToYNcm(GetN4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorNcm.X())<STRANGEUFLOW && fabs(vectorNcm.Y())<STRANGEUFLOW)
    return 0.;
  
  return (vectorNcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_phikcm_deg(EAxisConvention ac) const
{
  return GetKYsystem_phikcm(ac) * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_phiycm_deg(EAxisConvention ac) const
{
  return GetKYsystem_phiycm(ac) * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_phikcm_deg(EAxisConvention ac) const
{
  return GetKNsystem_phikcm(ac) * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_phincm_deg(EAxisConvention ac) const
{
  return GetKNsystem_phincm(ac) * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_phiycm_deg(EAxisConvention ac) const
{
  return GetYNsystem_phiycm(ac) * TMath::RadToDeg();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_phincm_deg(EAxisConvention ac) const
{
  return GetYNsystem_phincm(ac) * TMath::RadToDeg();
}

//_____________________________________________________________________
bool TKinematics2to3::operator==(const TKinematics2to3& rhs) const
{
  // Evaluate whether both objects currently refer to the same point in the
  // kinematical plane.

  // when the objects are the same, the answer is trivial
  if( this==&rhs ) return true;

  double areEqual = true;

  if( GetIsospin()!=rhs.GetIsospin() ) areEqual=false;

  if( fabs(GetStot())>1. ) {
    if( fabs((GetStot()-rhs.GetStot()))>1.e-6*fabs(GetStot()) )
      areEqual = false;
  } else {
    if( fabs(GetStot()-rhs.GetStot())>1.e-6 ) 
      areEqual = false;
  }

  if( fabs(GetSkn())>1. ) {
    if( fabs((GetSkn()-rhs.GetSkn()))>1.e-6*fabs(GetSkn()) )
      areEqual = false;
  } else {
    if( fabs((GetSkn()-rhs.GetSkn()))>1.e-6 )
      areEqual = false;
  }

  if( fabs(GetSyn())>1. ) {
    if( fabs((GetSyn()-rhs.GetSyn()))>1.e-6*fabs(GetSyn()) )
      areEqual = false;
  } else {
    if( fabs((GetSyn()-rhs.GetSyn()))>1.e-6 )
      areEqual = false;
  }

  if( fabs(GetTkg())>1. ) {
    if( fabs((GetTkg()-rhs.GetTkg()))>1.e-6*fabs(GetTkg()) )
      areEqual = false;
  } else {
    if( fabs((GetTkg()-rhs.GetTkg()))>1.e-6 )
      areEqual = false;
  }

  if( fabs(GetTyd())>1. ) {
    if( fabs((GetTyd()-rhs.GetTyd()))>1.e-6*fabs(GetTyd()) )
      areEqual = false;
  } else {
    if( fabs((GetTyd()-rhs.GetTyd()))>1.e-6 )
      areEqual = false;
  }

  if( fabs(GetQsquared())>1. ) {
    if( fabs((GetQsquared()-rhs.GetQsquared()))>1.e-6*fabs(GetQsquared()) )
      areEqual = false;
  } else {
    if( fabs((GetQsquared()-rhs.GetQsquared()))>1.e-6 )
      areEqual = false;
  }

  return areEqual;  
}

//_____________________________________________________________________
double TKinematics2to3::GetVarByGetter(const TString& variable) const
{
  // Return the value of the variable of type 'variable'. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.
  //
  // This function is used by TKinematics2to3::UpdateKinematics() to check the solution.

  double(TKinematics2to3::*ptrToVarGetter)() const; // pointer to TKinematics getters
  
  TString varName = variable;
  varName.ToLower();
  
  // Assign the correct Getter
  if( !strcmp(varName,"pk") ) {
    ptrToVarGetter = &TKinematics2to3::GetPk;
  }
  else if( !strcmp(varName,"py") ) {
    ptrToVarGetter = &TKinematics2to3::GetPy;
  }
  else if( !strcmp(varName,"pn") ) {
    ptrToVarGetter = &TKinematics2to3::GetPn;
  }
  else if( !strcmp(varName,"klab") ) {
    ptrToVarGetter = &TKinematics2to3::GetKlab;
  }
  else if( !strcmp(varName,"ek") ) {
    ptrToVarGetter = &TKinematics2to3::GetEk;
  }
  else if( !strcmp(varName,"ey") ) {
    ptrToVarGetter = &TKinematics2to3::GetEy;
  }
  else if( !strcmp(varName,"en") ) {
    ptrToVarGetter = &TKinematics2to3::GetEn;
  }
  else if( !strcmp(varName,"wlab") ) {
    ptrToVarGetter = &TKinematics2to3::GetWlab;
  }
  else if( !strcmp(varName,"thk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk_deg;
  }
  else if( !strcmp(varName,"thy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy_deg;
  }
  else if( !strcmp(varName,"thn_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn_deg;
  }
  else if( !strcmp(varName,"thk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk;
  }
  else if( !strcmp(varName,"thy") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy;
  }
  else if( !strcmp(varName,"thn") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn;
  }
  else if( !strcmp(varName,"costhk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthk;
  }
  else if( !strcmp(varName,"costhy") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthy;
  }
  else if( !strcmp(varName,"costhn") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthn;
  }
  else if( !strcmp(varName,"phik_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_deg_kS;
  }
  else if( !strcmp(varName,"phiy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_deg_kS;
  }
  else if( !strcmp(varName,"phin_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_deg_kS;
  }
  else if( !strcmp(varName,"phik") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_kS;
  }
  else if( !strcmp(varName,"phiy") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_kS;
  }
  else if( !strcmp(varName,"phin") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_kS;
  }
  else if( !strcmp(varName,"thky_deg") ||
	   !strcmp(varName,"thyk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThky_deg;
  }
  else if( !strcmp(varName,"thkn_deg") ||
	   !strcmp(varName,"thnk_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn_deg;
  }
  else if( !strcmp(varName,"thyn_deg") ||
	   !strcmp(varName,"thny_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn_deg;
  }
  else if( !strcmp(varName,"thky") ||
	   !strcmp(varName,"thyk" )) {
    ptrToVarGetter = &TKinematics2to3::GetThky;
  }
  else if( !strcmp(varName,"thkn") ||
	   !strcmp(varName,"thnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn;
  }
  else if( !strcmp(varName,"thyn") ||
	   !strcmp(varName,"thny") ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn;
  }
  else if( !strcmp(varName,"costhky") ||
	   !strcmp(varName,"costhyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthky;
  }
  else if( !strcmp(varName,"costhkn") ||
	   !strcmp(varName,"costhnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthkn;
  }
  else if( !strcmp(varName,"costhyn") ||
	   !strcmp(varName,"costhny") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthyn;
  }
  else if( !strcmp(varName,"qsquared") ) {
    ptrToVarGetter = &TKinematics2to3::GetQsquared;
  }
  else if( !strcmp(varName,"wtot") ) {
    ptrToVarGetter = &TKinematics2to3::GetWtot;
  }
  else if( !strcmp(varName,"wky") ||
	   !strcmp(varName,"wyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWky;
  }
  else if( !strcmp(varName,"wkn") ||
	   !strcmp(varName,"wnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWkn;
  }
  else if( !strcmp(varName,"wyn") ||
	   !strcmp(varName,"wny") ) {
    ptrToVarGetter = &TKinematics2to3::GetWyn;
  }
  else if( !strcmp(varName,"stot") ) {
    ptrToVarGetter = &TKinematics2to3::GetStot;
  }
  else if( !strcmp(varName,"sky") ||
	   !strcmp(varName,"syk") ) {
    ptrToVarGetter = &TKinematics2to3::GetSky;
  }
  else if( !strcmp(varName,"skn") ||
	   !strcmp(varName,"snk" )) {
    ptrToVarGetter = &TKinematics2to3::GetSkn;
  }
  else if( !strcmp(varName,"syn") ||
	   !strcmp(varName,"sny" )) {
    ptrToVarGetter = &TKinematics2to3::GetSyn;
  }
  else if( !strcmp(varName,"tkg") ||
	   !strcmp(varName,"tgk" )) {
    ptrToVarGetter = &TKinematics2to3::GetTkg;
  }
  else if( !strcmp(varName,"tyg") ||
	   !strcmp(varName,"tgy" )) {
    ptrToVarGetter = &TKinematics2to3::GetTyg;
  }
  else if( !strcmp(varName,"tng") ||
	   !strcmp(varName,"tgn" )) {
    ptrToVarGetter = &TKinematics2to3::GetTng;
  }
  else if( !strcmp(varName,"tkd") ||
	   !strcmp(varName,"tdk" )) {
    ptrToVarGetter = &TKinematics2to3::GetTkd;
  }
  else if( !strcmp(varName,"tyd") ||
	   !strcmp(varName,"tdy" )) {
    ptrToVarGetter = &TKinematics2to3::GetTyd;
  }
  else if( !strcmp(varName,"tnd") ||
	   !strcmp(varName,"tdn" )) {
    ptrToVarGetter = &TKinematics2to3::GetTnd;
  }
  else if( !strcmp(varName,"kypk") ||
	   !strcmp(varName,"ykpk" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Pk;
  }
  else if( !strcmp(varName,"kypy") ||
	   !strcmp(varName,"ykpy" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Py;
  }
  else if( !strcmp(varName,"kyek") ||
	   !strcmp(varName,"ykek" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ek;
  }
  else if( !strcmp(varName,"kyey") ||
	   !strcmp(varName,"ykey" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ey;
  }
  else if( !strcmp(varName,"kys") ||
	   !strcmp(varName,"yks" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_S;
  }
  else if( !strcmp(varName,"kyw") ||
	   !strcmp(varName,"ykw" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_W;
  }
  else if( !strcmp(varName,"kycosthkcm") ||
	   !strcmp(varName,"ykcosthkcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhkcm;
  }
  else if( !strcmp(varName,"kycosthycm") ||
	   !strcmp(varName,"ykcosthycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhycm;
  }
  else if( !strcmp(varName,"kythkcm") ||
	   !strcmp(varName,"ykthkcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thkcm;
  }
  else if( !strcmp(varName,"kythycm") ||
	   !strcmp(varName,"ykthycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thycm;
  }
  else if( !strcmp(varName,"kyphikcm") ||
	   !strcmp(varName,"ykphikcm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_kS;
  }
  else if( !strcmp(varName,"kyphiycm") ||
	   !strcmp(varName,"ykphiycm" )) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_kS;
  }
  else if( !strcmp(varName,"kyhpikcm_deg") ||
	   !strcmp(varName,"ykphikcm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_deg_kS;
  }
  else if( !strcmp(varName,"kyphiycm_deg") ||
	   !strcmp(varName,"ykphiycm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_deg_kS;
  }
  else if( !strcmp(varName,"knpk") ||
	   !strcmp(varName,"nkpk" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pk;
  }
  else if( !strcmp(varName,"knpn") ||
	   !strcmp(varName,"nkpn" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pn;
  }
  else if( !strcmp(varName,"knek") ||
	   !strcmp(varName,"nkek" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Ek;
  }
  else if( !strcmp(varName,"knen") ||
	   !strcmp(varName,"nken" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_En;
  }
  else if( !strcmp(varName,"kns") ||
	   !strcmp(varName,"nks" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_S;
  }
  else if( !strcmp(varName,"knw") ||
	   !strcmp(varName,"nkw" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_W;
  }
  else if( !strcmp(varName,"kncosthkcm") ||
	   !strcmp(varName,"nkcosthkcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhkcm;
  }
  else if( !strcmp(varName,"kncosthncm") ||
	   !strcmp(varName,"nkcosthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhncm;
  }
  else if( !strcmp(varName,"knthkcm") ||
	   !strcmp(varName,"nkthkcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thkcm;
  }
  else if( !strcmp(varName,"knthncm") ||
	   !strcmp(varName,"nkthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thncm;
  }
  else if( !strcmp(varName,"knphikcm") ||
	   !strcmp(varName,"nkphikcm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_kS;
  }
  else if( !strcmp(varName,"knphincm") ||
	   !strcmp(varName,"nkphincm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_kS;
  }
  else if( !strcmp(varName,"knphikcm_deg") ||
	   !strcmp(varName,"nkphikcm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_deg_kS;
  }
  else if( !strcmp(varName,"knphincm_deg") ||
	   !strcmp(varName,"nkphincm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_deg_kS;
  }
  else if( !strcmp(varName,"ynpy") ||
	   !strcmp(varName,"nypy" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Py;
  }
  else if( !strcmp(varName,"ynpn") ||
	   !strcmp(varName,"nypn" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Pn;
  }
  else if( !strcmp(varName,"yney") ||
	   !strcmp(varName,"nyey" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Ey;
  }
  else if( !strcmp(varName,"ynen") ||
	   !strcmp(varName,"nyen" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_En;
  }
  else if( !strcmp(varName,"yns") ||
	   !strcmp(varName,"nys" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_S;
  }
  else if( !strcmp(varName,"ynw") ||
	   !strcmp(varName,"nyw" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_W;
  }
  else if( !strcmp(varName,"yncosthycm") ||
	   !strcmp(varName,"nycosthycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhycm;
  }
  else if( !strcmp(varName,"yncosthncm") ||
	   !strcmp(varName,"nycosthncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhncm;
  }
  else if( !strcmp(varName,"ynthycm") ||
	   !strcmp(varName,"nythycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thycm;
  }
  else if( !strcmp(varName,"ynthncm") ||
	   !strcmp(varName,"nythncm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thncm;
  }
  else if( !strcmp(varName,"ynphiycm") ||
	   !strcmp(varName,"nyphiycm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_kS;
  }
  else if( !strcmp(varName,"ynphincm") ||
	   !strcmp(varName,"nyphincm" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_kS;
  }
  else if( !strcmp(varName,"ynphiycm_deg") ||
	   !strcmp(varName,"nyphiycm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_deg_kS;
  }
  else if( !strcmp(varName,"ynphincm_deg") ||
	   !strcmp(varName,"nyphincm_deg" ) ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_deg_kS;
  }
  else {
    std::cerr << "ERROR in  TKinematics2to3::GetVarByGetter(TString): "
	      << "Kinematic variable \"" << variable << "\" unknown!\n";
    exit(1);
  }

  return (this->*ptrToVarGetter)();
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::BoostToKYcm(const TLorentzVector& vector,
					    EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  return BoostToKYcm(ac) * vector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::BoostToKNcm(const TLorentzVector& vector,
					    EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  return BoostToKNcm(ac) * vector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::BoostToYNcm(const TLorentzVector& vector,
					    EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  return BoostToYNcm(ac) * vector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::BoostToGNcm(const TLorentzVector& vector,
					    EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  return BoostToGNcm(ac) * vector;
}
//_____________________________________________________________________
TLorentzVector TKinematics2to3::BoostFromKYcm(const TLorentzVector& vector,
					      EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostFromKYcm(ac) * vector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::BoostFromKNcm(const TLorentzVector& vector,
					      EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostFromKNcm(ac) * vector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::BoostFromYNcm(const TLorentzVector& vector,
					      EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostFromYNcm(ac) * vector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::BoostFromGNcm(const TLorentzVector& vector,
					      EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  return BoostFromGNcm(ac) * vector;
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::BoostToKYcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  TLorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TLorentzRotation rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  TVector3 boostVector = -1.0 * (rotateToKY*KYvector).BoostVector();
  TLorentzRotation boostCMKY = boostVector;

  return boostCMKY * rotateToKY;
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::BoostToKNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  TLorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TLorentzRotation rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  TVector3 boostVector = -1.0 * (rotateToKN*KNvector).BoostVector();
  TLorentzRotation boostCMKN = boostVector;

  return boostCMKN * rotateToKN;
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::BoostToYNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  TLorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TLorentzRotation rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  TVector3 boostVector = -1.0 * (rotateToYN*YNvector).BoostVector();
  TLorentzRotation boostCMYN = boostVector;

  return boostCMYN * rotateToYN;
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::BoostToGNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  TLorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TLorentzRotation rotateAxes;
  rotateAxes.RotateY(-TMath::Sign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  TLorentzVector KvectorCM = rotateAxes * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return rotateAxes * BoostToKYcm(ac);
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::BoostFromKYcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame.

  TLorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TLorentzRotation rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  TVector3 boostVector = -1.0 * (rotateToKY*KYvector).BoostVector();
  TLorentzRotation boostCMKY = boostVector;

  return  rotateToKY.Inverse() * boostCMKY.Inverse();
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::BoostFromKNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  TLorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TLorentzRotation rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  TVector3 boostVector = -1.0 * (rotateToKN*KNvector).BoostVector();
  TLorentzRotation boostCMKN = boostVector;

  return rotateToKN.Inverse() * boostCMKN.Inverse();
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::BoostFromYNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  TLorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TLorentzRotation rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  TVector3 boostVector = -1.0 * (rotateToYN*YNvector).BoostVector();
  TLorentzRotation boostCMYN = boostVector;

  return rotateToYN.Inverse() * boostCMYN.Inverse();
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::BoostFromGNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  TLorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TLorentzRotation rotateAxes;
  rotateAxes.RotateY(-TMath::Sign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  TLorentzVector KvectorCM = rotateAxes * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return BoostFromKYcm(ac) * rotateAxes.Inverse();
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostToKYcm(const TSpinor& spinor,
				     EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  return BoostSpinorToKYcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostToKNcm(const TSpinor& spinor,
				     EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  return BoostSpinorToKNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostToYNcm(const TSpinor& spinor,
				     EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  return BoostSpinorToYNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostToGNcm(const TSpinor& spinor,
				     EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  return BoostSpinorToGNcm(ac) * spinor;
}
//_____________________________________________________________________
TSpinor TKinematics2to3::BoostFromKYcm(const TSpinor& spinor,
				       EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostSpinorFromKYcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostFromKNcm(const TSpinor& spinor,
				       EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostSpinorFromKNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostFromYNcm(const TSpinor& spinor,
				       EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostSpinorFromYNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostFromGNcm(const TSpinor& spinor,
				       EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  return BoostSpinorFromGNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorToKYcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  TLorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  TVector3 boostVector 
    = -1.0 * (TLorentzRotation().RotateY(rotationAngle)*KYvector).BoostVector();
  TSpinor::LorentzRotation boostCMKY = boostVector;

  return boostCMKY * rotateToKY;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorToKNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  TLorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  TVector3 boostVector 
    = -1.0 * (TLorentzRotation().RotateY(rotationAngle)*KNvector).BoostVector();
  TSpinor::LorentzRotation boostCMKN = boostVector;

  return boostCMKN * rotateToKN;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorToYNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  TLorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  TVector3 boostVector 
    = -1.0 * (TLorentzRotation().RotateY(rotationAngle)*YNvector).BoostVector();
  TSpinor::LorentzRotation boostCMYN = boostVector;

  return boostCMYN * rotateToYN;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorToGNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  TLorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TSpinor::LorentzRotation rotateAxes;
  rotateAxes.RotateY(-TMath::Sign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  TLorentzVector KvectorCM 
    = TLorentzRotation().RotateY(-TMath::Sign(GvectorCM.Theta(),GvectorCM.X())) 
    * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return rotateAxes * BoostSpinorToKYcm(ac);
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorFromKYcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame.

  TLorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  TVector3 boostVector 
    = -1.0 * (TLorentzRotation().RotateY(rotationAngle)*KYvector).BoostVector();
  TSpinor::LorentzRotation boostCMKY = boostVector;

  return  rotateToKY.Inverse() * boostCMKY.Inverse();
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorFromKNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  TLorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  TVector3 boostVector 
    = -1.0 * (TLorentzRotation().RotateY(rotationAngle)*KNvector).BoostVector();
  TSpinor::LorentzRotation boostCMKN = boostVector;

  return rotateToKN.Inverse() * boostCMKN.Inverse();
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorFromYNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  TLorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  TVector3 boostVector 
    = -1.0 * (TLorentzRotation().RotateY(rotationAngle)*YNvector).BoostVector();
  TSpinor::LorentzRotation boostCMYN = boostVector;

  return rotateToYN.Inverse() * boostCMYN.Inverse();
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorFromGNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  TLorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TSpinor::LorentzRotation rotateAxes;
  rotateAxes.RotateY(-TMath::Sign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  TLorentzVector KvectorCM 
    =  TLorentzRotation().RotateY(-TMath::Sign(GvectorCM.Theta(),GvectorCM.X()))
    * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return BoostSpinorFromKYcm(ac) * rotateAxes.Inverse();
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionToKYcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  TLorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  TVector3 boostVector = -1.0 * (rotateToKY*KYvector).BoostVector();
  TLorentzQuaternion boostCMKY = boostVector;

  return boostCMKY * rotateToKY;
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionToKNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  TLorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  TVector3 boostVector = -1.0 * (rotateToKN*KNvector).BoostVector();
  TLorentzQuaternion boostCMKN = boostVector;

  return boostCMKN * rotateToKN;
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionToYNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  TLorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  TVector3 boostVector = -1.0 * (rotateToYN*YNvector).BoostVector();
  TLorentzQuaternion boostCMYN = boostVector;

  return boostCMYN * rotateToYN;
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionToGNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  TLorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TLorentzQuaternion rotateAxes;
  rotateAxes.RotateY(-TMath::Sign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  TLorentzVector KvectorCM = rotateAxes * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return rotateAxes * BoostQuaternionToKYcm(ac);
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionFromKYcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame.

  TLorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  TVector3 boostVector = -1.0 * (rotateToKY*KYvector).BoostVector();
  TLorentzQuaternion boostCMKY = boostVector;

  return  rotateToKY.Inverse() * boostCMKY.Inverse();
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionFromKNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  TLorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  TVector3 boostVector = -1.0 * (rotateToKN*KNvector).BoostVector();
  TLorentzQuaternion boostCMKN = boostVector;

  return rotateToKN.Inverse() * boostCMKN.Inverse();
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionFromYNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  TLorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*TMath::Sign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  TVector3 boostVector = -1.0 * (rotateToYN*YNvector).BoostVector();
  TLorentzQuaternion boostCMYN = boostVector;

  return rotateToYN.Inverse() * boostCMYN.Inverse();
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionFromGNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  TLorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TLorentzQuaternion rotateAxes;
  rotateAxes.RotateY(-TMath::Sign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  TLorentzVector KvectorCM = rotateAxes * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return BoostQuaternionFromKYcm(ac) * rotateAxes.Inverse();
}

//_____________________________________________________________________
void TKinematics2to3::Help()
{
  TMPI::Cout() << "I can't help you! "
	       << "In fact, can anyone?\n";
}

//_____________________________________________________________________
// FourVector<double> ToFourVector(const TLorentzVector& v)
// {
//   return FourVector<double>(v.E(),v.X(),v.Y(),v.Z());
// }
// 
// //_____________________________________________________________________
// TLorentzVector ToLorentzVector(const FourVector<double>& v)
// {
//   return TLorentzVector(v[1],v[2],v[3],v[0]);
// }

//_____________________________________________________________________
double TKinematics2to3::GetPhotonThreshold() const
{
  // Calculate the threshold photon lab energy for the current isospin
  // channel.
  return (( Power(fMk+fMy+fMn,2.) + fQsquared - fMd*fMd ) / 2. / fMd) + kUflow;
}

//_____________________________________________________________________
void TKinematics2to3::GetBoundsMomentumParticle1(double th1,
						 double& lower,
						 double& upper) const
{
  // For the current invariant mass 'w' and a given scattering angle 
  // 'th1' of the particle specified by the axis convention, the lab 
  // momentum |p_1| has an upper and lower bound.
  //
  // These bounds are calculated by imposing the constraint W_23>=m_2+m_3.
  // This leads to a quadratic equation for |p_1|, which we solve.

  double eg = fGvector->E();
  double stot = fMd*fMd - fQsquared + 2.*fMd*fGvector->E();
  
  double costh1 = Cos(th1);
  double m1=0., m23=0.;
  if( fAC==kK ) {
    m1 = fMk;
    m23 = fMy + fMn;
  } // particle1 = kaon
  else if( fAC==kY ) {
    m1 = fMy;
    m23 = fMk + fMn;
  } // particle1 = hyperon
  else {
    m1 = fMn;
    m23 = fMk + fMy;
  } // particle1 = nucleon

  // Define the quadratic equation...
  Polynomial eq(// |p_1|^2 coefficient
		(fQsquared+eg*eg)*costh1*costh1 - (eg+fMd)*(eg+fMd),
		// |p_1| coefficient
		TMath::Sqrt(fQsquared+eg*eg)*costh1*(stot+m1*m1-(m23)*(m23)),
		// constant coefficient
		(TMath::Power(stot+m1*m1-(m23)*(m23),2.)/4.
		 -m1*m1*(eg+fMd)*(eg+fMd)));

  // ...and solve
  std::vector<double> roots = eq.FindRealRoots();
  
  if( roots[0]>roots[1] ) {
    upper = TMath::Max(roots[0],0.);
    lower = TMath::Max(roots[1],0.);
  } else {
    upper = TMath::Max(roots[1],0.);
    lower = TMath::Max(roots[0],0.);
  }
}

//_____________________________________________________________________
double TKinematics2to3::GetLowerBoundMomentumParticle1(double th1) const
{
  // For the current invariant mass 'w' and a given scattering angle 
  // 'th1' of the particle specified by the axis convention, the lab 
  // momentum |p_1| has a lower bound.

  double lower,upper;
  GetBoundsMomentumParticle1(th1,lower,upper);
  
  return lower + kUflow;
}

//_____________________________________________________________________
double TKinematics2to3::GetUpperBoundMomentumParticle1(double th1) const
{
  // For the current invariant mass 'w' and a given scattering angle 
  // 'th1' of the particle specified by the axis convention, the lab 
  // momentum |p_1| has an upper bound.

  double lower,upper;
  GetBoundsMomentumParticle1(th1,lower,upper);
  
  return upper - kUflow;
}

//_____________________________________________________________________
double TKinematics2to3::GetMaximumScatterAngleParticle1() const
{
  // For the current invariant mass 'w' and lab momentum 'p_1' of the
  // particle specified by the axis convention, the lab scattering angle 
  // 'theta_1' has a maximum.
  //
  // This maximum is found by imposing the constraint W_23 >= m_2+m_3.
  
  // Find the 4-vector of particle 1 and sum masses of particles 2 and 3.
  TLorentzVector *p1;
  double m23;
  if( fAC==kK ) {
    p1 = fKvector;
    m23 = fMy + fMn;
  } // particle1 = kaon
  else if( fAC==kY ) {
    p1 = fYvector;
    m23 = fMk + fMn;
  } // particle1 = hyperon
  else {
    p1 = fNvector;
    m23 = fMk + fMy;
  } // particle1 = nucleon

  double stot = fMd*fMd - fQsquared + 2.*fMd*fGvector->E();

  double costh1Min = 
    ( m23*m23 + 2.*p1->E() * (fGvector->E() + fMd)
      - stot - p1->M2() )
    / 2. / p1->P() / fGvector->P();

  if( costh1Min > 1.-kUflow ) return 0.;

  return ACos(costh1Min) - kUflow;
}

//_____________________________________________________________________
// TString TKinematics2to3::HashName() const
// {
//   // Returns a unique string that specifies the current state of 
//   // the kinematics object.
// 
//   //Each kinematic point is fixed by the isospin channel, 
//   //the Mandelstam variables s_tot, s_kn, s_ny, t_gk, tdy and Q^2
//   double key1[7] = {GetIsospin(), GetQsquared(), 
// 		    GetWlab(), GetEk(),
// 		    GetCosthk(), GetYNsystem_costhycm(),
// 		    GetYNsystem_phiycm(kK)};
//   
//   // We create a hash key using the function in the cachetools folder
//   char hash1[sizeof(double)*7*2 + 1];
//   doublestochar(7, key1, hash1);
// 
//   return hash1;
// }

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::RotateToKac(EAxisConvention axisConvention) const
{
  // Object to rotate a 4vector defined in the LAB frame from the 'ac'
  // axis convention to the 'kK' axis convention.
  TLorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateToKac(fAC);

  else if( axisConvention == kY ) 
    rotation.RotateZ(-1.0*GetPhik(kY));

  else if( axisConvention == kN )
    rotation.RotateZ(-1.0*GetPhik(kN));

  // else( axisConvention == kK )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::RotateToYac(EAxisConvention axisConvention) const
{
  // Object to rotate a 4vector defined in the LAB frame from the 'ac'
  // axis convention to the 'kY' axis convention.
  TLorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateToYac(fAC);

  else if( axisConvention == kK ) 
    rotation.RotateZ(-1.0*GetPhiy(kK));

  else if( axisConvention == kN )
    rotation.RotateZ(-1.0*GetPhiy(kN));

  // else( axisConvention == kY )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
TLorentzRotation TKinematics2to3::RotateToNac(EAxisConvention axisConvention) const
{
  // Object to rotate a 4vector defined in the LAB frame from the 'ac'
  // axis convention to the 'kN' axis convention.
  TLorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateToNac(fAC);

  else if( axisConvention == kK ) 
    rotation.RotateZ(-1.0*GetPhin(kK));

  else if( axisConvention == kY )
    rotation.RotateZ(-1.0*GetPhin(kY));

  // else( axisConvention == kN )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::RotateSpinorToKac(EAxisConvention axisConvention) const
{
  // Object to rotate a Dirac spinor defined in the LAB frame from the 'ac'
  // axis convention to the 'kK' axis convention.
  TSpinor::LorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateSpinorToKac(fAC);

  else if( axisConvention == kY ) 
    rotation.RotateZ(-1.0*GetPhik(kY));

  else if( axisConvention == kN )
    rotation.RotateZ(-1.0*GetPhik(kN));

  // else( axisConvention == kK )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::RotateSpinorToYac(EAxisConvention axisConvention) const
{
  // Object to rotate a Dirac spinor defined in the LAB frame from the 'ac'
  // axis convention to the 'kY' axis convention.
  TSpinor::LorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateSpinorToYac(fAC);

  else if( axisConvention == kK ) 
    rotation.RotateZ(-1.0*GetPhiy(kK));

  else if( axisConvention == kN )
    rotation.RotateZ(-1.0*GetPhiy(kN));

  // else( axisConvention == kY )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::RotateSpinorToNac(EAxisConvention axisConvention) const
{
  // Object to rotate a Dirac spinor defined in the LAB frame from the 'ac'
  // axis convention to the 'kN' axis convention.
  TSpinor::LorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateSpinorToNac(fAC);

  else if( axisConvention == kK ) 
    rotation.RotateZ(-1.0*GetPhin(kK));

  else if( axisConvention == kY )
    rotation.RotateZ(-1.0*GetPhin(kY));

  // else( axisConvention == kN )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::RotateToKac(const TLorentzVector& vector,
					    EAxisConvention axisConvention) const
{
  // Rotate the 4vector 'vector' defined in the LAB frame from the 'ac'
  // axis convention to the 'kK' axis convention.  
  return RotateToKac(axisConvention) * vector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::RotateToYac(const TLorentzVector& vector,
					    EAxisConvention axisConvention) const
{
  // Rotate the 4vector 'vector' defined in the LAB frame from the 'ac'
  // axis convention to the 'kY' axis convention.  
  return RotateToYac(axisConvention) * vector;
}

//_____________________________________________________________________
TLorentzVector TKinematics2to3::RotateToNac(const TLorentzVector& vector,
					    EAxisConvention axisConvention) const
{
  // Rotate the 4vector 'vector' defined in the LAB frame from the 'ac'
  // axis convention to the 'kN' axis convention.  
  return RotateToNac(axisConvention) * vector;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::RotateToKac(const TSpinor& spinor,
				     EAxisConvention axisConvention) const
{
  // Rotate the Dirac spinor 'spinor' defined in the LAB frame from the 'ac'
  // axis convention to the 'kK' axis convention.  
  return RotateSpinorToKac(axisConvention) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::RotateToYac(const TSpinor& spinor,
				     EAxisConvention axisConvention) const
{
  // Rotate the Dirac spinor 'spinor' defined in the LAB frame from the 'ac'
  // axis convention to the 'kY' axis convention.  
  return RotateSpinorToYac(axisConvention) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::RotateToNac(const TSpinor& spinor,
				     EAxisConvention axisConvention) const
{
  // Rotate the Dirac spinor 'spinor' defined in the LAB frame from the 'ac'
  // axis convention to the 'kN' axis convention.  
  return RotateSpinorToNac(axisConvention) * spinor;
}
