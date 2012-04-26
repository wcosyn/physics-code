/*
 * TKinematics2to3WithLabAngles.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: August, 20 2009
 *
 */

///////////////////////////////////////////////////////////////////////////
/* Begin_Html
<center><h2>TKinematics2to3WithLabAngles</h2></center>
This is an extention to TKinematics2to3 which allows to calculate kinematics 
for g D -> K Y N.
<p>
All quantities are in MeV or rad.
<p>

The user needs to specify 6 variables to unambiguously fix the kinematics.
Obviously, he/she needs to made sure the variables of choice are independent. 
This class extents TKinematics2to3 with 1 extra situation
(extention of <a href="TKinematics2to3.html#indepvarCaseA">case A</a>):
<ul>
<li>The user provides info about the energy and the theta angle (either cos or angle) of the mother.
<li>A 3rd variable should either specify
<ul><li>the photon's energy
    <li>the invariant mass of the total system
    <li>the invariant mass of the children subsystem</ul>
<li>Finally, the user specifies the solid angle (theta+phi) of one of the children in the LAB frame.
<li>The virtualiy Q<sup>2</sup> of the photon is always one of the independent variables.
</ul> 

Such a set of independent variables can have two solutions.
The number of solutions can be checked with GetNumberOfSolutions() and the
separate solutions can be accessed with GetSolution(int).
<p>
End_Html */
///////////////////////////////////////////////////////////////////////////

#include "TMPI.h"
#include "TKinematics2to3WithLabAngles.h"
#include "constants.hpp"
#include <TLorentzRotation.h>
#include <TMath.h>
using TMath::Cos; using TMath::Sin;
using TMath::Abs; using TMath::Sqrt;
#include <Math/Polynomial.h>
using ROOT::Math::Polynomial;
#include <iostream>
using std::cout; using std::cerr;
#include <cassert>
#include <cstring>
#include <cstdio>
#include <vector>
using std::vector;

ClassImp(TKinematics2to3WithLabAngles)

//_____________________________________________________________________
const std::ofstream TKinematics2to3WithLabAngles::fgDevNull("/dev/null");

//_____________________________________________________________________
TKinematics2to3WithLabAngles::TKinematics2to3WithLabAngles(const char *name, const char *title,
							   int isospin, EAxisConvention ac, 
							   const char *format,
							   double var1, double var2, double var3,
							   double var4, double var5, double var6,
							   double Md,double Mn, double Mk,double My)
							  
  : TKinematics2to3(name,title,isospin,kK,"qsquared:wlab:pk:thk:yncosthycm:ynphiycm",
		    0.,0.,0.,0.,0.,0.,Md,Mn,Mk,My),
    fNrOfSolutions(0), fSolutions(new TKinematics2to3*[2])
{
  // Constructor
  for(int i=0; i<2; ++i) fSolutions[i] = 0;
  TKinematics2to3::SetAxisConvention(ac);
  ReadFormatString(format);
  TKinematics2to3::SetVar(var1,var2,var3,var4,var5,var6);
}

//_____________________________________________________________________
TKinematics2to3WithLabAngles::TKinematics2to3WithLabAngles(TRootIOCtor *rio)
  :  TKinematics2to3(rio), fNrOfSolutions(0), fSolutions(0)
{
  // ROOT I/O constructor
}

//_____________________________________________________________________
TKinematics2to3WithLabAngles::TKinematics2to3WithLabAngles(const TKinematics2to3WithLabAngles& rhs)
  : TKinematics2to3("","",1,kK,"qsquared:wlab:pk:thk:yncosthycm:ynphiycm",
		    0.,0.,0.,0.,0.,0.,rhs.fMd,rhs.fMn,rhs.fMk,rhs.fMy),
    fNrOfSolutions(rhs.fNrOfSolutions), fSolutions(new TKinematics2to3*[2])
{
  // Copy constructor
  TKinematics2to3::operator=(rhs);
  for(int i=0; i<fNrOfSolutions; ++i)
    fSolutions[i] = new TKinematics2to3(*rhs.fSolutions[i]);
  for(int i=fNrOfSolutions; i<2; ++i)
    fSolutions[i] = 0;
}

//_____________________________________________________________________
TKinematics2to3WithLabAngles& TKinematics2to3WithLabAngles::operator=(const TKinematics2to3WithLabAngles& rhs)
{
  // Assignment
  if( this!=&rhs ) {
    TKinematics2to3::operator=(rhs);
    fNrOfSolutions = rhs.fNrOfSolutions;
    for(int i=0; i<2; ++i) if(fSolutions[i]) delete fSolutions[i];
    for(int i=0; i<fNrOfSolutions; ++i)
      fSolutions[i] = new TKinematics2to3(*rhs.fSolutions[i]);
    for(int i=fNrOfSolutions; i<2; ++i)
      fSolutions[i] = 0;
  }
  return *this;
}

//_____________________________________________________________________
TKinematics2to3WithLabAngles::~TKinematics2to3WithLabAngles()
{
  // Destructor
  for(int i=0; i<2; ++i) if(fSolutions[i]) delete fSolutions[i];
  delete[] fSolutions;
}

//_____________________________________________________________________
void TKinematics2to3WithLabAngles::TurnOutputOff(bool off)
{
  // Turning the output off means redirecting cout and cerr to /dev/null.
  // Don't forget to turn the output back on at the end of the member function.
  static std::streambuf *devNullBuffer = fgDevNull.rdbuf();
  static std::streambuf *coutBuffer = std::cout.rdbuf();
  static std::streambuf *cerrBuffer = std::cerr.rdbuf();

  std::cout.flush();

  if( off ) {
    std::cout.rdbuf( devNullBuffer );
    std::cerr.rdbuf( devNullBuffer );
  }
  else {
    std::cout.rdbuf( coutBuffer );
    std::cerr.rdbuf( cerrBuffer );
  }
}

//_____________________________________________________________________
void TKinematics2to3WithLabAngles::ReadFormatString(const char *format)
{
  // Process the format string.
  TurnOutputOff();
  TKinematics2to3::ReadFormatString(format);
  TurnOutputOff(false);

  // Check whether the angles of particle 2 are really given in the lab
  if( GetSolutionType()[2] != 'L' ) {
    std::cerr << "ERROR in TKinematics2to3WithLabAngles::ReadFormatString(const char *): "
	      << "This class is reserved for format strings with the angles "
	      << "of particle 2 given in the LAB frame.\n";
    assert(1==0);
  }
}

//_____________________________________________________________________
double* TKinematics2to3WithLabAngles::GetVarArray(const TString& varName) const
{
  // Return an array of the variable of type 'variable'. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.
  //
  // The returned pointer has a size GetNumberOfSteps() and is owned by the user.
  TurnOutputOff();
  double *array = TKinematics2to3::GetVarArray(varName);
  TurnOutputOff(false);

  return array;
}

//_____________________________________________________________________
double* TKinematics2to3WithLabAngles::GetPhysicalVarArray(const TString& varName) const
{
  // Return an array of the variable of type 'variable'. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.
  //
  // The returned pointer has a size GetNumberOfPhysicalSteps() and is owned by the user.
  TurnOutputOff();
  double *array = TKinematics2to3::GetPhysicalVarArray(varName);
  TurnOutputOff(false);

  return array;
}

//_____________________________________________________________________
int TKinematics2to3WithLabAngles::GetNumberOfPhysicalSteps() const
{
  TurnOutputOff();
  int nrOfPhysicals = TKinematics2to3::GetNumberOfPhysicalSteps();
  TurnOutputOff(false);

  return nrOfPhysicals;
}

//_____________________________________________________________________
void TKinematics2to3WithLabAngles::UpdateKinematics()
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

  //    The user provides info about the energy (fE1Var)
  //    and the theta-angle (fTh1Var) of the particle that
  //    fixes the x-axis (fAC).
  //    One variable should either specify the photon energy
  //    or the invariant mass of the total system or the
  //    (2+3) subsystem (fEVar).
  //    In addition angular info needs to be given about a
  //    second final particle. Both theta (fTh2Var) and phi (fPhi2Var)
  //    give information in the LAB frame.

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
    fNrOfSolutions = 0;
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
      fNrOfSolutions = 1;
    }
    else {
      if( sol1 && sol2 ) {
	if( fabs(roots[0]-roots[1]) > STRANGEUFLOW*fabs(roots[0]))
	  // We pick solution arbitrarily for now. We will deal with this
	  // at the end of the function.
	  fNrOfSolutions = 2;
	  p2 = roots[0];
      }
      else {
	if( sol2 ) roots[0] = roots[1];
	p2 = roots[0];
	fNrOfSolutions = 1;
      }
    }
  } // solution found
      
  if(p2<0.) fIsPhysical = false;

  // set 4vector for particle 2
  vector2.SetXYZM(sinth2*Cos(*fPhi2Var)*p2,
		  sinth2*Sin(*fPhi2Var)*p2,
		  costh2*p2, m2);

  vector3 = vector23 - vector2;
 

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

  // In the end we need to initialize the solutions we obtained when we
  // determined the lab momentum of particle 2.
  // This is obviously not necessary when we didn't find a solution.
  for(int i=0; i<fNrOfSolutions; ++i) {
    
    // Determine the centre-of-mass angles
    double costh2cm=0., phi2cm=0.;
    
    // Build the vector in the lab and boost it!
    TLorentzVector vector2cm;
    vector2cm.SetXYZM(sinth2*Cos(*fPhi2Var)*roots[i],
		      sinth2*Sin(*fPhi2Var)*roots[i],
		      costh2*roots[i], m2);
    vector2cm = boostCM23 * rotateTo23 * vector2cm;

    costh2cm = vector2cm.CosTheta();
    phi2cm = vector2cm.Phi();

    // Determine the correct format string.
    char formatString[100];
    
    // the first 4 tokens in fFormat can be copied.
    int pos=0, countColon=0;
    while( countColon !=4 ) 
      if( fFormat[pos++] == ':' ) ++countColon;
    std::strncpy(formatString,fFormat,pos);
    formatString[pos] = '\0';
    
    // Add the 2 last tokens being 23costh2cm and 23phi2cm.
    char name23[3];
    if(fAC==kK) std::strcpy(name23,"yn");
    else if(fAC==kY) std::strcpy(name23,"kn");
    else if(fAC==kN) std::strcpy(name23,"ky");

    char name2[2];
    name2[0] = GetSolutionType()[1];
    name2[1] = '\0';

    // paste everything together
    std::sprintf(formatString,"%s%s%s%s%s%s%s%s%s",
		 formatString,name23,"costh",name2,"cm:",
		 name23,"phi",name2,"cm");
    
    // Create the solution
    if(fSolutions[i]) delete fSolutions[i];
    fSolutions[i] = new TKinematics2to3("","",fIsospin,fAC,
					formatString,
					GetVar(1), GetVar(2),
					GetVar(3), GetVar(4),
					costh2cm,phi2cm);

  } // end initializing fSolutions

//   for(int i=fNrOfSolutions; i<2; ++i) {
//     if(fSolutions[i]) {
//       delete fSolutions[i];
//       fSolutions[i] = 0;
//     }
//   } // clean up the remaining solutions
}

//_____________________________________________________________________
void TKinematics2to3WithLabAngles::Streamer(TBuffer &R__b)
{
  // Stream an object of class TKinematics2to3WithLabAngles.

  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
    TurnOutputOff();
    TKinematics2to3::Streamer(R__b);
    TurnOutputOff(false);
    fSolutions = new TKinematics2to3*[2];
    for(int i=0; i<2; ++i) fSolutions[i] = 0;
    UpdateKinematics();
    R__b.CheckByteCount(R__s, R__c, TKinematics2to3WithLabAngles::IsA());
  } else {
    R__c = R__b.WriteVersion(TKinematics2to3WithLabAngles::IsA(), kTRUE);
    TurnOutputOff();
    TKinematics2to3::Streamer(R__b);
    TurnOutputOff(false);
    R__b.SetByteCount(R__c, kTRUE);
  }
}

//_____________________________________________________________________
const TKinematics2to3& TKinematics2to3WithLabAngles::GetSolution(int solution) const
{
  // Retreive the solution. The number of solutions can be obtained from
  // GetNumberOfSolutions().
  // Solutions are numbered 0 and 1.
  if( solution >= fNrOfSolutions ) {
    std::cerr << "ERROR in TKinematics2to3WithLabAngles::GetSolution(int): "
	      << "index out-of-bounds.\n";
    assert(1==0);
  }
  return *fSolutions[solution];
}
