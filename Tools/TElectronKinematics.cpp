/*
 * TElectronKinematics.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: March,24 2010
 *
 */

#include "TElectronKinematics.h"
#include "TKinematics2to2.h"
#include "TKinematics2to3.h"
// #include <Structures.h>
#include "constants.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
using std::cerr;

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TElectronKinematics                                                   //
//                                                                       //
// In electroproduction reactions we separate the electron scattering
// vertex from the virtual photon - deuteron system.
// This class deals with the kinematics of the former. In conjunction
// with a TKinematics2to3 object the kinematics are completely fixed by
// * the incoming electron's beam energy
// * or the outgoing electron's scattering angle
// * or the virtual photon's transverse linear polarization.
//
// The virtual photon's transverse linear polarization is defined as
//BEGIN_LATEX #varepsilon = #left( 1+2#frac{|p_{#gamma}|^{2}}{Q^{2}}tan^{2}#frac{#theta_{e}}{2} #right)^{-1}  END_LATEX
// and can be obtained or set with GetEpsilon and SetEpsilon.
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TElectronKinematics)

//_____________________________________________________________________
TElectronKinematics::TElectronKinematics()
: fInput(kEpsilon), fBeamEnergy(0.), fEpsilon(0.), fScatterAngle(0.)
{
  // Protected Default Constructor
  //
  // Instances of TElectronKinematics should be created with
  // the named constructors.
}

//_____________________________________________________________________
TElectronKinematics::TElectronKinematics(TRootIOCtor* rio)
  : fInput(kEpsilon), fBeamEnergy(0.), fEpsilon(0.), fScatterAngle(0.)
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
TElectronKinematics::TElectronKinematics(const TElectronKinematics& rhs)
  : fInput(rhs.fInput), fBeamEnergy(rhs.fBeamEnergy), 
    fEpsilon(rhs.fEpsilon), fScatterAngle(rhs.fScatterAngle)
{
  // Copy Constructor
}

//_____________________________________________________________________
TElectronKinematics& TElectronKinematics::operator=(const TElectronKinematics& rhs)
{
  // Assignment
  if( this!=&rhs ) {
    fInput = rhs.fInput;
    fBeamEnergy = rhs.fBeamEnergy;
    fEpsilon = rhs.fEpsilon;
    fScatterAngle = rhs.fScatterAngle;
  }
  return *this;
}

//_____________________________________________________________________
TElectronKinematics* TElectronKinematics::CreateWithBeamEnergy(double beamEnergy)
{
  // Named constructor
  //
  // Fix the electron kinematics by specifying the energy
  // of the incoming electron beam in [MeV].
  //
  // The user is responsible for deleting the returned object.
  TElectronKinematics *kinematics = new TElectronKinematics();
  kinematics->SetBeamEnergy(beamEnergy);
  return kinematics;
}

//_____________________________________________________________________
TElectronKinematics* TElectronKinematics::CreateWithEpsilon(double epsilon)
{
  // Named constructor
  //
  // Fix the electron kinematics by specifying the 
  // virtual photon's transverse linear polarization.
  //
  // The user is responsible for deleting the returned object.
  TElectronKinematics *kinematics = new TElectronKinematics();
  kinematics->SetEpsilon(epsilon);
  return kinematics;
}

//_____________________________________________________________________
TElectronKinematics* TElectronKinematics::CreateWithCosScatterAngle(double scatterAngle)
{
  // Named constructor
  //
  // Fix the electron kinematics by specifying the cosine of the
  // electron's scattering angle in ]-1,1].
  //
  // The user is responsible for deleting the returned object.
  TElectronKinematics *kinematics = new TElectronKinematics();
  kinematics->SetCosScatterAngle(scatterAngle);
  return kinematics;
}

//_____________________________________________________________________
void TElectronKinematics::SetBeamEnergy(double beamEnergy)
{
  // Fix the electron kinematics by specifying the energy
  // of the incoming electron beam in [MeV].
  fInput = kBeamEnergy;
  fBeamEnergy = beamEnergy;
}

//_____________________________________________________________________
void TElectronKinematics::SetEpsilon(double epsilon)
{
  // Fix the electron kinematics by specifying the
  // virtual photon's transverse linear polarization.
  fInput = kEpsilon;
  fEpsilon = epsilon;
}

//_____________________________________________________________________
void TElectronKinematics::SetCosScatterAngle(double scatterAngle)
{
  // Fix the electron kinematics by specifying the cosine of the
  // electron's scattering angle in ]-1,1].
  fInput = kScatterAngle;
  fScatterAngle = scatterAngle;
}

//_____________________________________________________________________
double TElectronKinematics::GetBeamEnergy(const TKinematics2to2& tk) const
{
  // Returns the incoming electron's energy in [MeV]
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TElectronKinematics::GetBeamEnergy(TKinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fBeamEnergy;
}

double TElectronKinematics::GetBeamEnergy(const TKinematics2to3& tk) const
{
  // Returns the incoming electron's energy in [MeV]
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TElectronKinematics::GetBeamEnergy(TKinematics2to3&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fBeamEnergy;
}


//_____________________________________________________________________
double TElectronKinematics::GetEpsilon(const TKinematics2to2& tk) const
{
  // Returns the virtual photon's transverse linear polarization
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TElectronKinematics::GetEpsilon(TKinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fEpsilon;
}

double TElectronKinematics::GetEpsilon(const TKinematics2to3& tk) const
{
  // Returns the virtual photon's transverse linear polarization
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TElectronKinematics::GetEpsilon(TKinematics2to3&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fEpsilon;
}

//_____________________________________________________________________
double TElectronKinematics::GetCosScatterAngle(const TKinematics2to2& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TElectronKinematics::GetCosScatterAngle(TKinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fScatterAngle; 
}

double TElectronKinematics::GetCosScatterAngle(const TKinematics2to3& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TElectronKinematics::GetCosScatterAngle(TKinematics2to3&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fScatterAngle; 
}

//_____________________________________________________________________
double TElectronKinematics::GetTan2HalfAngle(const TKinematics2to2& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TElectronKinematics::GetCosScatterAngle(TKinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return ftan2HalfAngle; 
}


//_____________________________________________________________________
double TElectronKinematics::GetTan2HalfAngle(const TKinematics2to3& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TElectronKinematics::GetCosScatterAngle(TKinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return ftan2HalfAngle; 
}

//_____________________________________________________________________
bool TElectronKinematics::SolveKinematics(const TKinematics2to2& tk) const
{
  // Solve the electron scattering system kinematics and check whether 
  // the electron-system kinematics are compatible with those of the 
  // virtual photon - deuteron system.
  //
  // Returns 'false' when the kinematics are unphysical.

  double tmp;

  switch( fInput ) {
  case kBeamEnergy:
    // Solve the kinematics starting with the beam energy.
    if( fBeamEnergy <= tk.GetWlab() ) return false;
    
    // Find the scattering angle
    fScatterAngle = 1. - tk.GetQsquared()/2./fBeamEnergy/(fBeamEnergy-tk.GetWlab());
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Finally calculate epsilon
    if( fScatterAngle<= -1. ) fEpsilon = 0.;
    else if( tk.GetQsquared()<STRANGEUFLOW ) fEpsilon = 1.;
    else {
      ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
      fEpsilon = 1./(1.+2.*(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared())
		     *ftan2HalfAngle);
    }
    if( fEpsilon < 0. ) {
      if( fEpsilon < -STRANGEUFLOW ) return false;
      else fEpsilon = 0.;
    } else if( fEpsilon > 1. ) {
      if( (fEpsilon-1.) > STRANGEUFLOW ) return false;
      else fEpsilon = 1.;
    }
    break;

  case kEpsilon:
    // Solve the kinematics starting with epsilon.
    if( fEpsilon<0. || fEpsilon>1. ) return false;

    // Find the scattering angle
    if( tk.GetQsquared()<STRANGEUFLOW ) fScatterAngle = 1.;
    else {
      tmp = (1./fEpsilon-1.)/2./(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared());
      fScatterAngle = (1.-tmp)/(1.+tmp);
    }
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Calculate the beam energy
    if( fScatterAngle > 1. ) return false;
    else if( fScatterAngle == 1. ) fBeamEnergy = tk.GetWlab() + STRANGEUFLOW;
    else
      fBeamEnergy = tk.GetWlab()/2. + sqrt(tk.GetWlab()*tk.GetWlab()/4.
					   +tk.GetQsquared()/2./(1.-fScatterAngle));
    if( fBeamEnergy <= tk.GetWlab() ) return false;
    ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
    break;

  case kScatterAngle:
    // Solve the kinematics starting with the scattering angle
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Find the beam energy
    if( fScatterAngle == 1. ) fBeamEnergy = tk.GetWlab() + STRANGEUFLOW;
    else
      fBeamEnergy = tk.GetWlab()/2. + sqrt( tk.GetWlab()*tk.GetWlab()/4.
					    +tk.GetQsquared()/2./(1.-fScatterAngle) );
    if( fBeamEnergy <= tk.GetWlab() ) return false;

    // Calculate epsilon
    if( fScatterAngle<= -1. ) fEpsilon = 0.;
    else if( tk.GetQsquared()<STRANGEUFLOW ) fEpsilon = 1.;
    else {
      ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
      fEpsilon = 1./(1.+2.*(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared())
		     *ftan2HalfAngle);
    }
    if( fEpsilon < 0. ) {
      if( fEpsilon < -STRANGEUFLOW ) return false;
      else fEpsilon = 0.;
    } else if( fEpsilon > 1. ) {
      if( (fEpsilon-1.) > STRANGEUFLOW ) return false;
      else fEpsilon = 1.;
    }
    ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
    break;
    
  default:
    cerr << "ERROR in TElectronKinematics::SolveKinematics(TKinematics2to2&): "
	 << "Invalid input variable.\n";
    assert(1==0);
  }
  
  return true;
}

//_____________________________________________________________________
bool TElectronKinematics::SolveKinematics(const TKinematics2to3& tk) const
{
  // Solve the electron scattering system kinematics and check whether 
  // the electron-system kinematics are compatible with those of the 
  // virtual photon - deuteron system.
  //
  // Returns 'false' when the kinematics are unphysical.

  double tmp;

  switch( fInput ) {
  case kBeamEnergy:
    // Solve the kinematics starting with the beam energy.
    if( fBeamEnergy <= tk.GetWlab() ) return false;
    
    // Find the scattering angle
    fScatterAngle = 1. - tk.GetQsquared()/2./fBeamEnergy/(fBeamEnergy-tk.GetWlab());
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Finally calculate epsilon
    if( fScatterAngle<= -1. ) fEpsilon = 0.;
    else if( tk.GetQsquared()<STRANGEUFLOW ) fEpsilon = 1.;
    else {
      ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
      fEpsilon = 1./(1.+2.*(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared())
		     *ftan2HalfAngle);
    }
    if( fEpsilon < 0. ) {
      if( fEpsilon < -STRANGEUFLOW ) return false;
      else fEpsilon = 0.;
    } else if( fEpsilon > 1. ) {
      if( (fEpsilon-1.) > STRANGEUFLOW ) return false;
      else fEpsilon = 1.;
    }
    break;

  case kEpsilon:
    // Solve the kinematics starting with epsilon.
    if( fEpsilon<0. || fEpsilon>1. ) return false;

    // Find the scattering angle
    if( tk.GetQsquared()<STRANGEUFLOW ) fScatterAngle = 1.;
    else {
      tmp = (1./fEpsilon-1.)/2./(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared());
      fScatterAngle = (1.-tmp)/(1.+tmp);
    }
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Calculate the beam energy
    if( fScatterAngle > 1. ) return false;
    else if( fScatterAngle == 1. ) fBeamEnergy = tk.GetWlab() + STRANGEUFLOW;
    else
      fBeamEnergy = tk.GetWlab()/2. + sqrt(tk.GetWlab()*tk.GetWlab()/4.
					   +tk.GetQsquared()/2./(1.-fScatterAngle));
    if( fBeamEnergy <= tk.GetWlab() ) return false;
    ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
    break;

  case kScatterAngle:
    // Solve the kinematics starting with the scattering angle
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Find the beam energy
    if( fScatterAngle == 1. ) fBeamEnergy = tk.GetWlab() + STRANGEUFLOW;
    else
      fBeamEnergy = tk.GetWlab()/2. + sqrt( tk.GetWlab()*tk.GetWlab()/4.
					    +tk.GetQsquared()/2./(1.-fScatterAngle) );
    if( fBeamEnergy <= tk.GetWlab() ) return false;

    // Calculate epsilon
    if( fScatterAngle<= -1. ) fEpsilon = 0.;
    else if( tk.GetQsquared()<STRANGEUFLOW ) fEpsilon = 1.;
    else {
      ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
      fEpsilon = 1./(1.+2.*(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared())
		     *ftan2HalfAngle);
    }
    if( fEpsilon < 0. ) {
      if( fEpsilon < -STRANGEUFLOW ) return false;
      else fEpsilon = 0.;
    } else if( fEpsilon > 1. ) {
      if( (fEpsilon-1.) > STRANGEUFLOW ) return false;
      else fEpsilon = 1.;
    }
    ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
    break;
    
  default:
    cerr << "ERROR in TElectronKinematics::SolveKinematics(TKinematics2to2&): "
	 << "Invalid input variable.\n";
    assert(1==0);
  }
  
  return true;
}

void TElectronKinematics::GetLeptonVectors(const TKinematics2to2& kin, FourVector<double> &k_in, 
					 FourVector<double> &k_out){
  
  if( !SolveKinematics(kin) ) {
    cerr << "ERROR in TLeptonKinematics::GetLeptonVectors "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  
  double Eout=fBeamEnergy-kin.GetWlab();
  double pout=Eout;
  double kx=fBeamEnergy*pout*sqrt(1.-fScatterAngle*fScatterAngle)/kin.GetKlab();
  double kz=kin.GetWlab()*fBeamEnergy/kin.GetKlab() + (kin.GetQsquared())/(2.*kin.GetKlab());
  double kzprime=kin.GetWlab()*Eout/kin.GetKlab() + (-kin.GetQsquared())/(2.*kin.GetKlab());
  k_in[0]=fBeamEnergy; k_in[1]=k_out[1]=kx; k_in[2]=k_out[2]=0.;k_in[3]=kz;
  k_out[0]=Eout; k_out[3]=kzprime;
  
}

