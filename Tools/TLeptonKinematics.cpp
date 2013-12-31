/*
 * TLeptonKinematics.cpp
 *
 * Author: Wim Cosyn
 *
 * Work started on: Dec 12 2013
 *
 */

#include "TLeptonKinematics.h"
#include "TKinematics2to2.h"
#include "TKinematics2to3.h"
// #include <Structures.h>
#include "constants.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <gsl/gsl_poly.h>
#include <vector>
#include "Utilfunctions.hpp"

using namespace std;
using std::cerr;

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TLeptonKinematics                                                   //
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

ClassImp(TLeptonKinematics)

 const double TLeptonKinematics::masse= 0.510998910;
 const double TLeptonKinematics::massmu = 105.6583715;
 const double TLeptonKinematics::masstau = 1776.82;

//_____________________________________________________________________
TLeptonKinematics::TLeptonKinematics(Lepton type)
: fInput(kEpsilon), fBeamEnergy(0.), fScatterAngle(0.)
{
  // Protected Default Constructor
  //
  // Instances of TLeptonKinematics should be created with
  // the named constructors.
  switch(type){
    case(electron):
      mass=masse;
      break;
    case(muon):
      mass=massmu;
      break;
    case(tau):
      mass=masstau;
      break;
    default:
      cerr << "ERROR in TLeptonKinematics::constructor "
	<< "invalid lepton type.\n";
      assert(1==0);

  }
    
}

//_____________________________________________________________________
TLeptonKinematics::TLeptonKinematics(TRootIOCtor* rio)
  : fInput(kEpsilon), fBeamEnergy(0.), fScatterAngle(0.), mass(masse)
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
TLeptonKinematics::TLeptonKinematics(const TLeptonKinematics& rhs)
  : fInput(rhs.fInput), fBeamEnergy(rhs.fBeamEnergy), 
    fScatterAngle(rhs.fScatterAngle), mass(rhs.mass)
{
  // Copy Constructor
}

//_____________________________________________________________________
TLeptonKinematics& TLeptonKinematics::operator=(const TLeptonKinematics& rhs)
{
  // Assignment
  if( this!=&rhs ) {
    fInput = rhs.fInput;
    fBeamEnergy = rhs.fBeamEnergy;
//     fEpsilon = rhs.fEpsilon;
    fScatterAngle = rhs.fScatterAngle;
    mass = rhs.mass;
  }
  return *this;
}

//_____________________________________________________________________
TLeptonKinematics* TLeptonKinematics::CreateWithBeamEnergy(Lepton type, double beamEnergy)
{
  // Named constructor
  //
  // Fix the electron kinematics by specifying the energy
  // of the incoming electron beam in [MeV].
  //
  // The user is responsible for deleting the returned object.
  TLeptonKinematics *kinematics = new TLeptonKinematics(type);
  kinematics->SetBeamEnergy(beamEnergy);
  return kinematics;
}

//_____________________________________________________________________
// TLeptonKinematics* TLeptonKinematics::CreateWithEpsilon(Lepton type, double epsilon)
// {
//   // Named constructor
//   //
//   // Fix the electron kinematics by specifying the 
//   // virtual photon's transverse linear polarization.
//   //
//   // The user is responsible for deleting the returned object.
//   TLeptonKinematics *kinematics = new TLeptonKinematics(type);
//   kinematics->SetEpsilon(epsilon);
//   return kinematics;
// }

//_____________________________________________________________________
TLeptonKinematics* TLeptonKinematics::CreateWithCosScatterAngle(Lepton type, double scatterAngle)
{
  // Named constructor
  //
  // Fix the electron kinematics by specifying the cosine of the
  // electron's scattering angle in ]-1,1].
  //
  // The user is responsible for deleting the returned object.
  TLeptonKinematics *kinematics = new TLeptonKinematics(type);
  kinematics->SetCosScatterAngle(scatterAngle);
  return kinematics;
}

//_____________________________________________________________________
void TLeptonKinematics::SetBeamEnergy(double beamEnergy)
{
  // Fix the electron kinematics by specifying the energy
  // of the incoming electron beam in [MeV].
  fInput = kBeamEnergy;
  fBeamEnergy = beamEnergy;
}

//_____________________________________________________________________
// void TLeptonKinematics::SetEpsilon(double epsilon)
// {
//   // Fix the electron kinematics by specifying the
//   // virtual photon's transverse linear polarization.
//   fInput = kEpsilon;
//   fEpsilon = epsilon;
// }

//_____________________________________________________________________
void TLeptonKinematics::SetCosScatterAngle(double scatterAngle)
{
  // Fix the electron kinematics by specifying the cosine of the
  // electron's scattering angle in ]-1,1].
  fInput = kScatterAngle;
  fScatterAngle = scatterAngle;
}

//_____________________________________________________________________
double TLeptonKinematics::GetBeamEnergy(const TKinematics2to2& tk) const
{
  // Returns the incoming electron's energy in [MeV]
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TLeptonKinematics::GetBeamEnergy(TKinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fBeamEnergy;
}

double TLeptonKinematics::GetBeamEnergy(const TKinematics2to3& tk) const
{
  // Returns the incoming electron's energy in [MeV]
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TLeptonKinematics::GetBeamEnergy(TKinematics2to3&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fBeamEnergy;
}


//_____________________________________________________________________
// double TLeptonKinematics::GetEpsilon(const TKinematics2to2& tk) const
// {
//   // Returns the virtual photon's transverse linear polarization
//   // given the kinematics in the g*D system.
//   if( !SolveKinematics(tk) ) {
//     cerr << "ERROR in TLeptonKinematics::GetEpsilon(TKinematics2to2&) "
// 	 << "impossible kinematics.\n";
//     assert(1==0);
//   }
//   return fEpsilon;
// }
// 
// double TLeptonKinematics::GetEpsilon(const TKinematics2to3& tk) const
// {
//   // Returns the virtual photon's transverse linear polarization
//   // given the kinematics in the g*D system.
//   if( !SolveKinematics(tk) ) {
//     cerr << "ERROR in TLeptonKinematics::GetEpsilon(TKinematics2to3&) "
// 	 << "impossible kinematics.\n";
//     assert(1==0);
//   }
//   return fEpsilon;
// }

//_____________________________________________________________________
double TLeptonKinematics::GetCosScatterAngle(const TKinematics2to2& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TLeptonKinematics::GetCosScatterAngle(TKinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fScatterAngle; 
}

double TLeptonKinematics::GetCosScatterAngle(const TKinematics2to3& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in TLeptonKinematics::GetCosScatterAngle(TKinematics2to3&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fScatterAngle; 
}

// //_____________________________________________________________________
// double TLeptonKinematics::GetTan2HalfAngle(const TKinematics2to2& tk) const
// {
//   // Returns the cosine of the electron's scattering angle
//   // given the kinematics in the g*D system.
//   if( !SolveKinematics(tk) ) {
//     cerr << "ERROR in TLeptonKinematics::GetCosScatterAngle(TKinematics2to2&) "
// 	 << "impossible kinematics.\n";
//     assert(1==0);
//   }
//   return ftan2HalfAngle; 
// }
// 
// 
// //_____________________________________________________________________
// double TLeptonKinematics::GetTan2HalfAngle(const TKinematics2to3& tk) const
// {
//   // Returns the cosine of the electron's scattering angle
//   // given the kinematics in the g*D system.
//   if( !SolveKinematics(tk) ) {
//     cerr << "ERROR in TLeptonKinematics::GetCosScatterAngle(TKinematics2to2&) "
// 	 << "impossible kinematics.\n";
//     assert(1==0);
//   }
//   return ftan2HalfAngle; 
// }

//_____________________________________________________________________
bool TLeptonKinematics::SolveKinematics(const TKinematics2to2& tk) const
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
    fScatterAngle = (1. - (tk.GetQsquared()+mass*mass)/(2.*fBeamEnergy*(fBeamEnergy-tk.GetWlab())))
	    /sqrt(1-mass*mass/((fBeamEnergy-tk.GetWlab())*(fBeamEnergy-tk.GetWlab())));
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

//     // Finally calculate epsilon
//     if( fScatterAngle<= -1. ) fEpsilon = 0.;
//     else if( tk.GetQsquared()<STRANGEUFLOW ) fEpsilon = 1.;
//     else {
//       ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
//       fEpsilon = 1./(1.+2.*(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared())
// 		     *ftan2HalfAngle);
//     }
//     if( fEpsilon < 0. ) {
//       if( fEpsilon < -STRANGEUFLOW ) return false;
//       else fEpsilon = 0.;
//     } else if( fEpsilon > 1. ) {
//       if( (fEpsilon-1.) > STRANGEUFLOW ) return false;
//       else fEpsilon = 1.;
//     }
    break;

//   case kEpsilon:
//     // Solve the kinematics starting with epsilon.
//     if( fEpsilon<0. || fEpsilon>1. ) return false;
// 
//     // Find the scattering angle
//     if( tk.GetQsquared()<STRANGEUFLOW ) fScatterAngle = 1.;
//     else {
//       tmp = (1./fEpsilon-1.)/2./(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared());
//       fScatterAngle = (1.-tmp)/(1.+tmp);
//     }
//     if( fabs(fScatterAngle) > 1. ) {
//       if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
//      	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
//       else
// 	return false;
//     }
// 
//     // Calculate the beam energy
//     if( fScatterAngle > 1. ) return false;
//     else if( fScatterAngle == 1. ) fBeamEnergy = tk.GetWlab() + STRANGEUFLOW;
//     else
//       fBeamEnergy = tk.GetWlab()/2. + sqrt(tk.GetWlab()*tk.GetWlab()/4.
// 					   +tk.GetQsquared()/2./(1.-fScatterAngle));
//     if( fBeamEnergy <= tk.GetWlab() ) return false;
//     ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
//     break;

  case kScatterAngle:{
    // Solve the kinematics starting with the scattering angle
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Find the beam energy
//     if( fScatterAngle == 1. ) fBeamEnergy = tk.GetWlab() + STRANGEUFLOW;
    int count=0;
    if(abs(fScatterAngle)==1.){
      fBeamEnergy = (tk.GetWlab()+sqrt(tk.GetWlab()*tk.GetWlab()+(tk.GetQsquared())))*(tk.GetQsquared()+mass*mass)/(2.*tk.GetQsquared());
      double eout=fBeamEnergy-tk.GetWlab();
      double pout=sqrt(eout*eout-mass*mass);
      double a=(tk.GetQsquared()+mass*mass)/2.-fBeamEnergy*eout;
      double b=-fBeamEnergy*pout*fScatterAngle;
      if(SIGN(a)==SIGN(b)){
	count++;
// 	cout << tk.GetQsquared() << " " << -mass*mass+2.*fBeamEnergy*(eout-pout*fScatterAngle) << endl;
// 	cout << sqrt(pout*pout+fBeamEnergy*fBeamEnergy-2.*fScatterAngle*pout*fBeamEnergy) << " " << tk.GetKlab() << endl;      
      }
    }
    else if(fScatterAngle==0.){
      fBeamEnergy=(tk.GetWlab()+sqrt(tk.GetWlab()*tk.GetWlab()+2.*(tk.GetQsquared()+mass*mass)))/2.;
      count++;
    }
    else{
      //solve in GeV because otherwise the solver doesn't return the right solution for some reason...(overflows?)
      double sol[8];
      double coef[5];
      coef[0]=pow(tk.GetQsquared()+mass*mass,2.)/4.*1.E-12;
      coef[1]=(tk.GetQsquared()+mass*mass)*tk.GetWlab()*1.E-09;
      coef[2]=(-tk.GetQsquared()+(mass*mass-tk.GetWlab()*tk.GetWlab())*(fScatterAngle*fScatterAngle-1.))*1.E-06;
      coef[3]=2.*tk.GetWlab()*(fScatterAngle*fScatterAngle-1.)*1.E-03;
      coef[4]=1.-fScatterAngle*fScatterAngle;
      gsl_poly_complex_workspace * w
	  = gsl_poly_complex_workspace_alloc (5);
      
      int succes = gsl_poly_complex_solve (coef, 5, w, sol);
      if(succes!=0){      cerr << "ERROR in TLeptonKinematics::SolveKinematics(TKinematics2to2&): "
	  << "Gsl solver for lepton kinematics did not converge \n";
      assert(1==0);
      }


      gsl_poly_complex_workspace_free (w);
      for(int i=0;i<4;i++){
	if(abs(sol[i*2+1])<STRANGEUFLOW && sol[i*2] >=tk.GetWlab()*1.E-03){
	  double eout=sol[2*i]*1.E03-tk.GetWlab();
	  double pout=sqrt(eout*eout-mass*mass);
	  double a=(tk.GetQsquared()+mass*mass)/2.-sol[2*i]*1.E03*eout;
	  double b=-sol[2*i]*1.E03*pout*fScatterAngle;
//  	  cout << "sol " << i << " " << sol[i*2] << " " << sol[i*2+1] << " " << a << " " << b << endl;
	  if(SIGN(a)==SIGN(b)){
	    if(count==0){
	      count++;
	      fBeamEnergy=sol[2*i]*1.E03;	  
// 	      cout << tk.GetQsquared() << " " << -mass*mass+2.*fBeamEnergy*(eout-pout*fScatterAngle) << endl;
// 	      cout << sqrt(pout*pout+fBeamEnergy*fBeamEnergy-2.*fScatterAngle*pout*fBeamEnergy) << " " << tk.GetKlab() << endl;
	    }
	    else if(abs(sol[2*i]-fBeamEnergy*1.E-03)>1.E-06){
	      count++;
	      fBeamEnergy=sol[2*i]*1.E03;	  
// 	      cout << tk.GetQsquared() << " " << -mass*mass+2.*fBeamEnergy*(eout-pout*fScatterAngle) << endl;
// 	      cout << sqrt(pout*pout+fBeamEnergy*fBeamEnergy-2.*fScatterAngle*pout*fBeamEnergy) << " " << tk.GetKlab() << endl;
	    }
	  }
// 	  
	}
	
      }
    }

    if(count==0) return false;
    if(count>1){
    cerr << "ERROR in TLeptonKinematics::SolveKinematics(TKinematics2to2&): "
	<< "Multiple lepton solutions\n";
    assert(1==0);
    }
    if( fBeamEnergy <= tk.GetWlab() ) return false;

//     // Calculate epsilon
//     if( fScatterAngle<= -1. ) fEpsilon = 0.;
//     else if( tk.GetQsquared()<STRANGEUFLOW ) fEpsilon = 1.;
//     else {
//       ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
//       fEpsilon = 1./(1.+2.*(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared())
// 		     *ftan2HalfAngle);
//     }
//     if( fEpsilon < 0. ) {
//       if( fEpsilon < -STRANGEUFLOW ) return false;
//       else fEpsilon = 0.;
//     } else if( fEpsilon > 1. ) {
//       if( (fEpsilon-1.) > STRANGEUFLOW ) return false;
//       else fEpsilon = 1.;
//     }
//     ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
    break;
  }
  default:{
    cerr << "ERROR in TLeptonKinematics::SolveKinematics(TKinematics2to2&): "
	 << "Invalid input variable.\n";
    assert(1==0);
  }
  }
  return true;
}

//_____________________________________________________________________
bool TLeptonKinematics::SolveKinematics(const TKinematics2to3& tk) const
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
    fScatterAngle = (1. - (tk.GetQsquared()+mass*mass)/(2.*fBeamEnergy*(fBeamEnergy-tk.GetWlab())))
	    /sqrt(1-mass*mass/((fBeamEnergy-tk.GetWlab())*(fBeamEnergy-tk.GetWlab())));
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Finally calculate epsilon
//     if( fScatterAngle<= -1. ) fEpsilon = 0.;
//     else if( tk.GetQsquared()<STRANGEUFLOW ) fEpsilon = 1.;
//     else {
//       ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
//       fEpsilon = 1./(1.+2.*(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared())
// 		     *ftan2HalfAngle);
//     }
//     if( fEpsilon < 0. ) {
//       if( fEpsilon < -STRANGEUFLOW ) return false;
//       else fEpsilon = 0.;
//     } else if( fEpsilon > 1. ) {
//       if( (fEpsilon-1.) > STRANGEUFLOW ) return false;
//       else fEpsilon = 1.;
//     }
    break;

//   case kEpsilon:
//     // Solve the kinematics starting with epsilon.
//     if( fEpsilon<0. || fEpsilon>1. ) return false;
// 
//     // Find the scattering angle
//     if( tk.GetQsquared()<STRANGEUFLOW ) fScatterAngle = 1.;
//     else {
//       tmp = (1./fEpsilon-1.)/2./(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared());
//       fScatterAngle = (1.-tmp)/(1.+tmp);
//     }
//     if( fabs(fScatterAngle) > 1. ) {
//       if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
//      	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
//       else
// 	return false;
//     }
// 
//     // Calculate the beam energy
//     if( fScatterAngle > 1. ) return false;
//     else if( fScatterAngle == 1. ) fBeamEnergy = tk.GetWlab() + STRANGEUFLOW;
//     else
//       fBeamEnergy = tk.GetWlab()/2. + sqrt(tk.GetWlab()*tk.GetWlab()/4.
// 					   +tk.GetQsquared()/2./(1.-fScatterAngle));
//     if( fBeamEnergy <= tk.GetWlab() ) return false;
//     ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
//     break;

  case kScatterAngle:{
    // Solve the kinematics starting with the scattering angle
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < STRANGEUFLOW ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Find the beam energy
//     if( fScatterAngle == 1. ) fBeamEnergy = tk.GetWlab() + STRANGEUFLOW;
    int count=0;
    if(abs(fScatterAngle)==1.){
      fBeamEnergy = (tk.GetWlab()+sqrt(tk.GetWlab()*tk.GetWlab()+(tk.GetQsquared())))*(tk.GetQsquared()+mass*mass)/(2.*tk.GetQsquared());
      double eout=fBeamEnergy-tk.GetWlab();
      double pout=sqrt(eout*eout-mass*mass);
      double a=(tk.GetQsquared()+mass*mass)/2.-fBeamEnergy*eout;
      double b=-fBeamEnergy*pout*fScatterAngle;
      if(SIGN(a)==SIGN(b)){
	count++;
// 	cout << tk.GetQsquared() << " " << -mass*mass+2.*fBeamEnergy*(eout-pout*fScatterAngle) << endl;
// 	cout << sqrt(pout*pout+fBeamEnergy*fBeamEnergy-2.*fScatterAngle*pout*fBeamEnergy) << " " << tk.GetKlab() << endl;      
      }
    }
    else if(fScatterAngle==0.){
      fBeamEnergy=(tk.GetWlab()+sqrt(tk.GetWlab()*tk.GetWlab()+2.*(tk.GetQsquared()+mass*mass)))/2.;
      count++;
    }
    else{
      //solve in GeV because otherwise the solver doesn't return the right solution for some reason...(overflows?)
      double sol[8];
      double coef[5];
      coef[0]=pow(tk.GetQsquared()+mass*mass,2.)/4.*1.E-12;
      coef[1]=(tk.GetQsquared()+mass*mass)*tk.GetWlab()*1.E-09;
      coef[2]=(-tk.GetQsquared()+(mass*mass-tk.GetWlab()*tk.GetWlab())*(fScatterAngle*fScatterAngle-1.))*1.E-06;
      coef[3]=2.*tk.GetWlab()*(fScatterAngle*fScatterAngle-1.)*1.E-03;
      coef[4]=1.-fScatterAngle*fScatterAngle;
      gsl_poly_complex_workspace * w
	  = gsl_poly_complex_workspace_alloc (5);
      
      int succes = gsl_poly_complex_solve (coef, 5, w, sol);
      if(succes!=0){      cerr << "ERROR in TLeptonKinematics::SolveKinematics(TKinematics2to2&): "
	  << "Gsl solver for lepton kinematics did not converge \n";
      assert(1==0);
      }


      gsl_poly_complex_workspace_free (w);
      for(int i=0;i<4;i++){
	if(abs(sol[i*2+1])<STRANGEUFLOW && sol[i*2] >=tk.GetWlab()*1.E-03){
	  double eout=sol[2*i]*1.E03-tk.GetWlab();
	  double pout=sqrt(eout*eout-mass*mass);
	  double a=(tk.GetQsquared()+mass*mass)/2.-sol[2*i]*1.E03*eout;
	  double b=-sol[2*i]*1.E03*pout*fScatterAngle;
//  	  cout << "sol " << i << " " << sol[i*2] << " " << sol[i*2+1] << " " << a << " " << b << endl;
	  if(SIGN(a)==SIGN(b)){
	    if(count==0){
	      count++;
	      fBeamEnergy=sol[2*i]*1.E03;	  
// 	      cout << tk.GetQsquared() << " " << -mass*mass+2.*fBeamEnergy*(eout-pout*fScatterAngle) << endl;
// 	      cout << sqrt(pout*pout+fBeamEnergy*fBeamEnergy-2.*fScatterAngle*pout*fBeamEnergy) << " " << tk.GetKlab() << endl;
	    }
	    else if(abs(sol[2*i]-fBeamEnergy*1.E-03)>1.E-06){
	      count++;
	      fBeamEnergy=sol[2*i]*1.E03;	  
// 	      cout << tk.GetQsquared() << " " << -mass*mass+2.*fBeamEnergy*(eout-pout*fScatterAngle) << endl;
// 	      cout << sqrt(pout*pout+fBeamEnergy*fBeamEnergy-2.*fScatterAngle*pout*fBeamEnergy) << " " << tk.GetKlab() << endl;
	    }
	  }
// 	  
	}
	
      }
    }

    if(count==0) return false;
    if(count>1){
    cerr << "ERROR in TLeptonKinematics::SolveKinematics(TKinematics2to2&): "
	<< "Multiple lepton solutions\n";
    assert(1==0);
    }
    if( fBeamEnergy <= tk.GetWlab() ) return false;

//     // Calculate epsilon
//     if( fScatterAngle<= -1. ) fEpsilon = 0.;
//     else if( tk.GetQsquared()<STRANGEUFLOW ) fEpsilon = 1.;
//     else {
//       ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
//       fEpsilon = 1./(1.+2.*(1.+tk.GetWlab()*tk.GetWlab()/tk.GetQsquared())
// 		     *ftan2HalfAngle);
//     }
//     if( fEpsilon < 0. ) {
//       if( fEpsilon < -STRANGEUFLOW ) return false;
//       else fEpsilon = 0.;
//     } else if( fEpsilon > 1. ) {
//       if( (fEpsilon-1.) > STRANGEUFLOW ) return false;
//       else fEpsilon = 1.;
//     }
//     ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
    break;
  }
  default:{
    cerr << "ERROR in TLeptonKinematics::SolveKinematics(TKinematics2to2&): "
	 << "Invalid input variable.\n";
    assert(1==0);
  }
  }
  
  return true;
}

void TLeptonKinematics::GetLeptonVectors(const TKinematics2to2& kin, FourVector<double> &k_in, 
					 FourVector<double> &k_out){
  
  if( !SolveKinematics(kin) ) {
    cerr << "ERROR in TLeptonKinematics::GetLeptonVectors "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  
  double Eout=fBeamEnergy-kin.GetWlab();
  double pout=sqrt(Eout*Eout-mass*mass);
  double kx=fBeamEnergy*pout*sqrt(1.-fScatterAngle*fScatterAngle)/kin.GetKlab();
  double kz=kin.GetWlab()*fBeamEnergy/kin.GetKlab() + (mass*mass+kin.GetQsquared())/(2.*kin.GetKlab());
  double kzprime=kin.GetWlab()*Eout/kin.GetKlab() + (mass*mass-kin.GetQsquared())/(2.*kin.GetKlab());
  k_in[0]=fBeamEnergy; k_in[1]=k_out[1]=kx; k_in[2]=k_out[2]=0.;k_in[3]=kz;
  k_out[0]=Eout; k_out[3]=kzprime;
  
}
