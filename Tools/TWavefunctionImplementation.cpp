/*
 * TWavefunctionImplementation.cpp
 * 
 * Author: Pieter Vancraeyveld
 *         pieter.vancraeyveld@ugent.be
 *
 * Worked started February,5 2009
 */

////////////////////////////////////////////////////////////////////////
//                                                                    //
// TWavefunctionImplementation                                        //
//                                                                    //
// A wavefunction can be implemented in a number of different ways.   //
// TWavefunctionImplementation is the interface all wavefunction      //
// implementations should inherit from.                               //
//                                                                    //
// TWavefunctionImplementation is an abstract base class. It does     //
// not have an explicit constructor. We have defined a pure virtual   //
// copy constructor Clone() which all concrete derived classes should //
// implement.                                                         //
//                                                                    //
// All concrete wavefunction implementations should override the      //
// GetUr(double),GetWr(double),GetVSr(double),GetVTr(double),         //
// GetUp(double),GetWp(double),GetVSp(double),GetVTp(double)          //
// pure virtual member functions.                                     //
//                                                                    //
// We fix the following conventions concerning the units:             //
// * position in [fm]                                                 //
// * momentum in [MeV]                                                //
// * r-space wavefunctions in [fm^-1/2]                               //
// * p-space wavefunctions in [MeV^-3/2]                              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TWavefunctionImplementation.h"
#include <iostream>

using std::cout; using std::cerr; using std::endl;

ClassImp(TWavefunctionImplementation)

const double TWavefunctionImplementation::fgHbarc = 197.326968;

//_____________________________________________________________________
TWavefunctionImplementation::TWavefunctionImplementation(TRootIOCtor* rio)
  : TObject()
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
TWavefunctionImplementation::~TWavefunctionImplementation()
{
  // destructor
}

//_____________________________________________________________________
double TWavefunctionImplementation::Radial_p(int l, double p) const
{
  // Radial wavefunction in momentum space.
  // l = orbital angular momentum
  // p = relative momentum in [MeV]
  // Return value in units [MeV^-3/2]

  switch(l) {
  case 0: 
    return GetUp(p);
    break;
  case 1: 
    cerr << "WARNING in TWavefunctionImplementation::Radial_p(int,double): "
	 << "l=1 wavefunction is ambiguous.\n";
    break;
  case 2:
    return GetWp(p);
    break;
  default:
    cerr << "ERROR in TWavefunctionImplementation::Radial_p(int,double): "
	 << "l=" << l << " wavefunction unknown.\n";
  }
  
  return 0;
}

//_____________________________________________________________________
double TWavefunctionImplementation::Radial_r(int l, double r) const
{
  // Radial wavefunction in configuration space.
  // l = orbital angular momentum
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  switch(l) {
  case 0: 
    return GetUr(r);
    break;
  case 1: 
    cerr << "WARNING in TWavefunctionImplementation::Radial_r(int,double): "
	 << "l=1 wavefunction is ambiguous.\n";
    break;
  case 2:
    return GetWr(r);
    break;
  default:
    cerr << "ERROR in TWavefunctionImplementation::Radial_r(int,double): "
	 << "l=" << l << " wavefunction unknown.\n";
  }
  
  return 0;
}


