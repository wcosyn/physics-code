/*
 * TInterpolatingWavefunction.cpp
 *
 * Author: Pieter Vancraeyveld
 *         pieter.vancraeyveld@ugent.be
 *
 * Worked started February,16 2009
 *
 */

////////////////////////////////////////////////////////////////////////
//                                                                    //
// TInterpolatingWavefunction                                         //
//                                                                    //
// TInterpolatingWavefunction implements the abstract interface       //
// TWavefunctionImplementation.                                       //
//                                                                    //
// This class calculates (non-)relativistic wavefunctions by          //
// interpolating (linear) between a list of tabulated values.         //
//                                                                    //
// TWavefunctionImplementation fixes the following conventions        //
// concerning the units:                                              //
// * position in [fm]                                                 //
// * momentum in [MeV]                                                //
// * r-space wavefunctions in [fm^-1/2]                               //
// * p-space wavefunctions in [MeV^-3/2]                              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TInterpolatingWavefunction.h"
#include <iostream>
#include <TDeuteron.h>

using std::cout; using std::cerr; using std::endl;
using std::map; using std::pair;

ClassImp(TInterpolatingWavefunction)

//_____________________________________________________________________
TInterpolatingWavefunction::TInterpolatingWavefunction()
: TWavefunctionImplementation(0),
  fUr(0), fWr(0), fVTr(0), fVSr(0), fUp(0), fWp(0), fVTp(0), fVSp(0)
{
  // Constructor
}

//_____________________________________________________________________
TInterpolatingWavefunction::TInterpolatingWavefunction(const TInterpolatingWavefunction& rhs)
: TWavefunctionImplementation(0),
  fUr(0), fWr(0), fVTr(0), fVSr(0), fUp(0), fWp(0), fVTp(0), fVSp(0)
{
  // Copy constructor

  if(rhs.fUr) fUr = new map<double,double>(*rhs.fUr);
  if(rhs.fWr) fWr = new map<double,double>(*rhs.fWr);
  if(rhs.fVTr) fVTr = new map<double,double>(*rhs.fVTr);
  if(rhs.fVSr) fVSr = new map<double,double>(*rhs.fVSr);

  if(rhs.fUp) fUp = new map<double,double>(*rhs.fUp);
  if(rhs.fWp) fWp = new map<double,double>(*rhs.fWp);
  if(rhs.fVTp) fVTp = new map<double,double>(*rhs.fVTp);
  if(rhs.fVSp) fVSp = new map<double,double>(*rhs.fVSp);
}

//_____________________________________________________________________
TInterpolatingWavefunction& TInterpolatingWavefunction::operator=(const TInterpolatingWavefunction& rhs)
{
  // Assignment

  if(this!=&rhs) { // avoid self-assignment

    if(rhs.fUr) {
      if(fUr) *fUr = *rhs.fUr;
      else fUr = new map<double,double>(*rhs.fUr);
    } else {
      delete fUr;
      fUr = 0;
    }

    if(rhs.fWr) {
      if(fWr) *fWr = *rhs.fWr;
      else fWr = new map<double,double>(*rhs.fWr);
    } else {
      delete fWr;
      fWr = 0;
    }

    if(rhs.fVTr) {
      if(fVTr) *fVTr = *rhs.fVTr;
      else fVTr = new map<double,double>(*rhs.fVTr);
    } else {
      delete fVTr;
      fVTr = 0;
    }

    if(rhs.fVSr) {
      if(fVSr) *fVSr = *rhs.fVSr;
      else fVSr = new map<double,double>(*rhs.fVSr);
    } else {
      delete fVSr;
      fVSr = 0;
    }

    if(rhs.fUp) {
      if(fUp) *fUp = *rhs.fUp;
      else fUp = new map<double,double>(*rhs.fUp);
    } else {
      delete fUp;
      fUp = 0;
    }

    if(rhs.fWp) {
      if(fWp) *fWp = *rhs.fWp;
      else fWp = new map<double,double>(*rhs.fWp);
    } else {
      delete fWp;
      fWp = 0;
    }

    if(rhs.fVTp) {
      if(fVTp) *fVTp = *rhs.fVTp;
      else fVTp = new map<double,double>(*rhs.fVTp);
    } else {
      delete fVTp;
      fVTp = 0;
    }

    if(rhs.fVSp) {
      if(fVSp) *fVSp = *rhs.fVSp;
      else fVSp = new map<double,double>(*rhs.fVSp);
    } else {
      delete fVSp;
      fVSp = 0;
    }
  }

  return *this;
}

//_____________________________________________________________________
TInterpolatingWavefunction::~TInterpolatingWavefunction()
{
  // Destructor

  delete fUr;
  delete fWr;
  delete fVTr;
  delete fVSr;
  
  delete fUp;
  delete fWp;
  delete fVTp;
  delete fVSp;
}

//_____________________________________________________________________
void TInterpolatingWavefunction::RemoveUr()
{
  // Clear list of tabulated values
  if(fUr) fUr->clear();
}

//_____________________________________________________________________
void TInterpolatingWavefunction::RemoveWr()
{
  // Clear list of tabulated values
  if(fWr) fWr->clear();
}

//_____________________________________________________________________
void TInterpolatingWavefunction::RemoveVTr()
{
  // Clear list of tabulated values
  if(fVTr) fVTr->clear();
}

//_____________________________________________________________________
void TInterpolatingWavefunction::RemoveVSr()
{
  // Clear list of tabulated values
  if(fVSr) fVSr->clear();
}

//_____________________________________________________________________
void TInterpolatingWavefunction::RemoveUp()
{
  // Clear list of tabulated values
  if(fUp) fUp->clear();
}

//_____________________________________________________________________
void TInterpolatingWavefunction::RemoveWp()
{
  // Clear list of tabulated values
  if(fWp) fWp->clear();
}

//_____________________________________________________________________
void TInterpolatingWavefunction::RemoveVTp()
{
  // Clear list of tabulated values
  if(fVTp) fVTp->clear();
}

//_____________________________________________________________________
void TInterpolatingWavefunction::RemoveVSp()
{
  // Clear list of tabulated values
  if(fVSp) fVSp->clear();
}

//_____________________________________________________________________
void TInterpolatingWavefunction::SetUr(int n, double *r, double *list)
{
  // Add 'n' L=0 wavefunction values.
  // r in [fm]
  // list in [fm^-1/2]

  for(int i=0; i<n; ++i)
    AddUr(r[i],list[i]);  
}

//_____________________________________________________________________
void TInterpolatingWavefunction::SetWr(int n, double *r, double *list)
{
  // Add 'n' L=2 wavefunction values.
  // r in [fm]
  // list in [fm^-1/2]

  for(int i=0; i<n; ++i)
    AddWr(r[i],list[i]);  
}

//_____________________________________________________________________
void TInterpolatingWavefunction::SetVTr(int n, double *r, double *list)
{
  // Add 'n' L=1 triplet wavefunction values.
  // r in [fm]
  // list in [fm^-1/2]

  for(int i=0; i<n; ++i)
    AddVTr(r[i],list[i]);  
}

//_____________________________________________________________________
void TInterpolatingWavefunction::SetVSr(int n, double *r, double *list)
{
  // Add 'n' L=1 singlet wavefunction values.
  // r in [fm]
  // list in [fm^-1/2]

  for(int i=0; i<n; ++i)
    AddVSr(r[i],list[i]);  
}

//_____________________________________________________________________
void TInterpolatingWavefunction::SetUp(int n, double *p, double *list)
{
  // Add 'n' L=0 wavefunction values.
  // p in [MeV]
  // list in [MeV^-3/2]

  for(int i=0; i<n; ++i)
    AddUp(p[i],list[i]);  
}

//_____________________________________________________________________
void TInterpolatingWavefunction::SetWp(int n, double *p, double *list)
{
  // Add 'n' L=2 wavefunction values.
  // p in [MeV]
  // list in [MeV^-3/2]

  for(int i=0; i<n; ++i)
    AddWp(p[i],list[i]);  
}

//_____________________________________________________________________
void TInterpolatingWavefunction::SetVTp(int n, double *p, double *list)
{
  // Add 'n' L=1 triplet wavefunction values.
  // p in [MeV]
  // list in [MeV^-3/2]

  for(int i=0; i<n; ++i)
    AddVTp(p[i],list[i]);  
}

//_____________________________________________________________________
void TInterpolatingWavefunction::SetVSp(int n, double *p, double *list)
{
  // Add 'n' L=1 singlet wavefunction values.
  // p in [MeV]
  // list in [MeV^-3/2]

  for(int i=0; i<n; ++i)
    AddVSp(p[i],list[i]);  
}

//_____________________________________________________________________
void TInterpolatingWavefunction::AddUr(double r, double value)
{
  // Add L=0 wavefunction value.
  // r in [fm]
  // wavefunction in [fm^-1/2]

  if(!fUr) fUr = new map<double,double>;

  fUr->insert( pair<double,double>(r,value) );
}

//_____________________________________________________________________
void TInterpolatingWavefunction::AddWr(double r, double value)
{
  // Add L=2 wavefunction value.
  // r in [fm]
  // wavefunction in [fm^-1/2]

  if(!fWr) fWr = new map<double,double>;

  fWr->insert( pair<double,double>(r,value) );
}

//_____________________________________________________________________
void TInterpolatingWavefunction::AddVTr(double r, double value)
{
  // Add L=1 triplet wavefunction value.
  // r in [fm]
  // wavefunction in [fm^-1/2]

  if(!fVTr) fVTr = new map<double,double>;

  fVTr->insert( pair<double,double>(r,value) );
}

//_____________________________________________________________________
void TInterpolatingWavefunction::AddVSr(double r, double value)
{
  // Add L=1 singlet wavefunction value.
  // r in [fm]
  // wavefunction in [fm^-1/2]

  if(!fVSr) fVSr = new map<double,double>;

  fVSr->insert( pair<double,double>(r,value) );
}

//_____________________________________________________________________
void TInterpolatingWavefunction::AddUp(double p, double value)
{
  // Add L=0 wavefunction value.
  // p in [MeV]
  // wavefunction in [MeV^-3/2]

  if(!fUp) fUp = new map<double,double>;

  fUp->insert( pair<double,double>(p,value) );
}

//_____________________________________________________________________
void TInterpolatingWavefunction::AddWp(double p, double value)
{
  // Add L=2 wavefunction value.
  // p in [MeV]
  // wavefunction in [MeV^-3/2]

  if(!fWp) fWp = new map<double,double>;

  fWp->insert( pair<double,double>(p,value) );
}

//_____________________________________________________________________
void TInterpolatingWavefunction::AddVTp(double p, double value)
{
  // Add L=1 triplet wavefunction value.
  // p in [MeV]
  // wavefunction in [MeV^-3/2]

  if(!fVTp) fVTp = new map<double,double>;

  fVTp->insert( pair<double,double>(p,value) );
}

//_____________________________________________________________________
void TInterpolatingWavefunction::AddVSp(double p, double value)
{
  // Add L=1 singlet wavefunction value.
  // p in [MeV]
  // wavefunction in [MeV^-3/2]

  if(!fVSp) fVSp = new map<double,double>;

  fVSp->insert( pair<double,double>(p,value) );
}

//_____________________________________________________________________
double TInterpolatingWavefunction::GetUr(double r) const
{
  // L=0 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  return Interpolate(fUr,r);
}

//_____________________________________________________________________
double TInterpolatingWavefunction::GetWr(double r) const
{
  // L=2 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  return Interpolate(fWr,r);
}

//_____________________________________________________________________
double TInterpolatingWavefunction::GetVTr(double r) const
{
  // L=1 triplet wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  return Interpolate(fVTr,r);
}

//_____________________________________________________________________
double TInterpolatingWavefunction::GetVSr(double r) const
{
  // L=1 singlet wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  return Interpolate(fVSr,r);
}

//_____________________________________________________________________
double TInterpolatingWavefunction::GetUp(double p) const
{
  // L=0 wave in momentum space
  // p = relative position in [MeV]
  // Return value in units [MeV^-3/2]

  return Interpolate(fUp,p);
}

//_____________________________________________________________________
double TInterpolatingWavefunction::GetWp(double p) const
{
  // L=2 wave in momentum space
  // p = relative position in [MeV]
  // Return value in units [MeV^-3/2]

  return Interpolate(fWp,p);
}

//_____________________________________________________________________
double TInterpolatingWavefunction::GetVTp(double p) const
{
  // L=1 triplet wave in momentum space
  // p = relative position in [MeV]
  // Return value in units [MeV^-3/2]

  return Interpolate(fVTp,p);
}

//_____________________________________________________________________
double TInterpolatingWavefunction::GetVSp(double p) const
{
  // L=1 singlet wave in momentum space
  // p = relative position in [MeV]
  // Return value in units [MeV^-3/2]

  return Interpolate(fVSp,p);
}

double TInterpolatingWavefunction::GetUpoff(const TVector3& p) const{
 
  cerr << "Off-shell wave function not implemented in Interpolating wave function class" << endl;
  return 0.;
  
}

double TInterpolatingWavefunction::GetWpoff1(const TVector3& p) const{
 
  cerr << "Off-shell wave function not implemented in Interpolating wave function class" << endl;
  return 0.;
  
}

double TInterpolatingWavefunction::GetWpoff2(double pperp2) const{
 
  cerr << "Off-shell wave function not implemented in Interpolating wave function class" << endl;
  return 0.;
  
}


//_____________________________________________________________________
double TInterpolatingWavefunction::Interpolate(map<double,double>* list,
					       double q) const
{
  // Find the wavefunction value corresponding to 'q', by interpolating
  // between the known values found in the ordered 'list'
  if(std::isnan(q)) return 0.;
  if( !list || list->size()<2 ) {
    cerr << "WARNING in TInterpolatingWavefunction::Interpolate(map*,double): "
	 << "nothing to interpolate!\n";
    return 0;
  }

  // Find point v2=(q2,wf2) such that q < q2
  map<double,double>::iterator v2 = list->upper_bound(q);
  
  // check whether it is within bounds
  if( v2==list->end() ) {
    --v2; // move iterator one back
    if( v2->first!=q ) {
      cerr << q << " WARNING in TInterpolatingWavefunction::Interpolate(map*,double): "
	   << "requested point outside range of tabulated values.\n";
      return 0;
    }
  }
  
  if( v2==list->begin() ) {
    cerr << q << " WARNING in TInterpolatingWavefunction::Interpolate(map*,double): "
	 << "requested point outside range of tabulated values.\n";
    return 0;
  }
  
  // Find point v1=(q1,wf1) such that q1 <= q
  map<double,double>::iterator v1 = v2;
  --v1;

  // Now interpolation is children's play
  return v1->second + (q-v1->first) 
    * (v2->second-v1->second) / (v2->first-v1->first);
}

//______________________________________________________________________________
void TInterpolatingWavefunction::Streamer(TBuffer &R__b){}
// {
//   // Stream an object of class TInterpolatingWavefunction.

//   // We implement it ourselves because ROOT seems to have trouble
//   // with the pointers to map.

//    UInt_t R__s, R__c;
//    if (R__b.IsReading()) {
//       Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
//       TWavefunctionImplementation::Streamer(R__b);
//       {
//          delete fUr;
//          fUr = new map<double,double>;
//          map<double,double> &R__stl = *fUr;
//          R__stl.clear();
//          int R__i, R__n;
//          R__b >> R__n;
//          for (R__i = 0; R__i < R__n; R__i++) {
//             double R__t;
//             R__b >> R__t;
//             double R__t2;
//             R__b >> R__t2;
//             typedef double Value_t;
//             std::pair<Value_t const, double > R__t3(R__t,R__t2);
//             R__stl.insert(R__t3);
//          }
//       }
//       {
//          delete fWr;
//          fWr = new map<double,double>;
//          map<double,double> &R__stl = *fWr;
//          R__stl.clear();
//          int R__i, R__n;
//          R__b >> R__n;
//          for (R__i = 0; R__i < R__n; R__i++) {
//             double R__t;
//             R__b >> R__t;
//             double R__t2;
//             R__b >> R__t2;
//             typedef double Value_t;
//             std::pair<Value_t const, double > R__t3(R__t,R__t2);
//             R__stl.insert(R__t3);
//          }
//       }
//       {
//          delete fVTr;
//          fVTr = new map<double,double>;
//          map<double,double> &R__stl = *fVTr;
//          R__stl.clear();
//          int R__i, R__n;
//          R__b >> R__n;
//          for (R__i = 0; R__i < R__n; R__i++) {
//             double R__t;
//             R__b >> R__t;
//             double R__t2;
//             R__b >> R__t2;
//             typedef double Value_t;
//             std::pair<Value_t const, double > R__t3(R__t,R__t2);
//             R__stl.insert(R__t3);
//          }
//       }
//       {
//          delete fVSr;
//          fVSr = new map<double,double>;
//          map<double,double> &R__stl = *fVSr;
//          R__stl.clear();
//          int R__i, R__n;
//          R__b >> R__n;
//          for (R__i = 0; R__i < R__n; R__i++) {
//             double R__t;
//             R__b >> R__t;
//             double R__t2;
//             R__b >> R__t2;
//             typedef double Value_t;
//             std::pair<Value_t const, double > R__t3(R__t,R__t2);
//             R__stl.insert(R__t3);
//          }
//       }
//       {
//          delete fUp;
//          fUp = new map<double,double>;
//          map<double,double> &R__stl = *fUp;
//          R__stl.clear();
//          int R__i, R__n;
//          R__b >> R__n;
//          for (R__i = 0; R__i < R__n; R__i++) {
//             double R__t;
//             R__b >> R__t;
//             double R__t2;
//             R__b >> R__t2;
//             typedef double Value_t;
//             std::pair<Value_t const, double > R__t3(R__t,R__t2);
//             R__stl.insert(R__t3);
//          }
//       }
//       {
//          delete fWp;
//          fWp = new map<double,double>;
//          map<double,double> &R__stl = *fWp;
//          R__stl.clear();
//          int R__i, R__n;
//          R__b >> R__n;
//          for (R__i = 0; R__i < R__n; R__i++) {
//             double R__t;
//             R__b >> R__t;
//             double R__t2;
//             R__b >> R__t2;
//             typedef double Value_t;
//             std::pair<Value_t const, double > R__t3(R__t,R__t2);
//             R__stl.insert(R__t3);
//          }
//       }
//       {
//          delete fVTp;
//          fVTp = new map<double,double>;
//          map<double,double> &R__stl = *fVTp;
//          R__stl.clear();
//          int R__i, R__n;
//          R__b >> R__n;
//          for (R__i = 0; R__i < R__n; R__i++) {
//             double R__t;
//             R__b >> R__t;
//             double R__t2;
//             R__b >> R__t2;
//             typedef double Value_t;
//             std::pair<Value_t const, double > R__t3(R__t,R__t2);
//             R__stl.insert(R__t3);
//          }
//       }
//       {
//          delete fVSp;
//          fVSp = new map<double,double>;
//          map<double,double> &R__stl = *fVSp;
//          R__stl.clear();
//          int R__i, R__n;
//          R__b >> R__n;
//          for (R__i = 0; R__i < R__n; R__i++) {
//             double R__t;
//             R__b >> R__t;
//             double R__t2;
//             R__b >> R__t2;
//             typedef double Value_t;
//             std::pair<Value_t const, double > R__t3(R__t,R__t2);
//             R__stl.insert(R__t3);
//          }
//       }
//       R__b.CheckByteCount(R__s, R__c, TInterpolatingWavefunction::IsA());
//    } else {
//       R__c = R__b.WriteVersion(TInterpolatingWavefunction::IsA(), kTRUE);
//       TWavefunctionImplementation::Streamer(R__b);
//       {
//          map<double,double> &R__stl = *fUr;
//          int R__n=(&R__stl) ? int(R__stl.size()) : 0;
//          R__b << R__n;
//          if(R__n) {
//             map<double,double>::iterator R__k;
//             for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
//             R__b << ((*R__k).first );
//             R__b << ((*R__k).second);
//             }
//          }
//       }
//       {
//          map<double,double> &R__stl = *fWr;
//          int R__n=(&R__stl) ? int(R__stl.size()) : 0;
//          R__b << R__n;
//          if(R__n) {
//             map<double,double>::iterator R__k;
//             for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
//             R__b << ((*R__k).first );
//             R__b << ((*R__k).second);
//             }
//          }
//       }
//       {
//          map<double,double> &R__stl = *fVTr;
//          int R__n=(&R__stl) ? int(R__stl.size()) : 0;
//          R__b << R__n;
//          if(R__n) {
//             map<double,double>::iterator R__k;
//             for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
//             R__b << ((*R__k).first );
//             R__b << ((*R__k).second);
//             }
//          }
//       }
//       {
//          map<double,double> &R__stl = *fVSr;
//          int R__n=(&R__stl) ? int(R__stl.size()) : 0;
//          R__b << R__n;
//          if(R__n) {
//             map<double,double>::iterator R__k;
//             for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
//             R__b << ((*R__k).first );
//             R__b << ((*R__k).second);
//             }
//          }
//       }
//       {
//          map<double,double> &R__stl = *fUp;
//          int R__n=(&R__stl) ? int(R__stl.size()) : 0;
//          R__b << R__n;
//          if(R__n) {
//             map<double,double>::iterator R__k;
//             for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
//             R__b << ((*R__k).first );
//             R__b << ((*R__k).second);
//             }
//          }
//       }
//       {
//          map<double,double> &R__stl = *fWp;
//          int R__n=(&R__stl) ? int(R__stl.size()) : 0;
//          R__b << R__n;
//          if(R__n) {
//             map<double,double>::iterator R__k;
//             for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
//             R__b << ((*R__k).first );
//             R__b << ((*R__k).second);
//             }
//          }
//       }
//       {
//          map<double,double> &R__stl = *fVTp;
//          int R__n=(&R__stl) ? int(R__stl.size()) : 0;
//          R__b << R__n;
//          if(R__n) {
//             map<double,double>::iterator R__k;
//             for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
//             R__b << ((*R__k).first );
//             R__b << ((*R__k).second);
//             }
//          }
//       }
//       {
//          map<double,double> &R__stl = *fVSp;
//          int R__n=(&R__stl) ? int(R__stl.size()) : 0;
//          R__b << R__n;
//          if(R__n) {
//             map<double,double>::iterator R__k;
//             for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
//             R__b << ((*R__k).first );
//             R__b << ((*R__k).second);
//             }
//          }
//       }
//       R__b.SetByteCount(R__c, kTRUE);
//    }
// }


std::complex<double> TInterpolatingWavefunction::DeuteronPState(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const TVector3& p) const
{
  // Calculate the state of a deuteron in momentum space. 
  // deuteronPol is the deuteron polarization (times two)
  // nucleon1pol fixes the polarization of the first nucleon (the one that interacts with the photon). (times two)
  //nucleon2pol fixes spectator nucleon polarization (times two)
  
  // The Dnp vertex matrix element is
  // BEGIN_LATEX
  // #frac{1}{2E_{1}} #bar{u}(p_{2},m_{2})u(#dot{p}_{1},m_{1})#Gamma_{#mu}#phi^{#mu}_{D}(m_{D})
  // = #sum_{L} u_{L}(|#vec{p}_{1}|) #sum_{m_{L},m_{S}} <Lm_{L}Sm_{S}|1m_{D}><#frac{1}{2}m_{1}#frac{1}{2} m_{2}|1m_{S}> Y_{Lm_{L}}(#hat{p}_{1})
  // END_LATEX
  // with BEGIN_LATEX #dot{p}_{1} = ( #sqrt{#vec{p}_{1}^{2}+M_{1}^{2}},#vec{p}_{1}) END_LATEX the on-shell 4-momentum of particle 1.
  // The labelling particle 1 and 2 is arbitrary in the non-relativistic case (You can do the math).

  
  std::complex<double> state=0.;
// L=2, mL=mD-mS, S=1, mS=m1+m2
  for(int mS=-1; mS<=1; ++mS)
    state -=
      TDeuteron::Wavefunction::ClebschGordan(4,deuteronPol-mS*2,2,mS*2,2,deuteronPol)
      *TDeuteron::Wavefunction::ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,mS*2)
      *TDeuteron::Wavefunction::SphericalHarmonicCos(2,deuteronPol/2-mS,p.CosTheta(),p.Phi());
  state*=GetWp(p.Mag());
  // L=0, mL=0, S=1, mS=mD
  state += 
    TDeuteron::Wavefunction::ClebschGordan(0,0,2,deuteronPol,2,deuteronPol)
    *TDeuteron::Wavefunction::ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,deuteronPol)
    *TDeuteron::Wavefunction::SphericalHarmonicCos(0,0,p.CosTheta(),p.Phi())
    *GetUp(p.Mag());

  return state;
}
