/*
 * TDeuteron.h
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: February,5 2009
 *
 */

#ifndef TDEUTERON_H
#define TDEUTERON_H

#include "TWavefunctionImplementation.h"
//#include "TKinematics2to3.h"
#include <TNamed.h>
#include <TString.h>
#include "FourVector.h"
#include <complex>
#include "TSpinor.h"

// Forward declarations
//class TStrangePolarization;
class TLorentzVector;

class TDeuteron
{
  ClassDef(TDeuteron,0); // Deuteron namespace

 public:
  // Forward declare the polarization class
  class Polarization;

  ///////////////////////////////////////////////////////////////////////////
  //                                                                       //
  // TDeuteron::Wavefunction                                               //
  //                                                                       //
  // This class represents a generic deuteron wavefunction in strangecalc. //
  //                                                                       //
  ///////////////////////////////////////////////////////////////////////////

 public:
  class Wavefunction : public TNamed
  {
  private:
    Wavefunction(const TString& name="", const TString& title="");
  public:
    Wavefunction(TRootIOCtor*); // ROOT I/O Constructor
    Wavefunction(const Wavefunction&);
    Wavefunction& operator=(const Wavefunction&);
    virtual ~Wavefunction();

    // named constructors
    static Wavefunction *CreateCDBonnWavefunction(); // user owns Wavefunction
    static Wavefunction *CreateParisWavefunction();  // user owns Wavefunction
    static Wavefunction *CreateAV18Wavefunction(); // user owns Wavefunction
    static Wavefunction *CreateAV18bWavefunction(); // user owns Wavefunction
    static Wavefunction *CreateGrossIIbWavefunction(); // user owns Wavefunction
    static Wavefunction *CreateGrossL4Wavefunction(); // user owns Wavefunction
    static Wavefunction *CreateGrossWJC1Wavefunction(); // user owns Wavefunction
    static Wavefunction *CreateGrossWJC2Wavefunction(); // user owns Wavefunction
    static Wavefunction *CreateNijmegenWavefunction(); // user owns Wavefunction
    static Wavefunction *CreateTestWavefunction(); // user owns Wavefunction
    static Wavefunction *CreateWavefunction(const TString& name); // user owns Wavefunction 
  
    double GetUr(double r) const { return fImplementation->GetUr(r); } // r in [fm]
    double GetWr(double r) const { return fImplementation->GetWr(r); } // r in [fm]
    double GetVTr(double r) const { return fImplementation->GetVTr(r); } // r in [fm]
    double GetVSr(double r) const { return fImplementation->GetVSr(r); } // r in [fm]
    double GetUp(double p) const { return fImplementation->GetUp(p); } // p in [MeV]
    double GetWp(double p) const { return fImplementation->GetWp(p); } // p in [MeV]
    double GetUpoff(const TVector3& p) const { return fImplementation->GetUpoff(p); } // p in [MeV]
    double GetWpoff1(const TVector3& p) const { return fImplementation->GetWpoff1(p); } // p in [MeV]
    double GetWpoff2(double pperp2) const { return fImplementation->GetWpoff2(pperp2); } // p in [MeV]
    double GetVTp(double p) const { return fImplementation->GetVTp(p); } // p in [MeV]
    double GetVSp(double p) const { return fImplementation->GetVSp(p); } // p in [MeV]
    double Radial_r(int l, double r) const { return fImplementation->Radial_r(l,r); } // r in [fm]
    double Radial_p(int l, double p) const { return fImplementation->Radial_p(l,p); } // p in [MeV]
  
    bool   IsRelativistic() const { return fRelativistic; }

    FourVector<GammaStructure> Vertex(const FourVector<double>& p2,
				      const FourVector<double>& pD) const;
    FourVector<GammaStructure> Vertex(const TLorentzVector& p2,
				      const TLorentzVector& pD) const { return Vertex(ToFourVector(p2),ToFourVector(pD)); }
    std::complex<double> DeuteronPState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TVector3& p) const; 
    std::complex<double> DeuteronPState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const FourVector<double>& p1) const; 
    std::complex<double> DeuteronPState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TLorentzVector& p1) const; 
    std::complex<double> DeuteronPStateOff(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TVector3& p) const; 
    std::complex<double> DeuteronPStateOff(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const FourVector<double>& p1) const; 
    std::complex<double> DeuteronPStateOff(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TLorentzVector& p1) const; 
//     std::complex<double> DeuteronPState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const TVector3& p1) const; 
//     std::complex<double> DeuteronPState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const FourVector<double>& p1) const; 
//     std::complex<double> DeuteronPState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const TLorentzVector& p1) const;
    std::complex<double> DeuteronRState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TVector3& r) const; 
    std::complex<double> DeuteronRState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const FourVector<double>& r1) const; 
    std::complex<double> DeuteronRState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TLorentzVector& r1) const; 
//     std::complex<double> DeuteronRState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const TVector3& r1) const; 
//     std::complex<double> DeuteronRState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const FourVector<double>& r1) const; 
//     std::complex<double> DeuteronRState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const TLorentzVector& r1) const;
    static std::complex<double> SphericalHarmonic(unsigned int l, int m, double theta, double phi);
    static std::complex<double> SphericalHarmonicCos(unsigned int l, int m, double costheta, double phi);
    static double ClebschGordan(int two_j1,int two_m1,int two_j2,int two_m2,int two_j,int two_m);
    
  private:
    bool                         fRelativistic;   // relativistic wavefunction?
    TWavefunctionImplementation *fImplementation; //-> implementation of wavefunction

    static const double          kMd;             // deuteron mass

    ClassDef(Wavefunction,1); // Strangecalc wavefunction
    
  }; // class TDeuteron::Wavefunction


  ///////////////////////////////////////////////////////////////////////////
  //                                                                       //
  // TDeuteron::Polarization                                               //
  //                                                                       //
  // Polarization state of a deuteron                                      //
  //                                                                       //
  ///////////////////////////////////////////////////////////////////////////

  class Polarization : public TObject
  {
  public:   
    Polarization(double,double,int);
    Polarization(TRootIOCtor*); // ROOT I/O Constructor
    Polarization(const Polarization&);
    virtual ~Polarization();
    Polarization& operator=(const Polarization&);

    bool operator==(const Polarization&);
    operator int() const;
    //TString HashName() const;

    void SetTheta(double theta);
    void SetPhi(double phi);
    void SetState(int polState);

    Polarization& operator++();
    Polarization& operator--();

    FourVector<std::complex<double> > GetPolarizationVector() const;

  private:
    double fTheta;             // polar angle in rad
    double fPhi;               // azimuthal angle in rad
    int    fPolarizationState; // polarization state =-1, 0 or 1

    ClassDef(Polarization,1); // Polarization state of a deuteron

  }; // class TDeuteron::Polarization

}; // class TDeuteron

#endif
