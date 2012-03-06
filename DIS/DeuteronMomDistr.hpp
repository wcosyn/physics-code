#ifndef DEUTERONMOMDISTR_HPP
#define DEUTERONMOMDISTR_HPP

#include <TVector3.h>
#include <TDeuteron.h>
#include <TKinematics2to2.h>

using namespace std;

class DeuteronMomDistr{
  
public:
  DeuteronMomDistr(string name, double massi, int offshell, double sigma, double beta, double epsilon, double betaoff, double lambda);
  ~DeuteronMomDistr();
  double getMomDistrpw(TKinematics2to2 &kin, double phi) const;
  double getMomDistrfsi(TKinematics2to2 &kin, double phi);
  
private:
  TDeuteron::Wavefunction *wf;
  double sigma;  //MeV-2
  double beta;  //MeV-2
  double epsilon;
  double betaoff;
  double lambda;
  double massi;
  double massr;
  
  double przprime;
  double Wxprime2;
  int offshellset;
  
  void totdens_qt(const double qt, complex<double>* result, va_list ap);
  void totdens_qphi(const double qphi, complex<double>* result, va_list ap);
  void get_przprime(double pt, double Er, TKinematics2to2 *pkin);
  complex<double> scatter(double t);
  
};

#endif