#ifndef INCLUSIVECROSS_HPP
#define INCLUSIVECROSS_HPP

#include <string>
#include <TElectronKinematics.h>
#include <DeuteronStructure.hpp>
#include <TDeuteron.h>

using namespace std;

//forward declaration
class Resonance;


class InclusiveCross{
  
public:
  InclusiveCross(bool proton, string strucname, string wavename, TElectronKinematics &elec, int symm, int offshell, 
		 double sigmain=40.,double betain=8., double epsilonin=-0.5, double betaoffin=8., double lambdain=1.2);
  ~InclusiveCross();
  double calc_F2Dinc(double Q2,double x);
  void calc_F2DincFSI(double &fsi1, double &fsi2, double Q2,double x);
  //void calc_F2DincFSI2(double &fsi1, double &fsi2, double Q2,double x);
  
  
private:
  double sigma;  //MeV-2
  double beta;  //MeV-2
  double epsilon;
  double betaoff;
  double lambda;
  double massi;
  double massr;
  int symm;
  
  double przprime;
  double prz;
  double Wxprime2; //invariant mass of the X that is paired with the spectator that gets integrated first
  double otherWx2; //invariant mass of the other X that depends on the first ps integration
  int offshellset;

  TDeuteron::Wavefunction *wf;
  TElectronKinematics electron;
  DeuteronStructure structure;  

  void int_pr(double pr, double *result, va_list ap);
  void int_costheta_incl(double costheta, double *result, va_list ap);
  void int_pr_fsi(double pr, double *result, va_list ap);
  void int_costheta_incl_fsi(double costheta, double *result, va_list ap);
  void int_qt(double qt, double *results, va_list ap);
  void int_qphi(double qphi, double *results, va_list ap); 
  void get_prz(double pt, double Er, TKinematics2to2 *pkin, int first);
  complex<double> scatter(double t);

//   void int_pperp_fsi(double pperp, double *result, va_list ap);
//   void int_qt_bis(double qt, double *results, va_list ap);
//   void int_qphi_bis(double qphi, double *results, va_list ap); 
//   void get_przs(double pperp2, double pperp2other, double qvec, double nu, int first);
};


class InclusiveCrossRes{
  
public:
  InclusiveCrossRes(bool proton, string strucname, string wavename, TElectronKinematics &elec, int symm, int offshell, bool fixprop);
  ~InclusiveCrossRes();
  double calc_F2Dinc(double Q2,double x);
  void calc_F2DincFSI(double &fsi1, double &fsi2, double Q2,double x);
  //void calc_F2DincFSI2(double &fsi1, double &fsi2, double Q2,double x);
  void addResonance(Resonance &res);   
  vector<Resonance> getResonance_vec() const{return resonance_vec;}
  
private:
  double massi;
  double massr;
  int symm;
  bool fixprop;
  
  double przprime;
  double prz;
  double Wxprime2; //invariant mass of the X that is paired with the spectator that gets integrated first
  double otherWx2; //invariant mass of the other X that depends on the first ps integration
  int offshellset;

  TDeuteron::Wavefunction *wf;
  TElectronKinematics electron;
  DeuteronStructure structure;  

  vector<Resonance> resonance_vec;
  
  void int_pr(double pr, double *result, va_list ap);
  void int_costheta_incl(double costheta, double *result, va_list ap);
  void int_pr_fsi(double pr, double *result, va_list ap);
  void int_costheta_incl_fsi(double costheta, double *result, va_list ap);
  void int_qt(double qt, double *results, va_list ap);
  void int_qphi(double qphi, double *results, va_list ap); 
  void get_prz(double pt, double Er, TKinematics2to2 *pkin, int first, size_t res1, size_t res2);
  complex<double> scatter(double t, double Q2, size_t res1, size_t res2);
};


class Resonance{
public:
  Resonance(double m, complex<double> coeff, double sigma0=40., double sigmaslope=0.,
	    double betain=6.,double epsilonin=-0.5,double betaoffin=8., double lambdain=1.2);
  Resonance(const Resonance &copy);
  Resonance& operator=(const Resonance&);
  ~Resonance();
  double getSigma(double Q2) const;
  double getBeta(double Q2=0.) const{return beta;}
  double getBetaoff() const{return betaoff;}
  double getSigma0() const{return sigma0;}
  double getSigmaslope() const{return sigmaslope;}
  double getEpsilon() const{return epsilon;}
  double getLambda() const{return lambda;}
  double getMass2() const{return mass2;}
  double getMass() const{return mass;}
  complex<double> getCoeff() const{return coeff;}
  
private:
  double mass;
  double mass2;
  double sigma0;
  double sigmaslope;
  double beta;
  double betaoff;
  double epsilon;
  double lambda;
  complex<double> coeff;  
  
};
#endif