#ifndef INCLUSIVECROSSRES_HPP
#define INCLUSIVECROSSRES_HPP

#include <string>
#include <TElectronKinematics.h>
#include <DeuteronStructure.hpp>
#include <TDeuteron.h>

using namespace std;

//forward declaration
class Resonance;

class InclusiveCrossRes{
  
public:
  InclusiveCrossRes(bool proton, string strucname, string wavename, TElectronKinematics &elec, 
		    int symm, int offshell, bool fixprop, int t_choice);
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
  int t_choice;
  
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

#endif