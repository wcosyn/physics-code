#ifndef DOUBLENMODEL_HPP
#define DOUBLENMODEL_HPP

#include <FourVector.h>
#include <MeanFieldNucleusThick.hpp>
#include <GlauberGridThick.hpp>
#include <TKinematics2to3WithLabAngles.h>
#include <NucleonEMOperator.hpp>
#include <TVector3.h>
#include <pair.h>

#include <string>
#include <cstdarg>

using namespace std;

class DoubleNModel{
public:
  DoubleNModel(MeanFieldNucleusThick *pnucleus, bool setpw, bool setSRC, bool setCT, bool setcorr, int particletype1, int particletype2, string dir);
  ~DoubleNModel();
  void setSRC(int setSRC) {SRC=setSRC;}
  complex<double> getMatrixEl(const TKinematics2to3 &tk,int spin1, int spin2, int photonpol, 
			      int shellindex1, int shellindex2, int m1, int m2);
  complex<double> MC_helper(const double r1, const double costheta1, const double phi1, const double r2, const double costheta2, const double phi2) const;  
  bool getPW() const {return pw;};
  bool getSRC() const {return SRC;};
  bool getCT() const {return CT;};
  bool getCorr() const {return corr;};
  complex<double> getJ1contrup() const {return J1contrup;};
  complex<double> getJ1contrdown() const {return J1contrdown;};
  complex<double> getJ2contrup() const {return J2contrup;};
  complex<double> getJ2contrdown() const {return J2contrdown;};
  Pair * getPair1up() const {return pair1up;};
  Pair * getPair1down() const {return pair1down;};
  Pair * getPair2up() const {return pair2up;};
  Pair * getPair2down() const {return pair2down;};
  GlauberGridThick * getGridf1() const {return gridf1;};
  GlauberGridThick * getGridf2() const {return gridf2;};
  TVector3 getPf1() const {return pf1;};
  TVector3 getPf2() const {return pf2;};
  TVector3 getQvec3() const {return qvec3;};
  
  
  
private:
  bool pw;
  bool SRC;
  bool CT;
  bool corr; 
  int particletype1;
  int particletype2;
  MeanFieldNucleusThick *pnucl;
  NucleonEMOperator *J1;
  NucleonEMOperator *J2;
//   Matrix<1,4> barcontract;
//   Matrix<1,4> barcontract0up;
//   Matrix<1,4> barcontractminup;
//   Matrix<1,4> barcontractplusup;
//   Matrix<1,4> barcontract0down;
//   Matrix<1,4> barcontractmindown;
//   Matrix<1,4> barcontractplusdown;
  TVector3 pf1;
  TVector3 pf2;
  TVector3 qvec3;
  complex<double> J1contrup;
  complex<double> J1contrdown;
  complex<double> J2contrup;
  complex<double> J2contrdown;
  
  Pair *pair1up;
  Pair *pair1down;
  Pair *pair2up;
  Pair *pair2down;
  
  GlauberGridThick *gridf1;
  GlauberGridThick *gridf2;
  string homedir;
  
  

};

double matrixel_MC_real (double* x, size_t dim, void * p);
double matrixel_MC_imag (double* x, size_t dim, void * p);
void display_results (char *title, double result, double error);


#endif