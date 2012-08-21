#ifndef MODEL_HPP
#define MODEL_HPP

#include <FourVector.h>
#include <MeanFieldNucleusThick.hpp>
#include <AbstractFsiCTGrid.hpp>
#include <TKinematics2to2.h>
#include <NucleonEMOperator.hpp>
#include <TVector3.h>

#include <string>
#include <cstdarg>

using namespace std;

class Model{
public:
  Model(MeanFieldNucleusThick *pnucleus, int setSRC, int setthick, string dir);
  ~Model();
  void setSRC(int setSRC) {SRC=setSRC;}
  complex<double> getMatrixEl(TKinematics2to2 &tk,int spinout, int photonpol, int shellindex, int m, int CT, int pw);
  void getMatrixEl(TKinematics2to2 &tk, Matrix<2,3> &results, int shellindex, int m, int CT, int pw);
  complex<double> getFreeMatrixEl(TKinematics2to2 &tk, int spinin, int spinout, int photonpol);
  
private:
  int SRC;
  int thick;
  AbstractFsiCTGrid *grid;
  MeanFieldNucleusThick *pnucl;
  NucleonEMOperator *J;
  Matrix<1,4> barcontract;
  Matrix<1,4> barcontract0up;
  Matrix<1,4> barcontractminup;
  Matrix<1,4> barcontractplusup;
  Matrix<1,4> barcontract0down;
  Matrix<1,4> barcontractmindown;
  Matrix<1,4> barcontractplusdown;
  TVector3 pm;
  string homedir;
  int shell;
  int mm;
    /*! function that gets integrated over r, both fsi and fsi+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJR(const double r, complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi and fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJCosTheta(const double costheta, complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, both fsi and fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJPhi(const double phi, complex<double> *results, va_list ap);
    /*! function that gets integrated over r, both fsi and fsi+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJR12(const double r, complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi and fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJCosTheta12(const double costheta, complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, both fsi and fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intJPhi12(const double phi, complex<double> *results, va_list ap);
  
  

};

#endif