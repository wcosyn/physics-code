#ifndef RHODEUTERON_HPP
#define RHODEUTERON_HPP


#include <MeanFieldNucleusThick.hpp>
// #include <TKinematics2to2.h>
#include <DistMomDistrGrid.hpp>
#include <GlauberDecayGridThick.hpp>
#include <FastParticle.hpp>
#include <TDeuteron.h>
#include "DeuteronMomDistr.hpp"

#define NROFRES 9


class RhoDeuteron{
public:
  RhoDeuteron(const string wavefunction, const double p_max, const string dir);
  ~RhoDeuteron();
  void getCrosst(double *results, const double Ebeam, const double Q2, const double nu, const double t);
  void getCrossz(double *results, const double Ebeam, const double Q2, const double nu, const double z);
  void getCrosst_coh(double *results, const double Ebeam, const double Q2, const double nu, const double t);
  void getCrossz_coh(double *results, const double Ebeam, const double Q2, const double nu, const double z);
private:
  string homedir;  
  double pmax;
  DeuteronMomDistr deuteron;
  void totdens_qt(const double qt, complex<double>* result, va_list ap);
  void totdens_qphi(const double qphi, complex<double>* result, va_list ap);
  double calcCross_coh(TVector3 &pDvec);
  //   TKinematics2to2 *pkin;
  /*! function that gets integrated over pm, all different fsi outputs
   * \param pm [MeV] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intPmt(const double pm, double *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi and fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intCosThetat(const double costheta, double *results, va_list ap);
  /*! function that gets integrated over phi, both fsi and fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intPhit(const double phi, double *results, va_list ap);
  /*! function that gets integrated over pm, all different fsi outputs
   * \param pm [MeV] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intPmz(const double pm, double *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi and fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intCosThetaz(const double costheta, double *results, va_list ap);
  /*! function that gets integrated over phi, both fsi and fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intPhiz(const double phi, double *results, va_list ap);
  
  /*! function that gets integrated over pm, all different fsi outputs
   * \param pm [MeV] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intPmcoh(const double pm, complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi and fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intCosThetacoh(const double costheta, complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, both fsi and fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intPhicoh(const double phi, complex<double> *results, va_list ap);
  void getMomdistr(double *results, double Erho, double prho, double pzrho, 
			      double pxrho, TVector3 &pmvec, double nu, double qvec);
  double getfrontfactor(double qvec, double Erho, double prho, double pzrho, double s, double Q2, double t, double En);

};

#endif