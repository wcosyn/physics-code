#ifndef RHOTCROSS_HPP
#define RHOTCROSS_HPP


#include <MeanFieldNucleusThick.hpp>
// #include <TKinematics2to2.h>
#include <DistMomDistrGrid.hpp>
#include <GlauberDecayGridThick.hpp>
#include <FastParticle.hpp>

#define NROFRES 9


class RhoTCross{
public:
  RhoTCross(const int nucleus, const double pmax, const string dir);
  ~RhoTCross();
  void getCrosst(double *results, const double Ebeam, const double Q2, const double nu, const double t);
  void getCrossz(double *results, const double Ebeam, const double Q2, const double nu, const double z);
private:
  string homedir;  
  double pmax;
  MeanFieldNucleusThick nucleusthick;
  DistMomDistrGrid *pdistgrid;
  GlauberDecayGridThick *pfsigrid;
//   TKinematics2to2 *pkin;
  FastParticle *prho;
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
  
  void getMomdistr(double *results, double prho, double thetarho, double t, int shell, 
			    double pm, double pmcostheta, double pmphi);
  double getfrontfactor(double nu, double qvec, double Erho, double prho, double pzrho, double pxrho,
				 double s, double Q2, double mN, double t);
};

#endif