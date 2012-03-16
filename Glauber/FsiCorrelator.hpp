/*! \file FsiCorrelator.hpp 
 * \brief Contains declaration of class FsiCorrelator
 * \author Wim Cosyn
 * \date 16/08/2011
 * 
 */

#ifndef FSICORRELATOR_H
#define FSICORRELATOR_H

#include <string>

#include "MeanFieldNucleusThick.hpp"

/*! \brief A class for all things needed in correlated fsi/isi calculations
 * 
 * Contains a grid for the gamma function of a certain nucleus (normalization of conditional 2-b density). <BR>
 * Has the correlation function g(r). <BR>
 */
class FsiCorrelator{
public:
  /*! Constructor
   * \param inputnucleus pointer to an instance of MeanFieldNucleusThick
   * \param rgrid gridsize in r 
   * \param cthgrid gridsize in costheta
   * \param dir string that contains dir with all input
   */
  FsiCorrelator(MeanFieldNucleusThick *inputnucleus, const int rgrid, const int cthgrid,
    const string dir);
  ~FsiCorrelator();/*!< Destructor */
  
  const string getCorrFilename() const; /*!< returns filename for the gamma function grids */
  int getCorr_rgrid() const;/*!<  returns the gridsize in r*/
  int getCorr_cthgrid() const;/*!<  returns the gridsize in theta*/
  double getInvRstep() const;/*!<  returns the inverse of the stepsize of the grid in r*/
  double getInvCthstep() const;/*!<  returns the inverse of the stepsize of the grid in theta*/
  double getS_interp() const;/*!<  returns s, measure where the point for interpolation is between two indices in r*/
  double getT_interp() const;/*!<  returns t, measure where the point for interpolation is between two indices in theta*/
  double getComp_s_interp() const;/*!<  returns 1-s*/
  double getComp_t_interp() const;/*!<  returns 1-t*/
  int getRindex() const;/*!<  returns the index for r used in the 2d interpolation*/
  int getCthindex() const;/*!<  returns the index for theta used in the 2d interpolation*/
  void setRinterp(const double r);/*!< sets the rindex variable used in the 2d interpolation */
  void setCthinterp(double costheta); /*!< sets the thindex variable used in the 2d interpolation */
  double correlation(const double r) const; /*!< correlation function g(r) \param  r [fm] radial coordinate */
  
  /*!returns the value of the gamma grid for the total density at coordinate (r,theta) */
  double getCorrGridFull_interp(const double r, const double costheta);
  /*!returns the value of the gamma grid for the proton density at coordinate (r,theta) */
  double getCorrGridProton_interp(const double r, const double costheta);
  /*!returns the value of the gamma grid for the neutron density at coordinate (r,theta) */
  double getCorrGridNeutron_interp(const double r, const double costheta);
  /*!returns the value of the gamma grid for the total density at coordinate (theta), r has been set previously */
  double getCorrGridFull_interp(const double costheta);
  /*!returns the value of the gamma grid for the proton density at coordinate (theta), r has been set previously */
  double getCorrGridProton_interp(const double costheta);
  /*!returns the value of the gamma grid for the neutron density at coordinate (theta), r has been set previously */
  double getCorrGridNeutron_interp(const double costheta);
  /*!returns the value of the gamma grid for the total density at coordinate (r,theta) that has been set previously */
  double getCorrGridFull_interp() const;
  /*!returns the value of the gamma grid for the proton density at coordinate (r,theta) that has been set previously */
  double getCorrGridProton_interp() const;
  /*!returns the value of the gamma grid for the neutron density at coordinate (theta), r has been set previously */
  double getCorrGridNeutron_interp() const;
  /*!returns the value of the gamma grid for the proton or neutron density at coordinate (r,theta) 
   * \param r [fm] radial coordinate
   * \param costheta [rad] theta, spherical coordinates
   * \param proton selects proton (1) or neutron (0) density */
  double getCorrGrid_interp(const double r, const double costheta, bool proton);
  /*!returns the value of the gamma grid for the proton or neutron density at coordinate (theta), r has been set previously
   * \param theta [rad] theta, spherical coordinates
   * \param proton selects proton (1) or neutron (0) density */
  double getCorrGrid_interp(const double costheta, bool proton);
  /*!returns the value of the gamma grid for the proton or neutron density at coordinate (theta), r has been set previously
   \param proton selects proton (1) or neutron (0) density */
  double getCorrGrid_interp(bool proton) const;
  
  
  void printCorrGridFull() const; /*!< Print the gamma grid for the total density */
  void printCorrGridProton() const;/*!< Print the gamma grid for the proton density */
  void printCorrGridNeutron() const;/*!< Print the gamma grid for the neutron density */
  void printCorrGridAll() const;/*!< Print the gamma grid for the total,proton and neutron density */
  
private:
  string corrfilename; /*!< filename for the gamma grids */
  int corr_rgrid; /*!<  the gridsize in r*/
  int corr_cthgrid; /*!<  the gridsize in theta*/
  int corr_phigrid; /*!<  the gridsize in phi, not really used but somehow needed (it ain't broken so I keep it for now*/
  double invrstep; /*!<  the inverse of the stepsize of the grid in r*/
  double invcthstep;/*!<  the inverse of the stepsize of the grid in theta*/
  double s_interp; /*!<  s, measure where the point for interpolation is between two indices in r*/
  double t_interp; /*!<  t, measure where the point for interpolation is between two indices in theta*/
  double comp_s_interp; /*!<  1-s*/
  double comp_t_interp; /*!<  1-s*/
  int rindex; /*!<  the index for r used in the 3d interpolation*/
  int cthindex; /*!<  the index for theta used in the 3d interpolation*/
  double **corrgridfull; /*!< gamma grid for full density */
  double **corrgridproton; /*!< gamma grid for proton density*/
  double **corrgridneutron;  /*!< gamma grid for neutron density*/
  double **functionmatrix; /*!< grid that helps to construct corr grid*/
  double *x; /*!< grid that helps to construct corr grid*/
  
  MeanFieldNucleusThick *pnucleus; /*!<  the pointer to the MeanFieldNucleusThick instance */
  void setCorrFilename(const string dir); /*!< set the filename of the gamma grids */
  /*! Constructs the gamma grid for one of the densities 
   * \param density pointer to the density
   * \param corrgrid array with the gamma grid that will be filled
   */
  void constructCorrgrid(const double *density, double **corrgrid); 
  void constructMatrixSet(const double *density); /*!< construct matrix used to solve integration equation*/
  void newt(int *check);/*!< newton minimizationmethod for multi-dim */
  void ConstructJacobian(double **Jacobian);/*!< construct Jacobian for solving integration equation */
  double ReturnFunction(double *functionvalues);/*!< function that needs minimizing */
  double ReturnFunction();/*!< function that needs  minimizing */
  void corrmaintenance(int n, double *functionvalues, double **Jacobian, double *grad, double *xold, double *step); /*!< mem maintenance */
  void lnsrch(double *xold, double functionold, double *grad, double *step, double *functionnew,
	    double stpmax, int *check, double *functionvalues);/*!< linesearch minimalizer */
};

#endif
