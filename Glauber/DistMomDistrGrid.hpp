/*! \file DistMomDistrGrid.hpp 
 * \addtogroup Glauber
 * \brief Contains declaration of the class DistMomDistrGrid
 * \author Wim Cosyn
 * \date 16/04/2012
 * @{
 */

#ifndef DISTMOMDISTRGRID_H
#define DISTMOMDISTRGRID_H

#include <string>
#include <complex>
#include <TVector3.h>
#include <Matrix.h>

//forward declaration
class AbstractFsiCTGrid;

/*! \brief A class for a grid of a distorted momentum distribution
 * all momenta in MeV
 * all momentum distributions in fm^3!!!!!!!!!!!!!!!
 */
class DistMomDistrGrid{
public:
  /*! Constructor
   * \param shellindex shell index of the nucleus you want the mom distribution for. 
   * Starts from first proton shell, neutron shells follow
   * \param max_p maximum size of p (in MeV)
   * \param p_grid gridsize in p
   * \param cth_grid gridsize in cos theta
   * \param phi_grid gridsizein phi
   * \param pfsigrid pointer to a fsi/ct grid
   * \param dir string that contains dir with all input
   */
  DistMomDistrGrid(int shellindex, const double max_p, const int p_grid, const int cth_grid, const int phi_grid,
		   AbstractFsiCTGrid *pfsigrid, string dir);
  ~DistMomDistrGrid(); /*!< Destructor */
  
  int getShellindex() const{return shellindex;} /*!< shell level */
  const string getRho_Filename() const{return rho_filename;} /*!< returns filename for regular rho grid */
  const string getRhoCT_Filename() const{return rhoct_filename;} /*!< returns filename for rho+ct grid */
  const string getDir() const{return dir;}/*!< returns dir where all input/output is located */ 
  int getPgrid() const{return pgrid;}/*!<  returns the gridsize in r*/
  int getCthgrid() const{return cthgrid;}/*!<  returns the gridsize in theta*/
  int getPhigrid() const{return phigrid;}/*!<  returns the gridsize in phi*/
  double getMass() const{return mass;} /*!<  mass interacting nucleon*/
  double getPmax() const{return pmax;}/*!<  returns max p size of the grid*/
  double getInvPstep() const{return invpstep;}/*!<  returns the inverse of the stepsize of the grid in r*/
  double getInvCthstep() const{return invcthstep;}/*!<  returns the inverse of the stepsize of the grid in theta*/
  double getInvPhistep() const{return invphistep;}/*!<  returns the inverse of the stepsize of the grid in phi*/
  double getS_interp() const{return s_interp;}/*!<  returns s, measure where the point for interpolation is between two indices in r*/
  double getT_interp() const{return t_interp;}/*!<  returns t, measure where the point for interpolation is between two indices in theta*/
  double getU_interp() const{return u_interp;}/*!<  returns u, measure where the point for interpolation is between two indices in phi*/
  double getComp_s_interp() const{return comp_s_interp;}/*!<  returns 1-s*/
  double getComp_t_interp() const{return comp_t_interp;}/*!<  returns 1-t*/
  double getComp_u_interp() const{return comp_u_interp;}/*!<  returns 1-u*/
  int getPindex() const{return pindex;}/*!<  returns the index for r used in the 3d interpolation*/
  int getCthindex() const{return cthindex;}/*!<  returns the index for theta used in the 3d interpolation*/
  int getPhiindex() const{return phiindex;}/*!<  returns the index for phi used in the 3d interpolation */
  AbstractFsiCTGrid *getPfsigrid() const{return pfsigrid;}/*!<  returns the pointer to the MeanFieldNucleus instance */
  
  void printRho_grid(int gridindex);  /*!< Prints the rho grid for a certain situation*/
  void printRhopw_grid();  /*!< Prints the rho plane wave grid */

  void setPinterp(const double p); /*!< sets the rindex variable used in the 3d interpolation */
  void setCthinterp(const double costheta); /*!< sets the thindex variable used in the 3d interpolation */
  void setPhiinterp(const double phi); /*!< sets the phiindex variable used in the 3d interpolation */
  
  //interpolation functions, one without arguments is a pure virtual function
  /*!returns the value of the momentum distribution for a certain situation at coordinate (p,costheta,phi) */
  double getRhoGridFull_interpvec(int gridindex, const TVector3 &pvec); 
  /*!returns the value of the momentum distribution for a certain situation at coordinate (p,costheta,phi) */
  double getRhoGridFull_interp3(int gridindex, const double p, const double costheta, const double phi); 
  /*!returns the value of the momentum distribution for a certain situation at coordinate (costheta,phi), p has been set previously */
  double getRhoGridFull_interp2(int gridindex, const double costheta, const double phi);
  /*!returns the value of the momentum distribution for a certain situation at coordinate (phi), p&theta have been set previously */
  double getRhoGridFull_interp1(int gridindex, const double phi);
  /*!returns the value of the momentum distribution for a certain situation at coordinate (p,theta,phi) that has been set previously*/
  double getRhoGridFull_interp(int gridindex);
  double getRhopwGridFull_interp(double p);

   
  
   void fillGrids();  /*!< fills the momentum distribution grids that are used for interpolation */
   void updateGrids(AbstractFsiCTGrid *, int);  /*!< updates the momentum distribution grids that are used for interpolation */
  
  
  
private:
  int shellindex; /*!< shell level */
  string rho_filename; /*!< filename for regular rho grid */
  string rhoct_filename; /*!< filename for rho+ct grid */
  string dir;  /*!< dir where all input/output is located */ 
  double mass; /*!< mass interacting nucleon */ 
  double pmax; /*!< max size of p (MeV) */ 
  int level;  /*!< shell index level */ 
  int pgrid; /*!<  the gridsize in r*/
  int cthgrid; /*!<  the gridsize in theta*/
  int phigrid;/*!<  the gridsize in phi*/
  double invpstep; /*!<  the inverse of the stepsize of the grid in r*/
  double invcthstep;/*!<  the inverse of the stepsize of the grid in theta*/
  double invphistep;/*!<  the inverse of the stepsize of the grid in phi*/
  double s_interp; /*!<  s, measure where the point for interpolation is between two indices in r*/
  double t_interp; /*!<  t, measure where the point for interpolation is between two indices in theta*/
  double u_interp; /*!<  u, measure where the point for interpolation is between two indices in phi*/
  double comp_s_interp; /*!<  1-s*/
  double comp_t_interp; /*!<  1-t*/
  double comp_u_interp; /*!<  1-u*/
  int pindex; /*!<  the index for r used in the 3d interpolation*/
  int cthindex; ;/*!<  the index for theta used in the 3d interpolation*/
  int phiindex; /*!<  the index for phi used in the 3d interpolation */
  void setFilenames(string dir); /*!< set filenames of the grids \param dir dir where all input/output is located */ 
  bool filledgrid; /*!< denotes if the grid has been filled */
  bool filledctgrid; /*!< denotes if the ct grid has been filled */
  bool filledallgrid; /*!< denotes if allthe (possible) grids have been filled */
  TVector3 pvec_hit; /*!< coordinate where the mom distribution is calculated */
  double p_hit;  /*!< p coordinate of the momentum vector point, what the grid depens on */
  double costheta_hit;/*!< cos of theta coordinate of the momentum vector point, what the grid depens on */
  double sintheta_hit;/*!< sin of theta coordinate of the momentum vector point, what the grid depens on */
  double phi_hit;/*!< phi coordinate of the momentum vector point, what the grid depens on */
  double cosphi_hit;/*!< cos of phi coordinate of the momentum vector point, what the grid depens on */
  double sinphi_hit;/*!< sin of phi coordinate of the momentum vector point, what the grid depens on */
  Matrix<1,4> Upm_bar; /*!< holds a spinor that gest contracted with a nucleon MF spinor */
  
  double *rhopwgrid; /*!< plane-wave momentum distribution grid. only depends on norm of p */
  double ****rhogrid; /*!< non-ct momentum distribution grid */
  double ****rhoctgrid; /*!< ct momentum distribution grid */
    
  void constructpwGrid(); /*!< construct all grids */
  void constructAllGrids(); /*!< construct all grids */
  void readinRhoGrid(ifstream &infile); /*!< read in all grids */
  void writeoutRhoGrid(ofstream &outfile); /*!< write out all grids */
    
  void constructCtGrid();  /*!< construct only the fsi+ct grids */
  void readinRhoCtGrid(ifstream &infile); /*!< read in only the fsi+ct grids */
  void writeoutRhoCtGrid(ofstream &outfile); /*!< write out only the fsi+ct grids */

  AbstractFsiCTGrid *pfsigrid;/*!<  the pointer to the fsigrid instance */

  /*! function that gets integrated over r, both fsi and fsi+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intRhoR(const double r, complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi and fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intRhoCosTheta(const double costheta, complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, both fsi and fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intRhoPhi(const double phi, complex<double> *results, va_list ap);
  /*! function that gets integrated over r, only fsi+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intRhoRCT(const double r, complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), only fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intRhoCosThetaCT(const double costheta, complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, only fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi+Ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intRhoPhiCT(const double phi, complex<double> *results, va_list ap);
  
  /*! function that gets integrated over r, plane-wave thingy
   * \param r [fm] radial coordinate
   * \param result result: contains plane wave result
   * \param ap variable parameter list
   */
  void intRhoRpw(const double r, double *result, va_list ap);

  
  
  
};
/** @} */
#endif