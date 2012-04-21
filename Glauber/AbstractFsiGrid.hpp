/*! \file AbstractFsiGrid.hpp 
 * \brief Contains declaration of abstract class AbstractFsiGrid
 * \author Wim Cosyn
 * \date 16/08/2011
 * 
 */

#ifndef ABSTRACTFSIGRID_H
#define ABSTRACTFSIGRID_H

#include <string>
#include <complex>
#include <vector>
#include <TVector3.h>

#include "FastParticle.hpp"
#include "MeanFieldNucleus.hpp"

/*! \brief An abstract class for a ISI/FSI grid
 * 
 * Contains all necessary functions to operate a general fsi grid. Has a pointer to a nucleus for which the fsi grid is.  
 * A vector contains all the particles that undergo ISI or FSI.  Has interpolation and print functions for the grid.
 * 
 * Typically an object that is inherited from this class is operated as follows.<BR>
 * 1. Initialize object with constructor AbstractFsiGrid()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() <BR>
 * 4. Add particles that are knocked out from nucleus <BR>
 * 5. Interpolate grid for a certain point or print the grid or whatever...
 * 
 * 
 */
class AbstractFsiGrid{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in cos theta
   * \param phi_grid gridsizein phi
   * \param pnucl pointer to an instance of MeanFieldNucleus
   * \param dir string that contains dir with all input
   */
  AbstractFsiGrid(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleus *pnucl, string dir);
  virtual ~AbstractFsiGrid(); /*!< Destructor */
  
  const string getFsi_Filename() const; /*!< returns filename for regular fsi grid */
  const string getDir() const;/*!< returns dir where all input/output is located */ 
  int getRgrid() const;/*!<  returns the gridsize in r*/
  int getCthgrid() const;/*!<  returns the gridsize in theta*/
  int getPhigrid() const;/*!<  returns the gridsize in phi*/
  double getInvRstep() const;/*!<  returns the inverse of the stepsize of the grid in r*/
  double getInvCthstep() const;/*!<  returns the inverse of the stepsize of the grid in theta*/
  double getInvPhistep() const;/*!<  returns the inverse of the stepsize of the grid in phi*/
  double getS_interp() const;/*!<  returns s, measure where the point for interpolation is between two indices in r*/
  double getT_interp() const;/*!<  returns t, measure where the point for interpolation is between two indices in theta*/
  double getU_interp() const;/*!<  returns u, measure where the point for interpolation is between two indices in phi*/
  double getComp_s_interp() const;/*!<  returns 1-s*/
  double getComp_t_interp() const;/*!<  returns 1-t*/
  double getComp_u_interp() const;/*!<  returns 1-u*/
  int getRindex() const;/*!<  returns the index for r used in the 3d interpolation*/
  int getCthindex() const;/*!<  returns the index for theta used in the 3d interpolation*/
  int getPhiindex() const;/*!<  returns the index for phi used in the 3d interpolation */
  MeanFieldNucleus *getPnucleus() const;/*!<  returns the pointer to the MeanFieldNucleus instance */
  vector<FastParticle> & getParticles() {return particles;}/*!<  returns the vector with all particles subject to isi/fsi */
  vector<int> getKnockoutLevels() const;/*!< returns the vector with all the shell index values of knocked out nucleons */
  vector<int> getKnockoutM() const;/*!<  returns the vector with all the m values of knocked out nucleons */
  bool getAllinplane() const;/*!<  returns 1 if all the isi/fsi particles lie in the x-z plane */
  /*! add particle subject to isi/fsi, checks to see that the vector contains max one incoming particle
   *\param newparticle pointer to instance of FastParticle you want to add*/
  virtual void addParticle(FastParticle &newparticle); 
  void clearParticles(); /*!< Clear contents of isi/fsi particles vector*/
  void printParticles() const; /*!< Print contents of isi/fsi particles vector*/
  /*! add quantum numbers of particle that is knocked out from nucleus
   *\param level shellindex level
   *\param m mj quantum number (times two!!)  [-2j,2j]   
   */
  void addKnockout(const int level, const int m); 
  void clearKnockout(); /*!< clear the vectors of the knocked out particles*/
  void printKnockout() const; /*!< print info about the particles that are knocked out from the nucleus*/

  int getTotalProtonOut() const; /*!< get total protons subject to fsi/isi (isi counts negative)*/
  int getTotalNeutronOut() const; /*!< get total neutrons subject to fsi/isi (isi counts negative)*/
  int getProtonKnockout() const; /*!< get total protons that are knocked out from nucleus*/
  int getNeutronKnockout() const; /*!< get total neutrons that are knocked out from nucleus*/

  int getNumber_of_grids() const{return number_of_grids;}

  virtual void printFsi_grid()=0;  /*!< Prints the FSI grid for a certain situation, pure virtual function!! */
  virtual void print_grid(int gridindex)=0;  /*!< Prints the FSI grid for a certain situation, pure virtual function!! */

  void setRinterp(const double r); /*!< sets the rindex variable used in the 3d interpolation */
  void setCthinterp(const double costheta); /*!< sets the thindex variable used in the 3d interpolation */
  void setPhiinterp(const double phi); /*!< sets the phiindex variable used in the 3d interpolation */
  
  //interpolation functions, one without arguments is a pure virtual function
  /*!returns the value of the fsi grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiGridFull_interpvec(const TVector3 &rvec); 
  /*!returns the value of the fsi grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiGridFull_interp3(const double r, const double costheta, const double phi); 
  /*!returns the value of the fsi grid for a certain situation at coordinate (costheta,phi), r has been set previously */
  complex<double> getFsiGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiGridFull_interp1(const double phi);
  /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiGridFull_interp()=0;

  //interpolation functions, one without arguments is a pure virtual function
  /*!returns the value of the fsi grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiGridN_interpvec(const int grid, const TVector3 &rvec); 
  /*!returns the value of the fsi grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiGridN_interp3(const int grid, const double r, const double costheta, const double phi); 
  /*!returns the value of the fsi grid for a certain situation at coordinate (costheta,phi), r has been set previously */
  complex<double> getFsiGridN_interp2(const int grid, const double costheta, const double phi);
  /*!returns the value of the fsi grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiGridN_interp1(const int grid, const double phi);
  /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiGridN_interp(const int grid)=0;
  /*! 3d interpolation of a grid after a certain coordinate (r,costheta,phi) has been set! */
  virtual complex<double> getInterp(complex<double> ***grid);  //3d interpolation of grid
  
  
  virtual void fillGrids();  /*!< fills the fsi grids that are used for interpolation */
  
  
protected:
  virtual void setFilenames(string dir); /*!< set filenames of the grids \param dir dir where all input/output is located */ 
  bool filledgrid; /*!< denotes if the grid has been filled */
  bool filledallgrid; /*!< denotes if allthe (possible) grids have been filled */
  double r_hit;  /*!< r coordinate of the hard interaction point, what the grid depens on */
  //double theta_hit;/*!< theta coordinate of the hard interaction point, what the grid depens on */
  double costheta_hit;/*!< cos of theta coordinate of the hard interaction point, what the grid depens on */
  double sintheta_hit;/*!< sin of theta coordinate of the hard interaction point, what the grid depens on */
  double phi_hit;/*!< phi coordinate of the hard interaction point, what the grid depens on */
  double cosphi_hit;/*!< cos of phi coordinate of the hard interaction point, what the grid depens on */
  double sinphi_hit;/*!< sin of phi coordinate of the hard interaction point, what the grid depens on */
  string fsi_filename; /*!< filename for regular fsi grid */
  int number_of_grids; /*!<number of fsi grids */
  
private:
  string dir;  /*!< dir where all input/output is located */ 
  
  int rgrid; /*!<  the gridsize in r*/
  int cthgrid; /*!<  the gridsize in theta*/
  int phigrid;/*!<  the gridsize in phi*/
  double invrstep; /*!<  the inverse of the stepsize of the grid in r*/
  double invcthstep;/*!<  the inverse of the stepsize of the grid in theta*/
  double invphistep;/*!<  the inverse of the stepsize of the grid in phi*/
  double s_interp; /*!<  s, measure where the point for interpolation is between two indices in r*/
  double t_interp; /*!<  t, measure where the point for interpolation is between two indices in theta*/
  double u_interp; /*!<  u, measure where the point for interpolation is between two indices in phi*/
  double comp_s_interp; /*!<  1-s*/
  double comp_t_interp; /*!<  1-t*/
  double comp_u_interp; /*!<  1-u*/
  int rindex; /*!<  the index for r used in the 3d interpolation*/
  int cthindex; ;/*!<  the index for theta used in the 3d interpolation*/
  int phiindex; /*!<  the index for phi used in the 3d interpolation */
  bool allinplane; /*!<  1 if all the isi/fsi particles lie in the x-z plane */
  
  int totalprotonout; /*!<  total protons subject to fsi/isi (isi counts negative)*/
  int totalneutronout;  /*!<  total neutrons subject to fsi/isi (isi counts negative)*/
  int protonknockout;  /*!<  total protons that are knocked out from nucleus*/
  int neutronknockout; /*!<  total neutrons that are knocked out from nucleus*/
    
  vector<FastParticle> particles; /*!<  the vector with all particles subject to isi/fsi */
  vector<int> knockoutlevels; /*!< the vector with all the shell index values of knocked out nucleons */
  vector<int> knockoutm; /*!<  the vector with all the m values of knocked out nucleons */
  
  //pure virtual functions!!!
  virtual void constructAllGrids()=0; /*!< construct all grids */
  virtual void readinFsiGrid(ifstream &infile)=0; /*!< read in all grids */
  virtual void writeoutFsiGrid(ofstream &outfile)=0; /*!< write out all grids */
    
  MeanFieldNucleus *pnucleus;/*!<  the pointer to the MeanFieldNucleus instance */
  
  
};

#endif