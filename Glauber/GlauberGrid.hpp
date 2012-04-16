/*! \file GlauberGrid.hpp 
 * \brief Contains declaration of class GlauberGrid
 * \author Wim Cosyn
 * \date 16/08/2011
 * 
 */
#ifndef GLAUBERGRID_H
#define GLAUBERGRID_H

#include <cstdarg>

#include "AbstractFsiCTGrid.hpp"

/*! \brief A class for a RMSGA ISI/FSI grid, implementing the abstract AbstractFsiCTGrid class
 */
class GlauberGrid : public AbstractFsiCTGrid{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnucl pointer to an instance of MeanFieldNucleus
   * \param dir string that contains dir with all input
   */
  GlauberGrid(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleus *pnucl, string dir);
  virtual ~GlauberGrid();/*!< Destructor */
  virtual complex<double> getFsiGridFull_interp();   /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiCtGridFull_interp(); /*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiGridN_interp(int grid); /*!returns the value of  fsi grid number grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  
  virtual void printFsi_grid();/*!< Prints the FSI grid for a certain situation*/
  virtual void printFsi_ct_grid(); /*!< Prints the FSI+CT grid for a certain situation */
  virtual void print_grid(int gridindex);
  
protected:
  complex<double> *****fsi_grid; /*!< grid that contains the fsi factor */
  complex<double> *****fsi_ct_grid; /*!< grid that contains the fsi+ct factor*/
  int ** treshold; /*!< array that checks if calculated glauberphases are close to one, then doesn't compute them for larger r, to save computing time*/
    
  
private:

/*!< set filenames of the grids, includes "Gl" prefix 
   *\param dir dir where all input/output is located */ 
  virtual void setFilenames(string dir);  

  
  virtual void constructAllGrids(); /*!< construct both the fsi and the fsi+ct grids */
  virtual void constructCtGrid();  /*!< construct only the fsi+ct grids */
  /*! calculates the glauberphases for one gridpoint (both FSI and FSI+CT)
   * \param i grid index in r
   * \param j grid index in costheta
   * \param k grid index in phi
   */
  virtual void calcGlauberphasesBoth(const int i, const int j, const int k); 
  /*! calculates the glauberphases for one gridpoint (only FSI+CT)
   * \param i grid index in r
   * \param j grid index in costheta
   * \param k grid index in phi
   */
  virtual void calcGlauberphasesCt(const int i, const int j, const int k);
  virtual void readinFsiGrid(ifstream &infile); /*!< read in both the fsi and the fsi+ct grids */
  virtual void readinFsiCtGrid(ifstream &infile); /*!< read in only the fsi+ct grids */
  virtual void writeoutFsiGrid(ofstream &outfile); /*!< write out both the fsi and the fsi+ct grids */
  virtual void writeoutFsiCtGrid(ofstream &outfile); /*!< write out only the fsi+ct grids */
  
  /*! function that gets integrated over r, both fsi and fsi+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberR(const double r, complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi and fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberCosTheta(const double costheta, complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, both fsi and fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhi(const double phi, complex<double> *results, va_list ap);
  /*! function that gets integrated over r, only fsi+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberRCT(const double r, complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), only fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberCosThetaCT(const double costheta, complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, only fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi+Ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhiCT(const double phi, complex<double> *results, va_list ap);
  

};

#endif