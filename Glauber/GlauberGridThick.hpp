/*! \file GlauberGridThick.hpp 
 * \brief Contains declaration of class GlauberGridThick
 * \author Wim Cosyn
 * \date 16/08/2011
 * 
 */

#ifndef GLAUBERGRIDTHICK_H
#define GLAUBERGRIDTHICK_H

#include <cstdarg>

#include "AbstractFsiCTGridThick.hpp"

/*! \brief An class for a Glauber ISI/FSI grid using thickness approximation,implements AbstractFsiCTGridThick
 * 
 * Four grids are actually calculated.  fsi_grid contains both the regular fsi grid as well as the grid that includes SRC.  
 * Same thing for the fsi_ct_grid
 */
class GlauberGridThick : public AbstractFsiCTGridThick{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnuclthick pointer to an instance of MeanFieldNucleusThick
   * \param dir string that contains dir with all input
   */
  GlauberGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick *pnuclthick, 
		   string dir);/*!< Destructor */
  virtual ~GlauberGridThick();
  virtual complex<double> getFsiGridFull_interp(); /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiCtGridFull_interp();/*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiSrcGridFull_interp(); /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiSrcCtGridFull_interp();/*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/

  virtual void printFsi_grid();/*!< Prints the FSI grid for a certain situation*/
  virtual void printFsi_ct_grid(); /*!< Prints the FSI+CT grid for a certain situation */
  virtual void printFsi_src_grid();/*!< Prints the FSI+SRC grid for a certain situation*/
  virtual void printFsi_src_ct_grid();/*!< Prints the FSI+SRC+CT grid for a certain situation*/

  
protected:
  complex<double> *****fsi_grid; /*!< grid that contains the fsi factor, also has the fsi+src grid */
  complex<double> *****fsi_ct_grid; /*!< grid that contains the fsi+ct factor, also has the fsi+src+ct grid */
  int ** treshold; /*!< array that checks if calculated glauberphases are close to one, then doesn't compute them for larger r, to save computing time*/
  
  
private:

  
  virtual void setFilenames(string dir); /*!< set filenames of the grids \param dir dir where all input/output is located */ 

  
  virtual void constructAllGrids(); /*!< construct both the fsi(+src) and the fsi(+src)+ct grids */
  virtual void constructCtGrid(); /*!< construct only the fsi(+src)+ct grids */
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
  virtual void readinFsiGrid(ifstream &infile);  /*!< read in both the fsi(+src) and the fsi(+src)+ct grids */
  virtual void readinFsiCtGrid(ifstream &infile); /*!< read in only the fsi(+src)+ct grids */
  virtual void writeoutFsiGrid(ofstream &outfile); /*!< write out both the fsi(+src) and the fsi(+src)+ct grids */
  virtual void writeoutFsiCtGrid(ofstream &outfile); /*!< write out only the fsi(+src)+ct grids */
  
  /*! function that gets integrated over r, both fsi(+src) and fsi(+src)+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberR(const double r, complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi(+src) and fsi(+src)+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberCosTheta(const double costheta, complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, both fsi(+src) and fsi(+src)+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhi(const double phi, complex<double> *results, va_list ap);
  /*! function that gets integrated over r, only fsi(+src)+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberRCT(const double r, complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), only fsi(+src)+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberCosThetaCT(const double costheta, complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, only fsi(+src)+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhiCT(const double phi, complex<double> *results, va_list ap);
  

};

#endif