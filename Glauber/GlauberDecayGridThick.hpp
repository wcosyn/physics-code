/*! \file GlauberDecayGridThick.hpp 
 * \brief Contains declaration of class GlauberDecayGridThick
 * \author Wim Cosyn
 * \date 16/08/2011
 * \addtogroup Glauber
 * @{
 */

#ifndef GLAUBERDECAYGRIDTHICK_H
#define GLAUBERDECAYGRIDTHICK_H

#include <cstdarg>

#include "AbstractFsiCTDecayGridThick.hpp"

/*! \brief A class for a Glauber ISI/FSI grid using thickness approximation and including decay of fast particles,implements AbstractFsiCTDecayGridThick
 * 
 * Eight grids are actually calculated.  fsi_grid contains both the regular fsi grid as well as the grid that includes SRC.
 * Plus a one that includes Decay properties as well for each of those two.  
 * Same thing for the fsi_ct_grid. <BR>
 * 
 * Every function that takes a \param grid argument: <BR>
 * 0: RMSGA <BR>
 * 1: RMSGA+SRC <BR>
 * 2: RMSGA+Decay <BR>
 * 3: RMSGA+SRC+Decay <BR>
 * 4: RMSGA+CT <BR>
 * 5: RMSGA++CT+SRC <BR>
 * 6: RMSGA+CT+Decay <BR>
 * 7: RMSGA+CT+SRC+Decay aka "The works" <BR> 
 * 
 * 
 * Typically an object that is an instance from this class is operated as follows.<BR>
 * 1. Initialize object with constructor GlauberDecayGridThick()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() or updateGrids() <BR>
 * 4. Add particles that are knocked out from nucleus <BR>
 * 5. Interpolate grid for a certain point or print the grid or whatever...<BR>
 * 
 */
class GlauberDecayGridThick : public AbstractFsiCTDecayGridThick{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnuclthick pointer to an instance of MeanFieldNucleusThick
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   */
  GlauberDecayGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick *pnuclthick, 
		   double prec, int integrator, string dir);
  /*!< Destructor */
  virtual ~GlauberDecayGridThick();
  virtual complex<double> getFsiGridFull_interp(); /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiCtGridFull_interp();/*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiSrcGridFull_interp(); /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiSrcCtGridFull_interp();/*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiDecayGridFull_interp(); /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiCtDecayGridFull_interp();/*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiSrcDecayGridFull_interp(); /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual complex<double> getFsiSrcCtDecayGridFull_interp();/*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously
   * \param grid 0: RMSGA <BR>
  * 1: RMSGA+SRC <BR>
  * 2: RMSGA+DECAY <BR>
  * 3: RMSGA+SRC+DECAY <BR>
  * 4: RMSGA+CT <BR>
  * 5: RMSGA++CT+SRC <BR>
  * 6: RMSGA+CT+DECAY <BR>
  * 7: RMSGA+CT+SRC+DECAY AKA "THE WORKS" <BR> */
  virtual complex<double> getFsiGridN_interp(int grid); 

  virtual void printFsi_grid();/*!< Prints the FSI grid for a certain situation*/
  virtual void printFsi_ct_grid(); /*!< Prints the FSI+CT grid for a certain situation */
  virtual void printFsi_src_grid();/*!< Prints the FSI+SRC grid for a certain situation*/
  virtual void printFsi_src_ct_grid();/*!< Prints the FSI+SRC+CT grid for a certain situation*/
  virtual void printFsi_decay_grid();/*!< Prints the FSI grid for a certain situation*/
  virtual void printFsi_ct_decay_grid(); /*!< Prints the FSI+CT grid for a certain situation */
  virtual void printFsi_src_decay_grid();/*!< Prints the FSI+SRC grid for a certain situation*/
  virtual void printFsi_src_ct_decay_grid();/*!< Prints the FSI+SRC+CT grid for a certain situation*/
  /*! Prints the FSI grid for a certain situation
  * \param grid 0: RMSGA <BR>
  * 1: RMSGA+SRC <BR>
  * 2: RMSGA+DECAY <BR>
  * 3: RMSGA+SRC+DECAY <BR>
  * 4: RMSGA+CT <BR>
  * 5: RMSGA++CT+SRC <BR>
  * 6: RMSGA+CT+DECAY <BR>
  * 7: RMSGA+CT+SRC+DECAY AKA "THE WORKS" <BR> */
  virtual void print_grid(int grid);

  
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
/** @} */
#endif