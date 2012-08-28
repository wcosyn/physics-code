/*! \file AbstractFsiDecayGridThick.hpp 
 * \brief Contains declaration of abstract class AbstractFsiDecayGridThick
 * \author Wim Cosyn
 * \date 16/08/2011
 * \addtogroup Glauber
 * @{
 */

#ifndef ABSTRACTFSIDECAYGRIDTHICK_H
#define ABSTRACTFSIDECAYGRIDTHICK_H

#include <string>

#include "AbstractFsiGridThick.hpp"
#include "MeanFieldNucleusThick.hpp"
#include "FsiCorrelator.hpp"

/*! \brief An abstract class for a ISI/FSI grid using thickness approximation and including decay properties 
 * 
 * Contains all necessary functions to operate a general fsi grid. Has a pointer to a nucleus for which the fsi grid is.  
 * A vector contains all the particles that undergo ISI or FSI.  Has interpolation and print functions for the grid.
 * 
 * Typically an object that is inherited from this class is operated as follows.<BR>
 * 1. Initialize object with constructor AbstractFsiDecayGridThick()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() or updateGrids() <BR>
 * 4. Add particles that are knocked out from nucleus <BR>
 * 5. Interpolate grid for a certain point or print the grid or whatever...<BR>
 * 
 * 
 * 
 */
class AbstractFsiDecayGridThick : virtual public AbstractFsiGridThick{
  
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnuclthick pointer to an instance of MeanFieldNucleusThick
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   */
  AbstractFsiDecayGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick *pnuclthick, 
		       string dir);
  virtual ~AbstractFsiDecayGridThick();/*!< Destructor */
  


  virtual void printFsi_decay_grid()=0; /*!< Prints the FSI+decay grid for a certain situation, pure virtual function!! */
  virtual void printFsi_src_decay_grid()=0; /*!< Prints the FSI+SRC+decay grid for a certain situation, pure virtual function!! */

   /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiSrcDecayGridFull_interpvec(const TVector3 & rvec);
 /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiSrcDecayGridFull_interp3(const double r, const double costheta, const double phi);
  /*!returns the value of the fsi+src grid for a certain situation at coordinate (theta,phi), r has been set previously */
  complex<double> getFsiSrcDecayGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi+src grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiSrcDecayGridFull_interp1(const double phi);
  /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,costheta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiSrcDecayGridFull_interp()=0;
  
   /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiDecayGridFull_interpvec(const TVector3 & rvec);
 /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiDecayGridFull_interp3(const double r, const double costheta, const double phi);
  /*!returns the value of the fsi+src grid for a certain situation at coordinate (theta,phi), r has been set previously */
  complex<double> getFsiDecayGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi+src grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiDecayGridFull_interp1(const double phi);
  /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,costheta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiDecayGridFull_interp()=0;

protected:
   /*!< set filenames of the grids, includes "Thick" & "Decay" prefix
    * \param dir dir where all input/output is located */ 
  virtual void setFilenames(string dir); 
  
};
/** @} */
#endif