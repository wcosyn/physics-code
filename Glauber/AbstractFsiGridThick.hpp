/*! \file AbstractFsiGridThick.hpp 
 * \brief Contains declaration of abstract class AbstractFsiGridThick
 * \author Wim Cosyn
 * \date 16/08/2011
 * \addtogroup Glauber
 * @{
 */

#ifndef ABSTRACTFSIGRIDTHICK_H
#define ABSTRACTFSIGRIDTHICK_H

#include <string>

#include "AbstractFsiGrid.hpp"
#include "MeanFieldNucleusThick.hpp"
#include "FsiCorrelator.hpp"

/*! \brief An abstract class for a ISI/FSI grid using thickness approximation
 * 
 * Contains all necessary functions to operate a general fsi grid. Has a pointer to a nucleus for which the fsi grid is.  
 * A vector contains all the particles that undergo ISI or FSI.  Has interpolation and print functions for the grid.
 * 
 * Typically an object that is inherited from this class is operated as follows.<BR>
 * 1. Initialize object with constructor AbstractFsiGridThick()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() <BR>
 * 4. Interpolate grid for a certain point or print the grid or whatever...<BR>
 * 
 * 
 */
class AbstractFsiGridThick : virtual public AbstractFsiGrid{
  
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
  AbstractFsiGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick *pnuclthick, 
		       double prec, int integrator, string dir);
  virtual ~AbstractFsiGridThick();/*!< Destructor */
  
  const FsiCorrelator& getFsiCorrelator() const; /*!<  returns the pointer to the FsiCorrelator instance */
  FsiCorrelator& getFsiCorrelator(); /*!<  returns the pointer to the FsiCorrelator instance */

  MeanFieldNucleusThick *getPnucleusthick() const; /*!<  returns the pointer to the MeanFieldNucleusThick instance */

  virtual void printFsi_src_grid()=0; /*!< Prints the FSI+SRC grid for a certain situation, pure virtual function!! */

   /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiSrcGridFull_interpvec(const TVector3 & rvec);
 /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiSrcGridFull_interp3(const double r, const double costheta, const double phi);
  /*!returns the value of the fsi+src grid for a certain situation at coordinate (theta,phi), r has been set previously */
  complex<double> getFsiSrcGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi+src grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiSrcGridFull_interp1(const double phi);
  /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,costheta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiSrcGridFull_interp()=0;
  
protected:
   /*!< set filenames of the grids, includes "Thick" prefix
    * \param dir dir where all input/output is located */ 
  virtual void setFilenames(string dir); 
  
private:
  MeanFieldNucleusThick *pnucleusthick;/*!<  the pointer to the MeanFieldNucleusThick instance */
  FsiCorrelator fsicorrelator;/*!<  the pointer to the FsiCorrelator instance */
};
/** @} */
#endif