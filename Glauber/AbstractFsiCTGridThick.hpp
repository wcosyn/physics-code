/*! \file AbstractFsiCTGridThick.hpp 
 * \brief Contains declaration of abstract class AbstractFsiCTGridThick
 * \author Wim Cosyn
 * \date 18/08/2011
 * \addtogroup Glauber
 * @{
 */

#ifndef ABSTRACTFSICTGRIDTHICK_H
#define ABSTRACTFSICTGRIDTHICK_H

#include <string>

#include "AbstractFsiGridThick.hpp"
#include "AbstractFsiCTGrid.hpp"
#include "MeanFieldNucleusThick.hpp"
#include "FsiCorrelator.hpp"

/*! \brief An abstract class for a ISI/FSI grid using thickness approximation, including a CT grid and SRC grids
 * 
 * Typically an object that is inherited from this class is operated as follows.<BR>
 * 1. Initialize object with constructor AbstractFsiCTGridThick()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() or updateGrids() <BR>
 * 4. Add particles that are knocked out from nucleus <BR>
 * 5. Interpolate grid for a certain point or print the grid or whatever...<BR>
 * 

 */
class AbstractFsiCTGridThick : public AbstractFsiGridThick, public AbstractFsiCTGrid{
  
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnuclthick pointer to an instance of MeanFieldNucleusThick
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   */
  AbstractFsiCTGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick *pnuclthick, 
		         string dir);
  virtual ~AbstractFsiCTGridThick();/*!< Destructor */
  
  virtual void printFsi_src_ct_grid()=0; /*!< Prints the FSI+SRC+CT grid for a certain situation, pure virtual function!! */

  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiSrcCtGridFull_interpvec(const TVector3 &rvec);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiSrcCtGridFull_interp3(const double r, const double costheta, const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (costheta,phi), r has been set previously */
  complex<double> getFsiSrcCtGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiSrcCtGridFull_interp1(const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,costheta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiSrcCtGridFull_interp()=0;

  protected:
   /*!< set filenames of the grids, includes "Thick" prefix
    * \param dir dir where all input/output is located */ 
  virtual void setFilenames(string dir); 

};
/** @} */
#endif  