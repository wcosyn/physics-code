/*! \file AbstractFsiCTDecayGrid.hpp 
 * \brief Contains declaration of abstract class AbstractFsiCTDecayGrid
 * \author Wim Cosyn
 * \date 18/08/2011
 * \addtogroup Glauber
 * @{
 */

#ifndef ABSTRACTFSICTDECAYGRID_H
#define ABSTRACTFSICTDECAYGRID_H

#include <string>
#include <complex>
#include <vector>

#include "MeanFieldNucleus.hpp"
#include "AbstractFsiCTGrid.hpp"

/*! \brief An abstract class for a ISI/FSI with also a grid for CT ISI/FSI and/or Decay properties
 * 
 * Typically an object that is inherited from this class is operated as follows.<BR>
 * 1. Initialize object with constructor AbstractFsiCTDecayGrid()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() or updateGrids() <BR>
 * 4. Add particles that are knocked out from nucleus <BR>
 * 5. Interpolate grid for a certain point or print the grid or whatever...<BR>
 * 
 * 
 */
class AbstractFsiCTDecayGrid : public virtual AbstractFsiCTGrid{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnucl pointer to an instance of MeanFieldNucleus
   * \param prec precision you want in the integrations
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   */
  AbstractFsiCTDecayGrid(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleus *pnucl, 
			 double prec, string dir);
  virtual ~AbstractFsiCTDecayGrid(); /*!< Destructor */

  virtual void printFsi_decay_grid()=0; /*!< Prints the FSI+decay grid for a certain situation, pure virtual function!! */
  virtual void printFsi_ct_decay_grid()=0; /*!< Prints the FSI+ct+decay grid for a certain situation, pure virtual function!! */

  /*!returns the value of the fsi+ct+decay grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiCtDecayGridFull_interpvec(const TVector3 &rvec);
  /*!returns the value of the fsi+ct+decay grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiCtDecayGridFull_interp3(const double r, const double costheta, const double phi);
  /*!returns the value of the fsi+ct+decay grid for a certain situation at coordinate (costheta,phi), r has been set previously */
  complex<double> getFsiCtDecayGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi+ct+decay grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiCtDecayGridFull_interp1(const double phi);
  /*!returns the value of the fsi+ct+decay grid for a certain situation at coordinate (r,theta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiCtDecayGridFull_interp()=0;
  
  /*!returns the value of the fsi+decay grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiDecayGridFull_interpvec(const TVector3 &rvec);
  /*!returns the value of the fsi+decay grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiDecayGridFull_interp3(const double r, const double costheta, const double phi);
  /*!returns the value of the fsi+decay grid for a certain situation at coordinate (costheta,phi), r has been set previously */
  complex<double> getFsiDecayGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi+decay grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiDecayGridFull_interp1(const double phi);
  /*!returns the value of the fsi+decay grid for a certain situation at coordinate (r,theta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiDecayGridFull_interp()=0;

protected:
  virtual void setFilenames(string dir); /*!< set filenames of the grids \param dir dir where all input/output is located */ 


  
};
/** @} */
#endif
