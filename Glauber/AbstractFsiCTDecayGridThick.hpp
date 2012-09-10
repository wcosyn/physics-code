/*! \file AbstractFsiCTDecayGridThick.hpp 
 * \brief Contains declaration of abstract class AbstractFsiCTDecayGridThick
 * \author Wim Cosyn
 * \date 18/08/2011
 * \addtogroup Glauber
 * @{
 */

#ifndef ABSTRACTFSICTDECAYGRIDTHICK_H
#define ABSTRACTFSICTDECAYGRIDTHICK_H

#include <string>

#include "AbstractFsiDecayGridThick.hpp"
#include "AbstractFsiCTDecayGrid.hpp"
#include "MeanFieldNucleusThick.hpp"
#include "FsiCorrelator.hpp"

/*! \brief An abstract class for a ISI/FSI grid using thickness approximation, including a CT grid and Decay grids
 * 
 * Typically an object that is inherited from this class is operated as follows.<BR>
 * 1. Initialize object with constructor AbstractFsiCTDecayGridThick()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() or updateGrids() <BR>
 * 4. Add particles that are knocked out from nucleus <BR>
 * 5. Interpolate grid for a certain point or print the grid or whatever...<BR>
 * 
*/
class AbstractFsiCTDecayGridThick : public AbstractFsiDecayGridThick, public AbstractFsiCTDecayGrid{
  
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnuclthick pointer to an instance of MeanFieldNucleusThick
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param dir std::string that contains dir with all input
   */
  AbstractFsiCTDecayGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick *pnuclthick, 
		         double prec, int integrator, std::string dir);
  virtual ~AbstractFsiCTDecayGridThick();/*!< Destructor */
  
  virtual void printFsi_src_ct_grid()=0; /*!< Prints the FSI+SRC+CT grid for a certain situation, pure virtual function!! */
  virtual void printFsi_src_ct_decay_grid()=0; /*!< Prints the FSI+SRC+CT grid for a certain situation, pure virtual function!! */

  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,costheta,phi) */
  std::complex<double> getFsiSrcCtGridFull_interpvec(const TVector3 &rvec);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,costheta,phi) */
  std::complex<double> getFsiSrcCtGridFull_interp3(const double r, const double costheta, const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (costheta,phi), r has been set previously */
  std::complex<double> getFsiSrcCtGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (phi), r&theta have been set previously */
  std::complex<double> getFsiSrcCtGridFull_interp1(const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,costheta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual std::complex<double> getFsiSrcCtGridFull_interp()=0;

  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,costheta,phi) */
  std::complex<double> getFsiSrcCtDecayGridFull_interpvec(const TVector3 &rvec);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,costheta,phi) */
  std::complex<double> getFsiSrcCtDecayGridFull_interp3(const double r, const double costheta, const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (costheta,phi), r has been set previously */
  std::complex<double> getFsiSrcCtDecayGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (phi), r&theta have been set previously */
  std::complex<double> getFsiSrcCtDecayGridFull_interp1(const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,costheta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual std::complex<double> getFsiSrcCtDecayGridFull_interp()=0;

  protected:
   /*!< set filenames of the grids, includes "Thick" prefix
    * \param dir dir where all input/output is located */ 
  virtual void setFilenames(std::string dir); 

};
/** @} */
#endif  