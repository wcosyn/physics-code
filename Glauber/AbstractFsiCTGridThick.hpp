/*! \file AbstractFsiCTGridThick.hpp 
 * \brief Contains declaration of abstract class AbstractFsiCTGridThick
 * \author Wim Cosyn
 * \date 18/08/2011
 * 
 */

#ifndef ABSTRACTFSICTGRIDTHICK_H
#define ABSTRACTFSICTGRIDTHICK_H

#include <string>

#include "AbstractFsiGridThick.hpp"
#include "AbstractFsiCTGrid.hpp"
#include "MeanFieldNucleusThick.hpp"
#include "FsiCorrelator.hpp"

/*! \brief An abstract class for a ISI/FSI grid using thickness approximation, including a CT grid
 */
class AbstractFsiCTGridThick : public AbstractFsiGridThick, public AbstractFsiCTGrid{
  
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param th_grid gridsize in theta
   * \param phi_grid gridsizein phi
   * \param pnuclthick pointer to an instance of MeanFieldNucleusThick
   * \param dir string that contains dir with all input
   */
  AbstractFsiCTGridThick(const int r_grid, const int th_grid, const int phi_grid, MeanFieldNucleusThick *pnuclthick, 
		         string dir);
  virtual ~AbstractFsiCTGridThick();/*!< Destructor */
  
  virtual void printFsi_src_ct_grid()=0; /*!< Prints the FSI+SRC+CT grid for a certain situation, pure virtual function!! */

  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,theta,phi) */
  complex<double> getFsiSrcCtGridFull_interpvec(const TVector3 &rvec);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,theta,phi) */
  complex<double> getFsiSrcCtGridFull_interp3(const double r, const double theta, const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (theta,phi), r has been set previously */
  complex<double> getFsiSrcCtGridFull_interp2(const double theta, const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiSrcCtGridFull_interp1(const double phi);
  /*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiSrcCtGridFull_interp()=0;

  protected:
   /*!< set filenames of the grids, includes "Thick" prefix
    * \param dir dir where all input/output is located */ 
  virtual void setFilenames(string dir); 

};

#endif  