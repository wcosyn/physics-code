/*! \file AbstractFsiCTGrid.hpp 
 * \brief Contains declaration of abstract class AbstractFsiCTGrid
 * \author Wim Cosyn
 * \date 18/08/2011
 * 
 */

#ifndef ABSTRACTFSICTGRID_H
#define ABSTRACTFSICTGRID_H

#include <string>
#include <complex>
#include <vector>

#include "MeanFieldNucleus.hpp"
#include "AbstractFsiGrid.hpp"

/*! \brief An abstract class for a ISI/FSI with also a grid for CT ISI/FSI
 * 
 */
class AbstractFsiCTGrid : public virtual AbstractFsiGrid{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnucl pointer to an instance of MeanFieldNucleus
   * \param dir string that contains dir with all input
   */
  AbstractFsiCTGrid(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleus *pnucl, string dir);
  virtual ~AbstractFsiCTGrid(); /*!< Destructor */

  const string getFsi_Ct_Filename() const;/*!< returns filename for regular fsi+ct grid */

  virtual void printFsi_ct_grid()=0; /*!< Prints the FSI+CT grid for a certain situation, pure virtual function!! */

  /*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiCtGridFull_interpvec(const TVector3 &rvec);
  /*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,costheta,phi) */
  complex<double> getFsiCtGridFull_interp3(const double r, const double costheta, const double phi);
  /*!returns the value of the fsi+ct grid for a certain situation at coordinate (costheta,phi), r has been set previously */
  complex<double> getFsiCtGridFull_interp2(const double costheta, const double phi);
  /*!returns the value of the fsi+ct grid for a certain situation at coordinate (phi), r&theta have been set previously */
  complex<double> getFsiCtGridFull_interp1(const double phi);
  /*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual complex<double> getFsiCtGridFull_interp()=0;
  
  virtual void fillGrids();  /*!< fills the fsi grids that are used for interpolation */
  virtual void updateGrids();  /*!< updates the fsi grids that are used for interpolation */

protected:
  virtual void setFilenames(string dir); /*!< set filenames of the grids \param dir dir where all input/output is located */ 
  bool filledctgrid; /*!< denotes if the ct grid has been filled */
  string fsi_ct_filename; /*!< filename for fsi+ct grid */

  
private:

  virtual void constructCtGrid()=0;  /*!< construct only the fsi+ct grids */
  virtual void readinFsiCtGrid(ifstream &infile)=0; /*!< read in only the fsi+ct grids */
  virtual void writeoutFsiCtGrid(ofstream &outfile)=0; /*!< write out only the fsi+ct grids */
};

#endif
