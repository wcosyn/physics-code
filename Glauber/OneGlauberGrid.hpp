/*! \file OneGlauberGrid.hpp 
 * \brief Contains declaration of class GlauberGrid
 * \author Wim Cosyn
 * \date 16/08/2011
 * 
 */
#ifndef ONEGLAUBERGRID_H
#define ONEGLAUBERGRID_H

#include <cstdarg>

#include "GlauberGrid.hpp"

/*! \brief A class for a RMSGA ISI/FSI grid with one particle, exploiting the symmetry
 */
class OneGlauberGrid : public GlauberGrid{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param th_grid gridsize in theta
   * \param pnucl pointer to an instance of MeanFieldNucleus
   * \param dir string that contains dir with all input
   */
  OneGlauberGrid(const int r_grid, const int th_grid, MeanFieldNucleus *pnucl, string dir);
  virtual ~OneGlauberGrid();/*!< Destructor */
  /*! 2d interpolation of a grid after a certain coordinate (r,theta) has been set! */
  virtual complex<double> getInterp(complex<double> ***grid);  
  /*! add particle subject to isi/fsi, checks to see the vector only contains max one particle 
   *\param newparticle pointer to instance of FastParticle you want to add*/
  virtual void addParticle(FastParticle* newparticle);
  
  
private:
  /*! calculates the glauberphases for one gridpoint (both FSI and FSI+CT)
   * \param i grid index in r
   * \param j grid index in theta
   * \param k grid index in phi, irrelevant here (always 0)
   */
  virtual void calcGlauberphasesBoth(const int i, const int j, const int k); 
  /*! calculates the glauberphases for one gridpoint (only FSI+CT)
   * \param i grid index in r
   * \param j grid index in theta
   * \param k grid index in phi, irrelevant here (always 0)
   */
  virtual void calcGlauberphasesCt(const int i, const int j, const int k);
   /*! set filenames of the grids, includes "One" prefix 
    *\param dir dir where all input/output is located */ 
  virtual void setFilenames(string dir);
  
  
  /*! function that gets integrated over b, both fsi and fsi+ct grid output
   * \param b [fm] radial coordinate, cylindrical!!
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberB(const double b, double *results, va_list ap);
  /*! function that gets integrated over z, both fsi and fsi+ct grid output
   * \param z [fm] z coordinate, cylindrical!!
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberZ(const double z, double *results, va_list ap);
  /*! function that gets integrated over phi
   * \param phi [rad] phi coordinate, cylindrical!!!
   * \param result result: contains the glauberphases for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhi(const double phi, double *result, va_list ap);
  /*! function that gets integrated over b, both fsi and fsi+ct grid output
   * \param b [fm] radial coordinate, cylindrical!!
   * \param results result: contains the glauberphases (fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberBCT(const double b, double *results, va_list ap);
  /*! function that gets integrated over z, only fsi+ct grid output
   * \param z [fm] z coordinate, cylindrical!!
   * \param results result: contains the glauberphases (fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberZCT(const double z, double *results, va_list ap);
  
};

#endif