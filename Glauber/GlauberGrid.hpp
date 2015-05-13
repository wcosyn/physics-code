/*! \file GlauberGrid.hpp 
 * \brief Contains declaration of class GlauberGrid
 * \author Wim Cosyn
 * \date 16/08/2011
 * \addtogroup Glauber
 * @{
 */
#ifndef GLAUBERGRID_H
#define GLAUBERGRID_H

#include <cstdarg>

#include "AbstractFsiCTGrid.hpp"
#include <numint/numint.hpp>

/*! \brief A class for a RMSGA ISI/FSI grid, implementing the abstract AbstractFsiCTGrid class
 *
 * Typically an object that is an instance from this class is operated as follows.<BR>
 * 1. Initialize object with constructor GlauberGrid()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() or updateGrids() <BR>
 * 4. Add particles that are knocked out from nucleus  (using addKnockout()
)<BR>
 *  5. Interpolate grid for a certain point or print the grid or whatever... 
(use getFsiGridFull_interp3() for instance) <BR>
 * 
 * Every function that takes a \param grid argument: <BR>
 * 0: RMSGA <BR>
 * 1: RMSGA+CT <BR>
 */
class GlauberGrid : public AbstractFsiCTGrid{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnucl pointer to an instance of MeanFieldNucleus
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param dir std::string that contains dir with all input, should be the ./share subdir of the project!
   */
  GlauberGrid(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleus *pnucl, 
	      double prec, int integrator, std::string dir);
  virtual ~GlauberGrid();/*!< Destructor */
  virtual std::complex<double> getFsiGridFull_interp();   /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual std::complex<double> getFsiCtGridFull_interp(); /*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously
   * \param grid 0: RMSGA <BR>
  * 1: RMSGA+CT <BR>
  */
  virtual std::complex<double> getFsiGridN_interp(int grid); 
  
  virtual void printFsi_grid();/*!< Prints the FSI grid for a certain situation*/
  virtual void printFsi_ct_grid(); /*!< Prints the FSI+CT grid for a certain situation */
  /*! Prints a grid grid for a certain situation 
   * \param gridindex 0: RMSGA <BR>
  * 1: RMSGA+CT <BR>
  */  
  virtual void print_grid(int gridindex);
  
protected:
  std::complex<double> *****fsi_grid; /*!< grid that contains the fsi factor */
  std::complex<double> *****fsi_ct_grid; /*!< grid that contains the fsi+ct factor*/
  int ** treshold; /*!< array that checks if calculated glauberphases are close to one, then doesn't compute them for larger r, to save computing time*/
    
  
private:

/*! set filenames of the grids, includes "Gl" prefix 
   *\param dir dir where all input/output is located */ 
  virtual void setFilenames(std::string dir);  

  
  virtual void constructAllGrids(); /*!< construct both the fsi and the fsi+ct grids */
  virtual void constructCtGrid();  /*!< construct only the fsi+ct grids */
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
  virtual void readinFsiGrid(std::ifstream &infile); /*!< read in both the fsi and the fsi+ct grids */
  virtual void readinFsiCtGrid(std::ifstream &infile); /*!< read in only the fsi+ct grids */
  virtual void writeoutFsiGrid(std::ofstream &outfile); /*!< write out both the fsi and the fsi+ct grids */
  virtual void writeoutFsiCtGrid(std::ofstream &outfile); /*!< write out only the fsi+ct grids */
  
  /*! function that gets integrated over r, both fsi and fsi+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberR(const double r, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi and fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberCosTheta(const double costheta, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, both fsi and fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhi(const double phi, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over r, only fsi+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberRCT(const double r, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), only fsi+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberCosThetaCT(const double costheta, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, only fsi+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi+Ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhiCT(const double phi, std::complex<double> *results, va_list ap);
  
  /*! struct that is used for integrators */
  struct Ftor_one {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_one &p = * (Ftor_one *) param;
      p.f(ret,x[0],x[1],x[2],*p.grid,p.level,p.mm);
    }
    GlauberGrid *grid;/*!< pointer to the grid where the integration is performed */
    int level;/*!< knockout level */
    int mm; /*!< index in m_j */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param grid the grid instance
    * \param level knockout level
    * \param mm index in m_j
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, GlauberGrid & grid, int level, int mm);
  };
  /*! integrandum function (clean ones)*/
  static void klaas_one_bound(numint::vector_z &, double r, double costheta, double phi, GlauberGrid & grid, int level, int mm);
  /*! integrandum function (clean ones), only CT*/
  static void klaas_one_bound_ct(numint::vector_z &, double r, double costheta, double phi, GlauberGrid & grid, int level, int mm);

};




/** @} */
#endif