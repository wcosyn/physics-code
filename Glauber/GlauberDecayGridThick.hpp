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
#include <numint/numint.hpp>

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
 * 4. Add particles that are knocked out from nucleus  (using addKnockout()
)<BR>
 *  5. Interpolate grid for a certain point or print the grid or whatever... 
(use getFsiGridFull_interp3() for instance) <BR>
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
   * \param dir std::string that contains dir with all input, should be the ./share subdir of the project!
   */
  GlauberDecayGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick *pnuclthick, 
		   double prec, int integrator, std::string dir);
  /*!< Destructor */
  virtual ~GlauberDecayGridThick();
  virtual std::complex<double> getFsiGridFull_interp(); /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual std::complex<double> getFsiCtGridFull_interp();/*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual std::complex<double> getFsiSrcGridFull_interp(); /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual std::complex<double> getFsiSrcCtGridFull_interp();/*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual std::complex<double> getFsiDecayGridFull_interp(); /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual std::complex<double> getFsiCtDecayGridFull_interp();/*!returns the value of the fsi+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual std::complex<double> getFsiSrcDecayGridFull_interp(); /*!returns the value of the fsi+src grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual std::complex<double> getFsiSrcCtDecayGridFull_interp();/*!returns the value of the fsi+src+ct grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously
   * \param grid 0: RMSGA <BR>
  * 1: RMSGA+SRC <BR>
  * 2: RMSGA+DECAY <BR>
  * 3: RMSGA+SRC+DECAY <BR>
  * 4: RMSGA+CT <BR>
  * 5: RMSGA++CT+SRC <BR>
  * 6: RMSGA+CT+DECAY <BR>
  * 7: RMSGA+CT+SRC+DECAY AKA "THE WORKS" <BR> */
  virtual std::complex<double> getFsiGridN_interp(int grid); 
  
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
  std::complex<double> *****fsi_grid; /*!< grid that contains the fsi factor, also has the fsi+src grid */
  std::complex<double> *****fsi_ct_grid; /*!< grid that contains the fsi+ct factor, also has the fsi+src+ct grid */
  int ** treshold; /*!< array that checks if calculated glauberphases are close to one, then doesn't compute them for larger r, to save computing time*/
  
  
private:

  
  virtual void setFilenames(std::string dir); /*!< set filenames of the grids \param dir dir where all input/output is located */ 

  
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
  void intGlauberR(const double r, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi(+src) and fsi(+src)+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberCosTheta(const double costheta, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, both fsi(+src) and fsi(+src)+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhi(const double phi, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over r, only fsi(+src)+ct grid output
   * \param r [fm] radial coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberRCT(const double r, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over cos(theta), only fsi(+src)+ct grid output
   * \param costheta [] cos of theta coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberCosThetaCT(const double costheta, std::complex<double> *results, va_list ap);
  /*! function that gets integrated over phi, only fsi(+src)+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhiCT(const double phi, std::complex<double> *results, va_list ap);
  
  /*! struct that is used for integrators (dirty ones)*/
  struct Ftor_all {
    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_all &p = * (Ftor_all *) param;
      p.f(ret,x[0],x[1],x[2],*p.grid,p.proton);
    }
    GlauberDecayGridThick *grid; /*!< pointer to the grid where the integration is performed */
    int proton; /*!< proton glauberphase or not */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param grid the grid instance
    * \param proton proton glauberphase or not
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, GlauberDecayGridThick & grid, int proton);
  };
 /*! struct that is used for integrators (clean ones)*/
   struct Ftor_bound {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor_bound &p = * (Ftor_bound *) param;
      p.f(ret,x[0],x[1],x[2],*p.grid,p.proton);
    }
    GlauberDecayGridThick *grid;/*!< pointer to the grid where the integration is performed */
    int proton;/*!< proton glauberphase or not */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param grid the grid instance
    * \param proton proton glauberphase or not
    */
    void (*f)(numint::vector_d &, double x, double y, double z, GlauberDecayGridThick &, int proton);
  };

  /*! integrandum function  (dirty ones)*/
  static void klaas_int_all(numint::vector_z & ret, double r, double costheta, double phi, GlauberDecayGridThick & grid, int proton);
  /*! integrandum function (clean ones)*/
  static void klaas_int_bound(numint::vector_d &, double b, double z, double phi, GlauberDecayGridThick & grid, int proton);
  /*! integrandum function  (dirty ones), only CT*/
  static void klaas_int_all_ct(numint::vector_z & ret, double r, double costheta, double phi, GlauberDecayGridThick & grid, int proton);
  /*! integrandum function (clean ones), only CT*/
  static void klaas_int_bound_ct(numint::vector_d &, double b, double z, double phi, GlauberDecayGridThick & grid, int proton);
  /*! function that gets integrated over , both fsi(+src) and fsi(+src)+ct grid output
   * \param b [fm] radial coordinate in perp plane
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberb_bound(const double b, double *results, va_list ap);
  /*! function that gets integrated over cos(theta), both fsi(+src) and fsi(+src)+ct grid output
   * \param z [fm] z-axis (cylindrical)
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberz_bound(const double z, double *results, va_list ap);
  /*! function that gets integrated over phi, both fsi(+src) and fsi(+src)+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhi_bound(const double phi, double *results, va_list ap);
  /*! function that gets integrated over ,  only fsi(+src)+ct grid output
   * \param b [fm] radial coordinate in perp plane
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberb_bound_ct(const double b, double *results, va_list ap);
  /*! function that gets integrated over cos(theta),  only fsi(+src)+ct grid output
   * \param z [fm] z-axis (cylindrical)
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberz_bound_ct(const double z, double *results, va_list ap);
  /*! function that gets integrated over phi,  only fsi(+src)+ct grid output
   * \param phi [rad] phi coordinate
   * \param results result: contains the glauberphases (fsi(+src) and fsi(+src)+ct) for a gridpoint
   * \param ap variable parameter list
   */
  void intGlauberPhi_bound_ct(const double phi, double *results, va_list ap);

};
/** @} */
#endif
