/*! \file ROMEAGrid.hpp 
 * \brief Contains declaration of class ROMEAGrid
 * \author Wim Cosyn
 * \date 19/05/2014
 * \addtogroup Glauber
 * @{
 */
#ifndef ROMEAGRID_H
#define ROMEAGRID_H

#include <cstdarg>
#include <fstream>

#include "AbstractFsiGrid.hpp"
#include <numint/numint.hpp>

/*! \brief A class for a ROMEA ISI/FSI grid, implementing the abstract AbstractFsiGrid class
 *
 * We use optical potential from the most recent global of Cooper, Hama (Phys.Rev. C80 (2009) 034605) <BR>
 * Typically an object that is an instance from this class is operated as follows.<BR>
 * 1. Initialize object with constructor ROMEAGrid()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() [only nucleons here!] <BR>
 * 3. Add particles that are knocked out from nucleus  (using addKnockout() 
), DO THIS BEFORE STEP 4, otherwise you will be computing for the wrong nucleus...<BR>
 * 4. Call fillGrids() or updateGrids() <BR>
 *  5. Interpolate grid for a certain point or print the grid or whatever... 
(use getFsiGridFull_interp3() for instance) <BR>
 */
class ROMEAGrid : public AbstractFsiGrid{
public:
  /*! Different fits for the optical potential.  See paper Cooper et al. (Phys.Rev. C80 (2009) 034605) <BR>
   * FIT1:  (1) CA40-PB208 (P,P) 65-1040 MEV  FIT.1  <BR>
     FIT2: (2) CA40-PB208 (P,P) 65-1040 MEV  FIT.2  <BR>
     EDAIC12: (3) C12 (P,P)   29-1040 MEV (EDAI C12) <BR>
     EDAIO16: (4) O16 (P,P)   23-1040 MEV (EDAI O16) <BR>
     EDAICA40: (5) CA40 (P,P)  21-1040 MEV (EDAI CA40) <BR>
     EDAIZR90: (6) ZR90 (P,P)  22-800  MEV (EDAI ZR90) <BR>
     EDAIPB208: (7) PB208 (P,P) 21-1040 MEV (EDAI PB208) <BR>
     EDAD1: (8) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.1) <BR>
     EDAD2: (9) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.2) <BR>
    EDAD3: (10) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.3) <BR>
    UNDEMOCRATICFIT2: (11) HE4-PB208 (P,P) 21-1040 MEV (UNDEMOCRATIC-Fit2)  <BR>
    DEMOCRATICFIT1: (12) HE4-PB208 (P,P) 21-1040 MEV (DEMOCRATIC-Fit1)  <BR>
  */
  enum fit_type { FIT1=1,
		  FIT2=2,
		  EDAIC12=3,
		  EDAIO16=4,
		  EDAICA40=5,
		  EDAIZR90=6,
		  EDAIPB208=7,
		  EDAD1=8,
		  EDAD2=9,
		  EDAD3=10,
		  UNDEMOCRATICFIT2=11,
		  DEMOCRATICFIT1=12
  };
   /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnucl pointer to an instance of MeanFieldNucleus
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param opticalfit  Different fits for the optical potential.  See paper Cooper et al. (Phys.Rev. C80 (2009) 034605) <BR>
   * FIT1:  (1) CA40-PB208 (P,P) 65-1040 MEV  FIT.1  <BR>
     FIT2: (2) CA40-PB208 (P,P) 65-1040 MEV  FIT.2  <BR>
     EDAIC12: (3) C12 (P,P)   29-1040 MEV (EDAI C12) <BR>
     EDAIO16: (4) O16 (P,P)   23-1040 MEV (EDAI O16) <BR>
     EDAICA40: (5) CA40 (P,P)  21-1040 MEV (EDAI CA40) <BR>
     EDAIZR90: (6) ZR90 (P,P)  22-800  MEV (EDAI ZR90) <BR>
     EDAIPB208: (7) PB208 (P,P) 21-1040 MEV (EDAI PB208) <BR>
     EDADFIT1: (8) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.1) <BR>
     EDADFIT2: (9) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.2) <BR>
    EDADFIT3: (10) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.3) <BR>
    UNDEMOCRATICFIT2: (11) HE4-PB208 (P,P) 21-1040 MEV (UNDEMOCRATIC-Fit2)  <BR>
    DEMOCRATICFIT1: (12) HE4-PB208 (P,P) 21-1040 MEV (DEMOCRATIC-Fit1)  <BR>
   * \param dir std::string that contains dir with all input, should be the ./share subdir of the project!
   */
  ROMEAGrid(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleus *pnucl, 
	      double prec, int integrator, fit_type opticalfit, std::string dir);
  virtual ~ROMEAGrid();/*!< Destructor */
  virtual std::complex<double> getFsiGridFull_interp();   /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual std::complex<double> getFsiGridN_interp(const int grid); /*!returns the value of the fsi grid for a certain situation at coordinate (r,theta,phi) that has been set previously*/
  virtual void printFsi_grid();/*!< Prints the FSI grid for a certain situation of ejected particles*/
  virtual void print_grid(int gridindex);/*!< Prints the FSI grid for a certain situation of ejected particles*/
  std::complex<double> ** getOpticalpotential() {return opticalpotential;} /*!< returns pointer to the start of the opticalpotential array */
  double getOpticalstep(){return opticalstep;} /*!< [fm] stepsize in the grid opticalpotential */
  int getOpticalsize(){return opticalsize;} /*!< size in elements of the grid opticalpotential */
  virtual void constructAllGrids(); /*!< construct the optical fsi grid */
 
protected:
  std::complex<double> ***fsi_grid; /*!< grid that contains the fsi factor */
  int ** treshold; /*!< array that checks if calculated glauberphases are close to one, then doesn't compute them for larger r, to save computing time*/
  std::complex<double> **opticalpotential; /*!< grid that contains the central part of the optical potential */
  double opticalstep; /*!< [fm] stepsize in the grid opticalpotential */
  int opticalsize; /*!< size in elements of the grid opticalpotential */
  int fit; /*!< type of fit used for the optical potential */
  
private:

/*! set filenames of the grids, includes "ROMEA" prefix 
   *\param dir dir where all input/output is located */ 
  virtual void setFilenames(std::string dir);  

  
  /*! calculates the glauberphases for one gridpoint (both FSI and FSI+CT)
   * \param i grid index in r
   * \param j grid index in costheta
   * \param k grid index in phi
   */
  virtual void calcROMEAphase(const int i, const int j, const int k); 
  virtual void readinFsiGrid(std::ifstream &infile); /*!< read in both the fsi grid */
  virtual void writeoutFsiGrid(std::ofstream &outfile); /*!< write out both the fsi grid */
  
  
  /*! struct that is used for integrators */
  struct Ftor_one {

    /*! integrandum function */
    static void exec(const numint::array<double,1> &z, void *param, numint::vector_z &ret) {
      Ftor_one &p = * (Ftor_one *) param;
      p.f(ret,z[0],*p.grid, p.it);
    }
    ROMEAGrid *grid;/*!< pointer to the grid where the integration is performed */
    size_t it;/*!< iterator for particles array */
    /*! integrandum 
    * \param res results
    * \param z first integration variable
    * \param grid the grid instance
    * \param it index of knocket out particle
    */
    void (*f)(numint::vector_z & res, double z, ROMEAGrid & grid, size_t it);
  };
  /*! integrandum function (clean ones)*/
  static void klaas_one_bound(numint::vector_z &, double z, ROMEAGrid & grid, size_t it);

};




/** @} */
#endif