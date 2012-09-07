/*! \file DoubleNModel.hpp 
 * \brief Contains declaration of class Model, used to compute A(e,e'N) amplitudes
 * \author Wim Cosyn
 * \date 28/08/2012
 * 
 * \addtogroup Knockout
 * @{
 */

#ifndef DOUBLENMODEL_HPP
#define DOUBLENMODEL_HPP

#include <FourVector.h>
#include <MeanFieldNucleusThick.hpp>
#include <GlauberGridThick.hpp>
#include <TKinematics2to3WithLabAngles.h>
#include <NucleonEMOperator.hpp>
#include <TVector3.h>
#include <pair.h>

#include <string>
#include <cstdarg>

using namespace std;

/*! \brief A class used to compute A(e,e'NN) amplitudes, with possibility to include SRC in the wave funtions */
class DoubleNModel{
public:
    /*! Constructor
   * \param pnucleus pointer to a MF nucleus instance
   * \param setpw turn off FSI?
   * \param setSRC do you want SRC in your FSI?
   * \param setCT include CT in the FSI?
   * \param setcorr do you want to use correlated wave functions (from Maarten's code)
   * \param particletype1 type of first ejected nucleon
   * \param particletype2 type of second ejected nucleon
   * \param prec precision you want in the integrations
   * \param integrator which integrator (0:Wim's romberg fubini sort of thing, 1:Klaas thingy, 2:adaptive MIT thingy
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   */
  DoubleNModel(MeanFieldNucleusThick *pnucleus, bool setpw, bool setSRC, bool setCT, bool setcorr,
	       int particletype1, int particletype2, double prec, int integrator, string dir);
  ~DoubleNModel();  /*!< Destructor */
  void setSRC(int set_SRC) {SRC=set_SRC;}
  /*! compute matrix element 
   * \param tk contains 2to3 hadron kinematics
   * \param spin1 helicity of first final nucleon
   * \param spin2 helicity of second final nucleon
   * \param photonpol polarization of photon (-1,0,1)
   * \param shellindex1 shell of first ejected nucleon
   * \param shellindex2 shell of second ejected nucleon
   * \param m1 m_j of first ejected nucleon (times TWO!!!)
   * \param m2 m_j of second ejected nucleon (times TWO!!!)
   * \return matrix element [fm^3]
   */
  complex<double> getMatrixEl(const TKinematics2to3 &tk,int spin1, int spin2, int photonpol, 
			      int shellindex1, int shellindex2, int m1, int m2);
  
  /*! helper function for the MC integration, returns the function that gets integrated matrix element in r space
   * \param r1 [fm] spherical coordinates of first nucleon interaction
   * \param costheta1
   * \param phi1
   * \param r2 [fm] spherical coordinates of second nucleon
   * \param costheta2
   * \param phi2
   * \return matrix element in r-space [fm^-3]
   */
  complex<double> MC_helper(const double r1, const double costheta1, const double phi1, const double r2, const double costheta2, const double phi2) const;  
  bool getPW() const {return pw;}
  bool getSRC() const {return SRC;}
  bool getCT() const {return CT;}
  bool getCorr() const {return corr;}
  complex<double> getJ1contrup() const {return J1contrup;}
  complex<double> getJ1contrdown() const {return J1contrdown;}
  complex<double> getJ2contrup() const {return J2contrup;}
  complex<double> getJ2contrdown() const {return J2contrdown;}
  Pair * getPair1up() const {return pair1up;}
  Pair * getPair1down() const {return pair1down;}
  Pair * getPair2up() const {return pair2up;}
  Pair * getPair2down() const {return pair2down;}
  GlauberGridThick * getGridf1() const {return gridf1;}
  GlauberGridThick * getGridf2() const {return gridf2;}
  TVector3 getPf1() const {return pf1;}
  TVector3 getPf2() const {return pf2;}
  TVector3 getQvec3() const {return qvec3;}
  double getPrec() const {return prec;}
  
  
private:
  bool pw; /*!< plane-wave or not*/
  bool SRC;/*!< SRC in the FSI?*/
  bool CT;/*!< CT in the fsi?*/
  bool corr;  /*!< use correlated wave functions*/
  int particletype1; /*!< particletype of first nucleon*/
  int particletype2; /*!< particletype of second nucleon */
  double prec; /*!< precision in the integrations */
  int integrator; /*!< choice of integrator */
  MeanFieldNucleusThick *pnucl; /*!< pointer to nucleus */
  NucleonEMOperator *J1; /*!< nucleon FF for interaction with nucleon 1*/
  NucleonEMOperator *J2; /*!< nucleon FF for interaction with nucleon2*/
//   Matrix<1,4> barcontract;
//   Matrix<1,4> barcontract0up;
//   Matrix<1,4> barcontractminup;
//   Matrix<1,4> barcontractplusup;
//   Matrix<1,4> barcontract0down;
//   Matrix<1,4> barcontractmindown;
//   Matrix<1,4> barcontractplusdown;
  TVector3 pf1; /*!< final nucleon1 momentum*/
  TVector3 pf2; /*!< final nucleon2 momentum */
  TVector3 qvec3; /*!< virtual photon momentum */
  complex<double> J1contrup; /*!< intermediate contraction of spinor with current 4*4 operator*/
  complex<double> J1contrdown;/*!< intermediate contraction of spinor with current 4*4 operator*/
  complex<double> J2contrup;/*!< intermediate contraction of spinor with current 4*4 operator*/
  complex<double> J2contrdown;/*!< intermediate contraction of spinor with current 4*4 operator*/
  
  Pair *pair1up; /*!< Pair wave function for photon interacting with nucleon 1, intermediate spin up*/
  Pair *pair1down;/*!< Pair wave function for photon interacting with nucleon 1, intermediate spin down*/
  Pair *pair2up; /*!< Pair wave function for photon interacting with nucleon 2, intermediate spin up*/
  Pair *pair2down;/*!< Pair wave function for photon interacting with nucleon 2, intermediate spin down*/
  
  GlauberGridThick *gridf1; /*!< pointer to glaubergrid for first fast nucleon */
  GlauberGridThick *gridf2; /*!< pointer to glaubergrid for second fast nucleon */
  string homedir; /*!< share dir where all input is located */
  
  

};
/*! MC integration of the real part
 * \param x integration variables (r1, r2)
 * \param dim dimension of x[], 6
 * \param p pointer to instance of DoubleNModel
 * \return MC result
 */
double matrixel_MC_real (double* x, size_t dim, void * p);
/*! MC integration of the imag part
 * \param x integration variables (r1, r2)
 * \param dim dimension of x[], 6
 * \param p pointer to instance of DoubleNModel
 * \return MC result
 */
double matrixel_MC_imag (double* x, size_t dim, void * p);
/*! displays MC results with error
 */
void display_results (char *title, double result, double error);
/*! @} */

#endif