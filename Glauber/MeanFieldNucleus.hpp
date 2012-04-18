/*! \file MeanFieldNucleus.hpp 
 * \brief Contains MeanFieldNucleus class declaration
 * \author Wim Cosyn
 * \date 16/08/2011
 * 
 */


#ifndef MEANFIELDNUCLEUS_H
#define MEANFIELDNUCLEUS_H

#include <string>
#include <complex>

#define GRIDP 201 /*!< \def defines number of gridpoints for the Y_kappa grids*/


using namespace std;

#include "constants.hpp"

/*! \brief A Class for a nucleus with mean-field wave functions for the individual nucleons
 * 
 * Contains all useful variables for the nucleus (A,Z,N,excitation energy of the levels, etc)
 * Provides all quantum numbers for individual nucleons.  Functions to compute the wave function for a certain
 * nucleon at a certain coordinate.  
 */
class MeanFieldNucleus{
public:
  /*! \brief Constructor
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param dir string containing the dir were all input is located
   */
  MeanFieldNucleus(const int nucleus = 0, const string & dir = ".");
  
  ~MeanFieldNucleus();/*!< Destructor */
  
  const string getNucleusName() const; /*!< Returns a string with the name of the nucleus, f.i. "C" */
  int getA() const; /*!<  Returns the total number of nucleons */
  int getZ() const; /*!< Returns the total number of protons*/
  bool getOnlyOneProton() const; /*!< Returns bool that denotes if the final proton shell has an unpair number of protons*/
  int getFinalMProton() const;/*!< For not closed final proton shell, returns the m to which it is filled (arbitrary somewhat)*/
  bool getOnlyOneNeutron() const; /*!< Returns bool that denotes if the final proton shell has an unpair number of neutron*/
  int getFinalMNeutron() const; /*!< For not closed final neutron shell, returns the m to which it is filled (arbitrary somewhat)*/
  int getN() const;/*!< Returns the total number of neutrons*/
  int getPLevels() const;/*!< Returns total number of proton shells */
  int getNLevels() const;/*!< Returns total number of neutron shells */
  int getTotalLevels() const; /*!< Returns total number of shells*/
  const int* getN_array() const; /*!< Returns array with all the n quantum numbers of each shell*/
  const int* getKappas() const; /*!< Returns array with all the kappa quantum numbers of each shell*/
  const int* getL_array() const; /*!< Returns array with all the l quantum numbers of each shell*/
  const int* getLbar_array() const; /*!< Returns array with all the l_bar quantum numbers of each shell*/
  const int* getJ_array() const; /*!< Returns array with all the j quantum numbers of each shell*/
  double getMassA() const; /*!< Returns the nucleus mass*/
  double getMassA_min_1(int level) const; /*!< Returns the mass of A-1 (either Z-1 or N-1 depending on level)*/
  double getMassA_min_proton() const; /*!< Returns the mass of the A-1(Z-1,N) nucleus*/
  double getMassA_min_neutron() const; /*!< Returns the mass of the A-1(Z,N-1) nucleus*/
  double getMassA_min_pp() const {return massA_min_pp;}; /*!< Returns the mass of the A-2(Z-2,N) nucleus*/
  double getMassA_min_pn() const {return massA_min_pn;}; /*!< Returns the mass of the A-2(Z-1,N-1) nucleus*/
  double getMassA_min_nn() const {return massA_min_nn;}; /*!< Returns the mass of the A-2(Z,N-2) nucleus*/
  const double* getExcitation() const; /*!< Returns an array with all the excitation energies of the shells*/
  const string getInputfile() const; /*!< Returns a string with the location of the inputfile*/
  
  double getRange() const; /*!< Returns the range [fm] in r of the radial wave functions of the nucleons*/
  double getWF_r_step() const; /*!< Returns stepsize of the grid for the radial wave functions*/
  int getWF_r_lines() const; /*!< Returns the grid size of the radial wave functions grid*/
  double** getF() const; /*!< Returns the grids for F(r) */
  double** getG() const; /*!< Returns the grids for G(r) */
  double*** getYkappa() const; /*!< Returns the grids for Y_kappa^2 */
  double*** getYminkappa() const; /*!< Returns the grids for Y_{-kappa}^2 */
  
  /*! Computes the Dirac spinor of at (r,costheta,phi) for the nucleon from shell shellindex and orbital quantum number m
   * \param wave Dirac spinor of bound nucleon is stored here
   * \param shellindex shell of the nucleon
   * \param m j_z-proj quantum number of nucleon*2, so [-3 -1 1 3] for J=3/2 f.i.
   * \param r [fm] coordinate r, spherical
   * \param costheta cos{theta} of the coordinate, spherical
   * \param phi [rad] phi coordinate, spherical
   */
  void getWaveFunction(complex<double> *wave, const int shellindex, const int m, 
		  const double r, const double costheta, const double phi) const;
  /*! Returns value of F(r) for shellindex
   * \param shellindex shell of the nucleon
   * \param r [fm] coordinate r, spherical
   * \return F[shellindex](r)
   */
  double getWave_F(const int shellindex, const double r) const;
  /*! Returns value of B(r) for shellindex
   * \param shellindex shell of the nucleon
   * \param r [fm] coordinate r, spherical
   * \return B[shellindex](r)
   */
  double getWave_G(const int shellindex, const double r) const;
  /*! Returns value of Y^2_kappa(r) for shellindex
   * \param shellindex shell of the nucleon
   * \param m j_z-proj quantum number of nucleon, only positive due to symmetry!!! [0->1/2,1->3/2,etc]
   * \param costheta cos{theta} of the coordinate, spherical
   * \return Y^2_kappa[shellindex][m](r)
   */
  double getYkappacos(const int shellindex, const int m, const double costheta) const;
  /*! Returns value of Y^2_{-kappa}(r) for shellindex
   * \param shellindex shell of the nucleon
   * \param m j_z-proj quantum number of nucleon, only positive due to symmetry!!! [0->1/2,1->3/2,etc]
   * \param costheta cos{theta} of the coordinate, spherical
   * \return Y^2_{-kappa}[shellindex][m](r)
   */
  double getYminkappacos(const int shellindex, const int m, const double costheta) const;
  /*! Returns phi part of the wave function for J_z quantum number m and spin proj spin
   * \param m j_z-proj quantum number of nucleon*2, so [-3 -1 1 3] for J=3/2 f.i.
   * \param spin spin proj quantum number 1=up -1=down
   * \param phi [rad] phi coordinate, spherical
   * \return exp(I*0.5*phi*double(m-spin)); 
   */  
  complex<double> getWave_Phipart(const int m, const int spin, const double phi) const; 
  /*! Calculates the spherical harmonic Y_lm(x)
   * \param l orbital quantum number
   * \param m l_z
   * \param x argument, usually cos(theta)
   * \return Y_lm(x)
   */
  static double Spher_Harm(const int l, const int m, const double x);  //spherical harmonics without the phi part
  
  
private:
  string nucleusname;
  int A; /*!< Total nucleons*/
  int Z; /*!< Protons*/
  int N; /*!< Neutrons*/
  bool onlyoneproton; /*!< bool that denotes if the final proton shell has an unpair number of protons */
  int finalmproton; /*!< For not closed final proton shell, the m to which it is filled (arbitrary somewhat)*/
  bool onlyoneneutron; /*!<bool that denotes if the final proton shell has an unpair number of neutrons */
  int finalmneutron; /*!<For not closed final neutron shell, the m to which it is filled (arbitrary somewhat) */
  int plevels; /*!< number of proton shells*/
  int nlevels; /*!< number of neutron shells*/
  int totallevels; /*!< total number of shells*/
  int *n_array; /*!< array with the n quantum number for every shell*/
  int *kappas; /*!< array with the kappa quantum number for every shell*/
  int *l_array;/*!< array with the l quantum number for every shell*/
  int *lbar_array; /*!< array with the l_bar quantum number for every shell*/
  int *j_array; /*!< array with the j quantum number for every shell (2*j!!!!!, so integer)*/
  double massA; /*!< mass of the nucleus */
  double massA_min_proton; /*!< mass of A-1(Z-1)*/
  double massA_min_neutron; /*!< mass of A-1(Z)*/
  double massA_min_pp; /*!< mass of A-2(Z-2)*/
  double massA_min_pn; /*!< mass of A-2(Z-1)*/
  double massA_min_nn; /*!< mass of A-2(Z)*/
  double *excitation; /*!< array with excitation energy of every shell*/
  string inputfile; /*!< inputfile for the nucleus*/
  string rgridfile;/*!< file location of the radial wave funtions*/

  /*! Set the inputfile for the nucleus
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param dir string containing the dir were all input is located
   */
  void setInputfile(const int nucleus, const string &dir);//set inputfile
  /*! Set the inputfile for the radial wave function grids
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param dir string containing the dir were all input is located
   */
  void setRgridfile(const int nucleus, const string &dir);//set rgridfilefile
  /*! Set the string containing the name of the nucleus
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   */
  void setNucleusName(const int nucleus);
  
  double **F; /*!< Array with grid for F(r) of every shell*/
  double **G; /*!< Array with grid for G(r) of every shell*/
  double range; /*!< [fm]Range in r of the radial wave function grids*/
  double wf_r_step; /*!< stepsize of the radial wf grids*/
  int wf_r_lines; /*!< size of the radial wf grids*/
  
  /*! Read in the grids for F(r) and G(r) from file
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param dir string containing the dir were all input is located
   */
  void readRgrids(const int nucleus, const string &dir); 

  double ***Ykappa; /*!< Array with grid for Y_kappa(costheta) for every shell and m (only positive!!!,[0->1/2,1->3/2,etc]) */
  double ***Yminkappa; /*!< Array with grid for Y_{-kappa(costheta)} for every shell and m (only positive!!!,[0->1/2,1->3/2,etc]) */
  void constructThetaArray(); /*!< Construct the Y_kappa and Y_{-kappa} grids */
  
};
#endif