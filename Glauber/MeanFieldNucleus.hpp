/*! \file MeanFieldNucleus.hpp 
 * \brief Contains MeanFieldNucleus class declaration
 * \author Wim Cosyn
 * \date 16/08/2011
 * \addtogroup Glauber
 * @{
 */


#ifndef MEANFIELDNUCLEUS_H
#define MEANFIELDNUCLEUS_H

#include <string>
#include <complex>

#define GRIDP 201 /*!< \def defines number of gridpoints for the Y_kappa grids*/



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
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au], or use MeanFieldNucleus::[He|C|O|Fe|Pb|Al|Cu|Au], so you don't have to look up the numbering convention all the time. 
   * \param dir std::string containing the dir were all input is located
   */
  MeanFieldNucleus(const int nucleus, const std::string & dir = ".");
  MeanFieldNucleus(const MeanFieldNucleus&); /*!< copy constructor */
  ~MeanFieldNucleus();/*!< Destructor */
 
  /*! An enum so we don't have to look up which nucleus corresponds to which number all the time */
  enum Type { He, C, O, Fe, Pb, Al, Cu, Au };

  const std::string getNucleusName() const {return nucleusname;} /*!< Returns a std::string with the name of the nucleus, f.i. "C" */
  int getA() const {return A;} /*!<  Returns the total number of nucleons */
  int getZ() const {return Z;} /*!< Returns the total number of protons*/
  bool getOnlyOneProton() const {return onlyoneproton;} /*!< Returns bool that denotes if the final proton shell has an unpair number of protons*/
  int getFinalMProton() const { return finalmproton;}/*!< For not closed final proton shell, returns the m to which it is filled (arbitrary somewhat)*/
  bool getOnlyOneNeutron() const {return onlyoneneutron;} /*!< Returns bool that denotes if the final proton shell has an unpair number of neutron*/
  int getFinalMNeutron() const {return finalmneutron;} /*!< For not closed final neutron shell, returns the m to which it is filled (arbitrary somewhat)*/
  int getN() const {return N;}/*!< Returns the total number of neutrons*/
  int getPLevels() const { return plevels;}/*!< Returns total number of proton shells */
  int getNLevels() const {return nlevels;}/*!< Returns total number of neutron shells */
  int getTotalLevels() const {return totallevels;} /*!< Returns total number of shells*/
  const int* getN_array() const {return n_array;} /*!< Returns array with all the n quantum numbers of each shell*/
  const int* getKappas() const {return kappas;} /*!< Returns array with all the kappa quantum numbers of each shell*/
  const int* getL_array() const {return l_array;} /*!< Returns array with all the l quantum numbers of each shell*/
  const int* getLbar_array() const {return lbar_array;} /*!< Returns array with all the l_bar quantum numbers of each shell*/
  const int* getJ_array() const {return j_array;} /*!< Returns array with all the j quantum numbers of each shell, 2*j is returned so integer values!!!!*/
  double getMassA() const {return massA;} /*!< Returns the nucleus mass*/
  /*! Returns the mass of A-1 (either Z-1 or N-1 depending on level)
   * \param level shell level of minus 1 particle */
  double getMassA_min_1(int level) const {return level<getPLevels()?  getMassA_min_proton(): getMassA_min_neutron();} 
  double getMassA_min_proton() const {return massA_min_proton;} /*!< Returns the mass of the A-1(Z-1,N) nucleus*/
  double getMassA_min_neutron() const {return massA_min_neutron;} /*!< Returns the mass of the A-1(Z,N-1) nucleus*/
  double getMassA_min_pp() const {return massA_min_pp;} /*!< Returns the mass of the A-2(Z-2,N) nucleus*/
  double getMassA_min_pn() const {return massA_min_pn;} /*!< Returns the mass of the A-2(Z-1,N-1) nucleus*/
  double getMassA_min_nn() const {return massA_min_nn;} /*!< Returns the mass of the A-2(Z,N-2) nucleus*/
  const double* getExcitation() const {return excitation;} /*!< Returns an array with all the excitation energies of the shells*/
  const std::string getInputfile() const {return inputfile;} /*!< Returns a std::string with the location of the inputfile*/
  
  double getRange() const {return range;} /*!< Returns the range [fm] in r of the radial wave functions of the nucleons*/
  double getWF_r_step() const {return wf_r_step;} /*!< Returns stepsize of the grid for the radial wave functions*/
  int getWF_r_lines() const {return wf_r_lines;} /*!< Returns the grid size of the radial wave functions grid*/
  double** getF() const {return F;} /*!< Returns the grids for F(r) */
  double** getG() const {return G;} /*!< Returns the grids for G(r) */
  double*** getYkappa() const {return Ykappa;} /*!< Returns the grids for Y_kappa^2 */
  double*** getYminkappa() const {return Yminkappa;} /*!< Returns the grids for Y_{-kappa}^2 */
  
  /*! Computes the Dirac spinor of at (r,costheta,phi) for the nucleon from shell shellindex and orbital quantum number m
   * \param wave Dirac spinor of bound nucleon is stored here
   * \param shellindex shell of the nucleon
   * \param m j_z-proj quantum number of nucleon*2, so [-3 -1 1 3] for J=3/2 f.i.
   * \param r [fm] coordinate r, spherical
   * \param costheta cos{theta} of the coordinate, spherical
   * \param phi [rad] phi coordinate, spherical
   */
  void getWaveFunction(std::complex<double> *wave, const int shellindex, const int m, 
		  const double r, const double costheta, const double phi) const;
  /*! Returns value of F(r) for shellindex
   * \param shellindex shell of the nucleon
   * \param r [fm] coordinate r, spherical
   * \return [fm-{1/2}] F[shellindex](r)
   */
  double getWave_F(const int shellindex, const double r) const;
  /*! Returns value of F(r) for shellindex
   * \param shellindex shell of the nucleon
   * \param r [fm] coordinate r, spherical
   * \return [fm-{1/2}] G[shellindex](r) 
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
  std::complex<double> getWave_Phipart(const int m, const int spin, const double phi) const; 
  /*! Calculates the spherical harmonic Y_lm(x)
   * \param l orbital quantum number
   * \param m l_z
   * \param x argument, usually cos(theta)
   * \return Y_lm(x)
   */
  static double Spher_Harm(const int l, const int m, const double x);  //spherical harmonics without the phi part
  
  
private:
  std::string nucleusname; /*!< name of the nucleus*/
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
  std::string inputfile; /*!< inputfile for the nucleus*/
  std::string rgridfile;/*!< file location of the radial wave funtions*/

  /*! Set the inputfile for the nucleus
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param dir std::string containing the dir were all input is located
   */
  void setInputfile(const int nucleus, const std::string &dir);//set inputfile
  /*! Set the inputfile for the radial wave function grids
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param dir std::string containing the dir were all input is located
   */
  void setRgridfile(const int nucleus, const std::string &dir);//set rgridfilefile
  /*! Set the std::string containing the name of the nucleus
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   */
  void setNucleusName(const int nucleus); /*!< set the stirng with the name of the nucleus */
  
  double **F; /*!< Array with grid for F(r) of every shell*/
  double **G; /*!< Array with grid for G(r) of every shell*/
  double range; /*!< [fm]Range in r of the radial wave function grids*/
  double wf_r_step; /*!< stepsize of the radial wf grids*/
  int wf_r_lines; /*!< size of the radial wf grids*/
  
  /*! Read in the grids for F(r) and G(r) from file
   * \param nucleus an integer denoting which nucleus [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
   * \param dir std::string containing the dir were all input is located
   */
  void readRgrids(const int nucleus, const std::string &dir); 

  double ***Ykappa; /*!< Array with grid for Y_kappa(costheta) for every shell and m (only positive!!!,[0->1/2,1->3/2,etc]) */
  double ***Yminkappa; /*!< Array with grid for Y_{-kappa(costheta)} for every shell and m (only positive!!!,[0->1/2,1->3/2,etc]) */
  void constructThetaArray(); /*!< Construct the Y_kappa and Y_{-kappa} grids */
  
};
/** @} */
#endif
