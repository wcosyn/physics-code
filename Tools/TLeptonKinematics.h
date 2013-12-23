/*! \file TLeptonKinematics.h 
 * \brief Contains TLeptonKinematics class declaration
 * \author Wim Cosyn
 * \date 12/12/2013
 * \addtogroup MePhys
 * @{
 */

#ifndef TLEPTONKINEMATICS_H
#define TLEPTONKINEMATICS_H

#include "TObject.h"
#include <FourVector.h>

// Forward declarations
class TKinematics2to2;
class TKinematics2to3;

/*! \brief A Class that stores the lepton kinematics for a certain charged weak scattering process with incoming neutrino
 * In electroproduction reactions we separate the electron scattering
 * vertex from the virtual photon - deuteron system.
 * This class deals with the kinematics of the former. In conjunction
 * with a TKinematics2to3 object the kinematics are completely fixed by
 * the incoming electron's beam energy
 * or the outgoing electron's scattering angle
 * Objects should be created with the functions CreateWithBeamEnergy or CreateWithCosScatterAngle, 
 * these return a pointer to a newly created instance.  The user has to delete the Object after he is done with it!
 * Kinematics are solved when you pass a TKinematics2to2 or TKinematics2to3 object to the Getters.
 */

class TLeptonKinematics : public TObject
{
 public:
  TLeptonKinematics(TRootIOCtor*); /*!< ROOT I/O Constructor */
  TLeptonKinematics(const TLeptonKinematics&); /*!< copy constructor */
  TLeptonKinematics& operator=(const TLeptonKinematics&); /*!< assignment operator */
  virtual ~TLeptonKinematics() {} /*!< Destructor */

  virtual TLeptonKinematics* Clone(const char * ="") const { return new TLeptonKinematics(*this); } /*!< ROOT Cloner */

  enum Lepton { electron = 0,
		       muon = 1,
		       tau = 2 };

  /*! indirect constructorcreator 
   * \param lep lepton type
   * \param beam [MeV] incoming beam energy */
  static TLeptonKinematics* CreateWithBeamEnergy(Lepton lep, double beam); 
//   static TLeptonKinematics* CreateWithEpsilon(double);
  /*! indirect constructorcreator 
   * \param lep lepton type
   * \param costheta_l [-1,1] cosine of the lepton scattering angle */
  static TLeptonKinematics* CreateWithCosScatterAngle(Lepton lep, double costheta_l);

  void SetBeamEnergy(double beam ); /*!< \param beam [MeV] beam energy */
//   void SetEpsilon(double); 
  void SetCosScatterAngle(double costheta_l); /*!< \param costheta_l [-1,1] cosine of the lepton scattering angle */

  double GetBeamEnergy(const TKinematics2to2&) const; /*!< \return [MeV] energy of incoming beam */
//   double GetEpsilon(const TKinematics2to2&) const; 
  double GetCosScatterAngle(const TKinematics2to2&) const; /*!< \return costheta_l [-1,1] cosine of the lepton scattering angle */
//   double GetTan2HalfAngle(const TKinematics2to2&) const;
  double GetBeamEnergy(const TKinematics2to3&) const;
//   double GetEpsilon(const TKinematics2to3&) const; /*!< \return costheta_l [-1,1] cosine of the lepton scattering angle */
  double GetCosScatterAngle(const TKinematics2to3&) const;
  double GetLeptonMass() const{return mass;} /*!< [MeV] returns lepton mass */
//   double GetTan2HalfAngle(const TKinematics2to3&) const;
  static const double masse; /*!< [MeV] electronmass */
  static const double massmu; /*!< [MeV] muon mass */
  static const double masstau; /*!< [MeV] tau mass */
  
  /*! constructs the four vectors for the incoming and outgoing lepton
   * \param kin kinematics of the hadron piece
   * \param k_in [MeV] incoming lepton fourvector
   * \param k_out [MeV] outgoing lepton fourvector
   */
  void GetLeptonVectors(const TKinematics2to2& kin, FourVector<double> &k_in, FourVector<double> &k_out); 
  
 protected:
  /*! basic constructor
   * \param type outgoing lepton type
   */
  TLeptonKinematics(Lepton type); 
//   TLeptonKinematics(); /*!< default constructor */
  /*! Solves the kinematics 
   * \param kin contains all gauge boson + A -> B+C kinematics */
  bool SolveKinematics(const TKinematics2to2& kin) const; 
  /*! Solves the kinematics 
   * \param kin contains all gauge boson + A -> B+C+D kinematics */
  bool SolveKinematics(const TKinematics2to3& kin) const;
  
 private:
  /*! Controls input variable */
  enum InputVariable { kBeamEnergy = 0,
		       kEpsilon = 1,
		       kScatterAngle = 2 };

  InputVariable  fInput;        /*!< specifies the input variable */
  mutable double fBeamEnergy;   /*!< energy incoming electron beam in [MeV] */
//   mutable double fEpsilon;      /*!< virtual photon's transverse linear polarization */
  mutable double fScatterAngle; /*!< cos of electron scattering angle in ]-1,1] */
//   mutable double ftan2HalfAngle; /*!< tan squared of half of the lepton scattering angle */
  double mass; /*!< [MeV] mass of outgoing lepton */
  
  ClassDef(TLeptonKinematics,1); /*!< Root shizzle */

}; 
/** @} */
#endif 
