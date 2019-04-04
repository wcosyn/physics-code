#ifndef DEUT_GPD_T_SET_HPP
#define DEUT_GPD_T_SET_HPP

/**
 * @brief small class that holds deuteron (isosinglet) vector helicity amplitudes for a certain kinematics
 * 
 */
class Deut_GPD_T_set{

public:
/**
 * @brief constructor
 * //order of helicity amplitudes is ++,--,00,0+,-0,+0,0-,-+,+-
 */
Deut_GPD_T_set(const double amp_pp, const double amp_mm, const double amp_00, 
                const double amp_0p, const double amp_m0, const double amp_p0, const double amp_0m,
                const double amp_mp, const double amp_pm):
amp_pp(amp_pp),amp_mm(amp_mm),amp_00(amp_00),
amp_0p(amp_0p),amp_m0(amp_m0),amp_p0(amp_p0),amp_0m(amp_0m),
amp_mp(amp_mp),amp_pm(amp_pm)
{;}


Deut_GPD_T_set():
amp_pp(0.),amp_mm(0.),amp_00(0.),
amp_0p(0.),amp_m0(0.),amp_p0(0.),amp_0m(0.),
amp_mp(0.),amp_pm(0.)
{;}

~Deut_GPD_T_set(){;}


double getAmp_pp() const{return amp_pp;}
double getAmp_mm() const{return amp_mm;}
double getAmp_00() const{return amp_00;}
double getAmp_0p() const{return amp_0p;}
double getAmp_m0() const{return amp_m0;}
double getAmp_p0() const{return amp_p0;}
double getAmp_0m() const{return amp_0m;}
double getAmp_mp() const{return amp_mp;}
double getAmp_pm() const{return amp_pm;}

/**
 * @brief adding operator
 * 
 * @param rhs what we want to add
 * @return TransGPD_set sum of this + argument
 */
Deut_GPD_T_set operator+(const Deut_GPD_T_set& rhs) const;
/**
 * @brief right-multiply all elements with a scalar
 * 
 * @param sc what we multiply with
 * @return TransGPD_set result of the multiplication
 */
Deut_GPD_T_set operator*(const double sc) const;


private:

double amp_pp; ///< helicity amplitude ++ (outgoing,incoming)
double amp_mm; ///< helicity amplitude --
double amp_00; ///< helicity amplitude 00
double amp_0p; ///< helicity amplitude 0+
double amp_m0; ///< helicity amplitude -0
double amp_p0; ///< helicity amplitude +0
double amp_0m; ///< helicity amplitude 0-
double amp_mp; ///< helicity amplitude -+
double amp_pm; ///< helicity amplitude +-
};


#endif //TRNSGPD_SET_HPP