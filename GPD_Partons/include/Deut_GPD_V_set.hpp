#ifndef DEUT_GPD_V_SET_HPP
#define DEUT_GPD_V_SET_HPP

/**
 * @brief small class that holds deuteron (isosinglet) vector helicity amplitudes for a certain kinematics
 * 
 */
class Deut_GPD_V_set{

public:
/**
 * @brief constructor
 * 
 */
Deut_GPD_V_set(const double amp_pp, const double amp_00, const double amp_0p, const double amp_p0, const double amp_mp):
amp_pp(amp_pp),amp_00(amp_00),amp_0p(amp_0p),amp_p0(amp_p0),amp_mp(amp_mp)
{;}


Deut_GPD_V_set():
amp_pp(0.),amp_00(0.),amp_0p(0.),amp_p0(0.),amp_mp(0.)
{;}

~Deut_GPD_V_set(){;}


double getAmp_pp() const{return amp_pp;}
double getAmp_00() const{return amp_00;}
double getAmp_0p() const{return amp_0p;}
double getAmp_p0() const{return amp_p0;}
double getAmp_mp() const{return amp_mp;}
/**
 * @brief Get amplitude with certain index
 * 
 * @param index 0 is ++, 1 is 00, 2 is 0+, 3 is +0, 4 is -+
 * @return double deuteron helicity amplitude
 */
double getAmp(int index) const;
/**
 * @brief adding operator
 * 
 * @param rhs what we want to add
 * @return TransGPD_set sum of this + argument
 */
Deut_GPD_V_set operator+(const Deut_GPD_V_set& rhs) const;
/**
 * @brief right-multiply all elements with a scalar
 * 
 * @param sc what we multiply with
 * @return TransGPD_set result of the multiplication
 */
Deut_GPD_V_set operator*(const double sc) const;

Deut_GPD_V_set Reverse_xi() const;

private:

double amp_pp; ///< helicity amplitude ++ (outgoing,incoming)
double amp_00; ///< helicity amplitude 00
double amp_0p; ///< helicity amplitude 0+
double amp_p0; ///< helicity amplitude +0
double amp_mp; ///< helicity amplitude -+
};


#endif //TRNSGPD_SET_HPP