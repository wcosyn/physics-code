#ifndef TRANSGPD_SET_HPP
#define TRANSGPD_SET_HPP

/**
 * @brief small class that holds the nucleon transversity gpds for a certain kinematics and has methods to get certain combinations of them
 * 
 */
class TransGPD_set{

public:
/**
 * @brief constructor
 * 
 * @param HTu chiral odd GPD H_T for up quark
 * @param HTd chiral odd GPD H_T for down quark
 * @param EbarTu chiral odd GPD \bar{E}_T=2\tilde{H}_T+E_T for up quark
 * @param EbarTd chiral odd GPD \bar{E}_T=2\tilde{H}_T+E_T for down quark
 */
TransGPD_set(const double HTd, const double HTu, const double EbarTd, const double EbarTu):
HTd(HTd),HTu(HTu),EbarTd(EbarTd),EbarTu(EbarTu)
{}


TransGPD_set():
HTd(0.),HTu(0.),EbarTd(0.),EbarTu(0.)
{}

~TransGPD_set(){}

double getHT_singlet() const{return 0.5*(HTu+HTd);} ///< returns (H_T^u+H_T^d)/2
double getHT_vector() const{return 0.5*(HTu-HTd);} ///< returns (H_T^u-H_T^d)/2
/**
 * @brief returns isosinglet combination of \tilde{H}_T
 * 
 * @param model different implementations discussed in Phys. Rev. D 95, 094001 p3, top right column
 * [0]: \tilde{H}_T=0
 * [1]: \tilde{H}_T=H_T
 * [2]: \tilde{H}_T=-H_T
 * @return double getHtildeT_singlet 
 */
double getHtildeT_singlet(const int model) const;
/**
 * @brief returns isovector combination of \tilde{H}_T
 * 
 * @param model different implementations discussed in Phys. Rev. D 95, 094001 p3, top right column
 * [0]: \tilde{H}_T=0
 * [1]: \tilde{H}_T=H_T
 * [2]: \tilde{H}_T=-H_T
 * @return double getHtildeT_singlet 
 */
double getHtildeT_vector(const int model) const;
/**
 * @brief returns isosinglet combination of E_T
 * 
 * @param model different implementations discussed in Phys. Rev. D 95, 094001 p3, top right column
 * [0]: E_T=\bar{E}_T
 * [1]: E_T=\bar{E}_T-2H_T
 * [2]: E_T=\bar{E}_T+2H_T
 * @return double getET_singlet 
 */
double getET_singlet(const int model) const;
/**
 * @brief returns isovector combination of E_T
 * 
 * @param model different implementations discussed in Phys. Rev. D 95, 094001 p3, top right column
 * [0]: E_T=\bar{E}_T
 * [1]: E_T=\bar{E}_T-2H_T
 * [2]: E_T=\bar{E}_T+2H_T
 * @return double getET_singlet 
 */
double getET_vector(const int model) const;
/**
 * @brief returns \tilde{E}_T, put to zero here
 * 
 * @return 0.
 */
double getEtildeT_singlet() const{ return 0.;}
/**
 * @brief returns \tilde{E}_T, put to zero here
 * 
 * @return 0.
 */
double getEtildeT_vector() const{ return 0.;}

double getHTu() const{return HTu;}
double getHTd() const{return HTd;}
double getEbarTu() const{return EbarTu;}
double getEbarTd() const{return EbarTd;}

/**
 * @brief adding operator
 * 
 * @param rhs what we want to add
 * @return TransGPD_set sum of this + argument
 */
TransGPD_set operator+(const TransGPD_set& rhs) const;
/**
 * @brief right-multiply all elements with a scalar
 * 
 * @param sc what we multiply with
 * @return TransGPD_set result of the multiplication
 */
TransGPD_set operator*(const double sc) const;


private:

double HTd; ///< GPD H_T down quark
double HTu; ///< GPD H_T up quark
double EbarTd; ///< GPD \bar{E}_T=2\tilde{H}_T+E_T down quark
double EbarTu; ///< GPD \bar{E}_T=2\tilde{H}_T+E_T up quark

};


#endif //TRNSGPD_SET_HPP