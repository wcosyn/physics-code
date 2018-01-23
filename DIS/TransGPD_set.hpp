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

~TransGPD_set(){}

double getHT_singlet() const{return 0.5*(HTu+HTd);} ///< returns (H_T^u+H_T^d)/2
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
 * @brief returns \tilde{E}_T, put to zero here
 * 
 * @return double getEtildeT_singlet 
 */
double getEtildeT_singlet() const{ return 0.;}

double getHTu() const{return HTu;}
double getHTd() const{return HTd;}
double getEbarTu() const{return EbarTu;}
double getEbarTd() const{return EbarTd;}

private:

double HTu; ///< GPD H_T up quark
double HTd; ///< GPD H_T down quark
double EbarTu; ///< GPD \bar{E}_T=2\tilde{H}_T+E_T up quark
double EbarTd; ///< GPD \bar{E}_T=2\tilde{H}_T+E_T down quark

};


#endif //TRNSGPD_SET_HPP