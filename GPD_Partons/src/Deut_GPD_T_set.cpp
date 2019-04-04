#include "Deut_GPD_T_set.hpp"
#include <cassert>
#include <iostream>

using namespace std;



Deut_GPD_T_set Deut_GPD_T_set::operator+(const Deut_GPD_T_set& rhs) const{
    Deut_GPD_T_set result;
    result.amp_pp=this->amp_pp+rhs.amp_pp;
    result.amp_mm=this->amp_pp+rhs.amp_mm;
    result.amp_00=this->amp_00+rhs.amp_00;
    result.amp_0p=this->amp_0p+rhs.amp_0p;
    result.amp_m0=this->amp_p0+rhs.amp_m0;
    result.amp_p0=this->amp_0p+rhs.amp_p0;
    result.amp_0m=this->amp_p0+rhs.amp_0m;
    result.amp_mp=this->amp_mp+rhs.amp_mp;
    result.amp_pm=this->amp_mp+rhs.amp_pm;

    return result;
}

Deut_GPD_T_set Deut_GPD_T_set::operator*(const double sc) const{
    Deut_GPD_T_set result;
    result.amp_pp=this->amp_pp*sc;
    result.amp_mm=this->amp_mm*sc;
    result.amp_00=this->amp_00*sc;
    result.amp_0p=this->amp_0p*sc;
    result.amp_m0=this->amp_m0*sc;
    result.amp_p0=this->amp_p0*sc;
    result.amp_0m=this->amp_0m*sc;
    result.amp_mp=this->amp_mp*sc;
    result.amp_pm=this->amp_pm*sc;

    return result;
}
