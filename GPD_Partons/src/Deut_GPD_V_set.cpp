#include "Deut_GPD_V_set.hpp"
#include <cassert>
#include <iostream>

using namespace std;



Deut_GPD_V_set Deut_GPD_V_set::operator+(const Deut_GPD_V_set& rhs) const{
    Deut_GPD_V_set result;
    result.amp_pp=this->amp_pp+rhs.amp_pp;
    result.amp_00=this->amp_00+rhs.amp_00;
    result.amp_0p=this->amp_0p+rhs.amp_0p;
    result.amp_p0=this->amp_p0+rhs.amp_p0;
    result.amp_mp=this->amp_mp+rhs.amp_mp;

    return result;
}

Deut_GPD_V_set Deut_GPD_V_set::operator*(const double sc) const{
    Deut_GPD_V_set result;
    result.amp_pp=this->amp_pp*sc;
    result.amp_00=this->amp_00*sc;
    result.amp_0p=this->amp_0p*sc;
    result.amp_p0=this->amp_p0*sc;
    result.amp_mp=this->amp_mp*sc;

    return result;
}

Deut_GPD_V_set Deut_GPD_V_set::Reverse_xi() const{
    return Deut_GPD_V_set(this->amp_pp, this->amp_00, this->amp_p0, this->amp_0p, this->amp_mp);
}

double Deut_GPD_V_set::getAmp(int index) const{
    switch(index){
        case 0:
            return amp_pp;
            break;
        case 1:
            return amp_00;
            break;
        case 2:
            return amp_0p;
            break;
        case 3:
            return amp_p0;
            break;
        case 4:
            return amp_mp;
            break;
        default:
            cerr << "invalid index" << endl;
            return 0.;
            break;
    }

}