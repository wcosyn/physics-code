#include "TransGPD_set.hpp"
#include <cassert>
#include <iostream>

using namespace std;

double TransGPD_set::getHtildeT_singlet(const int model) const{

    switch(model){
        case 0:
            return 0.;
            break;
        case 1:
            return getHT_singlet();
            break;
        case 2:
            return -getHT_singlet();
            break;
        default:
            cerr << "invalid choice in TransGPD_set::getHtildeT_singlet" << endl;
            assert(1==0);
            return 0.;
            break;
    }

}

double TransGPD_set::getET_singlet(const int model) const{

    switch(model){
        case 0:
            return 0.5*(EbarTd+EbarTu);
            break;
        case 1:
            return 0.5*(EbarTd+EbarTu)-2.*getHT_singlet();
            break;
        case 2:
            return 0.5*(EbarTd+EbarTu)+2.*getHT_singlet();
            break;
        default:
            cerr << "invalid choice in TransGPD_set::getET_singlet" << endl;
            assert(1==0);
            return 0.;
            break;
    }
    
}


TransGPD_set TransGPD_set::operator+(const TransGPD_set& rhs) const{
    TransGPD_set result;
    result.HTu=this->HTu+rhs.HTu;
    result.HTd=this->HTd+rhs.HTd;
    result.EbarTu=this->EbarTu+rhs.EbarTu;
    result.EbarTd=this->EbarTd+rhs.EbarTd;

    return result;
}

TransGPD_set TransGPD_set::operator*(const double sc) const{
    TransGPD_set result;
    result.HTu=this->HTu*sc;
    result.HTd=this->HTd*sc;
    result.EbarTu=this->EbarTu*sc;
    result.EbarTd=this->EbarTd*sc;

    return result;
}
