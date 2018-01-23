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