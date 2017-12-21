#include <vector>
#include <constants.hpp>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>

/**
 * @brief conversion from helicity amplitudes to GPDs for spin 1 chiral odd quark GPDs.  We compute in a frame where phi=0
 * 
 * @param xi 
 * @param x 
 * @param t 
 * @param helamps 
 * @return std::vector<double> 
 */
std::vector<double> helamps_to_gpds(const double xi, const double x, const double t, const std::vector<double> & helamps){
    
    std::vector<double> gpds(9,0.);

    double D=sqrt((-4.*MASSn*MASSn*xi*xi/(1-xi*xi)-t)*(1-xi*xi))/(2.*MASSn);
    gpds.at(0)=2.*sqrt(2.)*xi/D/(1-xi*xi)*(helamps[0]-helamps[1])+2.*(1./(1.-xi)*helamps[3]+1./(1.+xi)*helamps[4]+2./(1+xi)*helamps[5]
                +2./(1-xi)*helamps[6]+2.*sqrt(2.)*D/(1.-xi*xi)*(helamps[7]-helamps[8]));
    gpds.at(1)=2.*sqrt(2.)*D/(1-xi*xi)*(helamps[0]-helamps[1])-2.*(1./(1.-xi)*helamps[3]-1./(1.+xi)*helamps[4]+2./(1+xi)*helamps[5]
                -2./(1-xi)*helamps[6]+2.*sqrt(2.)*xi/D/(1.-xi*xi)*(helamps[7]-helamps[8]));
    gpds.at(2)=1./(2.*D)*((1-xi)/(1+xi)*helamps[0]-(1+xi)/(1-xi)*helamps[1])+1./(sqrt(2.)*D*D)*((1-xi)/(1+xi)*helamps[5]-(1+xi)/(1-xi)*helamps[6])
                -2.*xi/(D*D*D*(1-xi*xi))*helamps[8];
    gpds.at(3)=1/D*(1./(1-xi)*helamps[0]+1./(1.-xi)*helamps[1])+1/(sqrt(2.)*D*D)*((1-xi)/(1+xi)*helamps[5]+(1+xi)/(1-xi)*helamps[6])+1./(2.*D)*helamps[7]
                -1./(2.*D*D*D)*(D*D-4.*xi/(1-xi*xi))*helamps[8];

    gpds.at(4)=-1./D*((1./2.+D/(1-xi))/(1+xi)/(1+xi)*helamps[0]+(1./2.+D/(1+xi))/(1-xi)/(1-xi)*helamps[1])+1./(D*(1.-xi*xi))*helamps[2]
                    +1/sqrt(2.)/(1-xi*xi)*(helamps[3]+helamps[4])-1./(sqrt(2.)*pow(1+xi,2.))*(1.-2.*xi/(D*D*(1-xi)))*helamps[5]
                    -1/(sqrt(2.)*pow(1-xi,2.))*(1.+2.*xi/(D*D*(1.+xi)))*helamps[6]-1./(2.*D*(1-xi*xi))*helamps[7]
                    +(D*D*(3.*xi+1-4.*xi*xi))/(2.*D*D*D*pow(1-xi*xi,2.))*helamps[8];
    
    gpds.at(5)=-1./(2.*D)*(helamps[0]+helamps[1]+helamps[7]+helamps[8]);
    gpds.at(6)=2.*(1-xi*xi)/(D*D*D)*helamps[8];
    gpds.at(7)=-1./D*(helamps[0]+helamps[1]+xi*(helamps[7]-helamps[8]));
    gpds.at(8)=-1./D*(helamps[0]+helamps[1]+(helamps[7]-helamps[8]));

    return gpds;
}