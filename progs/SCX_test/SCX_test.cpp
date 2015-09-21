#include "ClassGridThick_SCX.hpp"
#include "FastParticle.hpp"
#include "MeanFieldNucleusThick.hpp"

#ifndef SHAREDIR
#define SHAREDIR "/home/camel/Code/trunk/share"
#endif


/** exmaple density from the MeanFieldNucleusThick class 
 *  for actual SCX calculation you would want to use either
 *  getProtonDensity or getNeutronDensity.
 *  Always make sure you are using neutron densities if fastparticle is proton
 *  and proton densities if fastparticle is neut
 * **/
double meanFieldNucleusThickNeutronDensity(double r,void* param){
    MeanFieldNucleusThick* nuc = (MeanFieldNucleusThick*) param;
    if (r < nuc->getWF_r_step()) // going closer to 0 than this causes divergencies!
        r = nuc->getWF_r_step();
    else if ( r > nuc->getRange())
        return 0.; // outside RMAX density is for virtually zero
    return nuc->getNeutronDensity(r)/r/r; // note the division by r^2 to compensate r^2 already included in density function
}

/** custom density **/
double myDensity(double r,void* param){
    return ( r < 1. ) ? 1. : 0. ;
}

void testMeanFieldNucleus(){
    FastParticle fp(FastParticle::Particletype::P_CLASS_SCX,0,100.,0.,0.,0.,0.,SHAREDIR); // change elastic scattering nucleon (0,1) to SCX nucleon (8,9) 
    MeanFieldNucleusThick nuc(MeanFieldNucleus::C,SHAREDIR);
    char densName[256];
    sprintf(densName,"%s_neutronDensity",nuc.getNucleusName().c_str());
    ClassGridThick_SCX scx(meanFieldNucleusThickNeutronDensity,(void*)&nuc,densName,nuc.getRange(),&fp,100,100,SHAREDIR);
    scx.set_dens_fctr(6.);
    scx.constructGrid();
    for (double r=0.;r<nuc.getRange();r+=0.1){
        for (double th=0.;th<=M_PI+1e-10; th += 0.01*M_PI){ // +1e-10 for roundoff error preventing to go to exactly M_PI
            double phi=0.;
            double x[3] = {r,cos(th),phi};
            double res,err;
            scx.getInterp(x,res,err);
            std::cout << r*sin(th) << "\t" << r*cos(th) << "\t" << res << "\t" << err << std::endl;
        }
    }
}

void testCustomDensity(){
    FastParticle fp(FastParticle::Particletype::P_CLASS_SCX,0,1000.,0.,0.,0.,0.,SHAREDIR); // change elastic scattering nucleon (0,1) to SCX nucleon (8,9) 
    ClassGridThick_SCX scx(myDensity,NULL,"dummyDensity",2.0,&fp,100,120,SHAREDIR);
    scx.constructGrid();
    scx.printGrid();
}

int main(){
    testMeanFieldNucleus();
    return 0;
}
