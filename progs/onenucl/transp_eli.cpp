
#include <MeanFieldNucleusThick.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <Cross.hpp>
#include <FastParticle.hpp>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cassert>
#include <cstdio>

std::string getShareDir(){
    std::stringstream ss;
    ss << std::getenv("HOME") << "/Code/trunk/share";
    return ss.str();
}

double diffcross(double pm,MeanFieldNucleusThick& nuc,int level){
    /** SETTINGS **/
    double Q2    = 1.8*1e6; // GeV^2
    double xB    = 0.8;
    double omega = Q2/(2.*MASSP*xB);
    double q     = sqrt(Q2 + omega*omega);
    double pf    = 2000.; // momentum of final particle
    double Ebeam = 5014.; // HALLB beam energy,
    int    pw    = 0;     // pw (1) or fsi (0)?
    double abserr    = 1e-5; // I am guessing this is absolute error and not relative error...
    double screening = 0.;
    double src       = 0.;
    /** ========= **/

    double cos_q_pm = ( pf*pf - pm*pm - q*q) / ( 2.*pm*q) ;
    if( fabs( cos_q_pm ) > 1.){
        fprintf(stderr,"#[ERROR] cosine outside [-1,1] range: %8.3f for pm %8.3f\n",cos_q_pm,pm);
        return -999.;
    }

    TKinematics2to2 kin("kinematics_name","kinematics_title",nuc.getMassA(),nuc.getMassA_min_proton(),MASSP,"qsquared:wlab:pklab",Q2,omega,pm);
    TElectronKinematics* elec = TElectronKinematics::CreateWithBeamEnergy(Ebeam);
    Cross cross(*elec,&nuc,abserr,2,getShareDir(),screening,src);
    
    delete elec;
    return cross.getDiffCross(kin,1,0,0,0,pw,level,0.,1000000,1); // kinematics2to2,current [1,2,3], src, CT, pw, shellindex, phi,maxEval,lab
}

int main(){
    MeanFieldNucleusThick nuc(MeanFieldNucleus::C,getShareDir());
    printf("# Nucleus: %s, proton levels: %d, neutron levels: %d \n",nuc.getNucleusName().c_str(),nuc.getPLevels(),nuc.getNLevels());
    printf("# pm diffcross[levels]\n");
    for (double pm=201.;pm<300;pm+=1.){ // loop over missing momentum, pm
        printf("%8.3f",pm);
        for (int level=0;level<nuc.getTotalLevels(); level++){
            double c = diffcross(pm,nuc,level);
            printf("    %.6e",c);
        }
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
