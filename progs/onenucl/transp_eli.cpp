//program used to calculate some lead transparencies


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

double diffcross(MeanFieldNucleusThick& nuc,int level,double pm, double Q2, double omega,double Ebeam){
    int    current = 1;     // \in [1,2,3]
    int    thick   = 1;     // thick (1)
    int    src     = 0;     // src in fsi
    int    ct      = 0;     // color transparency
    int    pw      = 1;     // pw (1) or fsi (0)?

    double abserr    = 1e-5; // I am guessing this is absolute error and not relative error...
    bool user_sigma  = false;
    double screening = 0.;
    /** ========= **/
    
    
    TKinematics2to2 kin("kinematics_name","kinematics_title",nuc.getMassA(),nuc.getMassA_min_proton(),MASSP,"qsquared:wlab:pklab",Q2,omega,pm);
    TElectronKinematics* elec = TElectronKinematics::CreateWithBeamEnergy(Ebeam);
    Cross cross(*elec,&nuc,abserr,2,getShareDir(),user_sigma,screening);
    
    delete elec;
    return cross.getDiffCross(kin,current,thick,src,ct,pw,level,0.,25000,true, true); // kinematics2to2,current [1,2,3], src, CT, pw, shellindex, phi,maxEval,lab
}

void run(MeanFieldNucleusThick& nuc, double Q2,double pf){
    
    /** SETTINGS **/
    double xB    = 0.8;
    double omega = Q2/(2.*MASSP*xB);
    double q     = sqrt(Q2 + omega*omega);
    double Ebeam = 5014.; // HALLB beam energy,
    
    
    printf("#__PYHDR_START__\n");
    printf("# { \n#   \"nucleus\" : \"%s\",\n#   \"plevels\" : %d,\n#   \"nlevels\" : %d,\n#\n",nuc.getNucleusName().c_str(),nuc.getPLevels(),nuc.getNLevels());
    printf("#   \"Q2\"    : %.3e,\n#   \"pf\"    : %8.3f\n# }\n",Q2,pf);
    printf("#__PYHDR_STOP__\n");
    printf("# pm diffcross[sum protons] diffcross[sum neutrons] diffcross[levels]\n");
    for (double pm=0.;pm<600;pm+=10.){ // loop over missing momentum, pm

        double cos_q_pm = ( pf*pf - pm*pm - q*q) / ( 2.*pm*q) ;
        
        if( fabs( cos_q_pm ) > 1.){
            fprintf(stderr,"#[ERROR] cosine outside [-1,1] range: %8.3f for pm %8.3f\n",cos_q_pm,pm);
        } else {
            printf("%8.3f",pm);
            double clevels[nuc.getTotalLevels()];
	    memset(clevels,0.,nuc.getTotalLevels());
            double pleveltot = 0.;
            double nleveltot = 0.;
            for (int plevel=0;plevel<nuc.getPLevels();plevel++){ // iterate over proton shells
                clevels[plevel] = diffcross(nuc,plevel,pm,Q2,omega,Ebeam);
                pleveltot += clevels[plevel];
            }
            for (int nlevel=nuc.getPLevels();nlevel<nuc.getTotalLevels();nlevel++){ // iterate over neutron shells
                clevels[nlevel] = diffcross(nuc,nlevel,pm,Q2,omega,Ebeam);
                nleveltot += clevels[nlevel];
            }
            printf("   %.4e   %.4e  ",pleveltot,nleveltot);
            for (int level=0;level<nuc.getTotalLevels();level++){
                printf("   %.4e",clevels[level]);
            }
            printf("\n");
            fflush(stdout);
        }
    } 
}


void Q2pf_run(MeanFieldNucleusThick& nuc,int pi,int qi){
	double Q2l[3] = {1.2,2.45,3.70};
	double pfl[3] = {1400,1870,2340};
	run(nuc,Q2l[qi]*1e6,pfl[pi]);
}

// Q^2 [0,1,2] : [1.2,.2.45,3.7]
// pf  [0,1,2] : [1400,1870,2340]
int main(){
    MeanFieldNucleusThick nuc(MeanFieldNucleus::Pb,getShareDir());
    //run(nuc,1.8*1e6,2000.); // Q^2, pf, mean settings
    Q2pf_run(nuc,2,2);
    //run(nuc,1.2*1e6,1400.); // 00 setting Q^2, pf
    //run(nuc,2.45*1e6,1400.); // 01 setting q^2, pf
    //run(nuc,3.7*1e6,1400.); // 02 setting q^2, pf
    //run(nuc,1.2*1e6,1870.); // 10 setting q^2, pf
    //run(nuc,2.45*1e6,1870.); // 11 setting q^2, pf
    //run(nuc,3.70*1e6,1870.); // 12 setting q^2, pf
    //run(nuc,3.70*1e6,2340.); // 22 setting q^2, pf
    return 0;
}
