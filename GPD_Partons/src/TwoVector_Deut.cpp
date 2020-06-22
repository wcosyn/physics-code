#include "TwoVector_Deut.hpp"
#include <constants.hpp>
#include <mstwpdf.h>
#include <cassert>
#include <numint/numint.hpp>
#include <cmath>

using namespace std;

TwoVector_Deut::TwoVector_Deut(PARTONS::GPDService* pGPDService,
                                 PARTONS::GPDModule* pGPDModel, 
                                 PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule, string wfname, string pdfname, 
                                 unsigned int id):
pRunningAlphaStrongModule(pRunningAlphaStrongModule),deut_vector_grid(pGPDService,pGPDModel,wfname,id),
deut_tensor_grid(pdfname, wfname){
    ;
}

void TwoVector_Deut::integrandum_omega_general(numint::vector_z & result, const double u, const double z, double xi, double mandelstam_t, double scale, 
                                      double psq, double Qsq, int helampindex, TwoVector_Deut::Omega_pol omegapol,
                                      TwoVector_Deut::Photon_pol gamma,  int T_model, TwoVector_Deut& twovector){

    result = vector<complex <double> >(2,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) for(int i=0;i<2;i++) result[i]=0.;
    else{ 
        if(omegapol==TwoVector_Deut::komegaL){
            //only gpd u but isoscalar so d is equal to it (taken care of at end)
            result[0]= twovector.deut_vector_grid.getDeut_GPD_V_set(-xi*(2.*u-1),xi,mandelstam_t,scale,1, 50).getAmp(helampindex); 
            for(int i =1; i<2;i++) result[i] = result[0];
        }
        if(omegapol==TwoVector_Deut::komegaT){
            //only gpd u but isoscalar so d is equal to it (taken care of at end); grid contains i\sigma^+R matrix elements
            result[0]= twovector.deut_tensor_grid.getDeut_GPD_T_set(-xi*(2.*u-1),xi,mandelstam_t, scale, 1,T_model,50).getAmp(helampindex);
            for(int i =1; i<2;i++) result[i] = result[0];
        }

        // (Hq^u) * z^2 * barz^2 / u / baru * P [6z*barz,6u*baru from 2 DA taken into account]
        double integrand = (gamma == TwoVector_Deut::kgammaL?  36. 
                    *z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) )
                        : 36.  
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) ));
        for(int i=0;i<2;i++) result[i]*=integrand;                
            result[1]*=64./PI/PI/36./sqrt(z*(1-z)*u*(1-u));  //AdS DA
        
     }
    return;
}



void TwoVector_Deut::getCross_twovector(std::vector<double> & results, const double scale, const double xi, const double Q2, const double psq, 
                          TwoVector_Deut::Photon_pol gammapol, TwoVector_Deut::Omega_pol omegapol, int T_model){

    double f_lowermeson =  0.;  // meson decay constant for meson originating from lower part of the graph
    int helampmax= 0;
    if (omegapol == TwoVector_Deut::komegaT) {
        f_lowermeson=0.160;
        helampmax=8;
    }
    if (omegapol == TwoVector_Deut::komegaL) {
        f_lowermeson=0.197;
        helampmax=4;
    }
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi);//tmin [GeV^2]
//    mandelstam_t=floor(mandelstam_t*100.)/100.;
    mandelstam_t*=1.1;


    Ftor_2vector_general F;
    F.xi=xi;
    F.mandelstam_t=mandelstam_t;
    F.scale=scale;
    F.Qsq=Q2;
    F.psq=psq;
    F.omegapol = omegapol;
    F.gammapol = gammapol;
    F.T_model = T_model;
    F.pobj=this;
    
    numint::mdfunction<numint::vector_z,2> mdf;
    mdf.func = &Ftor_2vector_general::exec;
    mdf.param = &F;

    numint::array<double,2> lower = {{0.,0.}};
    numint::array<double,2> upper = {{0.5,1.}};  // symmetry in u exploited to simplify integrand!

    //gamma_T calculation
    F.f=integrandum_omega_general;
    unsigned count=0;
    results = vector<double>(2,0.);
    // integration over u and z
    for(int helamps=0;helamps<=helampmax;helamps++){
        F.helampindex=helamps;

        std::vector<complex<double> > integral(2,0.); //integrandum result array
        numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
        for(int index=0;index<2;index++){
            //cout << spinin << " " << spinout << " " << integral[index] << endl;
        
            //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection
            // sqrt(1-xi^2) compensated  in the helicity amplitudes
            integral[index] *= 16.*PI*PI*alpha_s*f_lowermeson*xi*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al.
            if(gammapol==TwoVector_Deut::kgammaL) integral[index] *= -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.
            if(gammapol==TwoVector_Deut::kgammaT) integral[index] *= PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.
            results[index]+=norm(integral[index])/256./pow(PI,3.)/xi/(1+xi)*((helampmax==4&&helamps!=1)?2:1); // Enberg et al Eq (12), corrected phase space factor!
            // if(index==0) cout << "new " << helamps << " " << 16.*PI*PI*alpha_s*f_lowermeson*xi*CF/Nc/psq/psq << " " << -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI)
            // << " " << norm(integral[index])/256./pow(PI,3.)/xi/(1+xi)*((helampmax==4&&helamps!=1)?2:1) << " " << ((helampmax==4&&helamps!=1)?2:1) << endl;
        }
    }

    //[Gev-2 -> nb conversion] + 1/3 avg initial spin + factor 2 because omega is ~ [(u+d)/sqrt(2)]^2 = 2u^2 for isoscalar
    //+ extra factor 1/2 from averaging over sin^2 theta in case of transverse vector meson (or from (RL+LR)/2 products)
    for(int index=0;index<2;index++) results[index]*=0.389379E06/3*2./(omegapol==TwoVector_Deut::komegaT? 2.:1.); 
    return;

}
