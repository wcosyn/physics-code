#include "TwoVector_Deut.hpp"
#include <constants.hpp>
#include <mstwpdf.h>
#include <cassert>
#include <numint/numint.hpp>
#include <cmath>

using namespace std;

TwoVector_Deut::TwoVector_Deut(PARTONS::GPDService* pGPDService,
                                 PARTONS::GPDModule* pGPDModel, 
                                 PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule, string wfname, string pdfname, unsigned int id):
pRunningAlphaStrongModule(pRunningAlphaStrongModule),deut_vector_grid(pGPDService,pGPDModel,wfname,id),
deut_tensor_grid(pdfname, wfname){
    ;
}

void TwoVector_Deut::integrandum_omega_general(numint::vector_z & result, double u, double z, double xi, double mandelstam_t, double scale, 
                                      double psq, double Qsq, int helampindex, TwoVector_Deut::Omega_pol omegapol,
                                      TwoVector_Deut::Photon_pol gamma,  TwoVector_Deut& twovector){

    result = vector<complex <double> >(6,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) for(int i=0;i<6;i++) result[i]=0.;
    else{ 
        if(omegapol==TwoVector_Deut::komegaL){
            result[0]= (twovector.deut_vector_grid.getDeut_GPD_V_set(-xi*(2.*u-1),xi,mandelstam_t,scale,1, 50).getAmp(helampindex)) 
            *sqrt(2.);  //why sqrt 2??
            for(int i =1; i<6;i++) result[i] = result[0];
        }
        if(omegapol==TwoVector_Deut::komegaT){
            result[0]= twovector.deut_tensor_grid.getDeut_GPD_T_set(-xi*(2.*u-1),xi,mandelstam_t, scale, 1,1,50).getAmp(helampindex)*sqrt(2.)*2.;  //check factors
            for(int i =1; i<6;i++) result[i] = result[0];
        }

        // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [6z*barz,6u*baru from 2 DA taken into account]
        double integrand = (gamma == TwoVector_Deut::kgammaL?  36. //sqrt factor from helamps accounted for in prefactor 
                    *z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) )
                        : 36.   //sqrt factor from helamps accounted for in prefactor 
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) ));
        for(int i=0;i<6;i++) result[i]*=integrand;                
        for(int i=0;i<3;i++){
            result[3+i]*=64./PI/PI/36./sqrt(z*(1-z)*u*(1-u));  //AdS DA
        }
     }
    return;
}


// integrandum for gamma_L amplitude
void TwoVector_Deut::integrandum_gammaL_rhoL(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq,
                                 int helampindex, TwoVector_Deut &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 

        // Amp^Isoscalar *sqrt(2) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        result[0]= (twovector.deut_vector_grid.getDeut_GPD_V_set(-xi*(2.*u-1),xi,mandelstam_t,scale,1, 50).getAmp(helampindex)) 
        *sqrt(2.)
                    *z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
    //cout << u << " " << z << " " << result[0] << endl;
    return;

}

// integrandum for gamma_T amplitude
void TwoVector_Deut::integrandum_gammaT_rhoL(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq, 
                                    int helampindex, TwoVector_Deut &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        

        // Amp^Isoscalar *sqrt(2) * (2z-1) * z * barz / u / baru * Q    [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards] [factor p_perp from nominator taken out]
        result[0]=(twovector.deut_vector_grid.getDeut_GPD_V_set(-xi*(2.*u-1),xi,mandelstam_t,scale,1,50).getAmp(helampindex))
        *sqrt(2.)
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
     return;

}


// integrandum for gamma_L amplitude
void TwoVector_Deut::integrandum_gammaL_rhoT(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq,
                                 int helampindex, TwoVector_Deut &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        // (Amp^isoscalar * sqrt(2) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        result[0]=twovector.deut_tensor_grid.getDeut_GPD_T_set(-xi*(2.*u-1),xi,mandelstam_t, scale, 1,1,50).getAmp(helampindex)*sqrt(2.)*2.* //helamp for T is w a factor 0.5
                    z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
    //cout << u << " " << z << " " << result[0] << endl;
    return;

}

// integrandum for gamma_T amplitude
void TwoVector_Deut::integrandum_gammaT_rhoT(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq, 
                                    int helampindex,  TwoVector_Deut &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 

        // (Amp^Isoscalar *sqrt(2) * (2z-1) * z * barz / u / baru * Q    [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards] [factor p_perp from nominator taken out]
        result[0]=twovector.deut_tensor_grid.getDeut_GPD_T_set(-xi*(2.*u-1),xi,mandelstam_t,scale, 1,1,50).getAmp(helampindex)*sqrt(2.)*2.* //helamp for T is w a factor 0.5
                    (2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
     return;

}


double TwoVector_Deut::getCross_gammaL_rhoL(const double scale, const double xi, const double Q2, const double psq){
    double frhoplus = 0.216; //omega decay constant [GeV]
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi);//tmin [GeV^2]
    //mandelstam_t=floor(mandelstam_t*100.)/100.;
    mandelstam_t*=1.1;

    Ftor_2vector F;
    F.xi=xi;
    F.mandelstam_t=mandelstam_t;
    F.scale=scale;
    F.Qsq=Q2;
    F.psq=psq;
    F.pobj=this;

    numint::mdfunction<numint::vector_d,2> mdf;
    mdf.func = &Ftor_2vector::exec;
    mdf.param = &F;

    numint::array<double,2> lower = {{0.,0.}};
    numint::array<double,2> upper = {{0.5,1.}};
    std::vector<double> integral(1,0.); //integrandum result array

    //gamma_L calculation
    F.f=integrandum_gammaL_rhoL;
    double result_L=0.;
    //cout << 16.*PI*PI*alpha_s*frhoplus*xi*CF/Nc/psq/psq << endl;        
    //cout << 2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
    for(int helampindex=0; helampindex<=4;helampindex++){
        F.helampindex=helampindex;
        unsigned count=0;
        // integration over u and z
        numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
        integral[0]*=36.; //prefactor DA (twice, squared)
        //cout << helampindex << " " << integral[0] << endl;
        
        //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection
        integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*CF/Nc/psq/psq; 
        //prefactors Eq (15) Enberg et al.
        integral[0] *= -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); 
        //cout << 1./256./pow(PI,3.)/xi/(1.+xi) << endl;
        result_L+=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.+xi)*(helampindex==1?1:2); // Enberg et al Eq (12) + symmetry factor helicity amplitudes, 00 one only once
        //cout << Qsq << " " << xi << " " << psq << " " << result_L*0.389379E06 << endl;  // factor to go from GeV-6 to nb GeV-4
        // cout << helampindex << " " << 16.*PI*PI*alpha_s*frhoplus*xi*CF/Nc/psq/psq << " " << 
        // -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << " " << pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.+xi)*(helampindex==1?1:2) << endl;
    }
    
    //cout << "final " << result_L*0.389379E06/3. << endl;
    //exit(1);
    return result_L*0.389379E06/3.*2.;  //[Gev-2 -> nb conversion] + avg initial spin + factor 2 because omega is ~ [(u+d)/sqrt(2)]^2 = 2u for isoscalar

}

double TwoVector_Deut::getCross_gammaT_rhoL(const double scale, const double xi, const double Q2, const double psq){
    double frhoplus = 0.216; //omega decay constant [GeV]
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi);//tmin [GeV^2]
    //mandelstam_t=floor(mandelstam_t*100.)/100.;
    mandelstam_t*=1.1;



    Ftor_2vector F;
    F.xi=xi;
    F.mandelstam_t=mandelstam_t;
    F.scale=scale;
    F.Qsq=Q2;
    F.psq=psq;
    F.pobj=this;

    numint::mdfunction<numint::vector_d,2> mdf;
    mdf.func = &Ftor_2vector::exec;
    mdf.param = &F;

    numint::array<double,2> lower = {{0.,0.}};
    numint::array<double,2> upper = {{0.5,1.}};
    std::vector<double> integral(1,0.); //integrandum result array

    //gamma_L calculation
    F.f=integrandum_gammaT_rhoL;
    double result_T=0.;
    for(int helamps=0;helamps<=4;helamps++){
        F.helampindex=helamps;
        unsigned count=0;
        // integration over u and z
        numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
        integral[0]*=36.; //prefactor DA (twice, squared)

        integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al.
        integral[0] *= PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.

        result_T+=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.+xi)*(helamps==1?1:2); // Enberg et al Eq (12) + helamp symmetry factor, 00 only once
    }
    return result_T*0.389379E06/3*2; //[Gev-2 -> nb conversion] + avg initial spin + factor 2 because omega is ~ [(u+d)/sqrt(2)]^2 = 2u for isoscalar

}

double TwoVector_Deut::getCross_gammaL_rhoT(const double scale, const double xi, const double Q2, const double psq){
    double frhoplus = 0.216; //rho+ decay constant [GeV]
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi);//tmin [GeV^2]
   //mandelstam_t=floor(mandelstam_t*100.)/100.;
    mandelstam_t*=1.1;

    Ftor_2vector F;
    F.xi=xi;
    F.mandelstam_t=mandelstam_t;
    F.scale=scale;
    F.Qsq=Q2;
    F.psq=psq;
    F.pobj=this;

    numint::mdfunction<numint::vector_d,2> mdf;
    mdf.func = &Ftor_2vector::exec;
    mdf.param = &F;

    numint::array<double,2> lower = {{0.,0.}};
    numint::array<double,2> upper = {{0.5,1.}};

    //gamma_L calculation
    F.f=integrandum_gammaL_rhoT;
    unsigned count=0;
    double result_L=0.;
    // integration over u and z
    for(int helamps=0;helamps<=8;helamps++){
        F.helampindex=helamps;
        std::vector<double> integral(1,0.); //integrandum result array
        numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
        
        integral[0]*=36.; //prefactor DA (twice, squared)
        
        //cout << helamps << " " << integral[0] << endl;
        
        //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection
        integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*CF/Nc/psq/psq; 
        //cout << 16.*PI*PI*alpha_s*frhoplus*xi*CF/Nc/psq/psq << endl;
        //prefactors Eq (15) Enberg et al.
        integral[0] *= -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); 
        //cout << -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << " " << frho0 << endl;
        //cout << 1./256./pow(PI,3.)/xi/(1.+xi) << endl;
        result_L+=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.+xi); // Enberg et al Eq (12)
        // cout << helamps << " " << 16.*PI*PI*alpha_s*frhoplus*xi*CF/Nc/psq/psq << " " << -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << " " <<
        // pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.+xi) << endl;

        
    }
    //cout << Qsq << " " << xi << " " << psq << " " << result_L*0.389379E06 << endl;  // factor to go from GeV-6 to nb GeV-4

    return result_L*0.389379E06/3.*2;  //[Gev-2 -> nb conversion] + avg initial spin+ factor 2 because omega is ~ [(u+d)/sqrt(2)]^2 = 2u for isoscalar

}


double TwoVector_Deut::getCross_gammaT_rhoT(const double scale, const double xi, const double Q2, const double psq){
    double frhoplus = 0.216; //rho+ decay constant [GeV]
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi);//tmin [GeV^2]
//    mandelstam_t=floor(mandelstam_t*100.)/100.;
    mandelstam_t*=1.1;


    Ftor_2vector F;
    F.xi=xi;
    F.mandelstam_t=mandelstam_t;
    F.scale=scale;
    F.Qsq=Q2;
    F.psq=psq;
    F.pobj=this;
    
    numint::mdfunction<numint::vector_d,2> mdf;
    mdf.func = &Ftor_2vector::exec;
    mdf.param = &F;

    numint::array<double,2> lower = {{0.,0.}};
    numint::array<double,2> upper = {{0.5,1.}};

    //gamma_T calculation
    F.f=integrandum_gammaT_rhoT;
    unsigned count=0;
    double result_T=0.;
    // integration over u and z
    for(int helamps=0;helamps<=8;helamps++){
        F.helampindex=helamps;
        std::vector<double> integral(1,0.); //integrandum result array
        numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
        integral[0]*=36.; //prefactor DA (twice, squared)
        //cout << helamps << " " << integral[0] << endl;

        integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al.
        integral[0] *= PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.
        // cout << PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
        result_T+=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.+xi); // Enberg et al Eq (12)

        
    }


    return result_T*0.389379E06/3.*2.; //[Gev-2 -> nb conversion] + avg initial spin+ factor 2 because omega is ~ [(u+d)/sqrt(2)]^2 = 2u for isoscalar

}



void TwoVector_Deut::getCross_twovector(std::vector<double> & results, const double scale, const double xi, const double Q2, const double psq, 
                          TwoVector_Deut::Photon_pol gammapol, TwoVector_Deut::Omega_pol omegapol){

    double f_lowermeson =  0.;  // meson decay constant for meson originating from lower part of the graph
    int helampmax= 0;
    if (omegapol == TwoVector_Deut::komegaT) {
        f_lowermeson=0.216;
        helampmax=8;
    }
    if (omegapol == TwoVector_Deut::komegaL) {
        f_lowermeson=0.216;
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
    F.pobj=this;
    
    numint::mdfunction<numint::vector_z,2> mdf;
    mdf.func = &Ftor_2vector_general::exec;
    mdf.param = &F;

    numint::array<double,2> lower = {{0.,0.}};
    numint::array<double,2> upper = {{0.5,1.}};  // symmetry in u exploited to simplify integrand!

    //gamma_T calculation
    F.f=integrandum_omega_general;
    unsigned count=0;
    results = vector<double>(12,0.);
    // integration over u and z
    for(int helamps=0;helamps<=helampmax;helamps++){
        F.helampindex=helamps;

        std::vector<complex<double> > integral(6,0.); //integrandum result array
        numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
        for(int index=0;index<6;index++){
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

    //[Gev-2 -> nb conversion] + 1/3 avg initial spin + factor 2 because omega is ~ [(u+d)/sqrt(2)]^2 = 2u for isoscalar
    for(int index=0;index<6;index++) results[index]*=0.389379E06/3*2.; 
    return;

}
