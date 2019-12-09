#include "TwoVector_Nucl.hpp"
#include <constants.hpp>
#include <mstwpdf.h>
#include <cassert>
#include <numint/numint.hpp>
#include <cmath>

#include "Deut_Conv_GPD_T.hpp"
using namespace std;

TwoVector_Nucl::TwoVector_Nucl(PARTONS::GPDService* pGPDService,
                                 PARTONS::GPDModule* pGPDModel, 
                                 PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule):
pGPDService(pGPDService),pGPDModel(pGPDModel),pRunningAlphaStrongModule(pRunningAlphaStrongModule),gpdTgrid("MSTW"){
    ;
}



// integrandum for gamma_L amplitude
void TwoVector_Nucl::integrandum_L(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq,
                                 int spinout, int spinin, TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        
        PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
            twovector.pGPDModel);
        // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        result[0]=(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution())
                    *z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
    //cout << u << " " << z << " " << result[0] << endl;
    return;

}

// integrandum for gamma_T amplitude
void TwoVector_Nucl::integrandum_T(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq, 
                                    int spinout, int spinint, TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        
        PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
            twovector.pGPDModel);

        // (Hq^u-Hq^d) * (2z-1) * z * barz / u / baru * Q    [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards] [factor p_perp from nominator taken out]
        result[0]=(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution())
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
     return;

}


// integrandum for gamma_L amplitude
void TwoVector_Nucl::integrandum_T_L(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq,
                                 int spinout, int spinin, TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        TransGPD_set gpds = twovector.gpdTgrid.getTransGPDSet(xi*(2.*u-1),xi,mandelstam_t*1.E06, scale);
        double gpdfactor = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,1,1,gpds).real();
        // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        result[0]=gpdfactor*
                    z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
    //cout << u << " " << z << " " << result[0] << endl;
    return;

}

// integrandum for gamma_T amplitude
void TwoVector_Nucl::integrandum_T_T(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq, 
                                    int spinout, int spinin,  TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        TransGPD_set gpds = twovector.gpdTgrid.getTransGPDSet(xi*(2.*u-1),xi,mandelstam_t*1.E06,scale);
        double gpdfactor = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.E06,mandelstam_t*1.E06,0.,1,1,gpds).real();  //only right?

        // (Hq^u-Hq^d) * (2z-1) * z * barz / u / baru * Q    [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards] [factor p_perp from nominator taken out]
        result[0]=gpdfactor
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
     return;

}


double TwoVector_Nucl::getCross_gammaL_rhoL(const double scale, const double xi, const double Q2, const double psq){
    double frhoplus = 0.198; //rho+ decay constant [GeV]
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);//tmin [GeV^2]

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
    numint::array<double,2> upper = {{1.,1.}};
    std::vector<double> integral(1,0.); //integrandum result array

    //gamma_L calculation
    F.f=integrandum_L;
    unsigned count=0;
    // integration over u and z
    //if(Q2>1.E-05) numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
    integral[0]*=36.; //prefactor DA (twice, squared)
    
    //cout << count << " " << integral[0] << endl;
    
    //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection
    integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)/(1+xi))*CF/Nc/psq/psq; 
    //cout << 16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)/(1+xi))*CF/Nc/psq/psq << endl;
    //prefactors Eq (15) Enberg et al.
    integral[0] *= -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); 
    //cout << 2.*PI*alpha_s*sqrt(Qsq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
    //cout << 1./256./pow(PI,3.)/xi/(1.-xi) << endl;
    double result_L=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.-xi); // Enberg et al Eq (12)
    //cout << Qsq << " " << xi << " " << psq << " " << result_L*0.389379E06 << endl;  // factor to go from GeV-6 to nb GeV-4

    return result_L*0.389379E06;  //[Gev-2 -> nb conversion]

}

double TwoVector_Nucl::getCross_gammaT_rhoL(const double scale, const double xi, const double Q2, const double psq){
    double frhoplus = 0.198; //rho+ decay constant [GeV]
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);//tmin [GeV^2]


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
    numint::array<double,2> upper = {{1.,1.}};
    std::vector<double> integral(1,0.); //integrandum result array

    //gamma_L calculation
    F.f=integrandum_T;
    unsigned count=0;
    // integration over u and z
    //numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
    integral[0]*=36.; //prefactor DA (twice, squared)
    cout << "int " << integral[0] << endl;

    integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)/(1+xi))*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al.
    cout << "pre1 " <<  16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)/(1+xi))*CF/Nc/psq/psq << " " << alpha_s << " " << frhoplus << " " << xi << " " << psq << endl;
    integral[0] *= PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.
    cout << "pre2 " << PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;

    double result_T=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.-xi); // Enberg et al Eq (12)
    cout << "pre3 " << 1./256./pow(PI,3.)/xi/(1.-xi) << endl;

    return result_T*0.389379E06; //[Gev-2 -> nb conversion]

}

double TwoVector_Nucl::getCross_gammaL_rhoT(const double scale, const double xi, const double Q2, const double psq){
    double frhoplus = 0.160; //rho+ decay constant [GeV]
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);//tmin [GeV^2]

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
    numint::array<double,2> upper = {{1.,1.}};

    //gamma_L calculation
    F.f=integrandum_T_L;
    unsigned count=0;
    double result_L=0.;
    // integration over u and z
    for(int spinin=-1;spinin<=1;spinin+=2){        
        for(int spinout=-1;spinout<=1;spinout+=2){        
            F.spinout=spinout;
            F.spinin=spinin;
            std::vector<double> integral(1,0.); //integrandum result array
            //if(Q2>1.E-05) numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
            
            integral[0]*=36.; //prefactor DA (twice, squared)
            //cout << spinin << " " << spinout << " " << integral[0] << endl;
            
            // cout << count << " " << integral[0] << endl;
            
            //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection
            integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)/(1+xi))*CF/Nc/psq/psq; 
            //cout << "pre1 " << 16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)/(1+xi))*CF/Nc/psq/psq << " " << alpha_s << " " << frhoplus << " " << xi << " " << psq <<  endl;
            //prefactors Eq (15) Enberg et al.
            integral[0] *= -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); 
            //cout << "pre2 " << -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
            //cout << "pre3 " << 1./256./pow(PI,3.)/xi/(1.-xi) << endl;
            result_L+=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.-xi); // Enberg et al Eq (12)

        }
    }
    //cout << Qsq << " " << xi << " " << psq << " " << result_L*0.389379E06 << endl;  // factor to go from GeV-6 to nb GeV-4

    return result_L*0.389379E06/2.;  //[Gev-2 -> nb conversion] + avg initial spin

}


double TwoVector_Nucl::getCross_gammaT_rhoT(const double scale, const double xi, const double Q2, const double psq){
    double frhoplus = 0.160; //rho+ decay constant [GeV]
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);//tmin [GeV^2]


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
    numint::array<double,2> upper = {{1.,1.}};

    //gamma_T calculation
    F.f=integrandum_T_T;
    unsigned count=0;
    double result_T=0.;
    // integration over u and z
    for(int spinin=-1;spinin<=1;spinin+=2){        
        for(int spinout=-1;spinout<=1;spinout+=2){        
            F.spinout=spinout;
            F.spinin=spinin;
            std::vector<double> integral(1,0.); //integrandum result array
            //numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
            integral[0]*=36.; //prefactor DA (twice, squared)
             cout << spinin << " " << spinout << " " << integral[0] << endl;

            integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)/(1+xi))*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al.
            cout << "pre1 " << 16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)/(1+xi))*CF/Nc/psq/psq << " " << frhoplus << endl;
            integral[0] *= PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.
            cout << "pre2 " << PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
            cout << "pre3 " << 1./256./pow(PI,3.)/xi/(1.-xi) << endl;
            // cout << PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
            result_T+=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.-xi); // Enberg et al Eq (12)

        }
    }


    return result_T*0.389379E06/2.; //[Gev-2 -> nb conversion] + avg initial spin

}
