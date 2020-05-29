#include "TwoVector_Nucl.hpp"
#include <constants.hpp>
#include <mstwpdf.h>
#include <cassert>
#include <numint/numint.hpp>
#include <cmath>

#include "Deut_Conv_GPD_T.hpp"
#include "Deut_Conv_GPD_V.hpp"
using namespace std;

TwoVector_Nucl::TwoVector_Nucl(PARTONS::GPDService* pGPDService,
                                 PARTONS::GPDModule* pGPDModel, 
                                 PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule):
pGPDService(pGPDService),pGPDModel(pGPDModel),pRunningAlphaStrongModule(pRunningAlphaStrongModule),gpdTgrid("MSTW"){
    ;
}



// integrandum for gamma_L amplitude
void TwoVector_Nucl::integrandum_rhoL_gammaL(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq,
                                 int spinout, int spinin, TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        
        PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
            twovector.pGPDModel);

        PARTONS::GPDKinematic gpdKinematic2(-xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult2 = twovector.pGPDService->computeGPDModel(gpdKinematic2,
            twovector.pGPDModel);

        // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        result[0]=(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                   +gpdResult2.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()
                    -gpdResult2.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution())
                    *z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
    //cout << u << " " << z << " " << result[0] << endl;
    return;

}

// integrandum for gamma_T amplitude
void TwoVector_Nucl::integrandum_rhoL_gammaT(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq, 
                                    int spinout, int spinint, TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        
        PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
            twovector.pGPDModel);

        PARTONS::GPDKinematic gpdKinematic2(-xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult2 = twovector.pGPDService->computeGPDModel(gpdKinematic2,
            twovector.pGPDModel);

        // (Hq^u-Hq^d) * (2z-1) * z * barz / u / baru * Q    [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards] [factor p_perp from nominator taken out]
        result[0]=(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    +gpdResult2.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()
                    -gpdResult2.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution())
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
    // cout << (2.*z-1)*z*(1.-z)/u/(1.-u)
    //                     *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) )<< " ";
    // u = 1.-u;
    // cout << (2.*z-1)*z*(1.-z)/u/(1.-u)
    //                     *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) ) << endl;                    
     return;

}

void TwoVector_Nucl::integrandum_rho_general(numint::vector_z & result, double u, double z, double xi, double mandelstam_t, double scale, 
                                      double psq, double Qsq, int spinout, int spinin, TwoVector_Nucl::Rho_pol rhopol, 
                                      TwoVector_Nucl::Photon_pol gamma, TwoVector_Nucl& twovector){

    result = vector<complex <double> >(12,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) for(int i=0;i<12;i++) result[i]=0.;
    else{ 
        if(rhopol==TwoVector_Nucl::krhoL){
            PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

            PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
                twovector.pGPDModel);

            double Hu = gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Hd = gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
            double Eu = gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Ed = gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
        
            PARTONS::GPDKinematic gpdKinematic2(-xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

            PARTONS::GPDResult gpdResult2 = twovector.pGPDService->computeGPDModel(gpdKinematic2,
                twovector.pGPDModel);

            double Hu2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Hd2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
            double Eu2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Ed2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

            //at t=tmin it will always be real
            complex<double> gpdfactor_HE_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2-Hd-Hd2,Eu+Eu2-Ed-Ed2,1)); 
            complex<double> gpdfactor_HE_is = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2+Hd+Hd2,Eu+Eu2+Ed+Ed2,1)); 
            complex<double> gpdfactor_H_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2-Hd-Hd2,0.,1)); 
            complex<double> gpdfactor_H_is = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2+Hd+Hd2,0.,1)); 

            result[0]= gpdfactor_HE_iv;  //H+E isovector
            result[1]= gpdfactor_H_iv;  //H isovector 
            result[3]= gpdfactor_HE_iv;  //H+E isovector 
            result[4]= gpdfactor_H_iv;  //H isovector 
            result[6]= gpdfactor_HE_is;  //H+E isoscalar 
            result[7]= gpdfactor_H_is;  //H isoscalar 
            result[9]= gpdfactor_HE_is;  //H+E isoscalar 
            result[10]= gpdfactor_H_is;  //H  isooscalar
        }
        if(rhopol==TwoVector_Nucl::kaxial){
            PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

            PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
                twovector.pGPDModel);

            double Hu = gpdResult.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Hd = gpdResult.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
            double Eu = gpdResult.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Ed = gpdResult.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

            PARTONS::GPDKinematic gpdKinematic2(-xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

            PARTONS::GPDResult gpdResult2 = twovector.pGPDService->computeGPDModel(gpdKinematic2,
                twovector.pGPDModel);

            double Hu2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Hd2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
            double Eu2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Ed2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

            //at t=tmin it will always be real
            complex<double> gpdfactor_HtEt_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2-Hd-Hd2,Eu+Eu2-Ed-Ed2,0)); 
            complex<double> gpdfactor_Ht_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2-Hd-Hd2,0.,0)); 
        
            result[0]= gpdfactor_HtEt_iv;  //H+E isovector
            result[1]= gpdfactor_Ht_iv;  //H isovector 
            result[3]= gpdfactor_HtEt_iv;  //H+E isovector 
            result[4]= gpdfactor_Ht_iv;  //H isovector 
        }
        if(rhopol==TwoVector_Nucl::krhoT){
            TransGPD_set gpds = twovector.gpdTgrid.getTransGPDSet(xi*(2.*u-1),xi,mandelstam_t*1.E06, scale);
            TransGPD_set gpds2 = twovector.gpdTgrid.getTransGPDSet(-xi*(2.*u-1),xi,mandelstam_t*1.E06, scale);
            //factor of 2 because of (u+/-d)/2 in helicity amplitudes!
            complex<double> gpd_model0_iv = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,0,1,0,gpds+gpds2)*2.;
            complex<double> gpd_model1_iv = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,1,1,0,gpds+gpds2)*2.;
            complex<double> gpd_model2_iv = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,2,1,0,gpds+gpds2)*2.;

            complex<double> gpd_model0_is = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,0,1,1,gpds+gpds2)*2.;
            complex<double> gpd_model1_is = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,1,1,1,gpds+gpds2)*2.;
            complex<double> gpd_model2_is = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,2,1,1,gpds+gpds2)*2.;

            // cout << xi*(2.*u-1) << " " << gpd_model0_iv.real() << " " << gpd_model1_iv.real() << " " << gpd_model2_iv.real()
            // << " " << gpd_model0_is.real() << " " << gpd_model1_is.real() << " " << gpd_model2_is.real() << endl;

            result[0]=result[3]= gpd_model0_iv;
            result[1]=result[4]= gpd_model1_iv;
            result[2]=result[5]= gpd_model2_iv;

            result[6]=result[9]= gpd_model0_is;
            result[7]=result[10]= gpd_model1_is;
            result[8]=result[11]= gpd_model2_is;

        }

        // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [6z*barz,6u*baru from 2 DA taken into account]
        double integrand = (gamma == TwoVector_Nucl::kgammaL?  36./sqrt(1-xi*xi) //sqrt factor from helamps accounted for in prefactor 
                    *z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) )
                        : 36./sqrt(1-xi*xi)   //sqrt factor from helamps accounted for in prefactor 
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) ));
        for(int i=0;i<12;i++) result[i]*=integrand;                
        for(int i=0;i<3;i++){
            result[3+i]*=64./PI/PI/36./sqrt(z*(1-z)*u*(1-u));  //AdS DA
            result[9+i]*=64./PI/PI/36./sqrt(z*(1-z)*u*(1-u));  //AdS DA
        }
     }
    return;
}



// rhoL gammaL but through helicity amplitudes
void TwoVector_Nucl::integrandum_rhoL_gammaL_helamps(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq,
                                 int spinout, int spinin, TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        
        PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
            twovector.pGPDModel);
        double H_isovec = gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

        double E_isovec = gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

        //at t=tmin it will always be real
        double gpdfactor = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,H_isovec,E_isovec,1)).real();
        // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        result[0]= gpdfactor/sqrt(1-xi*xi) //sqrt factor from helamps accounted for in prefactor 
                    *z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
    //cout << u << " " << z << " " << result[0] << endl;
    return ;

}

void TwoVector_Nucl::integrandum_rho_gammaT_general(numint::vector_z & result, double u, double z, double xi, double mandelstam_t, double scale, 
                                            double psq, double Qsq, int spinout, int spinin, TwoVector_Nucl::Rho_pol rhopol, TwoVector_Nucl& twovector){

    result = vector<complex <double> >(12,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) for(int i=0;i<12;i++) result[i]=0.;
    else{ 
        if(rhopol==TwoVector_Nucl::krhoL){
            PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

            PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
                twovector.pGPDModel);

            double Hu = gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Hd = gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
            double Eu = gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Ed = gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

            //at t=tmin it will always be real
            complex<double> gpdfactor_HE_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu-Hd,Eu-Ed,1)); 
            complex<double> gpdfactor_HE_is = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hd,Eu+Ed,1)); 
            complex<double> gpdfactor_H_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu-Hd,0.,1)); 
            complex<double> gpdfactor_H_is = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hd,0.,1)); 
        
            result[0]= gpdfactor_HE_iv;  //H+E isovector
            result[1]= gpdfactor_H_iv;  //H isovector 
            result[3]= gpdfactor_HE_iv;  //H+E isovector 
            result[4]= gpdfactor_H_iv;  //H isovector 
            result[6]= gpdfactor_HE_is;  //H+E isoscalar 
            result[7]= gpdfactor_H_is;  //H isoscalar 
            result[9]= gpdfactor_HE_is;  //H+E isoscalar 
            result[10]= gpdfactor_H_is;  //H  isooscalar
        }
        if(rhopol==TwoVector_Nucl::kaxial){
            PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

            PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
                twovector.pGPDModel);

            double Hu = gpdResult.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Hd = gpdResult.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
            double Eu = gpdResult.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Ed = gpdResult.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

            //at t=tmin it will always be real
            complex<double> gpdfactor_HtEt_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu-Hd,Eu-Ed,0)); 
            complex<double> gpdfactor_Ht_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu-Hd,0.,0)); 
        
            result[0]= gpdfactor_HtEt_iv;  //H+E isovector
            result[1]= gpdfactor_Ht_iv;  //H isovector 
            result[3]= gpdfactor_HtEt_iv;  //H+E isovector 
            result[4]= gpdfactor_Ht_iv;  //H isovector 
        }
        if(rhopol==TwoVector_Nucl::krhoT){
            TransGPD_set gpds = twovector.gpdTgrid.getTransGPDSet(xi*(2.*u-1),xi,mandelstam_t*1.E06, scale);
            //factor of 2 because helicity amplitudes have (u+/-d)/2
            complex<double> gpd_model0_iv = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,0,1,0,gpds)*2.;
            complex<double> gpd_model1_iv = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,1,1,0,gpds)*2.;
            complex<double> gpd_model2_iv = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,2,1,0,gpds)*2.;

            complex<double> gpd_model0_is = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,0,1,1,gpds)*2.;
            complex<double> gpd_model1_is = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,1,1,1,gpds)*2.;
            complex<double> gpd_model2_is = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,2,1,1,gpds)*2.;

            result[0]=result[3]= gpd_model0_iv;
            result[1]=result[4]= gpd_model1_iv;
            result[2]=result[5]= gpd_model2_iv;

            result[6]=result[9]= gpd_model0_is;
            result[7]=result[10]= gpd_model1_is;
            result[8]=result[11]= gpd_model2_is;

        }
        // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        double integrand = 36./sqrt(1-xi*xi)   //sqrt factor from helamps accounted for in prefactor 
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );
        for(int i=0;i<12;i++) result[i]*=integrand;                
        for(int i=0;i<3;i++){
            result[3+i]*=64./PI/PI/36./sqrt(z*(1-z)*u*(1-u));  //AdS DA
            result[9+i]*=64./PI/PI/36./sqrt(z*(1-z)*u*(1-u));  //AdS DA
        }
     }
    return;
                                            }


// integrandum for rhoL gamma_T amplitude, but using helicity amplitudes
void TwoVector_Nucl::integrandum_rhoL_gammaT_helamps(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq, 
                                    int spinout, int spinin, TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        
        PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult = twovector.pGPDService->computeGPDModel(gpdKinematic,
            twovector.pGPDModel);

        double H_isovec = gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

        double E_isovec = gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

        //at t=tmin it will always be real
        double gpdfactor = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,H_isovec,E_isovec,1).real());
        // (Hq^u-Hq^d) * (2z-1) * z * barz / u / baru * Q    [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards] [factor p_perp from nominator taken out]
        result[0]= gpdfactor/sqrt(1-xi*xi)   //sqrt factor from helamps accounted for in prefactor 
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
     return;

}

// integrandum for gamma_L amplitude
void TwoVector_Nucl::integrandum_rhoT_gammaL(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq,
                                 int spinout, int spinin, TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        TransGPD_set gpds = twovector.gpdTgrid.getTransGPDSet(xi*(2.*u-1),xi,mandelstam_t*1.E06, scale);
        double gpdfactor = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,0,1,0,gpds).real();
        // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        result[0]=gpdfactor*2./sqrt(1-xi*xi)* //sqrt factor is in prefactor of Enberg paper; factor 2 because isovector combination has /2 in it in the code
                    z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
    //cout << u << " " << z << " " << result[0] << endl;
    return;

}

// integrandum for gamma_T amplitude
void TwoVector_Nucl::integrandum_rhoT_gammaT(numint::vector_d &result, double u, double z, double xi, double mandelstam_t, double scale, double psq, double Qsq, 
                                    int spinout, int spinin,  TwoVector_Nucl &twovector){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        TransGPD_set gpds = twovector.gpdTgrid.getTransGPDSet(xi*(2.*u-1),xi,mandelstam_t*1.E06,scale);
        double gpdfactor = Deut_Conv_GPD_T::getGPD_odd_nucl(spinin,spinout,xi,mandelstam_t*1.E06,mandelstam_t*1.E06,0.,0,1,0,gpds).real();  //only right?

        // (Hq^u-Hq^d) * (2z-1) * z * barz / u / baru * Q    [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards] [factor p_perp from nominator taken out]
        result[0]=gpdfactor*2./sqrt(1-xi*xi) //sqrt factor is in prefactor of Enberg paper; factor 2 because isovector combination has /2 in it in the code
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
    numint::array<double,2> upper = {{0.5,1.}};
    std::vector<double> integral(1,0.); //integrandum result array

    //gamma_L calculation
    F.f=integrandum_rhoL_gammaL;
    unsigned count=0;
    // integration over u and z
    /* if(Q2>1.E-05) */ numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E04,integral,count,0);
    integral[0]*=36.; //prefactor DA (twice, squared)
    double kk=integral[0];
    //cout << "gammaL H " << " " << integral[0] << endl;
    
    //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection, corrected xi factor
    integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq; 
    //cout << "pre1 gammaL " << 16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq << endl;
    //prefactors Eq (15) Enberg et al.
    integral[0] *= -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); 
    //cout << "pre2 " << -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
    //cout << "pre3 " << 1./256./pow(PI,3.)/xi/(1+xi) << endl;
    double result_L=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1+xi); // Enberg et al Eq (12), corrected phase space factor!
    //cout << Qsq << " " << xi << " " << psq << " " << result_L*0.389379E06 << endl;  // factor to go from GeV-6 to nb GeV-4
    //cout << pow(16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq*2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI),2.)/256./pow(PI,3.)/xi/(1+xi)*0.389379E06 << " " << kk*kk << " " ;
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
    numint::array<double,2> upper = {{0.5,1.}};
    std::vector<double> integral(1,0.); //integrandum result array

    //gamma_L calculation
    F.f=integrandum_rhoL_gammaT;
    unsigned count=0;
    // integration over u and z
    numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E04,integral,count,0);
    integral[0]*=36.; //prefactor DA (twice, squared)
    double kk=integral[0];
    //cout << "gammaT H " << integral[0] << endl;

    integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al., corrected xi factor
    //cout << "pre1 gammaT " <<  16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq << endl;
    integral[0] *= PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.
    //cout << "pre2 " << PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;

    double result_T=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1+xi); // Enberg et al Eq (12), corrected phase space factor!
    //cout << "pre3 " << 1./256./pow(PI,3.)/xi/(1+xi) << endl;
    //cout << pow(16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq*PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI),2.)/256./pow(PI,3.)/xi/(1+xi)*0.389379E06 << " " << kk*kk << " " ;

    return result_T*0.389379E06; //[Gev-2 -> nb conversion]

}

double TwoVector_Nucl::getCross_gammaL_rhoT(const double scale, const double xi, const double Q2, const double psq, bool rhoT){
    double frhoplus = rhoT? 0.160: 0.198; //rho+ decay constant [GeV]
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
    if(rhoT) F.f=integrandum_rhoT_gammaL;
    else F.f=integrandum_rhoL_gammaL_helamps;
    unsigned count=0;
    double result_L=0.;
    // integration over u and z
    for(int spinin=-1;spinin<=1;spinin+=2){        
        for(int spinout=-1;spinout<=1;spinout+=2){        
            F.spinout=spinout;
            F.spinin=spinin;
            std::vector<double> integral(1,0.); //integrandum result array
            /* if(Q2>1.E-05) */ numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E04,integral,count,0);
            
            integral[0]*=36.; //prefactor DA (twice, squared)
            //cout << "gammaL " << spinin << " " << spinout << " " << integral[0] << endl;
            
            
            //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection
            // sqrt(1-xi^2) compensated in the helicity amplitudes
            integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq; 
            //cout << "pre1 helamps " << 16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq << endl;
            //prefactors Eq (15) Enberg et al.
            integral[0] *= -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); 
            //cout << "pre2 " << -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
            //cout << "pre3 " << 1./256./pow(PI,3.)/xi/(1+xi) << endl;
            result_L+=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1+xi); // Enberg et al Eq (12), corrected phase space factor!

        }
    }
    //cout << Qsq << " " << xi << " " << psq << " " << result_L*0.389379E06 << endl;  // factor to go from GeV-6 to nb GeV-4

    return result_L*0.389379E06/2.;  //[Gev-2 -> nb conversion] + avg initial spin

}


double TwoVector_Nucl::getCross_gammaT_rhoT(const double scale, const double xi, const double Q2, const double psq, bool rhoT){
    double frhoplus =  rhoT? 0.160: 0.198; //rho+ decay constant [GeV]
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
    if(rhoT) F.f=integrandum_rhoT_gammaT;
    else F.f=integrandum_rhoL_gammaT_helamps;
    unsigned count=0;
    double result_T=0.;
    // integration over u and z
    for(int spinin=-1;spinin<=1;spinin+=2){        
        for(int spinout=-1;spinout<=1;spinout+=2){        
            F.spinout=spinout;
            F.spinin=spinin;
            std::vector<double> integral(1,0.); //integrandum result array
            numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E04,integral,count,0);
            integral[0]*=36.; //prefactor DA (twice, squared)
            //cout << "gammaT " << spinin << " " << spinout << " " << integral[0] << endl;

            //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection
            // sqrt(1-xi^2) compensated  in the helicity amplitudes
            integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al.
            //cout << "pre1 helamps gammaT" << 16.*PI*PI*alpha_s*frhoplus*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq << endl;
            integral[0] *= PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.
            //cout << "pre2 " << PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
            //cout << "pre3 " << 1./256./pow(PI,3.)/xi/(1+xi) << endl;
            // cout << PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI) << endl;
            result_T+=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1+xi); // Enberg et al Eq (12), corrected phase space factor!

        }
    }


    return result_T*0.389379E06/2.; //[Gev-2 -> nb conversion] + avg initial spin

}


void TwoVector_Nucl::getCross_twovector(std::vector<double> & results, const double scale, const double xi, const double Q2, const double psq, 
                          TwoVector_Nucl::Photon_pol gammapol, TwoVector_Nucl::Rho_pol rhopol){


    double f_lowermeson =  0.;  // meson decay constant for meson originating from lower part of the graph
    if (rhopol == TwoVector_Nucl::krhoT) f_lowermeson=0.160;
    if (rhopol == TwoVector_Nucl::krhoL) f_lowermeson=0.198;
    if (rhopol == TwoVector_Nucl::kaxial) f_lowermeson=0.130;
    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    double mandelstam_t= -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);//tmin [GeV^2]


    Ftor_2vector_general F;
    F.xi=xi;
    F.mandelstam_t=mandelstam_t;
    F.scale=scale;
    F.Qsq=Q2;
    F.psq=psq;
    F.rhopol = rhopol;
    F.gammapol = gammapol;
    F.pobj=this;
    
    numint::mdfunction<numint::vector_z,2> mdf;
    mdf.func = &Ftor_2vector_general::exec;
    mdf.param = &F;

    numint::array<double,2> lower = {{0.,0.}};
    numint::array<double,2> upper = {{0.5,1.}};

    //gamma_T calculation
    F.f=integrandum_rho_general;
    unsigned count=0;
    results = vector<double>(12,0.);
    // integration over u and z
    for(int spinin=-1;spinin<=1;spinin+=2){        
        for(int spinout=-1;spinout<=1;spinout+=2){        
            F.spinout=spinout;
            F.spinin=spinin;

            std::vector<complex<double> > integral(12,0.); //integrandum result array
            numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E04,integral,count,0);
            for(int index=0;index<12;index++){
                //cout << spinin << " " << spinout << " " << integral[index] << endl;
            
                //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection
                // sqrt(1-xi^2) compensated  in the helicity amplitudes
                integral[index] *= 16.*PI*PI*alpha_s*f_lowermeson*xi*sqrt(1-xi*xi)*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al.
                if(gammapol==TwoVector_Nucl::kgammaL) integral[index] *= -2.*PI*alpha_s*sqrt(Q2)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI);; //prefactors Eq (17) Enberg et al.
                if(gammapol==TwoVector_Nucl::kgammaT) integral[index] *= PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.
                results[index]+=norm(integral[index])/256./pow(PI,3.)/xi/(1+xi); // Enberg et al Eq (12), corrected phase space factor!
            }
        }
    }

    //[Gev-2 -> nb conversion] + avg initial spin + extra factor of 1/2 for omega/rho0 wf (u+/-d)/sqrt(2) 
    //+ extra factor 1/2 from averaging over sin^2 theta in case of transverse vector meson
    for(int index=0;index<12;index++) results[index]*=0.389379E06/4./(rhopol==TwoVector_Nucl::krhoT? 2.:1.); 
    return;

}
