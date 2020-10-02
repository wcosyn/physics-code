#include "DiDVCS.hpp"
#include <constants.hpp>
//#include <mstwpdf.h>
#include <cassert>
#include <numint/numint.hpp>
#include <cmath>
#include <gsl/gsl_sf_dilog.h>

using namespace std;

DiDVCS::DiDVCS(PARTONS::GPDService* pGPDService,
                                 PARTONS::GPDModule* pGPDModel, 
                                 PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule):
pGPDService(pGPDService),pGPDModel(pGPDModel),pRunningAlphaStrongModule(pRunningAlphaStrongModule){
    ;
}




void DiDVCS::integrandum_DiDVCS(numint::vector_z & result, double x_over_xi,double xi, double mandelstam_t, double scale, 
                                      double Qsqin, double Qsqout, int spinout, int spinin,
                                      DiDVCS::Photon_pol gamma, DiDVCS& didvcs_obj){

    result = vector<complex <double> >(2,0.);
    //limits are well behaved
    if(x_over_xi==-1||x_over_xi==1) for(int i=0;i<2;i++) result[i]=0.;
    else{ 
        double x=x_over_xi*xi;
        PARTONS::GPDKinematic gpdKinematic(x,xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult = didvcs_obj.pGPDService->computeSingleKinematic(gpdKinematic,didvcs_obj.pGPDModel);

        double Hu = gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
        double Hd = gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
        double Eu = gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
        double Ed = gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
        // cout << x << " "<< Hu << " " << Hd << " " << Eu << " " << Ed << endl;
        // exit(1);
        // Hard part selects symmetric part in x!  (C^- for rho/omega; C^+ for pi )
        PARTONS::GPDKinematic gpdKinematic2(-x,xi,mandelstam_t, scale*scale, scale*scale);

        PARTONS::GPDResult gpdResult2 = didvcs_obj.pGPDService->computeSingleKinematic(gpdKinematic2,didvcs_obj.pGPDModel);

        double Hu2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
        double Hd2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
        double Eu2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
        double Ed2 = gpdResult2.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

        //at t=tmin it will always be real
        // complex<double> gpdfactor_HE_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2-Hd-Hd2,Eu+Eu2-Ed-Ed2,1)); 
        // complex<double> gpdfactor_HE_is = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2+Hd+Hd2,Eu+Eu2+Ed+Ed2,1)); 
        // complex<double> gpdfactor_H_iv = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2-Hd-Hd2,0.,1)); 
        // complex<double> gpdfactor_H_is = (Deut_Conv_GPD_V::getGPD_even_nucl(spinin, spinout,xi,mandelstam_t,mandelstam_t,0.,Hu+Hu2+Hd+Hd2,0.,1)); 

        double z = Qsqout/Qsqin*(xi*xi-x*x)/(4.*xi*xi);
        double sqt4zp1 = sqrt(4.*z+1);
        complex<double> integrand = ((gsl_sf_dilog(-2./(sqt4zp1-1))-gsl_sf_dilog(+2./(sqt4zp1+1)))*(-2.*z/sqt4zp1)
                    -2.-(log(z)-I_UNIT*PI)*(1.-4.*z/sqt4zp1*atanh(1./sqt4zp1)))*xi;  //the *xi at the end because I integrate over x/xi
        result[0]=integrand*(2.*(Hu+Hu2)-(Hd+Hd2))/3.;  //C-odd
        result[1]=integrand*(2.*(Eu+Eu2)-(Ed+Ed2))/3.;  //C-odd combination
     }
    return;
}





void DiDVCS::getCross_DiDVCS(std::vector<double> & results, const double scale, const double xi, const double Q2in, const double Q2out, 
                          const double mandelstam_t, DiDVCS::Photon_pol gammapol, const int maxintsteps){


    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    //double mandelstam_t= -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);//tmin [GeV^2]


    Ftor_DiDVCS_general F;
    F.xi=xi;
    F.mandelstam_t=mandelstam_t;
    F.scale=scale;
    F.Qsqin=Q2in;
    F.Qsqout=Q2out;
    F.gammapol = gammapol;
    F.pobj=this;
    
    numint::mdfunction<numint::vector_z,1> mdf;
    mdf.func = &Ftor_DiDVCS_general::exec;
    mdf.param = &F;

    numint::array<double,1> lower = {{0.}};
    numint::array<double,1> upper = {{1.}};  // symmetry in u exploited to simplify integrand!

    //gamma_T calculation
    F.f=integrandum_DiDVCS;
    unsigned count=0;
    results = vector<double>(2,0.);
    // integration over u and z
    // for(int spinin=-1;spinin<=1;spinin+=2){        
    //     for(int spinout=-1;spinout<=1;spinout+=2){        
    //         F.spinout=spinout;
    //         F.spinin=spinin;

            std::vector<complex<double> > integral(2,0.); //integrandum result array
            numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,maxintsteps,integral,count,0);
            cout << xi << " " << integral[0].real() << " " << integral[0].imag() << " " << integral[1].real() << " " << integral[1].imag() << endl;
            for(int index=0;index<2;index++){
                //cout << spinin << " " << spinout << " " << integral[index] << endl;
            
                //prefactors Eq (14) Enberg et al., not including s factor since it drops out in xsection
                // sqrt(1-xi^2) compensated  in the helicity amplitudes
                double K = pow(2,7.)*PI*PI*ALPHA*CF*CF*alpha_s*alpha_s/8.;
                integral[index] *= 3.*sqrt(2.)*K*frho0/Q2in/sqrt(Q2in*Q2out); //prefactors Eq (23) Pire et al.
//                results[index]+=norm(integral[index])/256./pow(PI,3.)/xi/(1+xi); // Enberg et al Eq (12), corrected phase space factor!
        //     }
        // }
    }

    //[Gev-2 -> nb conversion] + avg initial spin + extra factor of 1/2 for omega/rho0 wf (u+/-d)/sqrt(2) 
    //+ extra factor 1/2 from averaging over sin^2 theta in case of transverse vector meson (or from (RL+LR)/2 products)
    //for(int index=0;index<12;index++) results[index]*=0.389379E06/4./(rhopol==DiDVCS::krhoT? 2.:1.); 
    return;

}

void DiDVCS::getElectro_Cross_DiDVCS(std::vector<double> & results, std::vector<double> & resultsL, std::vector<double> & resultsT, 
                const double y, const double Q2){
    
    results= vector<double> (2,0.);
    
        
    double epsilon = (1-y)/(1-y+y*y/2.);
    for(int kk=0;kk<2;kk++) results[kk] = (resultsL[kk]+epsilon/2.*resultsT[kk])*ALPHA/2./PI*y/Q2/(1-epsilon);
    
}
