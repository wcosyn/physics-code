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
                                      DiDVCS& didvcs_obj){

    result = vector<complex <double> >(4,0.);
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
        // asymptotic DA result (see Eq. (23) Pire, Szym.,Wallon paper)
        complex<double> integrand = ((gsl_sf_dilog(-2./(sqt4zp1-1))-gsl_sf_dilog(+2./(sqt4zp1+1)))*(-2.*z/sqt4zp1)
                    -2.-(log(z)-I_UNIT*PI)*(1.-4.*z/sqt4zp1*atanh(1./sqt4zp1)))*xi;  //the *xi at the end because I integrate over x/xi
        // ads DA result, only imaginary part (delta function in k^2), evaulated through mathematica 
        complex<double> integrand2 = I_UNIT*PI*PI*(1.-2.*sqrt(z/(1.+4.*z)))*xi/6.*8./PI; //the *xi at the end because I integrate over x/xi
        result[0]=integrand*(2.*(Hu+Hu2)-(Hd+Hd2))/3.;  //CFF H
        result[1]=integrand*(2.*(Eu+Eu2)-(Ed+Ed2))/3.;  //CFF E
        result[2]=integrand2*(2.*(Hu+Hu2)-(Hd+Hd2))/3.;  //CFF H
        result[3]=integrand2*(2.*(Eu+Eu2)-(Ed+Ed2))/3.;  //CFF E
     }
    return;
}





void DiDVCS::getCross_DiDVCS(vector<double> & results, const double scale, const double t_rho, const double t_N, const double Q2in, const double Q2out, 
                          const double s2, const double s_eN, double y, const int maxintsteps, bool no_realpart){


    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = pRunningAlphaStrongModule->compute(scale*scale); //scale is 1 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;
    //double mandelstam_t= -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);//tmin [GeV^2]

    //kinematics
    double s_gN = s_eN*y;
    double s1=s_gN*Q2out/s2;
    double xi = (Q2in +s_gN*(Q2out/s2))/(2.*s_gN+Q2in -s_gN*(Q2out/s2));
//    cout << "xi =" << xi << endl;

    Ftor_DiDVCS_general F;
    F.xi=xi;
    F.mandelstam_t=t_N;
    F.scale=scale;
    F.Qsqin=Q2in;
    F.Qsqout=Q2out;
    F.pobj=this;
    
    numint::mdfunction<numint::vector_z,1> mdf;
    mdf.func = &Ftor_DiDVCS_general::exec;
    mdf.param = &F;

    numint::array<double,1> lower = {{0.}};
    numint::array<double,1> upper = {{1.}};  // symmetry in u exploited to simplify integrand!

    //gamma_T calculation
    F.f=integrandum_DiDVCS;
    unsigned count=0;
    results = vector<double>(4,0.);
    std::vector<complex<double> > integral(4,0.); //integrandum result array [dimensionless!]
    numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,maxintsteps,integral,count,0); 
    //cout << xi << " " << integral[0].real() << " " << integral[0].imag() << " " << integral[1].real() << " " << integral[1].imag() << endl;
    double K = pow(2,7.)*PI*PI*ALPHA*CF*CF*alpha_s*alpha_s/8.; //[dimensionless]
    if(no_realpart){
        integral[0]=complex<double> (0.,integral[0].imag());
        integral[1]=complex<double> (0.,integral[1].imag());
    }
    double CFF_asympt = pow(3.*sqrt(2.)*K*frho0/Q2in/sqrt(Q2in*Q2out),2.)*((1-xi*xi)*norm(integral[0])
                        -2.*xi*xi*((conj(integral[0])*integral[1]).real())-(xi*xi+t_N/4./MASSP_G/MASSP_G)*norm(integral[1])); //[GeV^-6]
    double CFF_ads = pow(3.*sqrt(2.)*K*frho0/Q2in/sqrt(Q2in*Q2out),2.)*((1-xi*xi)*norm(integral[2])
                        -2.*xi*xi*((conj(integral[2])*integral[3]).real())-(xi*xi+t_N/4./MASSP_G/MASSP_G)*norm(integral[3])); //[GeV^-6]
    results[0] = CFF_asympt*ALPHA/pow(1+xi,2.)/(48*pow(2.*PI,4.)*Q2out*(s2-t_rho)); //[GeV-10]
    results[2] = CFF_ads*ALPHA/pow(1+xi,2.)/(48*pow(2.*PI,4.)*Q2out*(s2-t_rho)); //[GeV-10]
    // cout << 6.*frho0*frho0*pow(ALPHA/Q2in,3.)*pow(alpha_s*CF,4.)*(1.-xi)/(1+xi)/pow(Q2out,2.)/(s2-t_rho)*0.389379E06 << " " << ((1-xi*xi)*norm(integral[0])
    //                     -2.*xi*((integral[0]*integral[1]).real())-(xi*xi+t_N/4./MASSP_G/MASSP_G)*norm(integral[1]))/(1.-xi*xi) << " ";
    //[Gev-2 -> nb conversion]
    results[0] *=0.389379E06; //[nb GeV-8]
    results[2] *=0.389379E06; //[nb GeV-8]

    double epsilon = (1-y)/(1-y+y*y/2.);
    results[1] = results[0]*ALPHA/2./PI*y/Q2in/(1-epsilon); //[nb GeV-10]
    results[3] = results[2]*ALPHA/2./PI*y/Q2in/(1-epsilon); //[nb GeV-10]
    return;

}
    
