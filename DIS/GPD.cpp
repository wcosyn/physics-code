#include <GPD.hpp>
#include <constants.hpp>

using namespace std;

vector< complex<double> > helamps_to_gpds(const double xi, const double x, const double t, const vector< complex <double> > & helamps){
    
    vector< complex<double> > gpds(9,0.);

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



void int_k3(numint::vector_z & res, double knorm, double kcosth, double kphi, TDeuteron::Wavefunction *wfref,
              double x, double xi, double t, int pold_in, int pold_out){
    res=numint::vector_z(1,0.);
    double Ek=sqrt(knorm*knorm+MASSn*MASSn);
    double kz=knorm*kcosth;
    double alpha1=1+kz/Ek;

    if((alpha1/2.*(1+xi)-fabs(x)-xi)<0.) return;
    if((alpha1/2.*(1+xi)-2.*xi)<0.) return;
    double t0=-4.*MASSn*MASSn*xi*xi/(1-xi*xi);
    double Delta_perp=sqrt((t0-t)*(1-xi*xi));

    double alphaprime=(alpha1*(1+xi)-4.*xi)/(1-xi);

    double xi_n=xi/(alpha1/2.*(1+xi)-xi);

    double x_n=x/(alpha1/2.*(1+xi)-xi);


    double ksinth=sqrt(1-kcosth*kcosth);

    double kxprime=knorm*ksinth*sin(kphi)+(1-alpha1/2.)/(1-xi)*Delta_perp;
    double kyprime=knorm*ksinth*cos(kphi);
    
    double kperpprime=sqrt(kxprime*kxprime+kyprime*kyprime);
    double Ekprime = sqrt((MASSn*MASSn+kperpprime*kperpprime)/(alphaprime*(2.-alphaprime)));
    double kzprime = (alphaprime -1)*Ekprime;
    TVector3 k_in(knorm*ksinth*cos(kphi),knorm*ksinth*sin(kphi),knorm*kcosth);
    TVector3 k_out(kxprime, kyprime, kzprime);

    complex<double> result=0.;
    for(int sigma2=-1;sigma2<=1;sigma2+=2){
        for(int sigma1in=-1; sigma1in<=1; sigma1in+=2){
            complex<double> wfin = wfref->DeuteronPState(2*pold_in, sigma2, sigma1in, k_in);
            for(int sigma1out=-1; sigma1out<=1; sigma1out+=2){
                complex<double> wfout=wfref->DeuteronPState(2*pold_out, sigma2, sigma1out, k_out);
                result+=wfin*conj(wfout)*getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0);
            }
        }
    }
    res[0]=result;
}

complex<double> getGPD_odd_nucl(int sigma_in, int sigma_out, double x, double xi, double t, double t0, double phi){

    if(sigma_in==sigma_out) return (sqrt(t0-t)/MASSn*getGPD_HtildeT(x,xi,t)
                            +(1-sigma_in*xi)*sqrt(t0-t)/(2.*MASSn)*getGPD_ET(x,xi,t)
                            +sigma_in*(1-sigma_in*xi)*sqrt(t0-t)/(2.*MASSn)*getGPD_EtildeT(x,xi,t))*exp(I_UNIT*phi);

    else return (sigma_in-1)/2*sqrt(1.-xi*xi)*getGPD_HT(x,xi,t)
                -sigma_in/2.*(t0-t)/MASSn*sqrt(1-xi*xi)*exp((sigma_in+1)*phi*I_UNIT)*getGPD_HtildeT(x,xi,t)
                +(sigma_in-1)*xi*xi/sqrt(1-xi*xi)*getGPD_ET(x,xi,t)
                -(sigma_in-1)/2.*xi/sqrt(1.-xi*xi)*getGPD_EtildeT(x,xi,t);

}


double getGPD_HT(double x,double xi, double t){
    return 0.;
}
double getGPD_HtildeT(double x,double xi, double t){
    return 0.;
}
double getGPD_ET(double x,double xi, double t){
    return 0.;
}
double getGPD_EtildeT(double x,double xi, double t){
    return 0.;
}