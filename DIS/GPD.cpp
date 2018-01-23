#include "GPD.hpp"
#include <constants.hpp>
#include <mstwpdf.h>
#include <cassert>
#include <numint/numint.hpp>

using namespace std;

GPD::GPD(const std::string &pdf_name){

    if(!pdf_name.compare("MSTW")){
        string file= string(HOMEDIR)+"/mstw_grids/mstw2008lo.00.dat";
        mstw = new c_mstwpdf(file,false,true);
    }
    else {
        mstw=NULL;
        cerr << "you have not chosen a valid pdf set in GPD::GPD()" << endl;
        assert(1==0);
        exit(1);
    }

  

}

GPD::~GPD(){
      if(mstw) delete mstw;
}


vector< complex<double> > GPD::helamps_to_gpds(const double xi, const double x, const double t, const vector< complex <double> > & helamps){
    
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



void GPD::int_k3(numint::vector_z & res, double knorm, double kcosth, double kphi, TDeuteron::Wavefunction *wfref,
              double x, double xi, double t, int pold_in, int pold_out, int model){
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

    double kxprime=knorm*ksinth*sin(kphi)+(1-alpha1/2.)/(1-xi)*Delta_perp;  //CHECK!!! k is not real momentum!
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
                result+=wfin*conj(wfout)*getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0,0.,model); ///CHECK PHI!!!  // change so all 4 pols are returned at once
            }
        }
    }
    res[0]=result*2.*knorm*knorm/Ek;
}

complex<double> GPD::getGPD_odd_nucl(int sigma_in, int sigma_out, double x, double xi, double t, double t0, double phi, int model){


    TransGPD_set gpd_nucl=getGK_param(x,xi,t);

    if(sigma_in==sigma_out) return (sqrt(t0-t)/MASSn*gpd_nucl.getHtildeT_singlet(model)
                            +(1-sigma_in*xi)*sqrt(t0-t)/(2.*MASSn)*gpd_nucl.getET_singlet(model)
                            +sigma_in*(1-sigma_in*xi)*sqrt(t0-t)/(2.*MASSn)*gpd_nucl.getEtildeT_singlet())*exp(I_UNIT*phi);

    else return (sigma_in-1)/2*sqrt(1.-xi*xi)*gpd_nucl.getHT_singlet()
                -sigma_in/2.*(t0-t)/MASSn*sqrt(1-xi*xi)*exp((sigma_in+1)*phi*I_UNIT)*gpd_nucl.getHtildeT_singlet(model)
                +(sigma_in-1)*xi*xi/sqrt(1-xi*xi)*gpd_nucl.getET_singlet(model)
                -(sigma_in-1)/2.*xi/sqrt(1.-xi*xi)*gpd_nucl.getEtildeT_singlet();

}


TransGPD_set GPD::getGK_param(const double x, const double xi, const double t){

    if(fabs(x)>1.) return TransGPD_set(0.,0.,0.,0.);
    GPD::Ftor_doubledistr F;
    F.x=x;
    F.xi=xi;
    F.t=t;
    F.gpdobject=this;
    numint::mdfunction<numint::vector_d,1> mdf;
    mdf.func = &Ftor_doubledistr::exec;
    mdf.param = &F;
    numint::vector_d ret(4,0.);
    F.f=GPD::DD_int_rho;

    double low=max(0.,(x-xi)/(1.-xi));
    double hi=min(1.,(x+xi)/(1.+xi));
    if(low>hi) return TransGPD_set(0.,0.,0.,0.);
    
    //  cout << "bounds " << x << " " << low << " " << hi << endl;
    if(low<1.E-09) low=1.E-09;
    if((1.-hi)<1.E-05) hi=1.-1.E-05;
    numint::array<double,1> lower = {{low}};
    numint::array<double,1> upper = {{hi}};

    unsigned count;
    numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,6E04,ret,count,0);

    return TransGPD_set(ret[0],ret[1],ret[2],ret[3]);

    
}

void GPD::DD_int_rho(numint::vector_d & res, double rho, double x, double xi, double t, GPD &gpdobject){

    res=numint::vector_d(4,0.);
    double eta=(x-rho)/xi;
    double temp=3./4.*((1.-rho)*(1.-rho)-eta*eta)/pow(1.-rho,3.)/xi;
    double HTdfront,HTufront, EbarTdfront, EbarTufront;
    gpdobject.getHTfront(rho,HTdfront,HTufront);
    gpdobject.getEbarTfront(rho,EbarTdfront, EbarTufront);
    
    double exp1=exp(-0.45*log(rho)*t*1.E-06);
    double exp2=exp1*exp(0.5*t*1.E-06);
    // cout << rho << " " << x<< " "<< xi << " " << t << " " << exp1 << " " << exp2 << " " << temp << " " << HTdfront << " " << HTufront << " " << EbarTdfront << " " << EbarTufront << endl;
    res[0]=exp1*temp*HTdfront;
    res[1]=exp1*temp*HTufront;
    res[2]=exp2*temp*EbarTdfront;
    res[3]=exp2*temp*EbarTufront;
    // cout << rho << " " << x << " " << xi << " " << t << " " << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << endl;    

}

void GPD::getHTfront(const double x, double &HTdfront, double &HTufront) const{
    mstw->update(x,2.);


    double d=mstw->parton(1,x,2.)/x;
    double dbar=mstw->parton(-1,x,2.)/x;
    double u=mstw->parton(2,x,2.)/x;
    double ubar=mstw->parton(-2,x,2.)/x;

    double Dd=-0.7*sqrt(x)*d;
    double Ddbar=-0.3*pow(x,0.4)*dbar;
    double Du=Dd/d/-0.7*u;
    double Dubar=Ddbar/dbar*ubar;


    HTdfront=-1.01*sqrt(x)*(1.-x)*(d-dbar+Dd-Ddbar);
    HTufront=0.78*sqrt(x)*(1.-x)*(u-ubar+Du-Dubar);
    return;
}

void GPD::getEbarTfront(const double x, double &EbarTdfront, double &EbarTufront) const{
    EbarTdfront=5.05*pow(x,-0.3)*pow(1.-x,5.);
    EbarTufront=EbarTdfront/(1.-x)/5.05*6.83;
    return;
}