#include "GPD.hpp"
#include <constants.hpp>
#include <mstwpdf.h>
#include <cassert>
#include <numint/numint.hpp>
#include <cmath> 

using namespace std;

GPD::GPD(const std::string &pdf_name, const std::string &wfname){

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

    wfref = TDeuteron::Wavefunction::CreateWavefunction(wfname);
    for(int i=0;i<=1000;i++){
        wf.AddUp(i,wfref->GetUp(i));
        wf.AddWp(i,wfref->GetWp(i));
    }

  

}

GPD::~GPD(){
      if(mstw) delete mstw;
      delete wfref;
}


vector< complex<double> > GPD::helamps_to_gpds(const double xi, const double t, const vector< complex <double> > & helamps){
    
    vector< complex<double> > gpds(9,0.);
    //symmetric frame, with momentum transfer in x,z plane, so phi=0
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



void GPD::int_k3(numint::vector_z & res, double alpha1, double kperp, double kphi, GPD &gpd,
              double x, double xi, double t, int pold_in, int pold_out, int model){
    res=numint::vector_z(1,0.);
    double Ek=sqrt((MASSn*MASSn+kperp*kperp)/alpha1/(2.-alpha1));
    double kz=sqrt(Ek*Ek-MASSn*MASSn-kperp*kperp);
    double knorm=sqrt(Ek*Ek-MASSn*MASSn);
    // cout << alpha1 << " " << kperp << " " << Ek << " " << knorm << endl;
    if(knorm>1.E03) {return;}
    
    double t0=-4.*MASSn*MASSn*xi*xi/(1-xi*xi);
    double Delta_perp=sqrt((t0-t)*(1-xi*xi)); //symmetric frame where \bar{P}^perp=0

    double alphaprime=(alpha1*(1+xi)-4.*xi)/(1-xi);

    double xi_n=xi/(alpha1/2.*(1+xi)-xi);

    double x_n=x/(alpha1/2.*(1+xi)-xi);
    // cout << "int " << x << " " << xi << " " << x_n << " " << xi_n << " " << xi/x << " " << xi_n/x_n << endl;
    double kcosth=kz/knorm;
    double ksinth=sqrt(1-kcosth*kcosth);
    double sinkphi,coskphi;
    sincos(kphi,&sinkphi,&coskphi);
    double kxprime=kperp*sinkphi+(1-alpha1/2.)/(1-xi)*Delta_perp;  //see (114) GPD note for kinematics
    double kyprime=kperp*coskphi;
    
    double kperpprime=sqrt(kxprime*kxprime+kyprime*kyprime);
    double Ekprime = sqrt((MASSn*MASSn+kperpprime*kperpprime)/(alphaprime*(2.-alphaprime)));
    double kzprime = (alphaprime -1)*Ekprime;
    TVector3 k_in(kperp*coskphi,kperp*sinkphi,kz);
    TVector3 k_out(kxprime, kyprime, kzprime);
    if(k_out.Mag2()>1.E06) {return;}

    double phin=atan2(2*xi_n*kyprime, 2*xi_n*kxprime+(1.-xi_n*(1.-alphaprime/2.)*Delta_perp)); //azimuthal angle of 2xi_n*\bar{p}_1+Delta
    //atan2(2.*xi_n*kyprime,2*xi_n*kperp*sinkphi+(1.+xi_n*(1.-alpha1/2.)*Delta_perp)) alternative expression

    
    complex<double> result=0.;
    for(int sigma2=-1;sigma2<=1;sigma2+=2){
        for(int sigma1in=-1; sigma1in<=1; sigma1in+=2){
            complex<double> wfin = gpd.getWf()->DeuteronPState(2*pold_in, sigma2, sigma1in, k_in);
            for(int sigma1out=-1; sigma1out<=1; sigma1out+=2){
                complex<double> wfout=gpd.getWf()->DeuteronPState(2*pold_out, sigma2, sigma1out, k_out);
                result+=wfin*conj(wfout)*gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0,phin,model); 
                // cout << wfin << " " << wfout << " " << gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0,phin,model) << " " << x_n << endl;
            }
        }
    }
    // cout << "int " << x << " " << result << " " << kperp << " " << alpha1 << endl;
    res[0]=result*2.*kperp/alpha1/(2.-alpha1);
}

complex<double> GPD::getGPD_odd_nucl(int sigma_in, int sigma_out, double x, double xi, double t, double t0, double phi, int model){


    TransGPD_set gpd_nucl=getGK_param(x,xi,t);
    double kinfac=sqrt(t0-t)/MASSn;
    if(sigma_in==sigma_out) return (kinfac*gpd_nucl.getHtildeT_singlet(model)
                            +(1-sigma_in*xi)*kinfac/2.*gpd_nucl.getET_singlet(model)
                            +sigma_in*(1-sigma_in*xi)*kinfac/2.*gpd_nucl.getEtildeT_singlet())*exp(I_UNIT*phi);

    else return -(sigma_in-1)*sqrt(1.-xi*xi)*gpd_nucl.getHT_singlet()
                -sigma_in/2.*kinfac*kinfac*sqrt(1-xi*xi)*exp((sigma_in+1)*phi*I_UNIT)*gpd_nucl.getHtildeT_singlet(model)
                +(sigma_in-1)*xi*xi/sqrt(1-xi*xi)*gpd_nucl.getET_singlet(model)
                -(sigma_in-1)*xi/sqrt(1.-xi*xi)*gpd_nucl.getEtildeT_singlet();

}


vector< complex<double> > GPD::gpd_conv(const double xi, const double x, const double t,const int model){
    if(fabs(x)>1.) return vector< complex<double> >(9,0.);
    GPD::Ftor_conv F;
    F.x=x;
    F.xi=xi;
    F.t=t;
    F.gpd=this;
    F.model=model;

    numint::mdfunction<numint::vector_z,3> mdf;
    mdf.func = &Ftor_conv::exec;
    mdf.param = &F;
    numint::vector_z out(1,0.);
    F.f=GPD::int_k3;
    double low=0.;
    if(fabs(x)>xi) low=2.*(fabs(x)+xi)/(1+xi);
    else low=4.*xi/(1+xi);
    cout << "low " << low << endl;
    numint::array<double,3> lower = {{low,0.,0.}};
    numint::array<double,3> upper = {{2.,1.E03,2.*PI}};

    unsigned count;
    vector< complex<double> > ret(9,0.);
    for(int i=0;i<9;i++){
        F.pold_in=i/3-1;
        F.pold_out=i%3-1;
        numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,6E04,out,count,0);
        ret[i]=out[0];
        cout << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << out[0] << " " << count << endl;
    }

    vector< complex<double> > result(9,0.);
    //order of helicity amplitudes is ++,--,00,0+,-0,+0,0-,-+,+-
    result[0]=ret[8];result[1]=ret[0];result[2]=ret[4];result[3]=ret[5];result[4]=ret[1];result[5]=ret[7];result[6]=ret[3];result[7]=ret[2];result[8]=ret[6];
    return result;

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
    // cout << "c " << ret[0] << " " << ret[1] << " " << ret[2] << " " << ret[3] << " " << count << endl;
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

