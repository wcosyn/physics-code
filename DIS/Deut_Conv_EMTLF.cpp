#include "Deut_Conv_EMTLF.hpp"
#include <constants.hpp>
#include <cassert>
#include <numint/numint.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <AuxFunction.hpp>

using namespace std;

Deut_Conv_EMTLF::Deut_Conv_EMTLF(const string &wfname){


    wfref = TDeuteron::Wavefunction::CreateWavefunction(wfname);
    for(int i=0;i<=1000;i++){
        wf.AddUp(i,wfref->GetUp(i));
        wf.AddWp(i,wfref->GetWp(i));
        wf.AddUr(i*0.02,wfref->GetUr(i*0.02));
        wf.AddWr(i*0.02,wfref->GetWr(i*0.02));        
    }
}

Deut_Conv_EMTLF::~Deut_Conv_EMTLF(){
      delete wfref;
}




void Deut_Conv_EMTLF::int_k3(numint::vector_z & res, double alpha1, double kperp, double kphi, Deut_Conv_EMTLF &emt,
              double t, int pold_in, int pold_out, double deltax, double GFF_A, double GFF_J, double GFF_D){
    res=numint::vector_z(3,0.);


    //All momenta/energies MeV!!!!
    double Ek=sqrt((MASSn*MASSn+kperp*kperp)/alpha1/(2.-alpha1)); 
    double kz=(alpha1-1.)*Ek;
    double knorm=sqrt(Ek*Ek-MASSn*MASSn);
    // cout << alpha1 << " " << kperp << " " << Ek << " " << knorm << endl;
    if(knorm>1.E03) {return;}

    // double Delta_perp=sqrt((t0-t)*(1-xi*xi)); //symmetric frame where \bar{P}^perp=0

    double alphaprime=alpha1;

    
    double sinkphi,coskphi;
    sincos(kphi,&sinkphi,&coskphi);
    double kxprime=kperp*coskphi+(1-alpha1/2.)*deltax;  //see Eq (114) GPD note for kinematics
    double kyprime=kperp*sinkphi;

    double kperpprime=sqrt(kxprime*kxprime+kyprime*kyprime);
    double Ekprime = sqrt((MASSn*MASSn+kperpprime*kperpprime)/(alphaprime*(2.-alphaprime)));
    double kzprime = (alphaprime -1)*Ekprime;
    TVector3 k_in(kperp*coskphi,kperp*sinkphi,kz);
    TVector3 k_out(kxprime, kyprime, kzprime);
    complex<double> p1R = (kperp*coskphi+kxprime)/2.+I_UNIT*(kperp*sinkphi+kyprime)/2.;
    complex<double> p1L = conj(p1R);
    if(k_out.Mag2()>1.E06) {return;}

    
    //atan2(2.*xi_n*kyprime,2*xi_n*kperp*sinkphi+(1.+xi_n*(1.-alpha1/2.))*deltax) alternative expression
    // cout << "angle " << 2*xi_n*kxprime+(1.-xi_n*(1.-alphaprime/2.)*deltax) << " " << 2*xi_n*kperp*sinkphi+(1.+xi_n*(1.-alpha1/2.)*deltax) << endl;
    // cout << "vec " << alpha1 << " " << kperp << " " << kphi << " " << alphaprime << " " << kperpprime << " " << atan2(kyprime,kxprime) << " " << xi_n << " " << x_n << " " << phin << endl;
    // cout << kperp*sinkphi << " " << kxprime << " " << (1-alphaprime/2.)/(1+xi)*deltax << " " << (1-alpha1/2.)/(1-xi)*deltax << " " << kperpprime*cos(atan2(kyprime,kxprime)) << endl;
    // cout << kxprime << " " << kperpprime*cos(atan2(kyprime,kxprime)) << " " << kyprime << " " << kperpprime*sin(atan2(kyprime,kxprime)) << endl;
    // exit(1);

    // double px=kperp*coskphi-alpha1/2.*deltax/2.;
    // double pxprime=kxprime+alphaprime/2.*deltax/2.;
    // double pbarx=0.5*(px+pxprime);
    // cout << "phi " << phin << " " << atan2(2.*xi_n*kyprime,2.*xi_n*pbarx+deltax) << " " << atan2(kyprime, kperp*coskphi+0.25*alpha1/xi*deltax) << endl;


    vector< complex<double> > result(3,0.);
    //Melosh rotation active nucleon
    vector<complex<double> > wf_in(4,0.), wf_out(4,0.);
    wf_in.at(0)=emt.getWf()->DeuteronPState(2*pold_in, -1, -1, k_in);
    wf_out.at(0)=emt.getWf()->DeuteronPState(2*pold_out, -1, -1, k_out);
    wf_in.at(1)=emt.getWf()->DeuteronPState(2*pold_in, -1, 1, k_in);
    wf_out.at(1)=emt.getWf()->DeuteronPState(2*pold_out, -1, 1, k_out);
    wf_in.at(2)=emt.getWf()->DeuteronPState(2*pold_in, 1, -1, k_in);
    wf_out.at(2)=emt.getWf()->DeuteronPState(2*pold_out, 1, -1, k_out);
    wf_in.at(3)=emt.getWf()->DeuteronPState(2*pold_in, 1, 1, k_in);
    wf_out.at(3)=emt.getWf()->DeuteronPState(2*pold_out, 1, 1, k_out);
    wf_in=lf_deut(Ek,k_in,wf_in);
    wf_out=lf_deut(Ekprime,k_out,wf_out);
    
    
    for(int sigma2=-1;sigma2<=1;sigma2+=2){  //spectator helicity
        for(int sigma1in=-1; sigma1in<=1; sigma1in+=2){ //initial active nucleon helicity
            for(int sigma1out=-1; sigma1out<=1; sigma1out+=2){ //final active nucleon helicity
                vector< complex<double>> temp = getGFF_nucl(sigma1in, sigma1out, alpha1, p1R, p1L, t, GFF_A, GFF_J, GFF_D);
                result[0]+=wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))*temp[0];
                result[1]+=wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))*temp[1];
                result[2]+=wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))*temp[2];
            }
        }
    }
    // cout << "int " << x << " " << result << " " << kperp << " " << alpha1 << endl;
    res[0]=result[0]*4.*kperp/alpha1/(2.-alpha1)/alphaprime*sqrt(Ek)*sqrt(Ekprime);    
    res[1]=result[1]*4.*kperp/alpha1/(2.-alpha1)/alphaprime*sqrt(Ek)*sqrt(Ekprime);    
    res[2]=result[2]*4.*kperp/alpha1/(2.-alpha1)/alphaprime*sqrt(Ek)*sqrt(Ekprime);    
    // cout << t0 << " " << alpha1 << " " << t0_n << " " << xi << " " << xi_n << " " << res[0].real() << endl;
    // cout << alpha1 << " " << kperp << " " << kphi << " " << alphaprime << " " << kperpprime << " " << atan2(kyprime,kxprime) << " " << xi_n << " " << x_n << " " << phin << " " << res[0].real() << " " << res[0].imag() << " "<<  abs(res[0]) << endl;
    // exit(1);
}

void Deut_Conv_EMTLF::int_k3_real(numint::vector_d & res, double alpha1, double kperp, double kphi, Deut_Conv_EMTLF &emt,
              double t, double deltax, double GFF_A, double GFF_J, double GFF_D){
    res=numint::vector_d(15,0.);


    //All momenta/energies MeV!!!!
    double Ek=sqrt((MASSn*MASSn+kperp*kperp)/alpha1/(2.-alpha1)); 
    double kz=(alpha1-1.)*Ek;
    double knorm=sqrt(Ek*Ek-MASSn*MASSn);
    // cout << alpha1 << " " << kperp << " " << Ek << " " << knorm << endl;
    if(knorm>1.E03) {return;}

    // double Delta_perp=sqrt((t0-t)*(1-xi*xi)); //symmetric frame where \bar{P}^perp=0

    double alphaprime=alpha1;

    
    double sinkphi,coskphi;
    sincos(kphi,&sinkphi,&coskphi);
    double kxprime=kperp*coskphi+(1-alpha1/2.)*deltax;  //see Eq (114) GPD note for kinematics
    double kyprime=kperp*sinkphi;

    double kperpprime=sqrt(kxprime*kxprime+kyprime*kyprime);
    double Ekprime = sqrt((MASSn*MASSn+kperpprime*kperpprime)/(alphaprime*(2.-alphaprime)));
    double kzprime = (alphaprime -1)*Ekprime;
    TVector3 k_in(kperp*coskphi,kperp*sinkphi,kz);
    TVector3 k_out(kxprime, kyprime, kzprime);
    complex<double> p1R = (kperp*coskphi+kxprime)/2.+I_UNIT*(kperp*sinkphi+kyprime)/2.;
    complex<double> p1L = conj(p1R);
    if(k_out.Mag2()>1.E06) {return;}

    for(int i=0;i<5;i++){
        int pold_in=(i+4)/3-1;
        int pold_out=(i+4)%3-1;
        //cout << i << " " << pold_in << " " << pold_out << " " << i*3 << endl;
        vector< double > result(3,0.);
        //Melosh rotation active nucleon
        vector<complex<double> > wf_in(4,0.), wf_out(4,0.);
        wf_in.at(0)=emt.getWf()->DeuteronPState(2*pold_in, -1, -1, k_in);
        wf_out.at(0)=emt.getWf()->DeuteronPState(2*pold_out, -1, -1, k_out);
        wf_in.at(1)=emt.getWf()->DeuteronPState(2*pold_in, -1, 1, k_in);
        wf_out.at(1)=emt.getWf()->DeuteronPState(2*pold_out, -1, 1, k_out);
        wf_in.at(2)=emt.getWf()->DeuteronPState(2*pold_in, 1, -1, k_in);
        wf_out.at(2)=emt.getWf()->DeuteronPState(2*pold_out, 1, -1, k_out);
        wf_in.at(3)=emt.getWf()->DeuteronPState(2*pold_in, 1, 1, k_in);
        wf_out.at(3)=emt.getWf()->DeuteronPState(2*pold_out, 1, 1, k_out);
        wf_in=lf_deut(Ek,k_in,wf_in);
        wf_out=lf_deut(Ekprime,k_out,wf_out);
        
        
        for(int sigma2=-1;sigma2<=1;sigma2+=2){  //spectator helicity
            for(int sigma1in=-1; sigma1in<=1; sigma1in+=2){ //initial active nucleon helicity
                for(int sigma1out=-1; sigma1out<=1; sigma1out+=2){ //final active nucleon helicity
                    vector< complex<double>> temp = getGFF_nucl(sigma1in, sigma1out, alpha1, p1R, p1L, t, GFF_A, GFF_J, GFF_D);
                    result[0]+=real(wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))*temp[0]);
                    result[1]+=real(wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))*temp[1]);
                    result[2]+=real(wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))*temp[2]);
                }
            }
        }
        // cout << "int " << x << " " << result << " " << kperp << " " << alpha1 << endl;
        res[i]=result[0]*4.*kperp/alpha1/(2.-alpha1)/alphaprime*sqrt(Ek)*sqrt(Ekprime);    
        res[5+i]=result[1]*4.*kperp/alpha1/(2.-alpha1)/alphaprime*sqrt(Ek)*sqrt(Ekprime);    
        res[10+i]=result[2]*4.*kperp/alpha1/(2.-alpha1)/alphaprime*sqrt(Ek)*sqrt(Ekprime);    
        // cout << t0 << " " << alpha1 << " " << t0_n << " " << xi << " " << xi_n << " " << res[0].real() << endl;
        // cout << alpha1 << " " << kperp << " " << kphi << " " << alphaprime << " " << kperpprime << " " << atan2(kyprime,kxprime) << " " << xi_n << " " << x_n << " " << phin << " " << res[0].real() << " " << res[0].imag() << " "<<  abs(res[0]) << endl;
        // exit(1);
    }
    //exit(1);
}

vector< complex<double> > Deut_Conv_EMTLF::getGFF_nucl(const int sigma_in, const int sigma_out, 
                        const double alpha_1, const complex<double> p1R, const complex<double> p1L, const double t,
                        const double GFF_A, const double GFF_J, const double GFF_D){

    vector< complex<double> > res(3,0.);

    
    // cout << "nucl " << xi << " " << gpd_nucl.getHtildeT_singlet(model) << " " << gpd_nucl.getET_singlet(model) << endl;
    if(sigma_in==sigma_out) {
        res[0] = alpha_1*alpha_1/4.*GFF_A;
        //res[1] = alpha_1/2.*(-2.*p1R/sqrt(-t)*GFF_A+sigma_in*GFF_J); //+R
        //res[1] = alpha_1/2.*(2.*p1L/sqrt(-t)*GFF_A+sigma_in*GFF_J); //+L
        res[1] = alpha_1/2.*(-(p1R-p1L)/sqrt(-t)*GFF_A+sigma_in*GFF_J); // +R +L average
        //res[2] = 4.*p1R*p1L*GFF_A+t*GFF_D+2.*sqrt(-t)*sigma_in*(p1R-p1L)*GFF_J; //RL
        res[2] = -4.*p1R*p1R*GFF_A+t*GFF_D+4.*sqrt(-t)*sigma_in*p1R*GFF_J;  //RR
        //res[2] = -4.*p1L*p1L*GFF_A+t*GFF_D-4.*sqrt(-t)*sigma_in*p1L*GFF_J; //LL
    }

    else{
        res[0] = (sigma_in==1?-1.:1.)*sqrt(-t)/(2.*MASSn)*pow(alpha_1/2.,2.)*(GFF_A-2.*GFF_J);
        //res[1] = (sigma_in==1?1.:-1.)*p1R/MASSn*alpha_1/2.*(GFF_A+(sigma_in==1?-2:(p1L/p1R-1.))*GFF_J); //+R
        //res[1] = (sigma_in==1?-1.:+1.)*p1L/MASSn*alpha_1/2.*(GFF_A+(sigma_in==-1?-2:(p1R/p1L-1.))*GFF_J); //+L
        res[1] = (sigma_in==1?1.:-1.)*p1R/MASSn*alpha_1/4.*(GFF_A+(sigma_in==1?-2:(p1L/p1R-1.))*GFF_J)
                + (sigma_in==1?-1.:+1.)*p1L/MASSn*alpha_1/4.*(GFF_A+(sigma_in==-1?-2:(p1R/p1L-1.))*GFF_J); //+L +R average
        //res[2] = (sigma_in==1?-1.:1.)*sqrt(-t)/(2.*MASSn)*(4.*p1R*p1L*GFF_A+t*GFF_D)-2.*(p1R-p1L)*sqrt(-t)/MASSn*GFF_J*(sigma_in==1?p1R:p1L); //RL
        res[2] = (sigma_in==1?-1.:1.)*sqrt(-t)/(2.*MASSn)*(-4.*p1R*p1R*GFF_A+t*GFF_D)-4.*sqrt(-t)*p1R/MASSn*GFF_J*(sigma_in==1?p1R:p1L); //RR
        //res[2] = (sigma_in==1?-1.:1.)*sqrt(-t)/(2.*MASSn)*(-4.*p1L*p1L*GFF_A+t*GFF_D)+4.*sqrt(-t)*p1L/MASSn*GFF_J*(sigma_in==1?p1R:p1L); //LL
    }
    return res;
}




vector< complex<double> > Deut_Conv_EMTLF::EMT_conv(const double t){
    double Delta_perp=sqrt(-t); //symmetric frame where \bar{P}^perp=0  //GeV

     double A = 1270.*1270/(1270*1270-t)*1430*1430/(1430*1430-t);
    double J = 0.5*1270*1270/(1270*1270-t)*1430*1430/(1430*1430-t);
    //double D = -2*1270*1270/(1270*1270-t)*1430*1430/(1430*1430-t)*800*800/(800*800-t); //used in original PRD
    double D = -1.275*pow(1.+t/(0.963*0.963),-2.) -1.30*pow(1.+t/(0.81*0.81),-2.);  //He/Zahed parametrization

    Deut_Conv_EMTLF::Ftor_conv F;
    F.t=t*1.E06; //[MeV^2]
    F.emt=this;
    F.GFF_A=A; 
    F.GFF_J=J;
    F.GFF_D=D;

    numint::mdfunction<numint::vector_z,3> mdf;
    mdf.func = &Ftor_conv::exec;
    mdf.param = &F;
    numint::vector_z out(3,0.);
    
    numint::array<double,3> lower = {{0.,0.,-PI}};
    numint::array<double,3> upper = {{2.,1.E03,PI}};

    double abserr=1.E-15;
    double relerr=1.E-03;

    unsigned count;
    vector< complex<double> > ret(9,0.);
    //we only need 5 Helicity amplitudes
    for(int i=0;i<9;i++){
        F.f=Deut_Conv_EMTLF::int_k3;
        F.pold_in=i/3-1;
        F.pold_out=i%3-1;
        F.deltax=Delta_perp*1.E03; //[MeV]

        int minEval=1E02;
        int maxEvalcuba=1E05;

        int maxEval=1E07;
        numint::cube_adaptive(mdf,lower,upper,abserr,relerr,5E02,maxEval,out,count,0);
        ret[i]=out[0];
        cout << F.pold_in << " " << F.pold_out << " " << out[0] << " " << out[1] << " " << out[2]*1.E-06 << " " << count << endl;
     }

    vector< complex<double> > result(5,0.);
    //order of helicity amplitudes is ++,00,0+,+0,-+
    result[0]=ret[8];result[1]=ret[4];result[2]=ret[7];result[3]=ret[5];result[4]=ret[6];
    return result;

}


vector< double > Deut_Conv_EMTLF::EMT_conv_real(const double t){
    double Delta_perp=sqrt(-t); //symmetric frame where \bar{P}^perp=0  //GeV

       cout << t << endl;

    double A = 1.270*1.270/(1.270*1.270-t)*1.430*1.430/(1.430*1.430-t);
    double J = 0.5*1.270*1.270/(1.270*1.270-t)*1.430*1.430/(1.430*1.430-t);
    double D = -2*1.270*1.270/(1.270*1.270-t)*1.430*1.430/(1.430*1.430-t)*.800*.800/(.800*.800-t); //used in original PRD
    //double D = -1.275*pow(1.+t/(0.963*0.963),-2.) -1.30*pow(1.+t/(0.81*0.81),-2.);  //He/Zahed parametrization

    Deut_Conv_EMTLF::Ftor_conv_real F;
    F.t=t*1.E06; //[MeV^2]
    F.emt=this;
    F.GFF_A=A; 
    F.GFF_J=J;
    F.GFF_D=D;

    numint::mdfunction<numint::vector_d,3> mdf;
    mdf.func = &Ftor_conv_real::exec;
    mdf.param = &F;
    numint::vector_d out(15,0.);
    
    numint::array<double,3> lower = {{0.,0.,-PI}};
    numint::array<double,3> upper = {{2.,1.E03,PI}};

    double abserr=1.E-15;
    double relerr=1.E-05;

    unsigned count;
    // vector< double > ret(15,0.);
    //we only need 5 Helicity amplitudes
    F.f=Deut_Conv_EMTLF::int_k3_real;
    F.deltax=Delta_perp*1.E03; //[MeV]

    int minEval=1E02;
    int maxEvalcuba=1E05;

    int maxEval=1E07;
    numint::cube_adaptive(mdf,lower,upper,abserr,relerr,5E02,maxEval,out,count,0);
    cout << t << " ";
    for(int i=0;i<15;i++) cout << out[i]*(i>9? 1.E-06:1.) << " ";
    cout << endl;
    //cout << F.pold_in << " " << F.pold_out << " " << out[0] << " " << out[1] << " " << out[2]*1.E-06 << endl;

    vector< double > result(15,0.);
    //order of helicity amplitudes is ++,00,0+,+0,-+
    //result[0]=ret[8];result[1]=ret[4];result[2]=ret[7];result[3]=ret[5];result[4]=ret[6];
    return result;

}


vector< complex<double> > Deut_Conv_EMTLF::lf_deut(const double Ek, const TVector3& k, const vector< complex<double> > &nonrelwf){
    vector< complex<double> > wf_out(4,0.);
    wf_out[0]=((k[2]+Ek+MASSn)*nonrelwf.at(0)+(-k[0]-I_UNIT*k[1])*nonrelwf.at(1))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    wf_out[1]=((k[2]+Ek+MASSn)*nonrelwf.at(1)+(k[0]-I_UNIT*k[1])*nonrelwf.at(0))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    wf_out[2]=((k[2]+Ek+MASSn)*nonrelwf.at(2)+(-k[0]-I_UNIT*k[1])*nonrelwf.at(3))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    wf_out[3]=((k[2]+Ek+MASSn)*nonrelwf.at(3)+(k[0]-I_UNIT*k[1])*nonrelwf.at(2))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    
    return wf_out;

}
