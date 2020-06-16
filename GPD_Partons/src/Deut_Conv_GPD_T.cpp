#include "Deut_Conv_GPD_T.hpp"
#include <constants.hpp>
#include <mstwpdf.h>
#include <cassert>
#include <numint/numint.hpp>

using namespace std;

Deut_Conv_GPD_T::Deut_Conv_GPD_T(const string &pdf_name, const string &wfname):
chiralodd_grid(pdf_name){


    wfref = TDeuteron::Wavefunction::CreateWavefunction(wfname);
    for(int i=0;i<=1000;i++){
        wf.AddUp(i,wfref->GetUp(i));
        wf.AddWp(i,wfref->GetWp(i));
    }
    grid = NULL;
    grid_set = false;
    ERBL_set = 0;
    grid_size=0;
    model_set=0;

}

Deut_Conv_GPD_T::~Deut_Conv_GPD_T(){
      delete wfref;
      if(grid != NULL) delete [] grid;
}


vector< complex<double> > Deut_Conv_GPD_T::helamps_to_gpds_T(const double xi, const double t, const vector< complex <double> > & helamps){

    vector< complex<double> > gpds(9,0.);
    //symmetric frame, with momentum transfer in x,z plane, so phi=0
    double D=sqrt((-4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi)-t)*(1-xi*xi))/(2.*MASSD_G);
    gpds.at(0)=2.*sqrt(2.)*xi/D/(1-xi*xi)*(helamps[0]-helamps[1])+2.*(1./(1.-xi)*helamps[3]+1./(1.+xi)*helamps[4])+2./(1+xi)*helamps[5]
                +2./(1-xi)*helamps[6]+2.*sqrt(2.)*D/(1.-xi*xi)*(helamps[7]-helamps[8]);
    gpds.at(1)=sqrt(2.)*(2.*D*D/(1-xi)-xi)/D/pow(1+xi,2.)*helamps[0]
                -sqrt(2.)*(2.*D*D/(1+xi)+xi)/D/pow(1-xi,2.)*helamps[1]
                +2.*sqrt(2.)*xi/(1-xi*xi)/D*helamps[2]
                -2./(1-xi*xi)*(helamps[3]-helamps[4])
                +2.*(1./pow(1+xi,2.)+2.*xi*xi/(D*D*(1-xi*xi)*(1+xi)))*helamps[5]
                -2.*(1./pow(1-xi,2.)+2.*xi*xi/(D*D*(1-xi*xi)*(1-xi)))*helamps[6]
                +sqrt(2.)*xi/(D*(1-xi*xi))*helamps[7]
                -sqrt(2)*xi/(D*pow(D*(1-xi*xi),2.))*(4.*xi*xi+D*D*(3.+xi*xi))*helamps[8];

    gpds.at(2)=-(1./(2.*D)*((1-xi)/(1+xi)*helamps[0]-(1+xi)/(1-xi)*helamps[1])+1./(sqrt(2.)*D*D)*((1-xi)/(1+xi)*helamps[5]-(1+xi)/(1-xi)*helamps[6])
                -2.*xi/(D*D*D*(1-xi*xi))*helamps[8]);
    gpds.at(3)=1/D*(1./(1+xi)*helamps[0]+1./(1.-xi)*helamps[1])+1/(sqrt(2.)*D*D)*((1-xi)/(1+xi)*helamps[5]+(1+xi)/(1-xi)*helamps[6])+1./(2.*D)*helamps[7]
                -1./(2.*D*D*D)*(D*D-4.*xi*xi/(1-xi*xi))*helamps[8];

    gpds.at(4)=sqrt(2.)*(-1./D*((1./2.+D*D/(1-xi))/(1+xi)/(1+xi)*helamps[0]+(1./2.+D*D/(1+xi))/(1-xi)/(1-xi)*helamps[1])+1./(D*(1.-xi*xi))*helamps[2]
                    +1/sqrt(2.)/(1-xi*xi)*(helamps[3]+helamps[4])-1./(sqrt(2.)*pow(1+xi,2.))*(1.-2.*xi/(D*D*(1-xi)))*helamps[5]
                    -1/(sqrt(2.)*pow(1-xi,2.))*(1.+2.*xi/(D*D*(1.+xi)))*helamps[6]-1./(2.*D*(1-xi*xi))*helamps[7]
                    -(D*D*(3.*xi*xi+1)+4.*xi*xi)/(2.*D*D*D*pow(1-xi*xi,2.))*helamps[8]);

    gpds.at(5)=-1./(2.*D)*(helamps[0]+helamps[1]+helamps[7]+helamps[8]);
    gpds.at(6)=2.*(1-xi*xi)/(D*D*D)*helamps[8];
    gpds.at(7)=-1./D*((1-xi)/(1+xi)*helamps[0]-(1+xi)/(1-xi)*helamps[1]
                    -xi*sqrt(2.)/D*((1-xi)/(1+xi)*helamps[5]+(1+xi)/(1-xi)*helamps[6])
                    -4.*xi*xi*xi/D/D/(1-xi*xi)*helamps[8]);
    gpds.at(8)=-1./D*(2.*xi*(helamps[0]/(1+xi)-helamps[1]/(1-xi))
                        +sqrt(2.)*xi/D*((1-xi)/(1+xi)*helamps[5]-(1+xi)/(1-xi)*helamps[6])
                        +helamps[7]-(4.*xi*xi/(D*D*(1-xi*xi))+1.)*helamps[8]);

    return gpds;
}

vector< complex<double> > Deut_Conv_GPD_T::gpds_to_helamps_T(const double xi, const double t, const vector< complex<double> > & gpds){

    vector< complex<double> > helamps(9,0.);
    //symmetric frame, with momentum transfer in x,z plane, so phi=0

    double D=sqrt((-4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi)-t)*(1-xi*xi))/(2.*MASSD_G);
    double t0=-4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi);
    // cout << (t0-t)*1.E-06 << " " << MASSD_G*1.E-03 << endl;
    helamps.at(0)=-D*(xi/(1.-xi)*(gpds[2]-gpds[3])+gpds[5]+(t0-t)/(8.*MASSD_G*MASSD_G)*gpds[6]+1./(2.*(1.-xi))*(gpds[7]-gpds[8]));

    helamps.at(1)=-D*(xi/(1.+xi)*(gpds[2]+gpds[3])+gpds[5]+(t0-t)/(8.*MASSD_G*MASSD_G)*gpds[6]-1./(2.*(1.+xi))*(gpds[7]+gpds[8]));

    helamps.at(2)=D*(-1./(2.*sqrt(2.))*(gpds[0]+xi*gpds[1])+2*xi/(1.-xi*xi)*gpds[2]+2.*D*D/(1.-xi*xi)*gpds[3]
                    +gpds[4]/sqrt(2.)+((t0-t)/(2.*MASSD_G*MASSD_G)-(1+xi*xi)/(1-xi*xi))*gpds[5]+(pow((t0-t)/(4.*MASSD_G*MASSD_G),2.)-pow(xi/(1-xi*xi),2.))*gpds[6]-xi/(1.-xi*xi)*gpds[7]
                    -(t0-t)/(4.*MASSD_G*MASSD_G)*gpds[8]);

    helamps.at(3)=1./sqrt(2.)*((1-xi)/(2.*sqrt(2.))*(gpds[0]-gpds[1])+(2*xi*xi/(1.-xi*xi)-D*D/(1-xi))*gpds[2]+(-2*xi*xi/(1.-xi*xi)+D*D*(3.*xi-1.)/(1-xi*xi))*gpds[3]
                    +1/sqrt(2.)*xi*(1-xi)*gpds[4]-(1-xi)*(t0-t)/(2.*MASSD_G*MASSD_G)*gpds[5]
                    -(1-xi)*(t0-t)/(4.*MASSD_G*MASSD_G)*((t0-t)/(4.*MASSD_G*MASSD_G)-xi/(1.-xi*xi))*gpds[6]+(-(1.+xi)*(t0-t)/(8.*MASSD_G*MASSD_G)+xi/(1-xi*xi))*gpds[7]
                    +((3.-xi)*(t0-t)/(8.*MASSD_G*MASSD_G)-xi/(1-xi*xi))*gpds[8]);                    

    helamps.at(4)=1./sqrt(2.)*((1+xi)/(2.*sqrt(2.))*(gpds[0]+gpds[1])-(2*xi*xi/(1.-xi*xi)-D*D/(1+xi))*gpds[2]+(-2*xi*xi/(1.-xi*xi)-D*D*(3.*xi+1.)/(1-xi*xi))*gpds[3]
                    -1/sqrt(2.)*xi*(1+xi)*gpds[4]-(1+xi)*(t0-t)/(2.*MASSD_G*MASSD_G)*gpds[5]
                    -(1+xi)*(t0-t)/(4.*MASSD_G*MASSD_G)*((t0-t)/(4.*MASSD_G*MASSD_G)+xi/(1.-xi*xi))*gpds[6]+((1.-xi)*(t0-t)/(8.*MASSD_G*MASSD_G)+xi/(1-xi*xi))*gpds[7]
                    +((3.+xi)*(t0-t)/(8.*MASSD_G*MASSD_G)+xi/(1-xi*xi))*gpds[8]);                    

    helamps.at(5)=D*D/sqrt(2.)/(1-xi)*((-gpds[2]+gpds[3])+2.*gpds[5]+((t0-t)/(4.*MASSD_G*MASSD_G)+xi/(1-xi*xi))*gpds[6]+1./2.*(gpds[7]-gpds[8]));

    helamps.at(6)=D*D/sqrt(2.)/(1+xi)*((gpds[2]+gpds[3])+2.*gpds[5]+((t0-t)/(4.*MASSD_G*MASSD_G)-xi/(1-xi*xi))*gpds[6]-1./2.*(gpds[7]+gpds[8]));

    helamps.at(7)=D*(2.*xi/(1-xi*xi)*(gpds[2]-xi*gpds[3])+(t0-t)/(8.*MASSD_G*MASSD_G)*gpds[6]+1./(1-xi*xi)*(xi*gpds[7]-gpds[8]));

    helamps.at(8)=D*(t0-t)/(8.*MASSD_G*MASSD_G)*gpds[6];
    return helamps;
}



std::complex<double > Deut_Conv_GPD_T::test(double x, double xi, double t, int pold_in, int pold_out, int model, bool right, double deltax){
    double alpha1=1.+xi;
    double kperp=-deltax/4.;
    double kphi=0.;
    double Ek=sqrt((MASSn*MASSn+kperp*kperp)/alpha1/(2.-alpha1));
    double kz=(alpha1-1.)*Ek;
    double knorm=sqrt(Ek*Ek-MASSn*MASSn);
    if(knorm>1.E03) {return 0.;}


    double alphaprime=(alpha1*(1+xi)-4.*xi)/(1-xi);

    double xi_n=xi/(alpha1/2.*(1+xi)-xi);
    double t0_n=-4.*MASSn*MASSn*xi_n*xi_n/(1-xi_n*xi_n); //[MeV^2]
    // double t0=-4.*MASSD*MASSD*xi*xi/(1-xi*xi);
    if(t>t0_n) return 0.;
    
    double x_n=x/(alpha1/2.*(1+xi)-xi);
    if(abs(x_n)>1.) return 0.;
    double sinkphi,coskphi;
    sincos(kphi,&sinkphi,&coskphi);
    double kxprime=kperp*coskphi+(1-alpha1/2.)/(1-xi)*deltax;  //see (114) GPD note for kinematics
    double kyprime=kperp*sinkphi;

    double kperpprime=sqrt(kxprime*kxprime+kyprime*kyprime);
    // cout << alpha1 << "  " << alphaprime << " " << kperp << " " << kxprime << endl;
    double Ekprime = sqrt((MASSn*MASSn+kperpprime*kperpprime)/(alphaprime*(2.-alphaprime)));
    double kzprime = (alphaprime -1)*Ekprime;
    TVector3 k_in(kperp*coskphi,kperp*sinkphi,kz);
    TVector3 k_out(kxprime, kyprime, kzprime);
    if(k_out.Mag2()>1.E06) {return 0.;}

    double phin=atan2(2*xi_n*kperp*sinkphi, 2*xi_n*(kperp*coskphi+alpha1/(4.*xi)*deltax));

    complex<double> result=0.;
    //Melosh rotation active nucleon
    vector<complex<double> > wf_in(4,0.), wf_out(4,0.);
    wf_in.at(0)=getWf()->DeuteronPState(2*pold_in, -1, -1, k_in);
    wf_out.at(0)=getWf()->DeuteronPState(2*pold_out, -1, -1, k_out);
    wf_in.at(1)=getWf()->DeuteronPState(2*pold_in, -1, 1, k_in);
    wf_out.at(1)=getWf()->DeuteronPState(2*pold_out, -1, 1, k_out);
    wf_in.at(2)=getWf()->DeuteronPState(2*pold_in, 1, -1, k_in);
    wf_out.at(2)=getWf()->DeuteronPState(2*pold_out, 1, -1, k_out);
    wf_in.at(3)=getWf()->DeuteronPState(2*pold_in, 1, 1, k_in);
    wf_out.at(3)=getWf()->DeuteronPState(2*pold_out, 1, 1, k_out);
    // wf_in=lf_deut(Ek,k_in,wf_in);
    // wf_out=lf_deut(Ekprime,k_out,wf_out);
    // cout << "in " << pold_in << " " << wf_in.at(0) << " " << wf_in.at(1) << " " << wf_in.at(2) << " " << wf_in.at(3) << endl;
    // cout << "out " << pold_out << " " << wf_out.at(0) << " " << wf_out.at(1) << " " << wf_out.at(2) << " " << wf_out.at(3) << endl;
    

    TransGPD_set gpd_nucl=chiralodd_grid.getTransGPDSet(x_n,xi_n,t,2.);
    // cout << x_n << " " << xi_n << " " << 2*xi/(1+xi*xi) << " " << gpd_nucl.getET_singlet(model) << endl;
    for(int sigma2=-1;sigma2<=1;sigma2+=2){
        for(int sigma1in=-1; sigma1in<=1; sigma1in+=2){
            for(int sigma1out=-1; sigma1out<=1; sigma1out+=2){
                result+=wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))
                             *getGPD_odd_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,model,right,1,gpd_nucl);
                // cout << pold_out << " " << pold_in << " " << sigma2 << " " << sigma1in << " " << sigma1out << " " << wf_in.at((sigma1in+1)/2+(sigma2+1)) << " " << 
                // wf_out.at((sigma1out+1)/2+(sigma2+1)) << " " << getGPD_odd_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,model,right,gpd_nucl) << " " <<  endl;

            }
        }
    }
    // cout << pold_out << " " << pold_in << " " << result << endl;
    return result*2.;

}

void Deut_Conv_GPD_T::int_k3(numint::vector_z & res, double alpha1, double kperp, double kphi, Deut_Conv_GPD_T &gpd,
              double x, double xi, double t, double scale, int pold_in, int pold_out, int model, bool right, double deltax){
    res=numint::vector_z(1,0.);

    double Ek=sqrt((MASSn*MASSn+kperp*kperp)/alpha1/(2.-alpha1));
    double kz=(alpha1-1.)*Ek;
    double knorm=sqrt(Ek*Ek-MASSn*MASSn);
    // cout << alpha1 << " " << kperp << " " << Ek << " " << knorm << endl;
    if(knorm>1.E03) {return;}

    // double Delta_perp=sqrt((t0-t)*(1-xi*xi)); //symmetric frame where \bar{P}^perp=0

    double alphaprime=(alpha1*(1+xi)-4.*xi)/(1-xi);

    double xi_n=xi/(alpha1/2.*(1+xi)-xi);
    double t0_n=-4.*MASSn*MASSn*xi_n*xi_n/(1-xi_n*xi_n); //[MeV^2]
    // double t0=-4.*MASSD*MASSD*xi*xi/(1-xi*xi);
    if(t>t0_n) return;
    
    double x_n=x/(alpha1/2.*(1+xi)-xi);
    // cout << "int " << x << " " << xi << " " << x_n << " " << xi_n << " " << xi/x << " " << xi_n/x_n << endl;
    double sinkphi,coskphi;
    sincos(kphi,&sinkphi,&coskphi);
    double kxprime=kperp*coskphi+(1-alpha1/2.)/(1-xi)*deltax;  //see Eq (114) GPD note for kinematics
    double kyprime=kperp*sinkphi;

    double kperpprime=sqrt(kxprime*kxprime+kyprime*kyprime);
    double Ekprime = sqrt((MASSn*MASSn+kperpprime*kperpprime)/(alphaprime*(2.-alphaprime)));
    double kzprime = (alphaprime -1)*Ekprime;
    TVector3 k_in(kperp*coskphi,kperp*sinkphi,kz);
    TVector3 k_out(kxprime, kyprime, kzprime);
    if(k_out.Mag2()>1.E06) {return;}

    //double phin=atan2(2*xi_n*kyprime, 2*xi_n*kxprime+(1.-xi_n*(1.-alphaprime/2.))*deltax); //azimuthal angle of 2xi_n*\bar{p}_1+Delta
    double phin=atan2(2*xi_n*kperp*sinkphi, 2*xi_n*(kperp*coskphi+alpha1/(4.*xi)*deltax)); //azimuthal angle of 2xi_n*\bar{p}_1+Delta
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


    complex<double> result=0.;
    //Melosh rotation active nucleon
    vector<complex<double> > wf_in(4,0.), wf_out(4,0.);
    wf_in.at(0)=gpd.getWf()->DeuteronPState(2*pold_in, -1, -1, k_in);
    wf_out.at(0)=gpd.getWf()->DeuteronPState(2*pold_out, -1, -1, k_out);
    wf_in.at(1)=gpd.getWf()->DeuteronPState(2*pold_in, -1, 1, k_in);
    wf_out.at(1)=gpd.getWf()->DeuteronPState(2*pold_out, -1, 1, k_out);
    wf_in.at(2)=gpd.getWf()->DeuteronPState(2*pold_in, 1, -1, k_in);
    wf_out.at(2)=gpd.getWf()->DeuteronPState(2*pold_out, 1, -1, k_out);
    wf_in.at(3)=gpd.getWf()->DeuteronPState(2*pold_in, 1, 1, k_in);
    wf_out.at(3)=gpd.getWf()->DeuteronPState(2*pold_out, 1, 1, k_out);
    wf_in=lf_deut(Ek,k_in,wf_in);
    wf_out=lf_deut(Ekprime,k_out,wf_out);
    

    TransGPD_set gpd_nucl=gpd.chiralodd_grid.getTransGPDSet(x_n,xi_n,t,scale);
    for(int sigma2=-1;sigma2<=1;sigma2+=2){  //spectator helicity
        for(int sigma1in=-1; sigma1in<=1; sigma1in+=2){ //initial active nucleon helicity
            // complex<double> wfin = gpd.getWf()->DeuteronPState(2*pold_in, sigma2, sigma1in, k_in);
            for(int sigma1out=-1; sigma1out<=1; sigma1out+=2){ //final active nucleon helicity
                // complex<double> wfout=gpd.getWf()->DeuteronPState(2*pold_out, sigma2, sigma1out, k_out);
                // result+=wfin*conj(wfout)
                //             *gpd.getGPD_odd_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,model,right,gpd_nucl);
                result+=wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))
                            *gpd.getGPD_odd_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,model,right,1,gpd_nucl);
                // cout << wf_in.at((sigma1in+1)/2+(sigma2+1)) << " " << gpd.getWf()->DeuteronPState(2*pold_in, sigma2, sigma1in, k_in) << endl;
                // cout << wf_out.at((sigma1out+1)/2+(sigma2+1)) << " " << gpd.getWf()->DeuteronPState(2*pold_out, sigma2, sigma1out, k_out) << endl;
                // cout << wfin << " " << wfout << " " << gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0,phin,model,1) << " " << x_n << endl;
                //symmetry tests nucleon level
                // cout << "N conj" << sigma1in << " " << sigma1out << " " << conj(gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0_n,phin,model,1))
                //     << " " << -gpd.getGPD_odd_nucl(sigma1out, sigma1in, x_n, -xi_n, t, t0_n,phin+PI,model,0) << endl;
                // cout << "N p" << sigma1in << " " << sigma1out << " " << gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0_n,phin,model,1)
                //     << " " << -gpd.getGPD_odd_nucl(-sigma1in, -sigma1out, x_n, xi_n, t, t0_n,-phin+PI,model,0) << endl;
                // cout << "N t" << sigma1in << " " << sigma1out << " " << gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0_n,phin,model,1)
                //     << " " << pow(-1.,(sigma1in-sigma1out)/2)*gpd.getGPD_odd_nucl(sigma1out, sigma1in, x_n, -xi_n, t, t0_n,-phin,model,0) << " " <<  pow(-1.,(sigma1in-sigma1out)/2) << endl;

            }
        }
    }
    // cout << "int " << x << " " << result << " " << kperp << " " << alpha1 << endl;
    res[0]=result*2.*kperp/alpha1/(2.-alpha1)/alphaprime*sqrt(Ek)*sqrt(Ekprime);
    // cout << t0 << " " << alpha1 << " " << t0_n << " " << xi << " " << xi_n << " " << res[0].real() << endl;
    // cout << alpha1 << " " << kperp << " " << kphi << " " << alphaprime << " " << kperpprime << " " << atan2(kyprime,kxprime) << " " << xi_n << " " << x_n << " " << phin << " " << res[0].real() << " " << res[0].imag() << endl;
    // exit(1);
}

void Deut_Conv_GPD_T::int_kprime3(numint::vector_z & res, double alphaprime, double kperpprime, double kphiprime, Deut_Conv_GPD_T &gpd,
              double x, double xi, double t, double scale, int pold_in, int pold_out, int model, bool right, double deltax){
    res=numint::vector_z(1,0.);
    // alphaprime=1.0016;
    // kperpprime=21.6618;
    // kphiprime=0.0873254;


    double Ekprime=sqrt((MASSn*MASSn+kperpprime*kperpprime)/alphaprime/(2.-alphaprime));
    double kzprime=(alphaprime-1.)*Ekprime;
    double knormprime=sqrt(Ekprime*Ekprime-MASSn*MASSn);
    // cout << alpha1 << " " << kperp << " " << Ek << " " << knorm << endl;
    if(knormprime>1.E03) {return;}

    // double Delta_perp=sqrt((t0-t)*(1-xi*xi)); //symmetric frame where \bar{P}^perp=0

    double alpha1=(alphaprime*(1-xi)+4.*xi)/(1+xi);

    double xi_n=xi/(alphaprime/2.*(1-xi)+xi);
    double t0_n=-4.*MASSn*MASSn*xi_n*xi_n/(1-xi_n*xi_n);  //[MeV]
    // double t0=-4.*MASSD*MASSD*xi*xi/(1-xi*xi);
    if(t>t0_n) return;

    double x_n=x/(alphaprime/2.*(1-xi)+xi);
    // cout << "int " << x << " " << xi << " " << x_n << " " << xi_n << " " << xi/x << " " << xi_n/x_n << endl;
    double sinkphiprime,coskphiprime;
    sincos(kphiprime,&sinkphiprime,&coskphiprime);
    double kx=kperpprime*coskphiprime-(1-alphaprime/2.)/(1+xi)*deltax;  //see (114) GPD note for kinematics
    double ky=kperpprime*sinkphiprime;

    double kperp=sqrt(kx*kx+ky*ky);
    double Ek = sqrt((MASSn*MASSn+kperp*kperp)/(alpha1*(2.-alpha1)));
    double kz = (alpha1 -1.)*Ek;
    TVector3 k_in(kx,ky,kz);
    TVector3 k_out(kperpprime*coskphiprime, kperpprime*sinkphiprime, kzprime);
    if(k_in.Mag2()>1.E06) {return;}

    //double phin=atan2(2*xi_n*kperpprime*sinkphiprime, 2*xi_n*kperpprime*coskphiprime+(1.-xi_n*(1.-alphaprime/2.))*deltax); //azimuthal angle of 2xi_n*\bar{p}_1+Delta
    double phin=atan2(2*xi_n*ky, 2*xi_n*(kx+alpha1/(4.*xi)*deltax));
    // cout << "vec " << alpha1 << " " << kperp << " " << atan2(ky,kx) << " " << alphaprime << " " << kperpprime << " " << kphiprime << " " << xi_n << " " << x_n << " " << phin << endl;
    // cout << kx << " " << kperpprime*sinkphiprime << " " << (1-alphaprime/2.)/(1+xi)*deltax << " " << (1-alpha1/2.)/(1-xi)*deltax << endl;
    // exit(1);
    //atan2(2.*xi_n*kyprime,2*xi_n*kperp*sinkphi+(1.+xi_n*(1.-alpha1/2.)*deltax)) alternative expression
    // cout << "angle " << 2*xi_n*kxprime+(1.-xi_n*(1.-alphaprime/2.)*deltax) << " " << 2*xi_n*kperp*sinkphi+(1.+xi_n*(1.-alpha1/2.)*deltax) << endl;

    complex<double> result=0.;
    //Melosh rotation active nucleon
    vector<complex<double> > wf_in(4,0.), wf_out(4,0.);
    wf_in.at(0)=gpd.getWf()->DeuteronPState(2*pold_in, -1, -1, k_in);
    wf_out.at(0)=gpd.getWf()->DeuteronPState(2*pold_out, -1, -1, k_out);
    wf_in.at(1)=gpd.getWf()->DeuteronPState(2*pold_in, -1, 1, k_in);
    wf_out.at(1)=gpd.getWf()->DeuteronPState(2*pold_out, -1, 1, k_out);
    wf_in.at(2)=gpd.getWf()->DeuteronPState(2*pold_in, 1, -1, k_in);
    wf_out.at(2)=gpd.getWf()->DeuteronPState(2*pold_out, 1, -1, k_out);
    wf_in.at(3)=gpd.getWf()->DeuteronPState(2*pold_in, 1, 1, k_in);
    wf_out.at(3)=gpd.getWf()->DeuteronPState(2*pold_out, 1, 1, k_out);
    wf_in=lf_deut(Ek,k_in,wf_in);
    wf_out=lf_deut(Ekprime,k_out,wf_out);
    

    TransGPD_set gpd_nucl=gpd.chiralodd_grid.getTransGPDSet(x_n,xi_n,t,scale);
    for(int sigma2=-1;sigma2<=1;sigma2+=2){
        for(int sigma1in=-1; sigma1in<=1; sigma1in+=2){
            //complex<double> wfin = gpd.getWf()->DeuteronPState(2*pold_in, sigma2, sigma1in, k_in);
            for(int sigma1out=-1; sigma1out<=1; sigma1out+=2){
                //complex<double> wfout=gpd.getWf()->DeuteronPState(2*pold_out, sigma2, sigma1out, k_out);
                result+=wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))
                            *gpd.getGPD_odd_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,model,right,1,gpd_nucl);

                // cout << wfin << " " << wfout << " " << gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0_n,phin,model,1) << " " << x_n << endl;
                // cout << "N conj" << sigma1in << " " << sigma1out << " " << conj(gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0_n,phin,model,1))
                //     << " " << -gpd.getGPD_odd_nucl(sigma1out, sigma1in, x_n, -xi_n, t, t0_n,phin+PI,model,0) << endl;
                // cout << "N p" << sigma1in << " " << sigma1out << " " << gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0_n,phin,model,1)
                //     << " " << -gpd.getGPD_odd_nucl(-sigma1in, -sigma1out, x_n, xi_n, t, t0_n,-phin+PI,model,0) << endl;
                // cout << "N t" << sigma1in << " " << sigma1out << " " << gpd.getGPD_odd_nucl(sigma1in, sigma1out, x_n, xi_n, t, t0_n,phin,model,1)
                //     << " " << pow(-1.,(sigma1in-sigma1out)/2)*gpd.getGPD_odd_nucl(sigma1out, sigma1in, x_n, -xi_n, t, t0_n,-phin,model,0) << " " <<  pow(-1.,(sigma1in-sigma1out)/2) << endl;

            }
        }
    }
    // cout << "int " << x << " " << result << " " << kperp << " " << alpha1 << endl;
    res[0]=result*2.*kperpprime/alphaprime/(2.-alphaprime)/alpha1*sqrt(Ek)*sqrt(Ekprime);

    // cout << alpha1 << " " << kperp << " " << atan2(ky,kx) << " " << alphaprime << " " << kperpprime << " " << kphiprime << " " << xi_n << " " << x_n << " " << phin << " " << res[0].real() << " " << res[0].imag() << endl;
    // exit(1);
}


complex<double> Deut_Conv_GPD_T::getGPD_odd_nucl(const int sigma_in, const int sigma_out, const double xi, const double t,
                                    const double t0, const double phi, const int model, const bool right, const bool singlet, const TransGPD_set &gpd_nucl){

    // cout << "nucl " << xi << " " << gpd_nucl.getHtildeT_singlet(model) << " " << gpd_nucl.getET_singlet(model) << endl;
    double kinfac=sqrt(t0-t)/MASSn;  // t in MeV^2 !! 
    if(singlet){
        if(right){
            if(sigma_in==sigma_out) return (kinfac*gpd_nucl.getHtildeT_singlet(model)
                                    +(1-sigma_in*xi)*kinfac/2.*gpd_nucl.getET_singlet(model)
                                    +sigma_in*(1-sigma_in*xi)*kinfac/2.*gpd_nucl.getEtildeT_singlet())*exp(I_UNIT*phi);

            else return -(sigma_in-1)*sqrt(1.-xi*xi)*gpd_nucl.getHT_singlet()
                        -sigma_in/2.*kinfac*kinfac*sqrt(1-xi*xi)*exp((sigma_in+1)*phi*I_UNIT)*gpd_nucl.getHtildeT_singlet(model)
                        +(sigma_in-1)*xi*xi/sqrt(1-xi*xi)*gpd_nucl.getET_singlet(model)
                        -(sigma_in-1)*xi/sqrt(1.-xi*xi)*gpd_nucl.getEtildeT_singlet();
        }
        else{
            if(sigma_in==sigma_out) return (kinfac*gpd_nucl.getHtildeT_singlet(model)
                                    +(1+sigma_in*xi)*kinfac/2.*gpd_nucl.getET_singlet(model)
                                    -sigma_in*(1+sigma_in*xi)*kinfac/2.*gpd_nucl.getEtildeT_singlet())*exp(-I_UNIT*phi);
            else return -(sigma_in+1)*sqrt(1.-xi*xi)*gpd_nucl.getHT_singlet()
                        -sigma_in/2.*kinfac*kinfac*sqrt(1-xi*xi)*exp((sigma_in-1)*phi*I_UNIT)*gpd_nucl.getHtildeT_singlet(model)
                        +(sigma_in+1)*xi*xi/sqrt(1-xi*xi)*gpd_nucl.getET_singlet(model)
                        -(sigma_in+1)*xi/sqrt(1.-xi*xi)*gpd_nucl.getEtildeT_singlet();
        }
    }
    else{
        if(right){
            if(sigma_in==sigma_out) return (kinfac*gpd_nucl.getHtildeT_vector(model)
                                    +(1-sigma_in*xi)*kinfac/2.*gpd_nucl.getET_vector(model)
                                    +sigma_in*(1-sigma_in*xi)*kinfac/2.*gpd_nucl.getEtildeT_vector())*exp(I_UNIT*phi);

            else return -(sigma_in-1)*sqrt(1.-xi*xi)*gpd_nucl.getHT_vector()
                        -sigma_in/2.*kinfac*kinfac*sqrt(1-xi*xi)*exp((sigma_in+1)*phi*I_UNIT)*gpd_nucl.getHtildeT_vector(model)
                        +(sigma_in-1)*xi*xi/sqrt(1-xi*xi)*gpd_nucl.getET_vector(model)
                        -(sigma_in-1)*xi/sqrt(1.-xi*xi)*gpd_nucl.getEtildeT_vector();
        }
        else{
            if(sigma_in==sigma_out) return (kinfac*gpd_nucl.getHtildeT_vector(model)
                                    +(1+sigma_in*xi)*kinfac/2.*gpd_nucl.getET_vector(model)
                                    -sigma_in*(1+sigma_in*xi)*kinfac/2.*gpd_nucl.getEtildeT_vector())*exp(-I_UNIT*phi);
            else return -(sigma_in+1)*sqrt(1.-xi*xi)*gpd_nucl.getHT_vector()
                        -sigma_in/2.*kinfac*kinfac*sqrt(1-xi*xi)*exp((sigma_in-1)*phi*I_UNIT)*gpd_nucl.getHtildeT_vector(model)
                        +(sigma_in+1)*xi*xi/sqrt(1-xi*xi)*gpd_nucl.getET_vector(model)
                        -(sigma_in+1)*xi/sqrt(1.-xi*xi)*gpd_nucl.getEtildeT_vector();
        }
    }
}


vector< complex<double> > Deut_Conv_GPD_T::gpd_conv(const double xi, const double x, const double t,const double scale, const int model){
    if(fabs(x)>1.) return vector< complex<double> >(9,0.);
    double t0=-4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi); // Gev^2
    double Delta_perp=sqrt((t0-t)*(1-xi*xi)); //symmetric frame where \bar{P}^perp=0  //GeV

    chiralodd_grid.getTransGPDSet(0.,0.,t*1.E06,scale); //just fill the grid here
    Deut_Conv_GPD_T::Ftor_conv F;
    F.x=x;
    F.t=t*1.E06; //[MeV^2]
    F.scale=scale;
    F.gpd=this;
    F.model=model;
    numint::mdfunction<numint::vector_z,3> mdf;
    mdf.func = &Ftor_conv::exec;
    mdf.param = &F;
    numint::vector_z out(1,0.);

    double low=0.,lowminus=0.,lowprime=0.,lowprimeminus=0.;
    //integration boundaries for alpha, originating from the Heavisides in the convolution (intermediate states w physical + momentum)
    if(fabs(x)>fabs(xi)) {
        low=2.*(fabs(x)+xi)/(1+xi);
        lowminus=2.*(fabs(x)-xi)/(1-xi);
        lowprime=2.*(fabs(x)-xi)/(1-xi);
        lowprimeminus=2.*(fabs(x)+xi)/(1+xi);
    }
    else{
        low=(xi>=0.? 4.*xi/(1+xi): 0.);
        lowminus=(xi>= 0.? 0.: -4.*xi/(1.-xi));
        lowprime=(xi>=0.? 0. : -4.*xi/(1.-xi));
        lowprimeminus=(xi>=0.? 4.*xi/(1+xi): 0.);
    }
    // cout << "low " << low << endl;
    numint::array<double,3> lower = {{low,0.,-PI}};
    numint::array<double,3> lowerminus = {{lowminus,0.,-PI}};
    numint::array<double,3> lowerprime = {{lowprime,0.,-PI}};
    numint::array<double,3> lowerprimeminus = {{lowprimeminus,0.,-PI}};
    numint::array<double,3> upper = {{2.,1.E03,PI}};

    double abserr=1.E-10;
    double relerr=1.E-05;

    unsigned count;
    vector< complex<double> > ret(9,0.);
    for(int i=0;i<9;i++){
        F.f=Deut_Conv_GPD_T::int_k3;
        F.pold_in=i/3-1;
        F.pold_out=i%3-1;
        F.xi=xi;
        F.right=0;
        F.deltax=Delta_perp*1.E03; //[MeV]

        string stf="statefile";
        int nregions,fail,countt;
        vector<complex<double> > err(9,0.);
        vector<complex<double> > prob(9,0.);
        int nvec=1;
        int flags=2;
        int seed=235;
        int minEval=1E02;
        int maxEvalcuba=1E06;

        int maxEval=1E04;
        numint::cube_adaptive(mdf,lower,upper,abserr,relerr,5E02,maxEval,out,count,0);
        //numint::cube_romb(mdf,lower,upper,abserr,relerr,out,count,0);
        //numint::vegas( mdf, lower,upper,nvec, relerr,abserr,flags,seed,minEval, maxEvalcuba,1000, 500, 1000, 0,(stf+"vegas").c_str(),countt,fail,out,err,prob );
        //numint::cuhre( mdf, lower,upper,nvec, relerr,abserr,flags,int(5E02), maxEvalcuba,11, stf.c_str() ,nregions,countt,fail,out,err,prob );
        ret[i]=out[0];
        // exit(1);
        // cout << "normal " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << out[0] << " " << count << endl;
        
                //for min conv model (also put deuteron wave function to only S-wave + no radial dep + H_Tnucleon=0!!!!)
        //  out[0]=ret[i]=test(F.x, F.xi, F.t, F.pold_in, F.pold_out, F.model, F.right, F.deltax);

        // //symmetry checks!!
        // cout << "normal " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << out[0] << " " << count << endl;
        // // F.f=Deut_Conv_GPD_T::int_kprime3;
        // // numint::cube_adaptive(mdf,lowerprime,upper,abserr,relerr,5E02,maxEval,out,count,0);
        // // //numint::cube_romb(mdf,lowerprime,upper,abserr,relerr,out,count,0);
        // // //numint::vegas( mdf, lower,upper,nvec, relerr,abserr,flags,seed,minEval, maxEvalcuba,1000, 500, 1000, 0,(stf+"vegas2").c_str(),countt,fail,out,err,prob );
        // // //numint::cuhre( mdf, lower,upper,nvec, relerr,abserr,flags,int(5E02), maxEvalcuba,11, stf.c_str() ,nregions,countt,fail,out,err,prob );
        // // cout << "normalprime " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << out[0] << " " << count << endl;
        // F.f=Deut_Conv_GPD_T::int_k3;
        // F.pold_in=i%3-1;
        // F.pold_out=i/3-1;
        // F.xi=-xi;
        // F.right=0;
        // F.deltax=-Delta_perp*1.E03;
        // out[0]=ret[i]=test(F.x, F.xi, F.t, F.pold_in, F.pold_out, F.model, F.right, F.deltax);
        // cout << "conj " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << -conj(out[0]) << " " << count << endl;
        // numint::cube_adaptive(mdf,lowerminus,upper,abserr,relerr,5E02,maxEval,out,count,0);
        // cout << "conj " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << -conj(out[0]) << " " << count << endl;
        // F.pold_in=-(i/3-1);
        // F.pold_out=-(i%3-1);
        // F.xi=xi;
        // F.right=0;
        // F.deltax=-Delta_perp*1.E03;
        // out[0]=ret[i]=test(F.x, F.xi, F.t, F.pold_in, F.pold_out, F.model, F.right, F.deltax);
        // cout << "P " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << -out[0] << " " << count << endl;
        // numint::cube_adaptive(mdf,lower,upper,abserr,relerr,5E02,maxEval,out,count,0);
        // cout << "P " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << -out[0] << " " << count << endl;
        // F.pold_in=i%3-1;
        // F.pold_out=i/3-1;
        // F.xi=-xi;
        // F.right=0;
        // F.deltax=Delta_perp*1.E03;
        // out[0]=ret[i]=test(F.x, F.xi, F.t, F.pold_in, F.pold_out, F.model, F.right, F.deltax);
        // cout << "T " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << pow(-1,F.pold_in-F.pold_out)*out[0] << " " << count << endl;
        // numint::cube_adaptive(mdf,lowerminus,upper,abserr,relerr,5E02,maxEval,out,count,0);
        // cout << "T " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << pow(-1,F.pold_in-F.pold_out)*out[0] << " " << count << endl;
    }

    vector< complex<double> > result(9,0.);
    //order of helicity amplitudes is ++,--,00,0+,-0,+0,0-,-+,+-
    result[0]=ret[8];result[1]=ret[0];result[2]=ret[4];result[3]=ret[7];result[4]=ret[3];result[5]=ret[5];result[6]=ret[1];result[7]=ret[6];result[8]=ret[2];
    return result;

}



vector< complex<double> > Deut_Conv_GPD_T::lf_deut(const double Ek, const TVector3& k, const vector< complex<double> > &nonrelwf){
    vector< complex<double> > wf_out(4,0.);
    wf_out[0]=((k[2]+Ek+MASSn)*nonrelwf.at(0)+(-k[0]-I_UNIT*k[1])*nonrelwf.at(1))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    wf_out[1]=((k[2]+Ek+MASSn)*nonrelwf.at(1)+(k[0]-I_UNIT*k[1])*nonrelwf.at(0))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    wf_out[2]=((k[2]+Ek+MASSn)*nonrelwf.at(2)+(-k[0]-I_UNIT*k[1])*nonrelwf.at(3))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    wf_out[3]=((k[2]+Ek+MASSn)*nonrelwf.at(3)+(k[0]-I_UNIT*k[1])*nonrelwf.at(2))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    
    return wf_out;

}

Deut_GPD_T_set Deut_Conv_GPD_T::getDeut_GPD_T_set(const double x, const double xi, const double t, const double scale, const bool ERBL, const int model,
            const int gridsize){
     //make a grid in x,xi since the integrals to compute the chiral odd gpds take some time, t is normally constant for a computation
    if(xi!=xi_grid||t!=t_grid||grid_set==false||ERBL!=ERBL_set||gridsize!=grid_size||model!=model_set){
        if(grid!=NULL) delete [] grid;
        grid = new Deut_GPD_T_set[gridsize+1];
        std::string filename = string(HOMEDIR)+"/gpd_deutgrids/T.xi"+to_string(xi)+".t"+to_string(t)+".mu"+to_string(scale)+".ERBL"+to_string(ERBL)+".model"+to_string(model)
        +".size"+to_string(gridsize);
        ifstream infile(filename.c_str());
        if(!infile.is_open()){
            cout << "constructing chiral odd deuteron helamps grid " << filename << endl;
            ofstream outfile(filename.c_str());

            for(int i=0;i<=gridsize;i++){
                double x=double(i)/gridsize*(ERBL? abs(xi): 1.)+(i==0? 1.E-04:0.);

                vector< complex<double> > result = gpd_conv(xi,x,t,scale,model);
                vector< complex<double> > resultmin = gpd_conv(xi,-x,t,scale, model);
                vector< complex<double> > total(9,0.);
                //factor of 2 because gpd_conv calculates helicity amplitudes, which differ by i\sigma^+R matrix elements by a factor of 2!
                for(int k=0; k<9; k++) total[k]=2.*(result[k]+resultmin[k]);
               grid[i]=Deut_GPD_T_set(total[0].real(),total[1].real(),total[2].real(),total[3].real(),
            total[4].real(),total[5].real(),total[6].real(),total[7].real(),total[8].real());
                vector < complex<double> > gpds=helamps_to_gpds_T(xi,t,total);
                outfile << x << " " << total[0].real() << " " << total[1].real() << " " << total[2].real() << " " << total[3].real() << " " << total[4].real() << 
                " " << total[5].real() << " " << total[6].real() << " " << total[7].real() << " " << total[8].real() << " " <<
                gpds[0].real() << " " << gpds[1].real() << " " << gpds[2].real() << " " << gpds[3].real() << " " << gpds[4].real() << 
                " " << gpds[5].real() << " " << gpds[6].real() << " " << gpds[7].real() << " " << gpds[8].real() << " " << endl;
            }
            outfile.close();
            
            cout << "construction done" << endl;
            grid_set=true;
            t_grid=t;
            xi_grid=xi;
            ERBL_set=ERBL;
            grid_size = gridsize;
            model_set=model;
        }
        else{
            cout << "Reading in deut T gpd grid " << filename << endl;
            string line;
            for(int i=0;i<=gridsize;i++){
                getline (infile,line);
                istringstream iss(line);
                vector<string> tokens{istream_iterator<string>{iss},
                      istream_iterator<string>{}};
                grid[i]=Deut_GPD_T_set(stod(tokens[1]),stod(tokens[2]),stod(tokens[3]),stod(tokens[4]),stod(tokens[5]),
                stod(tokens[6]),stod(tokens[7]),stod(tokens[8]),stod(tokens[9]));
            }
            infile.close();
            grid_set=true;
            t_grid=t;
            xi_grid=xi;
            ERBL_set=ERBL;
            grid_size = gridsize;
            model_set=model;
        }
   }
   //interpolation
    double index_i=0.;
    double frac_i=modf(x*gridsize/(ERBL? abs(xi): 1.),&index_i);
    //cout << x << " " << index_i << " "<< frac_i << endl;
    return grid[int(index_i)]*(1.-frac_i)+grid[int(index_i)+1]*(frac_i);
}
