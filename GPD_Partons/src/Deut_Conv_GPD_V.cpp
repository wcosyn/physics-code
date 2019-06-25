#include "Deut_Conv_GPD_V.hpp"
#include <constants.hpp>
#include <cassert>
#include <numint/numint.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

Deut_Conv_GPD_V::Deut_Conv_GPD_V(PARTONS::GPDService* pGPDService, PARTONS::GPDModule* pGPDModel, const string &wfname):
pGPDService(pGPDService),pGPDModel(pGPDModel),chiraleven_grid(pGPDService, pGPDModel),incH(1), incE(1)
{


    wfref = TDeuteron::Wavefunction::CreateWavefunction(wfname);
    for(int i=0;i<=1000;i++){
        wf.AddUp(i,wfref->GetUp(i));
        wf.AddWp(i,wfref->GetWp(i));
    }

}

Deut_Conv_GPD_V::~Deut_Conv_GPD_V(){
      delete wfref;
}



vector< complex<double> > Deut_Conv_GPD_V::helamps_to_gpds_V(const double xi, const double t, const vector< complex <double> > & helamps){

    vector< complex<double> > gpds(5,0.);
    //symmetric frame, with momentum transfer in x,z plane, so phi=0
    double t0=-4.*MASSD*MASSD*xi*xi/(1-xi*xi)*1.E-06;
    double D=(t0-t)/(4.*MASSD*MASSD)*1.E06;
    gpds.at(0)=1./3./pow(1.-xi*xi,2.)*(3.*pow(xi,4.)-7.*xi*xi+2-2.*D*(1-xi*xi))*helamps[0]+1./3./(1-xi*xi)*helamps[1]
        +sqrt(2./D*(1+xi)/(1-xi))/3./pow(1.-xi*xi,2.)*(D*(1.-xi*xi)+xi)*helamps[2]
        -sqrt(2./D*(1-xi)/(1+xi))/3./pow(1.-xi*xi,2.)*(D*(1.-xi*xi)-xi)*helamps[3]
        -1./3./D/pow(1-xi*xi,3.)*(2.*xi*xi+D*(3.*pow(xi,6.)-10.*pow(xi,4.)+9.*xi*xi-2))*helamps[4];

    gpds.at(1)=2./(1-xi*xi)*helamps[0]-sqrt(1./2./D*(1.+xi)/(1.-xi))/(1.-xi)*helamps[2]+sqrt(1./2./D*(1.-xi)/(1.+xi))/(1.+xi)*helamps[3]
        +2.*xi*xi/D/pow(1.-xi*xi,2.)*helamps[4];

    gpds.at(2)=-1/D*helamps[4];
    gpds.at(3)=-2.*xi/(1-xi*xi)*helamps[0]+sqrt(1./D/2.*(1.+xi)/(1.-xi))/(1.-xi)*helamps[2]+sqrt(1./D/2.*(1.-xi)/(1.+xi))/(1.+xi)*helamps[3]
    -2.*xi/D/pow(1.-xi*xi,2.)*helamps[4];


    gpds.at(4)=-1./pow(1.-xi*xi,2.)*(xi*xi+2.*D*(1.-xi*xi)+1)*helamps[0]+1./(1.-xi*xi)*helamps[1]
        +sqrt(2./D*(1+xi)/(1-xi))/pow(1.-xi*xi,2.)*(D*(1.-xi*xi)+xi)*helamps[2]
        -sqrt(2./D*(1-xi)/(1+xi))/pow(1.-xi*xi,2.)*(D*(1.-xi*xi)-xi)*helamps[3]
                -1./D/pow(1.-xi*xi,3.)*(2.*xi*xi+D*(1.-pow(xi,4.)))*helamps[4];


    return gpds;
}

vector< complex<double> > Deut_Conv_GPD_V::gpds_to_helamps_V(const double xi, const double t, const vector< complex<double> > & gpds){

    vector< complex<double> > helamps(5,0.);
    //symmetric frame, with momentum transfer in x,z plane, so phi=0

    double t0=-4.*MASSD*MASSD*xi*xi/(1-xi*xi)*1.E-06;
    double D=(t0-t)/(4.*MASSD*MASSD)*1.E06;
    // cout << (t0-t)*1.E-06 << " " << MASSD*1.E-03 << endl;
    helamps.at(0)=gpds[0]+D*gpds[2]-1./3.*gpds[4];

    helamps.at(1)=(-2.*D+(1+xi*xi)/(1-xi*xi))*gpds[0]+(2.*D-2.*xi*xi/(1-xi*xi))*gpds[1]
            -2.*(D*D*pow(1-xi*xi,2.)-xi*xi)/pow(1-xi*xi,2.)*gpds[2]+(2.*D*xi-2.*xi/(1-xi*xi))*gpds[3]
            +((1-xi*xi)-1./3.*(-2.*D+(1+xi*xi)/(1-xi*xi)))*gpds[4];

    helamps.at(2)=sqrt(2.*D*(1-xi*xi))/(1+xi)*(gpds[0]-(1.-xi)/2.*(gpds[1]-gpds[3])+(D*(1.-xi*xi)-xi)/(1.-xi*xi)*gpds[2]-1./3.*gpds[4]);

    helamps.at(3)=sqrt(2.*D*(1-xi*xi))/(1-xi)*(-gpds[0]+(1.+xi)/2.*(gpds[1]+gpds[3])-(D*(1.-xi*xi)+xi)/(1.-xi*xi)*gpds[2]+1./3.*gpds[4]);                    

    helamps.at(4)=-D*gpds[2];                    
    return helamps;
}



void Deut_Conv_GPD_V::int_k3(numint::vector_z & res, double alpha1, double kperp, double kphi, Deut_Conv_GPD_V &gpd,
              double x, double xi, double t, double scale, int pold_in, int pold_out, double deltax){
    res=numint::vector_z(1,0.);

    double Ek=sqrt((MASSn*MASSn+kperp*kperp)/alpha1/(2.-alpha1));
    double kz=(alpha1-1.)*Ek;
    double knorm=sqrt(Ek*Ek-MASSn*MASSn);
    // cout << alpha1 << " " << kperp << " " << Ek << " " << knorm << endl;
    if(knorm>1.E03) {return;}

    // double Delta_perp=sqrt((t0-t)*(1-xi*xi)); //symmetric frame where \bar{P}^perp=0

    double alphaprime=(alpha1*(1+xi)-4.*xi)/(1-xi);

    double xi_n=xi/(alpha1/2.*(1+xi)-xi);
    double t0_n=-4.*MASSn*MASSn*xi_n*xi_n/(1-xi_n*xi_n)*1.E-06;
    // double t0=-4.*MASSD*MASSD*xi*xi/(1-xi*xi);
    if(t>t0_n) return;
    
    double x_n=x/(alpha1/2.*(1+xi)-xi);
    if(abs(x_n)>=1.) return; //overflow
    //cout << "int " << x << " " << xi << " " << x_n << " " << xi_n << " " << xi/x << " " << xi_n/x_n << endl;
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
    
    // PARTONS::GPDKinematic gpdKinematic(x_n,abs(xi_n),t, 1., 1.);
    // PARTONS::GPDResult gpdResult = gpd.pGPDService->computeGPDModel(gpdKinematic,
    //         gpd.pGPDModel);
    // double H=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
    //     +gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
    // double E=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
    //     +gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
    vector<double> nuclgpd = gpd.chiraleven_grid.getVectorGPDSet(x_n,abs(xi_n),t,scale);
    double E=gpd.incE? nuclgpd[1] : 0.;
    double H= gpd.incH? nuclgpd[0] : 0.;
    if(isnan(E)||isnan(H)) return;

    for(int sigma2=-1;sigma2<=1;sigma2+=2){  //spectator helicity
        for(int sigma1in=-1; sigma1in<=1; sigma1in+=2){ //initial active nucleon helicity
            // complex<double> wfin = gpd.getWf()->DeuteronPState(2*pold_in, sigma2, sigma1in, k_in);
            for(int sigma1out=-1; sigma1out<=1; sigma1out+=2){ //final active nucleon helicity
                // complex<double> wfout=gpd.getWf()->DeuteronPState(2*pold_out, sigma2, sigma1out, k_out);
                // result+=wfin*conj(wfout)
                //             *gpd.getGPD_odd_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,model,right,gpd_nucl);
                result+=wf_in.at((sigma1in+1)/2+(sigma2+1))*conj(wf_out.at((sigma1out+1)/2+(sigma2+1)))
                            *gpd.getGPD_even_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,H,E);
                //if(isnan(real(result))) cout << result << " " << gpd.getGPD_even_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,H,E) << " " << H << " " << E << endl;
                //symmetry tests nucleon level (tested OK)
                // cout << "N conj " << sigma1in << " " << sigma1out << " " << conj(gpd.getGPD_even_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,H,E))
                //     << " " << gpd.getGPD_even_nucl(sigma1out, sigma1in, -xi_n, t, t0_n,phin+PI,H,E) << endl;
                // cout << "N p " << sigma1in << " " << sigma1out << " " << gpd.getGPD_even_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,H,E)
                //     << " " << gpd.getGPD_even_nucl(-sigma1in, -sigma1out, xi_n, t, t0_n,-phin+PI,H,E) << endl;
                // cout << "N t " << sigma1in << " " << sigma1out << " " << gpd.getGPD_even_nucl(sigma1in, sigma1out, xi_n, t, t0_n,phin,H,E)
                //     << " " << pow(-1.,(sigma1in-sigma1out)/2)*gpd.getGPD_even_nucl(sigma1out, sigma1in, -xi_n, t, t0_n,-phin,H,E) << " " <<  pow(-1.,(sigma1in-sigma1out)/2) << endl;

            }
        }
    }
    // cout << "int " << x << " " << result << " " << kperp << " " << alpha1 << endl;
    res[0]=result*2.*kperp/alpha1/(2.-alpha1)/alphaprime*sqrt(Ek)*sqrt(Ekprime);    // cout << t0 << " " << alpha1 << " " << t0_n << " " << xi << " " << xi_n << " " << res[0].real() << endl;
    // cout << alpha1 << " " << kperp << " " << kphi << " " << alphaprime << " " << kperpprime << " " << atan2(kyprime,kxprime) << " " << xi_n << " " << x_n << " " << phin << " " << res[0].real() << " " << res[0].imag() << endl;
    // exit(1);
}

complex<double> Deut_Conv_GPD_V::getGPD_even_nucl(const int sigma_in, const int sigma_out, const double xi, const double t,
                                    const double t0, const double phi,  const double gpd_H, const double gpd_E) const{

    // cout << "nucl " << xi << " " << gpd_nucl.getHtildeT_singlet(model) << " " << gpd_nucl.getET_singlet(model) << endl;
    double kinfac=sqrt(t0-t)/MASSn*1.E03;
    if(sigma_in==sigma_out) return sqrt(1-xi*xi)*gpd_H - xi*xi/sqrt(1-xi*xi)*gpd_E;

    else{
        if(sigma_in==-1) return -exp(-I_UNIT*phi)*kinfac/2.*gpd_E ;
        else return exp(I_UNIT*phi)*kinfac/2.*gpd_E ;
    }
}


vector< complex<double> > Deut_Conv_GPD_V::gpd_conv(const double xi, const double x, const double t, const double scale){
    if(fabs(x)>1.) return vector< complex<double> >(5,0.);
    double t0=-4.*MASSD*MASSD*xi*xi/(1-xi*xi)*1.E-06;
    double Delta_perp=sqrt((t0-t)*(1-xi*xi)); //symmetric frame where \bar{P}^perp=0

    
    Deut_Conv_GPD_V::Ftor_conv F;
    F.x=x;
    F.t=t;
    F.scale=scale;
    F.gpd=this;
    
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

    double abserr=1.E-15;
    double relerr=1.E-03;

    unsigned count;
    vector< complex<double> > ret(9,0.);
    //we only need 5 Helicity amplitudes
    for(int i=4;i<9;i++){
        F.f=Deut_Conv_GPD_V::int_k3;
        F.pold_in=i/3-1;
        F.pold_out=i%3-1;
        F.xi=xi;
        F.deltax=Delta_perp;

        string stf="statefile";
        int nregions,fail,countt;
        vector<complex<double> > err(9,0.);
        vector<complex<double> > prob(9,0.);
        int nvec=1;
        int flags=2;
        int seed=235;
        int minEval=1E02;
        int maxEvalcuba=1E06;

        int maxEval=1E06;
        numint::cube_adaptive(mdf,lower,upper,abserr,relerr,5E02,maxEval,out,count,0);
        ret[i]=out[0];

        //symmetry checks!! [seems ok, but quite large differences for the small helicity amps]
        // cout << "normal " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << out[0] << " " << count << endl;
        // F.pold_in=i%3-1;
        // F.pold_out=i/3-1;
        // F.xi=-xi;
        // F.deltax=-Delta_perp;
        // numint::cube_adaptive(mdf,lowerminus,upper,abserr,relerr,5E02,maxEval,out,count,0);
        // cout << "conj " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << conj(out[0]) << " " << count << endl;
        // F.pold_in=-(i/3-1);
        // F.pold_out=-(i%3-1);
        // F.xi=xi;
        // F.deltax=-Delta_perp;
        // numint::cube_adaptive(mdf,lower,upper,abserr,relerr,5E02,maxEval,out,count,0);
        // cout << "P " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << out[0] << " " << count << endl;
        // F.pold_in=i%3-1;
        // F.pold_out=i/3-1;
        // F.xi=-xi;
        // F.deltax=Delta_perp;
        // numint::cube_adaptive(mdf,lowerminus,upper,abserr,relerr,5E02,maxEval,out,count,0);
        // cout << "T " << x << " " << i << " " << F.pold_in << " " << F.pold_out << " " << pow(-1,F.pold_in-F.pold_out)*out[0] << " " << count << endl;

     }

    vector< complex<double> > result(5,0.);
    //order of helicity amplitudes is ++,00,0+,+0,-+
    result[0]=ret[8];result[1]=ret[4];result[2]=ret[7];result[3]=ret[5];result[4]=ret[6];
    return result;

}



vector< complex<double> > Deut_Conv_GPD_V::lf_deut(const double Ek, const TVector3& k, const vector< complex<double> > &nonrelwf){
    vector< complex<double> > wf_out(4,0.);
    wf_out[0]=((k[2]+Ek+MASSn)*nonrelwf.at(0)+(-k[0]-I_UNIT*k[1])*nonrelwf.at(1))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    wf_out[1]=((k[2]+Ek+MASSn)*nonrelwf.at(1)+(k[0]-I_UNIT*k[1])*nonrelwf.at(0))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    wf_out[2]=((k[2]+Ek+MASSn)*nonrelwf.at(2)+(-k[0]-I_UNIT*k[1])*nonrelwf.at(3))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    wf_out[3]=((k[2]+Ek+MASSn)*nonrelwf.at(3)+(k[0]-I_UNIT*k[1])*nonrelwf.at(2))/sqrt(2.*(Ek+MASSn)*(Ek+k[2]));
    
    return wf_out;

}


Deut_GPD_V_set Deut_Conv_GPD_V::getDeut_GPD_V_set(const double x, const double xi, const double t, const double scale, const bool ERBL){
    //make a grid in x,xi since the integrals to compute the chiral odd gpds take some time, t is normally constant for a computation
    if(xi!=xi_grid||t!=t_grid||grid_set==false||ERBL!=ERBL_set){
        std::string filename = string(HOMEDIR)+"/gpd_deutgrids/V_grid.xi"+to_string(xi)+".t"+to_string(t)+".mu"+to_string(scale)+".EBRL"+to_string(ERBL);
        ifstream infile(filename.c_str());
        if(!infile.is_open()){
            cout << "constructing chiral even deuteron helamps grid " << filename << endl;
            ofstream outfile(filename.c_str());

            for(int i=0;i<=100;i++){
                double x=0.02*(i-50)*(ERBL? abs(xi): 1.)+(i==50? 1.E-04:0.);
                vector< complex<double> > result = gpd_conv(xi,x,t, scale);
                vector< complex<double> > gpd = helamps_to_gpds_V(xi,t,result);
                grid[i]=Deut_GPD_V_set(result[0].real(),result[1].real(),result[2].real(),result[3].real(),result[4].real());
                outfile << x << " " << result[0].real() << " " << result[1].real() << " " << result[2].real() << " " << result[3].real() << " " << result[4].real()
                <<" " << gpd[0].real() << " " << gpd[1].real() << " " << gpd[2].real() << " " << gpd[3].real() << " " << gpd[4].real() << " " << endl;
            }
            outfile.close();
            
            cout << "construction done" << endl;
            grid_set=true;
            t_grid=t;
            xi_grid=xi;
            ERBL_set=ERBL;
        }
        else{
            cout << "Reading in grid " << filename << endl;
            string line;
            for(int i=0;i<=100;i++){
                getline (infile,line);
                istringstream iss(line);
                vector<string> tokens{istream_iterator<string>{iss},
                      istream_iterator<string>{}};
                grid[i]=Deut_GPD_V_set(stod(tokens[1]),stod(tokens[2]),stod(tokens[3]),stod(tokens[4]),stod(tokens[5]));
            }
            infile.close();
            grid_set=true;
            t_grid=t;
            xi_grid=xi;
            ERBL_set=ERBL;
        }
   }
   //interpolation
    double index_i=0.;
    double frac_i=modf(x*50/(ERBL? abs(xi): 1.)+50,&index_i);

    return grid[int(index_i)]*(1.-frac_i)+grid[int(index_i)+1]*(frac_i);
}
