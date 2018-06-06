#include<OldDeuteron.hpp>
#include<TDeuteron.h>
#include<TInterpolatingWavefunction.h>
#include<TVector.h>


#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;


#include <iostream>
#include <cmath>
#include "constants.hpp"
#include <numint/numint.hpp>
#include <numint/numint2Cuba.hpp>

void int_3d(numint::vector_z &, double p, double costheta, double phi, TInterpolatingWavefunction &wf, TVector3 &pvec);



int main(int argc, char *argv[]){
    OldDeuteron wf1("AV18");
    OldDeuteron wf2("AV18b");
    OldDeuteron wf3("CDBonn");
    OldDeuteron wf4("Paris");
    
    double sum1=0.,sum2=0.,sum3=0.,sum4=0.;

    //normalization test
    for(int i=0; i<=1000; i++){
        sum1+=i*i*(pow(wf1.U(i),2.)+pow(wf1.W(i),2.));
        sum2+=i*i*(pow(wf2.U(i),2.)+pow(wf2.W(i),2.));
        sum3+=i*i*(pow(wf3.U(i),2.)+pow(wf3.W(i),2.));
        sum4+=i*i*(pow(wf4.U(i),2.)+pow(wf4.W(i),2.));
    }

    cout << "normalization OldDeuteron " << sum1*4.*PI << " " << sum2*4.*PI << " " << sum3*4.*PI << " " << sum4*4.*PI << endl;
    

    TDeuteron::Wavefunction *wfref = TDeuteron::Wavefunction::CreateWavefunction("Paris");
    TInterpolatingWavefunction wf;
    double r_grid=20.E-03;
    for(int i=0;i<=1000;i++){
        wf.AddUp(i,wfref->GetUp(i));
        wf.AddWp(i,wfref->GetWp(i));
        wf.AddUr(i*r_grid,wfref->GetUr(i*r_grid));
        wf.AddWr(i*r_grid,wfref->GetWr(i*r_grid));
        
    }

    for(int dspin=-2;dspin<=2;dspin+=2){
        for(int pspin=-1;pspin<=1;pspin+=2){
            for(int nspin=-1;nspin<=1;nspin+=2){
                cout << dspin << " " << pspin << " " << nspin << " " << wf1.deuteronwf(dspin,1,pspin,nspin,100.,1.,2.) 
                    << " " << wf.DeuteronPState(dspin,nspin,pspin,TVector3(100.*sin(1)*cos(2),100.*sin(1)*sin(2),100.*cos(1))) << endl;
            }
        }
    }

    //CG-test!
    cout << "Testing Clebsch Gordan" << endl;
    cout << TDeuteron::Wavefunction::ClebschGordan(4,4,1,1,5,5) << " " << 1. << endl;

    double pnorm=0.,rnorm=0.,quad=0.;
    for(int i=0;i<=100;i++){
        rnorm+=pow(wf.GetUr(i*r_grid*10),2.)+pow(wf.GetWr(i*r_grid*10),2.);
        pnorm+=i*i*10*10*(pow(wf.GetUp(i*10),2.)+pow(wf.GetWp(i*10),2.));
        quad+=pow(i*r_grid*10,2.)*(sqrt(8.)*wf.GetUr(i*r_grid*10)*wf.GetWr(i*r_grid*10)-pow(wf.GetWr(i*r_grid*10),2.));
        
    } 
    cout << "R-space norm " << rnorm*r_grid*10 << endl;
    cout << "P-space norm " << pnorm*10 << endl;
    cout << "quadrupole moment " << quad*r_grid/20.*10 << endl;

    //FT radial functions
    double testsum=0.,rad=1.;
    for(int i=1;i<=1000;i++){
        testsum+=i*i*sin(i*rad*INVHBARC)/(i*rad*INVHBARC)*wf.GetUp(i);
    }
    cout << testsum*sqrt(2./PI)*pow(INVHBARC,3./2.) << " " << wf.GetUr(rad)/rad << endl;
    double testsum2=0.;
    for(int i=1;i<=1000;i++){
        double xx=i*rad*INVHBARC;
        testsum2+=i*i*((3/(xx*xx)-1)*sin(xx)/xx-3.*cos(xx)/(xx*xx))*wf.GetWp(i);
    }
    cout << testsum2*sqrt(2./PI)*pow(INVHBARC,3./2.) << " " << wf.GetWr(rad)/rad << endl;

    struct Ftor {  //Carbon

        static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
            Ftor &p = * (Ftor *) param;
            p.f(ret,x[0],x[1],x[2],*p.wf,p.pvec);
        }
        TInterpolatingWavefunction *wf;
        TVector3 pvec;
        void (*f)(numint::vector_z &, double p, double costheta, double phi, TInterpolatingWavefunction &wf, TVector3 &pvec);

    };

    Ftor F;
    F.wf=&wf;
    F.pvec=TVector3(50.,40.,30.);
    numint::mdfunction<numint::vector_z,3> mdfH;
    mdfH.func = &Ftor::exec;
    mdfH.param = &F;

    numint::array<double,3> lowerH = {{0.,-1.,0.}};  
    numint::array<double,3> upperH = {{1000.*r_grid,1.,2.*PI}};
    
    F.f=int_3d;
    vector<complex<double> > adapnorm(13,0.); 
    unsigned count=0;
    numint::cube_adaptive(mdfH,lowerH,upperH,1.E-30,1.E-03,2E02,4E04,adapnorm,count,0); 


    cout << "3dnorm " << adapnorm[0] << endl;

    //3D FT wf, OK!!
    cout << endl << endl << "Checking FT of wf " << endl;
    for(int dspin=-2;dspin<=2;dspin+=2){
        for(int pspin=-1;pspin<=1;pspin+=2){
            for(int nspin=-1;nspin<=1;nspin+=2){
                cout << adapnorm[1+(dspin+2)/2*4+(pspin+1)+(nspin+1)/2] << " " << wf.DeuteronPState(dspin,nspin,pspin,TVector3(50.,40.,30.)) << endl;
            }
        }
    }

    delete wfref;

}


void int_3d(numint::vector_z &res, double r, double costheta, double phi, TInterpolatingWavefunction &wf, TVector3 &pvec){

    res=numint::vector_z(13,0.);
    double d3dnorm=0.;
    double cosphi=cos(phi);
    double sinphi=sin(phi);
    double sintheta =sqrt(1.-costheta*costheta);
    complex<double> expf=exp(-I_UNIT*(pvec[0]*r*sintheta*cosphi+pvec[1]*r*sintheta*sinphi+pvec[2]*r*costheta)*INVHBARC);
    for(int dspin=-2;dspin<=2;dspin+=2){
        for(int pspin=-1;pspin<=1;pspin+=2){
            for(int nspin=-1;nspin<=1;nspin+=2){
                d3dnorm+=norm(wf.DeuteronRState(dspin,nspin,pspin,TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
                res[1+(dspin+2)/2*4+(pspin+1)+(nspin+1)/2]=expf
                    *r*wf.DeuteronRState(dspin,nspin,pspin,TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta))*pow(2.*PI*HBARC,-3./2.);
            }
        }
    }
    res[0]=d3dnorm/3.;
}