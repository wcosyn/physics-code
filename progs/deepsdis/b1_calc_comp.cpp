//calculate b1 from the tagged structure functions

//calculate Azz from tagged structure functions, accounting for virtual photon angle etc.

//has some explicit checks between structure functions b_i and the cross section structure functions (L,T etc)
//also comparison with the EPW definitions etc.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <NuclStructure.hpp>
#include <NucleonStructure.hpp>

#include <TDeuteron.h>
#include <numint/numint.hpp>


void k_int(numint::vector_d & res, double knorm, double kcosth, TDeuteron::Wavefunction *wfref,
          double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq, NucleonStructure &strp);


int main(int argc, char *argv[])
{
  double Q2 = atof(argv[2])*1.E06; 
  string wf = argv[3];
  double Ein= atof(argv[5])*1.E03; //beam energy in GeV
  string nuclstruc = argv[4];
  for(int i=1;i<=100;i++){
    double x=1./100.*i;
    
    TDeuteron::Wavefunction *wfref;
    wfref = TDeuteron::Wavefunction::CreateWavefunction(wf);
  //   for(int i=1;i<100;i++){
  //     double x=0.01*i;
  //     NuclStructure s1(1,Q2,x,0,"SLAC");
  //     NuclStructure s2(1,Q2,x,0,"HMRS");
  //     cout << x << " " << s1.getF1() << " " << s2.getF1() << " " << s1.getF2() << " " << s2.getF2() << endl;
  //   }
  //   exit(1);
    NucleonStructure strp(nuclstruc);
    
    
    
    struct Ftor_b1 {

      /*! integrandum function */
      static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
        Ftor_b1 &p = * (Ftor_b1 *) param;
        p.f(ret,x[0],x[1], p.wfref,p.Q2,p.x,p.nu,p.qvec,p.gamma,p.eps,p.thetaq,p.strp);
      }
      TDeuteron::Wavefunction *wfref;
      double Q2;
      double x;
      double nu;
      double qvec;
      double gamma;
      double eps;
      double thetaq;
      NucleonStructure strp;
      
      
      
      void (*f)(numint::vector_d & res, double knorm, double kcosth, TDeuteron::Wavefunction *wfref,
                double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq, NucleonStructure &strp);
    };
    
    double nu=Q2/(2.*MASSD*x);
    double qvec=sqrt(Q2+nu*nu);
    double gamma=sqrt(Q2)/nu;
  //     double gamma=0.;
    double Eout=Ein-nu;
    double y=nu/Ein;
//     if(y<1.){
      double eps=(1-y-gamma*gamma*y*y/4.)/(1-y+y*y/2+gamma*gamma*y*y/4.);
      double thetae=asin(sqrt(Q2/4/Ein/Eout))*2;
      double thetaq=acos((Ein*Ein+qvec*qvec-Eout*Eout)/2/Ein/qvec);

      
      
      
      
      NucleonStructure stp(nuclstruc);

      numint::array<double,2> lower = {{0.,-1.}};
      numint::array<double,2> upper = {{600.,1.}};
      Ftor_b1 F;
      F.wfref=wfref;
      F.Q2=Q2;
      F.x=x;
      F.nu=nu;
      F.qvec=qvec;
      F.gamma=gamma;
      F.eps=eps;
      F.thetaq=thetaq;
      F.strp=strp;
      
      double F1p=0.,F2p=0.,F1n=0.,F2n=0.;
        strp.getF_xQ(F1p,F2p,1,x*MASSD/MASSn,Q2);  
        strp.getF_xQ(F1n,F2n,0,x*MASSD/MASSn,Q2);  

      numint::mdfunction<numint::vector_d,2> mdf;
      mdf.func = &Ftor_b1::exec;
      mdf.param = &F;
      numint::vector_d ret(12,0.);
      F.f=k_int;
      int res=90;
      unsigned count=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E05,ret,count,0);
      
      double R=NucleonStructure::getr1998(x,Q2);
      double kap=1+gamma*gamma;

      double term1=(1.+3.*cos(2*thetaq))*(2.*(1-eps)*kap-gamma*gamma/3./kap+eps*(kap*kap+kap+1)*2/3/kap);
      double term2=sin(2.*thetaq)*sqrt(2.*eps*(1+eps))*gamma*(kap+2)/kap;
      double term3=(1-2.*cos(2.*thetaq))*eps*gamma*gamma/kap;

      double factor=1./8./(1+eps*R)*(term1+term2+term3);

      double A=-1./2./(1.+gamma*gamma)*(ret[4]+ret[7]);
      double D=-2.*ret[7]/gamma/gamma;
      double C=(-2.*ret[6]/gamma-2.*D)/(1.+gamma*gamma);
      double B=(2.*(1+gamma*gamma)*A-(1.+gamma*gamma)*C-D-ret[5])/pow(1.+gamma*gamma,2.);

    double b1=A;
    double b2=x*(B+C+D);
    double b3=x*(B+C-2.*D)/3.;
    double b4=x*(B-2.*C+D)/3.;
    
    
    //q | x_N |  thetaq | b1 | b1_EPW from F_UT/F_TT ratio | b1 from azz hermes way | b1 from azz improved approx way | azz | F_T/2 | sum F1*2 | extra mult factor from azz to b1 | 
    //F_TT/angle corrected nominateor | F_UT / F_UT+epsF_UL | eps | gamma | b1-4
    cout << setprecision(4) <<  " " << qvec << " " << 2.*x << " "  << thetaq*RADTODEGR << " " << sqrt(2./3.)*ret[10]/ret[11] << " " << ret[0] << " " << sqrt(2./3.)*ret[4]/ret[8]*2.*(F1p+F1n)*-3./2. << " " << sqrt(2./3.)*ret[10]/ret[11]*2.*(F1p+F1n)*-3./2. 
      << " " << sqrt(2./3.)*ret[10]/ret[11]/factor*2.*(F1p+F1n)*-3./2. << " " <<ret[8]/2. << " " << 2.*(F1p+F1n) << " " 
      << factor << " " << ret[4]/ret[10] << " " << ret[8]/ret[11] << " "<< eps << " "<< b1 << " " << b2 << " " << b3 << " " << b4
      << " " << b2/2./x/b1 << " " << gamma << endl;
      
      
      //x_N | b1 | b1_EPW from F_UT/F_TT ratio | b1 from azz hermes way | b1 from azz improved approx way | azz | F_T/2 | sum F1*2 | extra mult factor from azz to b1 | 
      // F_TL | F_TT | F_Tcosphi | F_Tcos2phi | azz nom | eps | thetaq
  //     cout << setprecision(4) << 2.*x << " " << ret[0] << " " << sqrt(2./3.)*ret[4]/ret[8]*2.*(F1p+F1n)*-3./2. << " " << sqrt(2./3.)*ret[10]/ret[11]*2.*(F1p+F1n)*-3./2. 
  //       << " " << sqrt(2./3.)*ret[10]/ret[11]/factor*2.*(F1p+F1n)*-3./2. << " " << sqrt(2./3.)*ret[10]/ret[11] << " " << ret[8]/2. << " " << 2.*(F1p+F1n) << " " 
  //       << factor << " " << ret[4]/ret[10] << " " << ret[8]/ret[11] << " "<< eps << " " << thetaq*RADTODEGR << " " << ret[12] << " " << ret[13] << " " << ret[14] << " " << ret[15] << " " << ret[13]/2./x/ret[12] << " " << gamma << endl;
      cout << Q2 << " " << 2.*x  << " " << b1 << " " << b2 << " " << b3 << " " << b4 << endl;
        
      //checks
  //     cout <<   1/x*sqrt(2./3.)*(2*(1+gamma*gamma)*x*b1-pow(1+gamma*gamma,2.)*(1./3.*b2+b3+b4)-(1+gamma*gamma)*(1./3.*b2-b4)-(1./3.*b2-b3)) << ret[5] << endl;
  //     cout << -1/x*sqrt(2./3.)*(2*(1+gamma*gamma)*x*b1-gamma*gamma*(1./6.*b2-1./2.*b3)) << ret[4] << endl;
  //     cout << -gamma/2./x*sqrt(2./3.)*((1+gamma*gamma)*(1./3.*b2-b4)+(2./3.*b2-2.*b3)) << ret[6] << endl;
  //     cout << -1/x*sqrt(2./3.)*gamma*gamma*(1./6.*b2-1./2.*b3) << ret[7] << endl;
        //ret[0] b1 total
      //ret[1] b1 sd
      //ret[2] b1 dd
      //ret[3] b1 total check
      //ret[4-7] tensor struc functions
      //ret[8-9] FL,FT
      //ret[10-11] azz nom, denom
      
//     }    
  }
  
}

void k_int(numint::vector_d & res, double knorm, double costh, TDeuteron::Wavefunction *wfref,
	      double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq,
              NucleonStructure &strp){
    
  
  res=numint::vector_d(12,0.);
 /* res[0]=res[1]=res[2]=pow(wfref->GetWp(knorm),2.)+pow(wfref->GetUp(knorm),2.)*knorm*knorm;
  cout << knorm << " " << costh << endl;
  return;
 */ 
  double dens_tensor_tot = wfref->GetUp(knorm)*wfref->GetWp(knorm)/sqrt(2.)+pow(wfref->GetWp(knorm),2.)/4.;  // UW/sqrt(2)+W^2/4
  double dens_tensor_SD  = wfref->GetUp(knorm)*wfref->GetWp(knorm)/sqrt(2.); // UW/sqrt(2)
  double dens_tensor_DD = pow(wfref->GetWp(knorm),2.)/4.; // W*W/4  
  double dens_U = pow(wfref->GetWp(knorm),2.)+pow(wfref->GetUp(knorm),2.);

  double sinth2=1.-costh*costh;
  
  
  double Ek=sqrt(MASSn*MASSn+knorm*knorm);
  
  double alpha_i=1-knorm*costh/Ek;
  double alpha_s=2.-alpha_i;
  
  double ps_min = MASSD*alpha_s/2.;
  double ps_plus = (MASSn*MASSn+knorm*knorm*sinth2)/ps_min;
  double Es=(ps_min+ps_plus)/2.;
  double ps_z=(-ps_min+ps_plus)/2.;
  double ps_norm=sqrt(ps_z*ps_z+knorm*knorm*sinth2);
  double costh_s= ps_z/ps_norm;
  double Einoff=MASSD-Es;
//   alpha_i=(Es+ps_z)/MASSD*2.;
  
  //Wx^2> massi^2
  double lowerlimit= -(-MASSN*MASSN-Q2+pow(Einoff,2.)-ps_norm*ps_norm+2.*(Einoff)*nu)/(2.*qvec*ps_norm);
  if(lowerlimit>costh_s) {res[0]=res[1]=res[2]=0.;return;}
  if(1.<lowerlimit){ res[0]=res[1]=res[2]=0.; return;}
//   double higherlimit= -(-MASSN*MASSN-Q2+pow(Einoff,2.)-ps_norm*ps_norm+2.*Einoff*nu)/(2.*qvec*ps_norm);
//   if(higherlimit<costh_s) {res[0]=res[1]=res[2]=0.;return;}
//   if(-1.>higherlimit){ res[0]=res[1]=res[2]=0.; return;}
  
  double piq=(Einoff*nu+ps_z*qvec);  //vec{pi}=-vec{pr} in PW!
  double nutilde=piq/MASSn;
  double xtilde=Q2/(2*piq);
  double F1p=0.,F2p=0.,F1n=0.,F2n=0.;
  strp.getF_xQ(F1p,F2p,1,xtilde,Q2);  
  strp.getF_xQ(F1n,F2n,0,xtilde,Q2);  
//   F1p=F1n=0.;
  double factor=(2.*(F1p+F1n)+knorm*knorm*sinth2/piq*(F2p+F2n))*(6*costh*costh-2.)+knorm*knorm/piq*(F2p+F2n)*sinth2*sinth2;

  double F_U_T = (2.*(F1p+F1n)+knorm*knorm*sinth2/piq*(F2p+F2n))*dens_U;
  double F_U_L = (-2.*(F1p+F1n)+2.*(F2p+F2n)*pow(MASSD*sqrt(1+gamma*gamma)/gamma-(Es*qvec-nu*ps_z)/sqrt(Q2),2.)/piq)*dens_U;
  double F_tensor_T = -3./2.*(2.*(F1p+F1n)+knorm*knorm*sinth2/piq*(F2p+F2n))*dens_tensor_tot*(6.*costh*costh-2.);
  double F_tensor_L = -3./2.*(-2.*(F1p+F1n)+2.*(F2p+F2n)*pow(MASSD*sqrt(1+gamma*gamma)/gamma-(Es*qvec-nu*ps_z)/sqrt(Q2),2.)/piq)*dens_tensor_tot*(6.*costh*costh-2.);
  double F_tensor_cosphi = 6.*(2.*knorm*sqrt(sinth2)/piq*(MASSD*sqrt(1+gamma*gamma)/gamma-(Es*qvec-nu*ps_z)/sqrt(Q2))*(F2p+F2n))*dens_tensor_tot*costh*sqrt(sinth2);
  double F_tensor_cos2phi = -3./2.*(knorm*knorm*sinth2/piq*(F2p+F2n))*dens_tensor_tot*sinth2;
  
 
  //last factor is normalization for lightcone 1/alpha_i converted to k momentum (k is spectator momentum here, be careful with signs!)
  factor*=3./4./(1+gamma*gamma)*knorm*knorm/alpha_i /*/2./(1.-(Es-knorm*costh)/MASSD)*/; 

  res[0]=dens_tensor_tot*factor;
  res[1]=dens_tensor_SD*factor;
  res[2]=dens_tensor_DD *factor;
  res[3]=-1./2./(1.+gamma*gamma)*knorm*knorm/alpha_i*(F_tensor_T+F_tensor_cos2phi); //b1
    
  res[4]=F_tensor_T*knorm*knorm/alpha_i;
  res[5]=F_tensor_L*knorm*knorm/alpha_i;
  res[6]=F_tensor_cosphi*knorm*knorm/alpha_i;
  res[7]=F_tensor_cos2phi*knorm*knorm/alpha_i;
  
  res[8]=F_U_T*knorm*knorm/alpha_i;
  res[9]=F_U_L*knorm*knorm/alpha_i;
  
  res[10]=(0.25+0.75*cos(2.*thetaq))*(res[4]+eps*res[5])+0.75*sin(2*thetaq)*sqrt(2.*eps*(1.+eps))*res[6]+0.75*(1.-cos(2.*thetaq))*eps*res[7];
  res[11]=res[8]+eps*res[9];
  
  
    
    
}
