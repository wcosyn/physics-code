//calculate b1 from the tagged structure functions

//calculate Azz from tagged structure functions, accounting for virtual photon angle etc.

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

void y_int(numint::vector_d & res, double y, double pperp, TDeuteron::Wavefunction *wfref,
          double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq, NucleonStructure &strp);

void p_int(numint::vector_d & res, double pnorm, double costh, TDeuteron::Wavefunction *wfref,
          double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq, NucleonStructure &strp);

int main(int argc, char *argv[])
{
  double Q2 = atof(argv[1])*1.E06; 
  string wf = argv[2];
  double Ein= atof(argv[4])*1.E03; //beam energy in GeV
  string nuclstruc = argv[3];
  
  TDeuteron::Wavefunction *wfref;
  wfref = TDeuteron::Wavefunction::CreateWavefunction(wf);

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
  
  for(int i=1;i<200;i++){
    double x=0.005*i;
    double nu=Q2/(2.*MASSD*x);
    double qvec=sqrt(Q2+nu*nu);
    double gamma=sqrt(Q2)/nu;
//     double gamma=0.;
    double Eout=Ein-nu;
    double y=nu/Ein;
    if(/*y<1.*/1){
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
      if(x<0.5){
        strp.getF_xQ(F1p,F2p,1,x*MASSD/MASSn,Q2);  
        strp.getF_xQ(F1n,F2n,0,x*MASSD/MASSn,Q2);  
      }
      numint::mdfunction<numint::vector_d,2> mdf;
      mdf.func = &Ftor_b1::exec;
      mdf.param = &F;
      numint::vector_d ret(12,0.);
      F.f=k_int;
      int res=90;
      unsigned count=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E05,ret,count,0);
      cout << 2.*x << " " << ret[0] << " " <<   ret[1] << " " << ret[2] << " " << ret[3] << " "<< ret[8]/2. << " " << 2.*(F1p+F1n) << " " << 1./(1+gamma*gamma) <</*" " << sqrt(2./3.)*ret[4]/ret[8]*2.*(F1p+F1n)*-3./2. << " " << sqrt(2./3.)*ret[10]/ret[11]*2.*(F1p+F1n)*-3./2. << " " << sqrt(2./3.)*ret[10]/ret[11] <<*/ endl;
      //ret[0] b1 total
      //ret[1] b1 sd
      //ret[2] b1 dd
      //ret[3] b1 total check
      //ret[4-7] tensor struc functions
      //ret[8-9] FUUT,FUUL
      //ret[10-11] azz nom, denom
      

      
//pint and yint part (shunzo's formulas)
      
//       numint::array<double,2> lower = {{2.*x,0.}};
//       numint::array<double,2> upper = {{2.,2.E03}};
// 
//             double uplimit=sqrt((MASSn-2.22)*2.*(MASSD-MASSn+2.22));


//       double uplimit=sqrt(MASSD*MASSD-MASSn*MASSn);
//       
//       numint::array<double,2> lower = {{0.,-1.}};
//       numint::array<double,2> upper = {{uplimit,1.}};
//       Ftor_b1 F;
//       F.wfref=wfref;
//       F.Q2=Q2;
//       F.x=x;
//       F.nu=nu;
//       F.qvec=qvec;
//       F.gamma=gamma;
//       F.eps=eps;
//       F.thetaq=thetaq;
//       F.strp=strp;
//       
//       double F1p=0.,F2p=0.,F1n=0.,F2n=0.;
// //       strp.getF_xQ(F1p,F2p,1,x*MASSD/MASSn,Q2);  
// //       strp.getF_xQ(F1n,F2n,0,x*MASSD/MASSn,Q2);  
// 
//       numint::mdfunction<numint::vector_d,2> mdf;
//       mdf.func = &Ftor_b1::exec;
//       mdf.param = &F;
//       numint::vector_d ret(3,0.);
//       F.f=p_int;
//       int res=90;
//       unsigned count=0;
//   //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
//       res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E05,ret,count,0);
//       cout << 2.*x << " " << ret[0] << " " <<   ret[1] << " " << ret[2] << endl;
//       
      
    }    
    
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
//   xtilde=2.*x/alpha_i;
  double F1p=0.,F2p=0.,F1n=0.,F2n=0.;
  if(xtilde<1.){
    strp.getF_xQ(F1p,F2p,1,xtilde,Q2);  
    strp.getF_xQ(F1n,F2n,0,xtilde,Q2);  
  }
//   F1p=F1n=0.;
  double factor=(2.*(F1p+F1n)+knorm*knorm*sinth2/piq*(F2p+F2n))*(6*costh*costh-2.)+knorm*knorm/piq*(F2p+F2n)*sinth2*sinth2;

  double F_U_T = (2.*(F1p+F1n)+knorm*knorm*sinth2/piq*(F2p+F2n))*dens_U;
  double F_U_L = (-2.*(F1p+F1n)+2.*(F2p+F2n)*pow(MASSD*sqrt(1+gamma*gamma)/gamma-(Es*qvec-nu*ps_z)/sqrt(Q2),2.)/piq)*dens_U;
  double F_tensor_T = -sqrt(3./2.)*(2.*(F1p+F1n)+knorm*knorm*sinth2/piq*(F2p+F2n))*dens_tensor_tot*(6.*costh*costh-2.);
  double F_tensor_L = -sqrt(3./2.)*(-2.*(F1p+F1n)+2.*(F2p+F2n)*pow(MASSD*sqrt(1+gamma*gamma)/gamma-(Es*qvec-nu*ps_z)/sqrt(Q2),2.)/piq)*dens_tensor_tot*(6.*costh*costh-2.);
  double F_tensor_cosphi = 2.*sqrt(6.)*(2.*knorm*sqrt(sinth2)/piq*(MASSD*sqrt(1+gamma*gamma)/gamma-(Es*qvec-nu*ps_z)/sqrt(Q2))*(F2p+F2n))*dens_tensor_tot*costh*sqrt(sinth2);
  double F_tensor_cos2phi = -sqrt(3./2.)*(knorm*knorm*sinth2/piq*(F2p+F2n))*dens_tensor_tot*sinth2;
  
  //last factor is normalization for lightcone 1/alpha_i converted to k momentum (k is spectator momentum here, be careful with signs!)
//   alpha_i=1.;
//    gamma=0.;
  factor*=3./4.*knorm*knorm/alpha_i/(1+gamma*gamma); // /(1.-(Ek+knorm*costh)/MASSD); 

  res[0]=dens_tensor_tot*factor;  // without 1+ gamma^2 factor!
  res[1]=dens_tensor_SD*factor;
  res[2]=dens_tensor_DD *factor;
  res[3]=-sqrt(3./8.)/(1.+gamma*gamma)*knorm*knorm/alpha_i*(F_tensor_T+F_tensor_cos2phi); //b1
    
  res[4]=F_tensor_T*knorm*knorm/alpha_i;
  res[5]=F_tensor_L*knorm*knorm/alpha_i;
  res[6]=F_tensor_cosphi*knorm*knorm/alpha_i;
  res[7]=F_tensor_cos2phi*knorm*knorm/alpha_i;
  
  res[8]=F_U_T*knorm*knorm/alpha_i;
  res[9]=F_U_L*knorm*knorm/alpha_i;
  
  res[10]=(0.25+0.75*cos(thetaq))*(res[4]+eps*res[5])+0.75*sin(2*thetaq)*sqrt(2.*eps*(1.+eps))*res[6]+0.75*(1.-cos(2.*thetaq))*eps*res[7];
  res[11]=res[8]+eps*res[9];
  
    
    
}



void y_int(numint::vector_d & res, double y, double pt, TDeuteron::Wavefunction *wfref,
	      double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq,
              NucleonStructure &strp){
  
  res=numint::vector_d(3,0.);
  double F1p=0.,F2p=0.,F1n=0.,F2n=0.;
  double xtilde=2.*x/y;
  if(xtilde<1.){
    strp.getF_xQ(F1p,F2p,1,xtilde,Q2);  
    strp.getF_xQ(F1n,F2n,0,xtilde,Q2);  
  }
  double result=0.5*(F1p+F1n)/y;
  double nu_q=1./sqrt(1+gamma*gamma);
  double massamin1=MASSD-MASSn+2.22;
  double xx1=MASSn*y-(1+gamma*gamma)*(massamin1);
  double xx2=MASSn*MASSn*y*y+(1+gamma*gamma)*(pt*pt-2.*massamin1*(MASSn-2.22));
  if((xx1*xx1-xx2)<0.){res=numint::vector_d(3,0.);return;}
  double p0=xx1+sqrt(xx1*xx1-xx2);
  double pz=(p0-MASSn*y)*nu_q;
  double pnorm=sqrt(pz*pz+pt*pt);
  double cc1=3.*pz*pz/pnorm/pnorm-1.;
  double a0=3./4./sqrt(2.)/PI;
  double a2=3./16./PI;
  result*=2.*PI*y/(1/(MASSn*nu_q)+pz/MASSn/massamin1);
  result*=pt*cc1;
  
  res[0]=result*(a0*wfref->GetUp(pnorm)*wfref->GetWp(pnorm)+a2*pow(wfref->GetWp(pnorm),2.));
  res[1]=result*(a0*wfref->GetUp(pnorm)*wfref->GetWp(pnorm));
  res[2]=result*(a2*pow(wfref->GetWp(pnorm),2.));
 }
  
void p_int(numint::vector_d & res, double pnorm, double costh, TDeuteron::Wavefunction *wfref,
	      double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq,
              NucleonStructure &strp){
  
  res=numint::vector_d(3,0.);
  double F1p=0.,F2p=0.,F1n=0.,F2n=0.;
  double pz=pnorm*costh;
  double pt=pnorm*sqrt(1.-costh*costh);
  
  double massamin1=MASSD-MASSn+2.22;
  double p0=MASSn-2.22-pnorm*pnorm/2./(massamin1);
//   double p0=MASSD-sqrt(MASSn*MASSn+pnorm*pnorm);
  double y=(p0*nu-pz*qvec)/MASSn/nu;
//   cout << pt << " " << pz << " " << p0 << " " << y << " " << 2.*x/y << endl;
  if(y<2.*x) return;
  double xtilde=2.*x/y;
  if(xtilde<1.){
    strp.getF_xQ(F1p,F2p,1,xtilde,Q2);  
    strp.getF_xQ(F1n,F2n,0,xtilde,Q2);  
  }
  double result=0.5*(F1p+F1n)/y;
  result*=2.*PI*y;
  double cc1=3.*costh*costh-1.;
  double a0=3./4./sqrt(2.)/PI;
  double a2=3./16./PI;
//   double alpha_i=(p0-pz)/MASSD*2.;
//   result/=alpha_i;
  result*=pnorm*pnorm*cc1;
  
  res[0]=result*(a0*wfref->GetUp(pnorm)*wfref->GetWp(pnorm)+a2*pow(wfref->GetWp(pnorm),2.));
  res[1]=result*(a0*wfref->GetUp(pnorm)*wfref->GetWp(pnorm));
  res[2]=result*(a2*pow(wfref->GetWp(pnorm),2.));
 }
    