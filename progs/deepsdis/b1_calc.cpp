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
#include <Utilfunctions.hpp>


#include <TDeuteron.h>
#include <numint/numint.hpp>


/*
 * @brief  declarations of the functions that we will later use as integrands, k_int integrates over the momentum of the nucleons inside the deuteron
 * 
 * @param[out] res results of the integrations 
      //ret[0] b1 total
      //ret[1] b1 sd part (of deuteron wf)
      //ret[2] b1 dd part (of deuteron wf)
      //ret[3] b1 total check
      //ret[4-7] tensor struc functions F_T that allow to calculate b1-4.
      //ret[8-9] FUUT,FUUL -> unpolarized structure functions
      //ret[10-11] azz asymmetry nominator, denominator
 * @param knorm [MeV] magnitude nucleon momenta
 * @param kcosth  [] cos(theta) of the spherical coordinate of the nucleon momentum
 * @param wfref deuteron wave function object
 * @param Q2 [MeV^2] photon Q^2
 * @param x  [] Bjorken x
 * @param nu [MeV] photon energy
 * @param qvec [MeV] photon momentum
 * @param gamma [] invariant that controls target mass corrections ~M^2/Q^2
 * @param eps [] invariant used in lepton kinematics
 * @param thetaq [rad] angle between beam and photon 3momenta in target lab grame
 * @param strp object that contains nucleon structure functions
 */
void k_int(numint::vector_d & res, double knorm, double kcosth, TDeuteron::Wavefunction *wfref,
          double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq, NucleonStructure &strp);

void y_int(numint::vector_d & res, double y, double pperp, TDeuteron::Wavefunction *wfref,
          double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq, NucleonStructure &strp);

void p_int(numint::vector_d & res, double pnorm, double costh, TDeuteron::Wavefunction *wfref,
          double Q2, double x, double nu, double qvec, double gamma, double eps, double thetaq, NucleonStructure &strp);

int main(int argc, char *argv[])
{

  string arg_names[5]={"exec name", 
                            "Q^2 [GeV^2]", 
                            "wave function",
                            "nucleon structure function parametrization",
                            "beam energy [GeV]"};

  std::cout << "Compiled from file: " << __FILE__ << std::endl;
  Bookkeep(argc,argv,arg_names);  



  double Q2 = atof(argv[1])*1.E06; // input in GeV^2, converted to MeV^2
  string wf = argv[2]; // deuteron wf parametrization
  double Ein= atof(argv[4])*1.E03; //beam energy in GeV, converted to MeV
  string nuclstruc = argv[3]; // Nucleon structure function parametrization
  
  // create deuteron wf object
  TDeuteron::Wavefunction *wfref;
  wfref = TDeuteron::Wavefunction::CreateWavefunction(wf);

  // create nucleon structure function object (gives you F2n, F2p)
  NucleonStructure strp(nuclstruc);
  
  
  // structure that we need to perform the numerical integration, contains the parameters and function we integrate over [abstract definition]
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
  
  // for loop over Bjorken x values
  for(int i=1;i<200;i++){
    double x=0.005*i; //Bjorken x, invariant for deuteron [0-2]
    double nu=Q2/(2.*MASSD*x); //photon energy rest frame [MeV]
    double qvec=sqrt(Q2+nu*nu); //photon momentum norm rest frame [MeV]
    double gamma=sqrt(Q2)/nu; // gamma, see paper (invariant)
//     double gamma=0.;
    double Eout=Ein-nu; //scattered electron energy [MeV]
    double y=nu/Ein; // inelasticity [0-1] (invariant)
    if(/*y<1.*/1){
      double eps=(1-y-gamma*gamma*y*y/4.)/(1-y+y*y/2+gamma*gamma*y*y/4.); //epsilon, see paper (invariant)
      double thetae=asin(sqrt(Q2/4/Ein/Eout))*2; //angle between beam and scattered electron 3momenta [radians]
      double thetaq=acos((Ein*Ein+qvec*qvec-Eout*Eout)/2/Ein/qvec);  //angle between beam and virutal photon 3momenta [radians]

      //cout << 2.*x << " " << thetaq*RADTODEGR << " " << (0.25+0.75*cos(thetaq)) << " " << 0.75*sin(2*thetaq) << " " << 0.75*(1.-cos(2.*thetaq)) << " " ;
      //cout << eps << " " << sqrt(2.*eps*(1.+eps)) << " ";
      cout << 2. *x << " ";
      double thetaq1 = -acos(-1.0/3.0)/2.0;
      NucleonStructure stp(nuclstruc);  //create instance nucleon structure 

      numint::array<double,2> lower = {{0.,-1.}}; //first index momentum norm [MeV], second is costheta nucleon
      numint::array<double,2> upper = {{600.,1.}};
      //Initialize an instance of the previously defined structure used in the numerical integration
      // and then assign all the variables it contains
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
      
      //nucleon structure functions at that x,Q^2 for reference (with the deuteron results)
      double F1p=0.,F2p=0.,F1n=0.,F2n=0.;
      if(x<0.5){
        strp.getF_xQ(F1p,F2p,1,x*MASSD/MASSn,Q2);  
        strp.getF_xQ(F1n,F2n,0,x*MASSD/MASSn,Q2);  
      }
      //initialize 2D integrand
      numint::mdfunction<numint::vector_d,2> mdf;
      mdf.func = &Ftor_b1::exec;
      mdf.param = &F;
      //we will generate 12 results
      numint::vector_d ret(12,0.);
      F.f=k_int;  //k_int is used as integrand
      int res=90;
      unsigned count=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E05,ret,count,0); //integration carried out
      //cout << 2.*x << " " << ret[0] << " " <<   ret[1] << " " << ret[2] << " " << ret[3] << " "<< ret[8]/2. << " " << 2.*(F1p+F1n) << " " << 1./(1+gamma*gamma) <</*" " << sqrt(2./3.)*ret[4]/ret[8]*2.*(F1p+F1n)*-3./2. << " " << sqrt(2./3.)*ret[10]/ret[11]*2.*(F1p+F1n)*-3./2. << " " << sqrt(2./3.)*ret[10]/ret[11] <<*/;
      //cout << ret[0] << " " << ret[3] << " "<< ret[8]/2. << " " << 2.*(F1p+F1n) << " "  << ret[4] << " " << ret[5] << " " << ret[6] << " " << ret[7]  << " ";
      

      double A=-1./2./(1.+gamma*gamma)*(ret[4]+ret[7]);
      double D=-2.*ret[7]/gamma/gamma;
      double C=(-2.*ret[6]/gamma-2.*D)/(1.+gamma*gamma);
      double B=(2.*(1+gamma*gamma)*A-(1.+gamma*gamma)*C-D-ret[5])/pow(1.+gamma*gamma,2.);

      double b1=A;
      double b2=x*(B+C+D);
      double b3=x*(B+C-2.*D)/3.;
      double b4=x*(B-2.*C+D)/3.;

      //b1 actual
      //cout << ret[0];
      //cout << " " << ret[3]<< " " << b1<< " " << b2<< " " << endl;
      
      double numerator = (0.25+0.75*cos(2.*thetaq))*(ret[4]+eps*ret[5])+0.75*sin(2*thetaq)*sqrt(2.*eps*(1.+eps))*ret[6]+0.75*(1.-cos(2.*thetaq))*eps*ret[7];
      //numerator for theta = -arccos(-1/3)
      double numerator1 = (0.25+0.75*cos(2.*thetaq1))*(ret[4]+eps*ret[5])+0.75*sin(2*thetaq1)*sqrt(2.*eps*(1.+eps))*ret[6]+0.75*(1.-cos(2.*thetaq1))*eps*ret[7];

      //new structure functions

      //double FuTlll= -2./3. * pow(gamma,2.) *(1.+pow(gamma,2.));
      //double FuTllt= -(2.+5./3.*pow(gamma,2.));
      //double Fcoslt = -gamma * (1.+pow(gamma,2.)/3.);
      //double Fcostt = -pow(gamma,2.)/3.;



     //coefficients for structure functions


      double FuTlll= -2./3.*pow(gamma,4.);  
      double FuTllt= -(2.+5./3.*pow(gamma,2.));
      double Fcoslt = -gamma * (1.+pow(gamma,2.)/3.);
      double Fcostt = -pow(gamma,2.)/3.;

      double Tll = (3.0 *cos(2*thetaq)+1.)/12.0;

      double Tll1 = (3.0 *cos(2*0)+1.)/12.0;

      double tltcos= (sin(2*thetaq))/4.0;

      double tltcos1 = (sin(2*0))/4.0;

      double tttcos = (1.-cos(2.*thetaq))/4.0;

      double tttcos1 = (1.-cos(2.*0))/4.0;
      //Alan
      double numeratorTest = 2.*(Tll*(FuTllt+eps*FuTlll)+tltcos*(pow(2.*eps*(1.+eps),0.5)*Fcoslt)+tttcos*eps*Fcostt);

      double numeratorTest1 = 2.*(Tll1*(FuTllt+eps*FuTlll)+tltcos1*(pow(2.*eps*(1.+eps),0.5)*Fcoslt)+tttcos1*eps*Fcostt);

      double denomTest = ret[8]*(1.+pow(gamma,2.0));
      //brandon
      double numeratorNew = (1./4.+3./4.* cos(2.*thetaq))*((2.0*(1.0+pow(gamma,2.0))*(eps-1))+2.0*(1.0/3.0*(pow(gamma,2.0)/2.0)-eps*(1.0+pow(gamma,2.0))-eps))-3.0/4.0*sin(2*thetaq)*(pow((2*eps*(1+eps)),1.0/2.0))*(1.0+pow(gamma,2.0)*gamma/6.0+1./3.*gamma)-3.0/4.*(1-cos(2*thetaq))/6.*eps*pow(gamma,2.0);
      double numeratorNew_0= (1./4.+3./4.* cos(2.*0))*((2.0*(1.0+pow(gamma,2.0))*(eps-1))+2.0*(1.0/3.0*(pow(gamma,2.0)/2.0)-eps*(1.0+pow(gamma,2.0))-eps))-3.0/4.0*sin(2*0)*(pow((2*eps*(1+eps)),1.0/2.0))*(1.0+pow(gamma,2.0)*gamma/6.0+1./3.*gamma)-3.0/4.*(1-cos(2*0))/6.*eps*pow(gamma,2.0);


      double numerator_0 = (ret[4]+eps*ret[5]);
      double denominator = ret[8]+eps*ret[9];
      double Azz = 2./3.*numerator/denominator;
      double Azz1 = 2./3. * numerator1/denominator;

      double b1_estimateTestE = Azz * denomTest/numeratorTest;
      double b1_estimateTest = Azz * denomTest/numeratorTest1;
      double b1_estimateNew= Azz* denomTest/numeratorNew;
      double Azz_0 = 2./3.*numerator_0/denominator;
      double b1_estimateTestE0 = Azz_0 * denominator/numeratorTest;
      double b1_estimate_New0= Azz_0*denomTest/numeratorNew_0;
      double b1_estimate = -3/2.*Azz*ret[8]/2.;
      double b1_estimate1 = -3/2.*Azz1*ret[8]/2.;
      double b1_estimateTest0= Azz_0 * denomTest/numeratorTest1;
      double b1_estimate_0 = -3/2.*Azz_0*ret[8]/2.;


  
      //tensor structure functions and coefficients (polarized along electron direction)
      /*
      double TF1=2./3.* (0.25+0.75*cos(2.*thetaq))*(ret[4]);
      double TF2=2./3.* (0.25+0.75*cos(2.*thetaq))*(eps*ret[5]);
      double TF3= 2./3.* (0.75*sin(2*thetaq)*sqrt(2.*eps*(1.+eps))*ret[6]);
      double TF4=2./3.* (0.75*(1.-cos(2.*thetaq))*eps*ret[7]);

      //tensor structure functions and coefficients (polarized along virtual photon direction)
      double TF1_0=2./3.* (0.25+0.75*cos(2.*0))*(ret[4]);
      double TF2_0=2./3.* (0.25+0.75*cos(2.*0))*(eps*ret[5]);
      double TF3_0= 2./3.* (0.75*sin(2*0)*sqrt(2.*eps*(1.+eps))*ret[6]);
      double TF4_0=2./3.* (0.75*(1.-cos(2.*0))*eps*ret[7]);
      */

    //tensor structure functions and coefficients (polarized along electron direction)
      double TF1= (0.25+0.75*cos(2.*thetaq))*(ret[4]);
      double TF2= (0.25+0.75*cos(2.*thetaq))*(eps*ret[5]);
      double TF3=  (0.75*sin(2*thetaq)*sqrt(2.*eps*(1.+eps))*ret[6]);
      double TF4= (0.75*(1.-cos(2.*thetaq))*eps*ret[7]);

      //tensor structure functions and coefficients (polarized along virtual photon direction)
      double TF1_0= (0.25+0.75*cos(2.*0))*(ret[4]);
      double TF2_0= (0.25+0.75*cos(2.*0))*(eps*ret[5]);
      double TF3_0= (0.75*sin(2*0)*sqrt(2.*eps*(1.+eps))*ret[6]);
      double TF4_0= (0.75*(1.-cos(2.*0))*eps*ret[7]);

      double b1C= 1./3. * (1+pow(gamma,2.0))*(-1.+eps)*(1+3*cos(2*thetaq))*b1;
      double b2C= 1./(36*x)*(pow(gamma,2.0)-(6.+9.*pow(gamma,2.0)+2*pow(gamma,4.0))*eps-3.*(6*eps+2*pow(gamma,4.0)*eps+pow(gamma,2.0)*(-1.+5.*eps))*cos(2.*thetaq)+3.*pow(2.,0.5)*pow(gamma,2.0)*(3.+pow(gamma,2.)*pow(eps*(1+eps),0.5)*sin(2.*thetaq)))*b2;
      double b3C= pow(gamma,2.0)/(12.*x)*(1.+eps+2.*pow(gamma,2.)*eps+3.*(1.+5.*eps+2.*pow(gamma,2.)*eps)*cos(2.*thetaq)+6.*pow(2.,0.5)*pow(eps*(1+eps),0.5)*sin(2*thetaq))*b3;
      double b4C= pow(gamma,2.0)/(12.*x)*(2.*(1.+pow(gamma,2.0))*(eps+3.*eps*cos(2*thetaq)+3.*pow(2,.5)*pow(eps*(1+eps),0.5)*cos(thetaq)*sin(thetaq)))*b4;


      double b1C0= 1./3. * (1+pow(gamma,2.0))*(-1.+eps)*(1+3*cos(2*0))*b1;
      double b2C0= 1./(36*x)*(pow(gamma,2.0)-(6.+9.*pow(gamma,2.0)+2*pow(gamma,4.0))*eps-3.*(6*eps+2*pow(gamma,4.0)*eps+pow(gamma,2.0)*(-1.+5.*eps))*cos(2.*0)+3.*pow(2.,0.5)*pow(gamma,2.0)*(3.+pow(gamma,2.)*pow(eps*(1+eps),0.5)*sin(2.*0)))*b2;
      double b3C0= pow(gamma,2.0)/(12.*x)*(1.+eps+2.*pow(gamma,2.)*eps+3.*(1.+5.*eps+2.*pow(gamma,2.)*eps)*cos(2.*0)+6.*pow(2.,0.5)*pow(eps*(1+eps),0.5)*sin(2*0))*b3;
      double b4C0= pow(gamma,2.0)/(12.*x)*(2.*(1.+pow(gamma,2.0))*(eps+3.*eps*cos(2*0)+3.*pow(2,.5)*pow(eps*(1+eps),0.5)*cos(0)*sin(0)))*b4;



    
    
     cout<< numerator <<" " << numerator_0<<" " << TF1<<" " <<TF2<<" " <<TF3<<" " <<TF4<<" " <<TF1_0 <<" " <<TF2_0<<" " << TF3_0<<" " << TF4_0 ;
     cout<<" "<<b1C<<" " <<b2C<<" " <<b3C<<" " <<b4C<<" " <<b1C0<<" " <<b2C0<<" " <<b3C0<<" " <<b4C0<<endl;

//     cout<<gamma<<endl; 


      //cout<< Azz <<" " << Azz_0<<" " << ret[4]<<" " <<ret[5]<<" " <<ret[6]<<" " <<ret[7]<<" " <<b1 <<" " <<b2<<" " << b3<<" " << b4<<endl;

      // explanation on what is what
      // original 2 estimates                        test on different direction             Brandons                                          new approximations (test= theta=0)     (E= theta is finite)
      //cout <<  b1_estimate << " " << b1_estimate_0 << " " << b1_estimate1 <<" "<< b1_estimateNew << " "<< b1_estimate_New0 << " " << b1_estimateTest << " "<< b1_estimateTest0<< " " << b1_estimateTestE << " "<< numeratorTest<< " " << numeratorTest1  << " "<< thetaq<< " "<< eps<< " "<< gamma<< endl;

      //cout<< -3./2.*ret[8]/2. << " "<< denomTest/numeratorTest << " " << denomTest/numeratorTest1 << " "<< denomTest/numeratorNew<<  " "<< denomTest/numeratorNew_0<<endl;

      //cout<< " "<< Azz << " " << Azz_0 << endl;
      // cout<<" "<< Azz  << endl;
      //integral result output
      //ret[0] b1 total
      //ret[1] b1 sd
      //ret[2] b1 dd
      //ret[3] b1 total check
      //ret[4-7] tensor struc functions
      //ret[8-9] FUUT,FUUL
      //ret[10-11] azz num, denom
      

      
//pint and yint part (shunzo's formulas in the b1 PRD paper), used for comparison
      
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

  double F_U_T = (2.*(F1p+F1n)+knorm*knorm*sinth2/piq*(F2p+F2n))*dens_U;  //second term does not survive bjorken limit
  double F_U_L = (-2.*(F1p+F1n)+2.*(F2p+F2n)*pow(MASSD*sqrt(1+gamma*gamma)/gamma-(Es*qvec-nu*ps_z)/sqrt(Q2),2.)/piq)*dens_U;
  double F_tensor_T = -3./2.*(2.*(F1p+F1n)+knorm*knorm*sinth2/piq*(F2p+F2n))*dens_tensor_tot*(6.*costh*costh-2.);
  double F_tensor_L = -3./2.*(-2.*(F1p+F1n)+2.*(F2p+F2n)*pow(MASSD*sqrt(1+gamma*gamma)/gamma-(Es*qvec-nu*ps_z)/sqrt(Q2),2.)/piq)*dens_tensor_tot*(6.*costh*costh-2.);
  double F_tensor_cosphi = 6*(2.*knorm*sqrt(sinth2)/piq*(MASSD*sqrt(1+gamma*gamma)/gamma-(Es*qvec-nu*ps_z)/sqrt(Q2))*(F2p+F2n))*dens_tensor_tot*costh*sqrt(sinth2);
  double F_tensor_cos2phi = -3./2.*(knorm*knorm*sinth2/piq*(F2p+F2n))*dens_tensor_tot*sinth2;
  
  //last factor is normalization for lightcone 1/alpha_i converted to k momentum (k is spectator momentum here, be careful with signs!)
//   alpha_i=1.;
//    gamma=0.;
  factor*=3./4.*knorm*knorm/alpha_i/(1+gamma*gamma); // /(1.-(Ek+knorm*costh)/MASSD); 

  res[0]=dens_tensor_tot*factor;  // without 1+ gamma^2 factor!
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
    