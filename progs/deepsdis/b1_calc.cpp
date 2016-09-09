#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <constants.hpp>
#include <NuclStructure.hpp>
#include <TDeuteron.h>
#include <numint/numint.hpp>


void k_int(numint::vector_d & res, double knorm, double kcosth, TDeuteron::Wavefunction *wfref,
	      string nuclstruc, double Q2, double x, double nu, double qvec, double gamma);


int main(int argc, char *argv[])
{
  double Q2 = atoi(argv[1])*1.E06; // parse from argv or something
  string wf = argv[2];
  string nuclstruc = argv[3];
  
  TDeuteron::Wavefunction *wfref;
  wfref = TDeuteron::Wavefunction::CreateWavefunction(wf);
//   for(int i=1;i<100;i++){
//     double x=0.01*i;
//     NuclStructure s1(1,Q2,x,0,"SLAC");
//     NuclStructure s2(1,Q2,x,0,"HMRS");
//     cout << x << " " << s1.getF1() << " " << s2.getF1() << " " << s1.getF2() << " " << s2.getF2() << endl;
//   }
//   exit(1);
    
  
  
  
  struct Ftor_b1 {

    /*! integrandum function */
    static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
      Ftor_b1 &p = * (Ftor_b1 *) param;
      p.f(ret,x[0],x[1], p.wfref,p.nuclstruc,p.Q2,p.x,p.nu,p.qvec,p.gamma);
    }
    TDeuteron::Wavefunction *wfref;
    string nuclstruc;
    double Q2;
    double x;
    double nu;
    double qvec;
    double gamma;
    
    
    void (*f)(numint::vector_d & res, double knorm, double kcosth, TDeuteron::Wavefunction *wfref,
	      string nuclstruc, double Q2, double x, double nu, double qvec, double gamma);
  };
  
  for(int i=1;i<80;i++){
    double x=0.0125*i;
    double nu=Q2/(2.*MASSD*x);
    double qvec=sqrt(Q2+nu*nu);
    double gamma=sqrt(Q2)/nu;
//     double gamma=0.;
    
    numint::array<double,2> lower = {{0.,-1.}};
    numint::array<double,2> upper = {{600.,1.}};
    Ftor_b1 F;
    F.wfref=wfref;
    F.nuclstruc=nuclstruc;
    F.Q2=Q2;
    F.x=x;
    F.nu=nu;
    F.qvec=qvec;
    F.gamma=gamma;
    
    numint::mdfunction<numint::vector_d,2> mdf;
    mdf.func = &Ftor_b1::exec;
    mdf.param = &F;
    numint::vector_d ret(3,0.);
    F.f=k_int;
    int res=90;
    unsigned count=0;
//     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
    res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E05,ret,count,0);
    cout << 2.*x << " " << ret[0] << " " << ret[1] << " " << ret[2] << endl;

    
  }
}

void k_int(numint::vector_d & res, double knorm, double costh, TDeuteron::Wavefunction *wfref,
	      string nuclstruc, double Q2, double x, double nu, double qvec, double gamma){
    

  res=numint::vector_d(3,0.);
 /* res[0]=res[1]=res[2]=pow(wfref->GetWp(knorm),2.)+pow(wfref->GetUp(knorm),2.)*knorm*knorm;
  cout << knorm << " " << costh << endl;
  return;
 */ 
  res[0]=wfref->GetUp(knorm)*wfref->GetWp(knorm)/sqrt(2.)+pow(wfref->GetWp(knorm),2.)/4.;  // UW/sqrt(2)+W^2/4
  res[1]=wfref->GetUp(knorm)*wfref->GetWp(knorm)/sqrt(2.); // UW/sqrt(2)
  res[2]=pow(wfref->GetWp(knorm),2.)/4.; // W*W/4  
 
  double sinth2=1.-costh*costh;
  double Es=sqrt(MASSn*MASSn+knorm*knorm);
  double Einoff=MASSD-Es;
  //Wx^2> massi^2
  double lowerlimit= -(-MASSn*MASSn-Q2+pow(MASSD-Es,2.)-knorm*knorm+2.*(MASSD-Es)*Q2/(2.*MASSD*x))/(2.*qvec*knorm);
  if(lowerlimit>costh) {res[0]=res[1]=res[2]=0.;return;}
  if(1.<lowerlimit){ res[0]=res[1]=res[2]=0.; return;}
  
  double piq=(Einoff*nu+knorm*costh*qvec);  //vec{pi}=-vec{pr} in PW!
  double nutilde=piq/MASSn;
  double mi_off = Einoff*Einoff-knorm*knorm; //effective mass off-shell nucleon SQUARED
  double xtilde=Q2/(2*piq);
  double nuoffshell=(mi_off-MASSn*MASSn+2.*piq)/(2.*MASSn); //(m_i+qoffshell)^2=(p_i+q)^2
  double xoffshell=Q2/(2.*MASSn*nuoffshell); //xoffshell consistent with Q^2,m_i,W
  double W_sq=-Q2+2.*piq+mi_off;
//   cout << x << " " << knorm << " " << costh << " " << W_sq << " " << lowerlimit << " " << -(-W_sq-Q2+pow(MASSD-Es,2.)-knorm*knorm+2.*(MASSD-Es)*Q2/(2.*MASSD*x))/(2.*qvec*knorm) << endl;
//   if(W_sq<0.){ res[0]=res[1]=res[2]=0.; return;}
  NuclStructure str_p(1,Q2,xtilde,W_sq,nuclstruc);
  NuclStructure str_n(0,Q2,xtilde,W_sq,nuclstruc);
//   cout << knorm << " " << costh << " " <<  x << " " << xoffshell << " " << sqrt(W_sq) << " " << Q2 << " " << nuoffshell << " " << sqrt(mi_off) << " " << Einoff << endl;
  double F1p,F2p,F1n,F2n;
  if(!(str_p.getName().compare("SLAC"))||!(str_p.getName().compare("CTEQ"))){
    F2p=str_p.getF2();
    F2n=str_n.getF2();
    double R=NuclStructure::getr1998(x,Q2);
    F1p=F2p/2./xoffshell/(R+1)*(1+gamma*gamma);
    F1n=F2n/2./xoffshell/(R+1)*(1+gamma*gamma);
  }
  else {str_p.getF(F1p,F2p); str_n.getF(F1n,F2n);}
//   F2p=F2n=0.;
  double factor=(2.*(F1p+F1n)+knorm*knorm*sinth2/piq*(F2p+F2n))*(6*costh*costh-2.)+knorm*knorm/piq*(F2p+F2n)*sinth2*sinth2;
  
  factor*=3./4./(1.+gamma*gamma)*knorm*knorm/(1.+knorm*costh/Es) /*/2./(1.-(Es-knorm*costh)/MASSD)*/;
  
  res[0]*=factor;
  res[1]*=factor;
  res[2]*=factor;
//   if(std::isnan(res[0]))  cout << knorm<< " " << costh << " " << xoffshell << " " << W_sq << " " << lowerlimit << " " << endl;
//   cout << knorm << " " << costh << " " << F2p << " " << F2n << endl;
    
    
    
    
}
