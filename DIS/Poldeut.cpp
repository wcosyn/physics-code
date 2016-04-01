#include "Poldeut.hpp"
#include <FourVector.h>
#include <NuclStructure.hpp>
#include <TVector3.h>
#include <cassert>
#include <TDeuteron.h>
#include <constants.hpp>

using namespace std;

Poldeut::Poldeut(string wf_name, string struc_name):
strucname(struc_name),
wfname(wf_name),
s_wave(1),
d_wave(1),
melosh(1)
{
  wfref = TDeuteron::Wavefunction::CreateWavefunction(wfname);
  for(int i=0;i<=1000;i++){
    wf.AddUp(i,wfref->GetUp(i));
    wf.AddWp(i,wfref->GetWp(i));
  }

  
  
}


Poldeut::~Poldeut(){
  delete wfref;
}

Poldeut::Poldeut(const Poldeut& rhs){
  strucname=rhs.strucname;
  wfname=rhs.wfname;
  wfref = TDeuteron::Wavefunction::CreateWavefunction(wfname);
  wf=rhs.wf;
  s_wave=rhs.s_wave;
  d_wave=rhs.d_wave;
  melosh=rhs.melosh;

}


Poldeut& Poldeut::operator=(const Poldeut& rhs){
  if(this!=&rhs) { // avoid self-assignment
    strucname=rhs.strucname;
    wfname=rhs.wfname;
    wfref = TDeuteron::Wavefunction::CreateWavefunction(wfname);
    wf=rhs.wf;
    s_wave=rhs.s_wave;
    d_wave=rhs.d_wave;
    melosh=rhs.melosh;
  }
  return *this;

}


void Poldeut::calc_Double_Asymm(double x, double Q2, double y, double alpha_s, double pt, bool proton, double &F_LSL, double &F_U){
  
  double pd_plus =MASSD;
  FourVector<double> pd_mu(MASSD,0.,0.,0.);
  double ps_plus = alpha_s*pd_plus/2.;
  double ps_min = (MASSn*MASSn+pt*pt)/ps_plus;
  double Es = (ps_plus + ps_min)*0.5;
  double psz =(ps_plus - ps_min)*0.5;
  double psnorm  = sqrt(psz*psz + pt*pt);
  FourVector<double> ps_mu(Es,pt,0.,psz);
  FourVector<double> pi_mu = pd_mu - ps_mu;
  double pi_plus=pd_plus-ps_plus;
  double pi_perp = -pt;
  
  double nu = 2.*Q2/MASSD/x;
  double qvec = sqrt(Q2+nu*nu);
  FourVector<double> q_mu(nu,0.,0.,qvec);
  FourVector<double> el_mu(qvec,0.,0.,nu);
  el_mu*=1./sqrt(Q2); //longitudinal basis vector
  FourVector<double> et_mu(0.,1.,0.,0.);
  
  double gamma_d = sqrt(Q2)*MASSD/(pd_mu*q_mu);
  double epsilon=(1.-y-pow(gamma_d*y/2.,2.))/(1.-y+y*y/2.+pow(gamma_d*y/2.,2.));
  
  
  double alpha_i = 2.-alpha_s;
  double k_perp=-pt;  //minus sign!!
  double Ek = sqrt((MASSn*MASSn+k_perp*k_perp)/alpha_i/alpha_s);
  double kfz=(alpha_i-1.)*Ek;
  double knorm = sqrt(kfz*kfz+k_perp*k_perp);
  double thetak = atan2(k_perp,kfz);
  double k_plus = Ek+kfz;
//   cout << kfz << " " << k_perp << " " << thetak << endl;
  
  double rhou, rho_l_x, rho_l_z;
  getDensities(knorm,thetak,rhou,rho_l_z,rho_l_x);
  
  if(melosh) Melosh_rot(k_perp, k_plus, Ek, rho_l_z, rho_l_x);
  
  double x_nucl = Q2/2./(pi_mu*q_mu);
//   cout << pi_mu << " " << q_mu << endl;
//   cout << x_nucl << endl;
  NuclStructure nucl(proton,Q2,x_nucl,0,strucname);
  double F1,F2,g1;
  g1=nucl.getG1_grsv2000(proton,x,Q2);
  nucl.getF(F1,F2);
  double gamma_n = sqrt(Q2)*MASSn/(pi_mu*q_mu);
//   cout << gamma_n << endl;
  F1=(1+gamma_n*gamma_n)/(1.+0.18)/(2.*x_nucl)*F2;
//   cout << F1 << " " << F2 << " " << g1 << " " << gamma_n << endl;
  
  FourVector<double> snz(0.5*(pi_plus+(pi_perp*pi_perp-MASSn*MASSn)/pi_plus),pi_perp,0.,0.5*(pi_plus-(pi_perp*pi_perp-MASSn*MASSn)/pi_plus));
  snz*=1./MASSn;
  FourVector<double> snx(pi_perp/pi_plus,1.,0.,-pi_perp/pi_plus);
  
  double factor=MASSD*sqrt(1.+1./gamma_d/gamma_d)-pi_mu*el_mu;
  F_U = (epsilon*(-2.*F1+factor*factor*2.*F2/(pi_mu*q_mu))+(2.*F1+pt*pt/(pi_mu*q_mu)*F2))*rhou/3.*4*pow(2.*PI,3.);
  F_LSL = gamma_n*g1*((el_mu*snz)*rho_l_z+(el_mu*snx)*rho_l_x)*4*pow(2.*PI,3.);
//   cout << rhou << " " << rho_l_x << " " << rho_l_z << endl;
  
  cout << x << " " << Q2 << " " << x_nucl << " " << alpha_s << " " << pt << " " << psnorm << " " << knorm << " " << k_perp << " " << F_U << " " << F_LSL << " " << F_LSL/F_U << endl;
  
  
}
  
void Poldeut::Melosh_rot(double k_perp, double k_plus, double Ek, double &rho_l_z, double &rho_l_x){
  double nom=k_plus*(Ek+MASSn);
  
  double new_rho_z= (1.-k_perp*k_perp/nom)*rho_l_z +k_perp*(k_plus+MASSn)/nom*rho_l_x;
  double new_rho_x=(1.-k_perp*k_perp/nom)*rho_l_x -k_perp*(k_plus+MASSn)/nom*rho_l_z;
//   cout << new_rho_z << " " << rho_l_z << " " << new_rho_x << " " << rho_l_x << " " << (1.-k_perp*k_perp/nom) << " " << k_perp*(k_plus+MASSn)/nom << endl;
  rho_l_z=new_rho_z;
  rho_l_x=new_rho_x;
  
  return;
}


void Poldeut::getDensities(double knorm, double thetak, double &rhou, double &rho_l_z, double &rho_l_x){

  double U=getU(knorm);
  double W=getW(knorm);
  double Ek=sqrt(MASSn*MASSn+knorm*knorm);
  double alphai=Ek+knorm*cos(thetak);
  double factor=Ek/(alphai*alphai*4.*PI);
  rhou=3.*(U*U+W*W)  * factor;
  rho_l_z = 2.*(U*U+U*W/(2.*sqrt(2.))*(3.*cos(2.*thetak)+1.)+W*W/4.*(3.*cos(2.*thetak)-1.)) * factor;
  rho_l_x = 6./4.*(sqrt(2.)*U*W+W*W)*sin(2.*thetak) * factor;
//   cout << sin(2.*thetak) << endl;
  return;
}
