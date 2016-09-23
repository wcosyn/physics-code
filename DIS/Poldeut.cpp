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
  FourVector<double> q_mu(nu,0.,0.,-qvec);
  FourVector<double> el_mu(qvec,0.,0.,-nu);
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
  
  double rhou, rho_l_x, rho_l_z, rho_tensor_u;
  getDensities(knorm,thetak,rhou,rho_l_z,rho_l_x,rho_tensor_u);
  
  if(melosh) Melosh_rot(k_perp, k_plus, Ek, rho_l_z, rho_l_x);
  
  double x_nucl = Q2/2./(pi_mu*q_mu);
//   cout << pi_mu << " " << q_mu << endl;
//    cout << x_nucl << " " << pi_mu endl;
  NuclStructure nucl(proton,Q2,x_nucl,0,strucname);
  double F1,F2,g1,g1pg2;
  g1=NuclStructure::getG1_grsv2000(proton,x,Q2);
  g1pg2=NuclStructure::getG1plusg2_grsv2000(proton,x,Q2);
//   cout << g1 << " " << g1pg2 << " " << g1pg2-g1 << endl;
  nucl.getF(F1,F2);
  double gamma_n = sqrt(Q2)*MASSn/(pi_mu*q_mu);
//   cout << gamma_n << endl;
  F1=(1+gamma_n*gamma_n)/(1.+0.18)/(2.*x_nucl)*F2;
//   cout << F1 << " " << F2 << " " << g1 << " " << gamma_n << endl;
  
  FourVector<double> snz(0.5*(pi_plus+(pi_perp*pi_perp-MASSn*MASSn)/pi_plus),pi_perp,0.,0.5*(pi_plus-(pi_perp*pi_perp-MASSn*MASSn)/pi_plus));
  snz*=1./MASSn;
  FourVector<double> snx(pi_perp/pi_plus,1.,0.,-pi_perp/pi_plus);
  
  double factor=MASSD*sqrt(1.+1./gamma_d/gamma_d)-ps_mu*el_mu;
  double ff=(epsilon*(-2.*F1+factor*factor*2.*F2/(pi_mu*q_mu))+(2.*F1+pt*pt/(pi_mu*q_mu)*F2))/F2*Ek/alpha_i/alpha_i*4.*pow(2.*PI,3.);
  F_U = (epsilon*(-2.*F1+factor*factor*2.*F2/(pi_mu*q_mu))+(2.*F1+pt*pt/(pi_mu*q_mu)*F2))*rhou;
  F_LSL = gamma_n*(((el_mu*snz)*g1pg2-q_mu*snz/(pi_mu*q_mu)*factor*(g1pg2-g1))*rho_l_z+((el_mu*snx)*g1pg2-q_mu*snx/(pi_mu*q_mu)*factor*(g1pg2-g1))*rho_l_x);
  double F_UTLL = -(epsilon*(-2.*F1+factor*factor*2.*F2/(pi_mu*q_mu))+(2.*F1+pt*pt/(pi_mu*q_mu)*F2))*rho_tensor_u;
//   cout << (el_mu*snz)*rho_l_z << " " << (el_mu*snx)*rho_l_x << " " << snz << " " << snx << endl;
//   cout << rhou << " " << rho_l_x << " " << rho_l_z << endl;
//   cout << "dd " << rhou/3./Ek*alpha_i*alpha_i*pow((MASSn*MASSn-pi_mu*pi_mu)/wfref->getResidu(),2.) << " " << (MASSn*MASSn-pi_mu*pi_mu)/wfref->getResidu()<< " " << Ek << " " << alpha_i << " " << rhou << endl;
  cout << x << " " << Q2 << " " << x_nucl << " " << alpha_s << " " << pt << " " << psnorm << " " << knorm << " " << k_perp << " " << F_U << " " << F_LSL << " " << F_LSL/F_U << " " << (MASSn*MASSn-pi_mu*pi_mu)*1.E-06 << " " << F_U*pow((MASSn*MASSn-pi_mu*pi_mu)/wfref->getResidu(),2.)/ff << " " << F2 << " " << F_UTLL/F_U/sqrt(3./2.) << endl;
  
  
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

//densities include C factor defined in the paper
void Poldeut::getDensities(double knorm, double thetak, double &rhou, double &rho_l_z, double &rho_l_x, double &rho_tensor_u){

  
//   TVector3 kk(200.*sin(1.)*cos(2.),200.*sin(1.)*sin(2.),200.*cos(1.));
//   double U=getU(200.);
//   double W=getW(200.);
//   cout << wf.DeuteronPState(2,1,1,kk)*sqrt(4.*PI) << " " << U-(3.*cos(1.)*cos(1.)-1)/sqrt(8.)*W << endl;
//   cout << wf.DeuteronPState(2,-1,1,kk)*sqrt(4.*PI) << " " << wf.DeuteronPState(2,1,-1,kk)*sqrt(4.*PI) << " " << -3.*cos(1.)*sin(1.)*exp(I_UNIT*2.)/sqrt(8.)*W << endl;
//   cout << wf.DeuteronPState(2,-1,-1,kk)*sqrt(4.*PI) << " " << -3.*sin(1.)*sin(1.)*exp(I_UNIT*4.)/sqrt(8.)*W << endl;
//   cout << wf.DeuteronPState(0,1,1,kk)*sqrt(4.*PI) << " " << -3./2.*W*cos(1.)*sin(1.)*exp(-I_UNIT*2.) << endl;
//   cout << wf.DeuteronPState(0,1,-1,kk)*sqrt(4.*PI) << " "<<  wf.DeuteronPState(0,-1,1,kk)*sqrt(4.*PI) << " " << U/sqrt(2.)+W/2.*(3.*cos(1.)*cos(1.)-1) << endl;
//   cout << wf.DeuteronPState(0,-1,-1,kk)*sqrt(4.*PI) << " " << 3./2.*exp(I_UNIT*2.)*W*sin(1.)*cos(1.) << endl;
//   cout << wf.DeuteronPState(-2,1,1,kk)*sqrt(4.*PI) << " " << -3/sqrt(8)*W*sin(1.)*sin(1.)*exp(-I_UNIT*4.) << endl;
//   cout << wf.DeuteronPState(-2,1,-1,kk)*sqrt(4.*PI) << " " << wf.DeuteronPState(-2,-1,1,kk)*sqrt(4.*PI) << " " << 3/sqrt(8)*W*cos(1.)*sin(1.)*exp(-I_UNIT*2.) << endl;
//   cout << wf.DeuteronPState(-2,-1,-1,kk)*sqrt(4.*PI) << " " << U-W/sqrt(8)*(3.*cos(1.)*cos(1.)-1.) << endl;
// //   exit(1);
  
  double U=getU(knorm);
  double W=getW(knorm);
  
//   TVector3 k(knorm*sin(thetak),0,knorm*cos(thetak));
//   double dens0=0.,densplus=0.,densmin=0.;
//   for(int i1=-1;i1<=1;i1+=2){
//     for(int i2=-1;i2<=1;i2+=2){
//       dens0+=norm(wf.DeuteronPState(0,i1,i2,k));
//       densplus+=norm(wf.DeuteronPState(2,i1,i2,k));
//       densmin+=norm(wf.DeuteronPState(-2,i1,i2,k));
//     }
//   }
//   cout << dens0 << " " << densplus << " " << densmin << endl;
    

  
  double Ek=sqrt(MASSn*MASSn+knorm*knorm);
  double t= pow(MASSD-Ek,2.)-knorm*knorm;
  double alphai=(Ek+knorm*cos(thetak))/MASSn;
  double factor=Ek/(alphai*alphai*4.*PI)*4.*pow(2.*PI,3.);
  rhou=(U*U+W*W)  * factor;
  rho_l_z = 2.*(U*U-U*W/(2.*sqrt(2.))*(3.*cos(2.*thetak)+1.)+W*W/4.*(3.*cos(2.*thetak)-1.)) * factor;
  rho_l_x = 3./2.*(-sqrt(2.)*U*W+W*W)*sin(2.*thetak) * factor;
  rho_tensor_u = sqrt(3./2.)*(U*W/sqrt(2.)+W*W/4.)*(3.*cos(2.*thetak)+1)*factor;
//   cout << rhou/factor/4./PI << " " << dens0+densplus+densmin << endl;
//   cout << rho_tensor_u/factor/4./PI/sqrt(3./2.)*3. << " " << densplus+densmin-2.*dens0 << " " << rho_tensor_u/factor/4./PI/sqrt(3./2.)*3./(densplus+densmin-2.*dens0) << endl;
//   cout << alphai << " " << U*U << " " << U*W << " " << W*W << endl;
//   cout << 3.*rho_tensor_u/rhou/sqrt(3./2.) << endl;
//   exit(1);
  return;
}
