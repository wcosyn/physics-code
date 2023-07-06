#include "Poldeut.hpp"
#include <FourVector.h>
#include <NuclStructure.hpp>
#include <NucleonStructure.hpp>
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
  double pi_min=(MASSn*MASSn+pt*pt)/pi_plus;
  FourVector<double> pi_mu_onshell((pi_plus + pi_min)*0.5,-pt,0.,(pi_plus - pi_min)*0.5);
  
  double nu = 2.*Q2/MASSD/x;
  double qvec = sqrt(Q2+nu*nu);
  FourVector<double> q_mu(nu,0.,0.,-qvec);
  FourVector<double> q_onshell = q_mu+pi_mu-pi_mu_onshell;
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
  
  double x_nucl = Q2/2./(pi_mu_onshell*q_onshell);
//   cout << pi_mu << " " << q_mu << endl;
//    cout << x_nucl << " " << pi_mu endl;
  NuclStructure nucl(proton,Q2,x_nucl,0,strucname);
  double F1,F2,g1,g1pg2;
  g1=NuclStructure::getG1_grsv2000(proton,x,Q2);
  g1pg2=NuclStructure::getG1plusg2_grsv2000(proton,x,Q2);
//   cout << g1 << " " << g1pg2 << " " << g1pg2-g1 << endl;
  nucl.getF(F1,F2);
  double gamma_n = sqrt(Q2)*MASSn/(pi_mu_onshell*q_onshell);
//   cout << gamma_n << endl;
  F1=(1+gamma_n*gamma_n)/(1.+0.18)/(2.*x_nucl)*F2;
//   cout << F1 << " " << F2 << " " << g1 << " " << gamma_n << endl;
  
  FourVector<double> snz(0.5*(pi_plus+(pi_perp*pi_perp-MASSn*MASSn)/pi_plus),pi_perp,0.,0.5*(pi_plus-(pi_perp*pi_perp-MASSn*MASSn)/pi_plus));
  snz*=1./MASSn;
  FourVector<double> snx(pi_perp/pi_plus,1.,0.,-pi_perp/pi_plus);
  
  double factor=MASSD*sqrt(1.+1./gamma_d/gamma_d)-ps_mu*el_mu;
  double ff=(epsilon*(-2.*F1+factor*factor*2.*F2/(pi_mu_onshell*q_onshell))+(2.*F1+pt*pt/(pi_mu_onshell*q_onshell)*F2))/F2*Ek/alpha_i/alpha_i*4.*pow(2.*PI,3.);
  F_U = (epsilon*(-2.*F1+factor*factor*2.*F2/(pi_mu_onshell*q_onshell))+(2.*F1+pt*pt/(pi_mu_onshell*q_onshell)*F2))*rhou;
  F_LSL = gamma_n*(((el_mu*snz)*g1pg2-q_onshell*snz/(pi_mu_onshell*q_onshell)*factor*(g1pg2-g1))*rho_l_z+((el_mu*snx)*g1pg2-q_onshell*snx/(pi_mu_onshell*q_onshell)*factor*(g1pg2-g1))*rho_l_x);
  double F_UTLL = -(epsilon*(-2.*F1+factor*factor*2.*F2/(pi_mu_onshell*q_onshell))+(2.*F1+pt*pt/(pi_mu_onshell*q_onshell)*F2))*rho_tensor_u;
//   cout << (el_mu*snz)*rho_l_z << " " << (el_mu*snx)*rho_l_x << " " << snz << " " << snx << endl;
//   cout << rhou << " " << rho_l_x << " " << rho_l_z << endl;
//   cout << "dd " << rhou/3./Ek*alpha_i*alpha_i*pow((MASSn*MASSn-pi_mu*pi_mu)/wfref->getResidu(),2.) << " " << (MASSn*MASSn-pi_mu*pi_mu)/wfref->getResidu()<< " " << Ek << " " << alpha_i << " " << rhou << endl;
  cout << x << " " << Q2 << " " << x_nucl << " " << alpha_s << " " << pt << " " << psnorm << " " << knorm << " " << k_perp << " " << F_U << " " << F_LSL << " " << F_LSL/(F_U+1/sqrt(6.)*F_UTLL) << " " << (MASSn*MASSn-pi_mu*pi_mu)*1.E-06 << " " << F_U*pow((MASSn*MASSn-pi_mu*pi_mu)/wfref->getResidu(),2.)/ff << " " << F2 << " " << F_UTLL/F_U/sqrt(3./2.) << endl;
  
  
}


double Poldeut::calc_dsigma_unpol(const double x, const double Q2, const double s, const double alpha_p, const double tprime, bool proton){
  double pt2=1./4.*alpha_p*(2.-alpha_p)*MASSD*MASSD-MASSn*MASSn-alpha_p*tprime/2.;

  double k_perp=-sqrt(pt2);  //minus sign!!
  double Ek = sqrt((MASSn*MASSn+k_perp*k_perp)/alpha_p/(2.-alpha_p));
  double kfz=(alpha_p-1.)*Ek;
  double knorm = sqrt(kfz*kfz+k_perp*k_perp);
  double thetak = atan2(k_perp,kfz);
  double k_plus = Ek+kfz;

  
  double y=Q2/x*2./(s-MASSD*MASSD);
  double gamma_d = x*x*MASSD*MASSD/(Q2);
  double epsilon=(1.-y-pow(gamma_d*y/2.,2.))/(1.-y+y*y/2.+pow(gamma_d*y/2.,2.));
  double gamma_n = pow(2.*x/(2.-alpha_p)*MASSn,2.)/Q2;
  //cout << y << " " << gamma_d << " " << gamma_n << " " << epsilon << endl; 

  double U=getU(knorm);
  double W=getW(knorm);
  NucleonStructure nucl(strucname);
  double F1=0.,F2=0;
  nucl.getF_xQ(F1,F2,proton,x/(2.-alpha_p),Q2);

  double dens=Ek*(U*U+W*W)/(2.-alpha_p)/4./PI;
  
  
  double F1d = 2.*pow(2.*PI,3.)*2.*dens/(2.-alpha_p)*F1;
  double F2d = 2.*pow(2.*PI,3.)*dens*F2;

  double strucunpol2 = 2.*F1d+epsilon*((1+gamma_d*gamma_d)*F2d/x*2.-2.*F1d);

  double strucunpol = 2.*pow(2.*PI,3.)*2.*dens/(2.-alpha_p)*(2.*F1+epsilon*((1+gamma_n*gamma_n)/x*(2.-alpha_p)*F2-2.*F1));
  cout << sqrt(pt2) << " " << dens << " " << strucunpol << " " << strucunpol2 << " " << ALPHA*ALPHA/4./pow(2.*PI,2.)*y*y/Q2/Q2/(1.-epsilon) << " "; //Mev^-2, Mev^-2, Mev^-4
  double dsigma = strucunpol*ALPHA*ALPHA/4./pow(2.*PI,2.)*y*y/Q2/Q2/(1.-epsilon);  //MeV-6
  return dsigma*HBARC*HBARC*1.E19;  //nb/GeV^-4

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


void Poldeut::Tensor_Compare_nonrel(double &lfratio, double &nrratio, double alpha_p, double pt){
  double Ek_sq = (MASSn*MASSn+pt*pt)/alpha_p/(2.-alpha_p);
  double knorm = sqrt(Ek_sq-MASSn*MASSn);
  double U_lf = getU(knorm);
  double W_lf = getW(knorm);

  double p_p_plus = alpha_p/2.*MASSD;
  double p_p_minus = (MASSn*MASSn+pt*pt)/p_p_plus;
  double p_z = (p_p_plus - p_p_minus)/2.;
  double pnorm = sqrt(p_z*p_z+pt*pt);

  double U_nr = getU(pnorm);
  double W_nr = getW(pnorm);

  lfratio = -((1-3./2.*pt*pt/knorm/knorm)*(2*U_lf+W_lf/sqrt(2.))*W_lf/sqrt(2.))/(U_lf*U_lf + W_lf*W_lf) ;
  lfratio = U_lf;

  nrratio = -((1-3./2.*pt*pt/pnorm/pnorm/(1+p_z/MASSn))*(2*U_nr+W_nr/sqrt(2.))*W_nr/sqrt(2.))/(U_nr*U_nr + W_nr*W_nr) ;
  nrratio = W_lf;

}



void Poldeut::getLFdistributions(double alpha_p, double pt, double &S_L, double &S_T, double &deltaS_L, double &deltaS_T, double &deltaT_S_L, double &deltaT_S_T){

  double E = sqrt((pt*pt+MASSn*MASSn)/(alpha_p*(2.-alpha_p))); //MeV
  double kz = E*(alpha_p-1.);
  double knorm = sqrt(pt*pt+kz*kz);
  double costheta = kz/knorm;
  double sintheta = pt/knorm;
  if(knorm>1.E03) { S_L=S_T=1.;deltaS_L=deltaS_T=deltaT_S_L=deltaT_S_T=std::numeric_limits<double>::quiet_NaN() ; return;}

  double f0 = getU(knorm)*sqrt(E/(4.*PI));  //MeV-1
  double f2 = getW(knorm)*sqrt(E/(4.*PI)); //MeV-1

  double S_unpol = (f0*f0 + f2*f2)/(2-alpha_p);
  double S_tensor_L = -1./(2.-alpha_p)*(2.*f0+f2/sqrt(2.))*f2/sqrt(2.)*(1.-3./2.*sintheta*sintheta); 
  double S_tensor_T = -1./(2.-alpha_p)*(2.*f0+f2/sqrt(2.))*f2/sqrt(2.)*(1.-3./2.*costheta*costheta); 

  S_L = S_unpol+S_tensor_L;
  S_T = S_unpol+S_tensor_T;

  double C0L = 1-(E+kz)*pt*pt/(E+MASSn)/(MASSn*MASSn+pt*pt);
  double C2L = 1-(E+2.*MASSn)*(E+kz)*pt*pt/(MASSn*MASSn+pt*pt)/(knorm*knorm);

  deltaS_L = 1/(2.-alpha_p)*(f0-f2/sqrt(2.))*(C0L*f0-C2L*f2/sqrt(2.));

  double C0T = -(E+kz)*(E-kz+MASSn)*pt/(E+MASSn)/(MASSn*MASSn+pt*pt);
  double C2T = (E+kz)*(-knorm*knorm+kz*(E+2.*MASSn))*pt/(knorm*knorm*(MASSn*MASSn+pt*pt));

  deltaS_T = 1/(2.-alpha_p)*(f0-f2/sqrt(2.))*(C0T*f0-C2T*f2/sqrt(2.));

  double D0L = (E+kz)*(E-kz+MASSn)*pt/(E+MASSn)/(MASSn*MASSn+pt*pt); 
  double D2L = -(E+kz)*(-knorm*knorm+(E+MASSn/2.)*kz)*pt/(knorm*knorm*(MASSn*MASSn+pt*pt));
  
  deltaT_S_L = 1/(2.-alpha_p)*(f0-f2/sqrt(2.))*(D0L*f0+D2L*sqrt(2.)*f2);

  double D0T = 1.-(E+kz)*pt*pt/(E+MASSn)/(MASSn*MASSn+pt*pt);
  double D2T = 1.-(E+MASSn/2)*(E+kz)*pt*pt/(knorm*knorm*(MASSn*MASSn+pt*pt));
  
  deltaT_S_T = 1/(2.-alpha_p)*(f0-f2/sqrt(2.))*(D0T*f0+D2T*sqrt(2.)*f2);

  return;
}