#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std; 

#include "FastParticle.hpp"
#include <constants.hpp>
#include <Utilfunctions.hpp>

using namespace std;

FastParticle::FastParticle(const int type, const int inc, const double momentum,
			   const double ptheta, const double pphi, const double hard_scale, const double ggamma, const std::string dir)
:particletype(type),
incoming(inc), 
p(momentum),
theta(ptheta),
costheta(cos(ptheta)),
phi(pphi),
hardscale(hard_scale),
sigma_decay_p(0.),
sigma_decay_n(0.),
userset(0),
Gamma(ggamma){
  //cout << "Initializing FastParticle object: ";
  ex=sin(theta)*cos(phi);
  ey=sin(theta)*sin(phi);
  ez=cos(theta);
  switch(particletype){
    case(0):
      //cout << "proton" << endl;
      setGlauberParameters(p,sigmap,beta2p,epsilonp,sigman,beta2n,epsilonn);
      nkt_sq=0.35*0.35*9;
      lc=2*p*HBARC/1./1e06;
      mass = MASSP;     
      break;
    case(1):
      //cout <<"neutron" << endl;
      setGlauberParameters(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp);
      nkt_sq=0.35*0.35*9;
      lc=2*p*HBARC/1./1e06;
      mass = MASSN;
      break;
    case(2):
      //cout <<"pi+" << endl;
//       setPionGlauberData(p,sigmap,beta2p,epsilonp,sigman,beta2n,epsilonn,dir);
//interpolating static arrays is waaaay faster of course.
      sigmap = interpolate(sigmap_array,log10(p),0.01,151,250.);
      beta2p = interpolate(beta2p_array,log10(p),0.01,151,250.);
      epsilonp = interpolate(epsp_array,log10(p),0.01,151,250.);
      sigman = interpolate(sigman_array,log10(p),0.01,151,250.);
      beta2n = interpolate(beta2n_array,log10(p),0.01,151,250.);
      epsilonn = interpolate(epsn_array,log10(p),0.01,151,250.);
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/0.7/1e06;
      mass = MASSPI;
      break;
    case(3):
      //cout <<"pi-" << endl;
//       setPionGlauberData(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp,dir);
      sigmap = interpolate(sigmap_array,log10(p),0.01,151,250.);
      beta2p = interpolate(beta2p_array,log10(p),0.01,151,250.);
      epsilonp = interpolate(epsp_array,log10(p),0.01,151,250.);
      sigman = interpolate(sigman_array,log10(p),0.01,151,250.);
      beta2n = interpolate(beta2n_array,log10(p),0.01,151,250.);
      epsilonn = interpolate(epsn_array,log10(p),0.01,151,250.);
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/0.7/1e06;
      mass = MASSPI;
      break;
    case(4):
      //cout <<"rho0" << endl;
//       setPionGlauberData(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp,dir);
      sigman=sigmap = (interpolate(sigman_array,log10(p),0.01,151,250.)+interpolate(sigmap_array,log10(p),0.01,151,250.))*0.5;
      beta2n=beta2p = (interpolate(beta2n_array,log10(p),0.01,151,250.)+interpolate(beta2p_array,log10(p),0.01,151,250.))*0.5;
      epsilonn = epsilonp = (interpolate(epsn_array,log10(p),0.01,151,250.)+interpolate(epsp_array,log10(p),0.01,151,250.))*0.5;
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/0.7/1e06;
      mass = MASSRHO;
      sigma_decay_n=sigma_decay_p=sigman+sigmap;
      break;
    case(5):
      //cout <<"rho0 more CT" << endl;
//       setPionGlauberData(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp,dir);
      sigman=sigmap = (interpolate(sigman_array,log10(p),0.01,151,250.)+interpolate(sigmap_array,log10(p),0.01,151,250.))*0.5;
      beta2n=beta2p = (interpolate(beta2n_array,log10(p),0.01,151,250.)+interpolate(beta2p_array,log10(p),0.01,151,250.))*0.5;
      epsilonn = epsilonp = (interpolate(epsn_array,log10(p),0.01,151,250.)+interpolate(epsp_array,log10(p),0.01,151,250.))*0.5;
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/0.5/1e06;
      mass = MASSRHO;
      sigma_decay_n=sigma_decay_p=sigman+sigmap;
      break;
    case(6):
      //cout <<"rho0 less CT" << endl;
//       setPionGlauberData(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp,dir);
      sigman=sigmap = (interpolate(sigman_array,log10(p),0.01,151,250.)+interpolate(sigmap_array,log10(p),0.01,151,250.))*0.5;
      beta2n=beta2p = (interpolate(beta2n_array,log10(p),0.01,151,250.)+interpolate(beta2p_array,log10(p),0.01,151,250.))*0.5;
      epsilonn = epsilonp = (interpolate(epsn_array,log10(p),0.01,151,250.)+interpolate(epsp_array,log10(p),0.01,151,250.))*0.5;
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/1./1e06;
      mass = MASSRHO;
      sigma_decay_n=sigma_decay_p=sigman+sigmap;
      break;
    case(7):
      //cout <<"double pion" << endl;
//       setPionGlauberData(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp,dir);
      sigman=sigmap = (interpolate(sigman_array,log10(p),0.01,151,250.)+interpolate(sigmap_array,log10(p),0.01,151,250.));
      beta2n=beta2p = (interpolate(beta2n_array,log10(p),0.01,151,250.)+interpolate(beta2p_array,log10(p),0.01,151,250.))*0.5;
      epsilonn = epsilonp = (interpolate(epsn_array,log10(p),0.01,151,250.)+interpolate(epsp_array,log10(p),0.01,151,250.))*0.5;
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/1./1e06;
      mass = MASSPI*2.;
//       sigma_decay_n=sigma_decay_p=sigman+sigmap;
      break;
    default:
      cerr << "Particle type is not yet supported!!: " << particletype << endl;
      exit(1);
  }
  E = sqrt(mass*mass+p*p);
  decay_dil = Gamma*sqrt(1.-p*p/(E*E));
  //if(hard_scale<nkt_sq) cout << "Warning: hard scale is too low for color transparency effects for this particle, just so you know!" << endl;
  
}



FastParticle::FastParticle(const int type, const int inc, const TVector3 &pvec, 
			   const double hard_scale, const double Gamma, const std::string dir)
:particletype(type),
incoming(inc), 
p(pvec.Mag()),
theta(pvec.Theta()),
costheta(pvec.CosTheta()),
phi(pvec.Phi()),
hardscale(hard_scale),
sigma_decay_p(0.),
sigma_decay_n(0.),
userset(0){
  //cout << "Initializing FastParticle object: ";
  ex=sin(theta)*cos(phi);
  ey=sin(theta)*sin(phi);
  ez=cos(theta);
  switch(particletype){
    case(0):
      //cout << "proton" << endl;
      setGlauberParameters(p,sigmap,beta2p,epsilonp,sigman,beta2n,epsilonn);
      nkt_sq=0.35*0.35*9;
      lc=2*p*HBARC/1./1e06;
      mass = MASSP;     
      break;
    case(1):
      //cout <<"neutron" << endl;
      setGlauberParameters(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp);
      nkt_sq=0.35*0.35*9;
      lc=2*p*HBARC/1./1e06;
      mass = MASSN;
      break;
    case(2):
      //cout <<"pi+" << endl;
//       setPionGlauberData(p,sigmap,beta2p,epsilonp,sigman,beta2n,epsilonn,dir);
      sigmap = interpolate(sigmap_array,log10(p),0.01,151,250.);
      beta2p = interpolate(beta2p_array,log10(p),0.01,151,250.);
      epsilonp = interpolate(epsp_array,log10(p),0.01,151,250.);
      sigman = interpolate(sigman_array,log10(p),0.01,151,250.);
      beta2n = interpolate(beta2n_array,log10(p),0.01,151,250.);
      epsilonn = interpolate(epsn_array,log10(p),0.01,151,250.);
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/0.7/1e06;
      mass = MASSPI;
      break;
    case(3):
      //cout <<"pi-" << endl;
//       setPionGlauberData(p,sigmap,beta2p,epsilonp,sigman,beta2n,epsilonn,dir);
      sigmap = interpolate(sigmap_array,log10(p),0.01,151,250.);
      beta2p = interpolate(beta2p_array,log10(p),0.01,151,250.);
      epsilonp = interpolate(epsp_array,log10(p),0.01,151,250.);
      sigman = interpolate(sigman_array,log10(p),0.01,151,250.);
      beta2n = interpolate(beta2n_array,log10(p),0.01,151,250.);
      epsilonn = interpolate(epsn_array,log10(p),0.01,151,250.);
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/0.7/1e06;
      mass = MASSPI;
      break;
    case(4):
      //cout <<"rho0" << endl;
//       setPionGlauberData(p,sigmap,beta2p,epsilonp,sigman,beta2n,epsilonn,dir);
      sigman=sigmap = (interpolate(sigman_array,log10(p),0.01,151,250.)+interpolate(sigmap_array,log10(p),0.01,151,250.))*0.5;
      beta2n=beta2p = (interpolate(beta2n_array,log10(p),0.01,151,250.)+interpolate(beta2p_array,log10(p),0.01,151,250.))*0.5;
      epsilonn = epsilonp = (interpolate(epsn_array,log10(p),0.01,151,250.)+interpolate(epsp_array,log10(p),0.01,151,250.))*0.5;
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/0.7/1e06;
      mass = MASSRHO;
      sigma_decay_n=sigma_decay_p=sigman+sigmap;
      break;
    case(5):
      //cout <<"rho0 more CT" << endl;
//       setPionGlauberData(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp,dir);
      sigman=sigmap = (interpolate(sigman_array,log10(p),0.01,151,250.)+interpolate(sigmap_array,log10(p),0.01,151,250.))*0.5;
      beta2n=beta2p = (interpolate(beta2n_array,log10(p),0.01,151,250.)+interpolate(beta2p_array,log10(p),0.01,151,250.))*0.5;
      epsilonn = epsilonp = (interpolate(epsn_array,log10(p),0.01,151,250.)+interpolate(epsp_array,log10(p),0.01,151,250.))*0.5;
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/0.5/1e06;
      mass = MASSRHO;
      sigma_decay_n=sigma_decay_p=sigman+sigmap;
      break;
    case(6):
      //cout <<"rho0 less CT" << endl;
//       setPionGlauberData(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp,dir);
      sigman=sigmap = (interpolate(sigman_array,log10(p),0.01,151,250.)+interpolate(sigmap_array,log10(p),0.01,151,250.))*0.5;
      beta2n=beta2p = (interpolate(beta2n_array,log10(p),0.01,151,250.)+interpolate(beta2p_array,log10(p),0.01,151,250.))*0.5;
      epsilonn = epsilonp = (interpolate(epsn_array,log10(p),0.01,151,250.)+interpolate(epsp_array,log10(p),0.01,151,250.))*0.5;
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/1./1e06;
      mass = MASSRHO;
      sigma_decay_n=sigma_decay_p=sigman+sigmap;
      break;
    case(7):
      //cout <<"double pion" << endl;
//       setPionGlauberData(p,sigman,beta2n,epsilonn,sigmap,beta2p,epsilonp,dir);
      sigman=sigmap = (interpolate(sigman_array,log10(p),0.01,151,250.)+interpolate(sigmap_array,log10(p),0.01,151,250.));
      beta2n=beta2p = (interpolate(beta2n_array,log10(p),0.01,151,250.)+interpolate(beta2p_array,log10(p),0.01,151,250.))*0.5;
      epsilonn = epsilonp = (interpolate(epsn_array,log10(p),0.01,151,250.)+interpolate(epsp_array,log10(p),0.01,151,250.))*0.5;
      nkt_sq=0.35*0.35*4;
      lc=2*p*HBARC/1./1e06;
      mass = MASSPI*2.;
//       sigma_decay_n=sigma_decay_p=sigman+sigmap;
      break;
    default:
      cerr << "Particle type is not yet supported!!: " << particletype << endl;
      exit(1);
  }
  E = sqrt(mass*mass+p*p);
  decay_dil = Gamma/sqrt(1.-p*p/E*E);
  //if(hard_scale<nkt_sq) cout << "Warning: hard scale is too low for color transparency effects for this particle, just so you know!" << endl;
  
}

FastParticle::FastParticle(const FastParticle &Copy):
particletype(Copy.getParticletype()),
incoming(Copy.getIncoming()),
p(Copy.getP()),
theta(Copy.getTheta()),
costheta(Copy.getCosTheta()),
phi(Copy.getPhi()),
ex(Copy.getEx()),
ey(Copy.getEy()),
ez(Copy.getEz()),
hitz(Copy.getHitz()),
hitbnorm(Copy.getHitbnorm()),
beta2p(Copy.getBeta2p()),
sigmap(Copy.getSigmap()),
epsilonp(Copy.getEpsilonp()),
beta2n(Copy.getBeta2n()),
sigman(Copy.getSigman()),
epsilonn(Copy.getEpsilonn()),
hardscale(Copy.getHardScale()),
nkt_sq(Copy.getNkt_sq()),
lc(Copy.getLc()),
mass(Copy.getMass()),
E(Copy.getE()),
decay_dil(Copy.getDecay_dil()),
sigma_decay_n(Copy.getSigma_decay_n()),
sigma_decay_p(Copy.getSigma_decay_p()),
userset(Copy.userset)
{
  for(int i=0;i<3;i++) hitb[i]=Copy.getHitbvec()[i];
}


FastParticle& FastParticle::operator=(const FastParticle& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
    particletype = rhs.getParticletype();
    incoming = rhs.getIncoming();
    p = rhs.getP();
    theta = rhs.getTheta();
    costheta = rhs.getCosTheta();
    phi = rhs.getPhi();
    ex = rhs.getEx();
    ey = rhs.getEy();
    ez = rhs.getEz();
    hitz = rhs.getHitz();
    hitbnorm = rhs.getHitbnorm();
    beta2p = rhs.getBeta2p();
    sigmap = rhs.getSigmap();
    epsilonp = rhs.getEpsilonp();
    beta2n = rhs.getBeta2n();
    sigman = rhs.getSigman();
    epsilonn = rhs.getEpsilonn();
    hardscale = rhs.getHardScale();
    nkt_sq = rhs.getNkt_sq();
    lc = rhs.getLc();
    mass = rhs.getMass();
    E = rhs.getE();
    decay_dil = rhs.getDecay_dil();
    sigma_decay_n = rhs.getSigma_decay_n();
    sigma_decay_p = rhs.getSigma_decay_p();
    userset = rhs.userset;
    for(int i=0;i<3;i++) hitb[i]=rhs.getHitbvec()[i];
     
     
  }

  return *this;
}



FastParticle::~FastParticle(){
  //cout << "Deleting Particle object..." << endl;
}

//set coord of hard interaction point 
void FastParticle::setHitcoord(double r, double costheta, double sintheta, double cosphi, double sinphi){
  hitz=r*costheta*getEz()+r*sintheta*cosphi*getEx()+r*sintheta*sinphi*getEy();
  hitb[0]=r*sintheta*cosphi-hitz*getEx();
  hitb[1]=r*sintheta*sinphi-hitz*getEy();
  hitb[2]=r*costheta-hitz*getEz();
  hitbnorm=sqrt(hitb[0]*hitb[0]+hitb[1]*hitb[1]+hitb[2]*hitb[2]);
}
  

//return z coord of a point along momentem
double  FastParticle::calcZ(double r, double costheta, double sintheta, double cosphi, double sinphi){
  return r*costheta*getEz()+r*sintheta*cosphi*getEx()+r*sintheta*sinphi*getEy();
}  

//get sigma for a particle, make distinction between scatt w proton or neutron
double FastParticle::getSigma_decay(int level, MeanFieldNucleus *pnucleus) const{
  return level<pnucleus->getPLevels()? getSigma_decay_p() : getSigma_decay_n();
}

//get sigma for a particle, make distinction between scatt w proton or neutron
double FastParticle::getSigma(int level, MeanFieldNucleus *pnucleus) const{
  return level<pnucleus->getPLevels()? getSigmap() : getSigman();
}

//get epsilon for a particle, make distinction between scatt w proton or neutron
double FastParticle::getEpsilon(int level, MeanFieldNucleus *pnucleus) const{
  return level<pnucleus->getPLevels()? getEpsilonp() : getEpsilonn();
}

//get betasq for a particle, make distinction between scatt w proton or neutron
double FastParticle::getBetasq(int level, MeanFieldNucleus *pnucleus) const{
  return level<pnucleus->getPLevels()? getBeta2p() : getBeta2n();
}


complex<double> FastParticle::getScatterfront(int level, MeanFieldNucleus *pnucleus) const{
  return getSigma(level,pnucleus)
			  *(1.-I_UNIT*getEpsilon(level,pnucleus))
			  /(4.*PI*getBetasq(level,pnucleus));
}

complex<double> FastParticle::getScatterfront(bool proton) const{
  return getSigma(proton)
			  *(1.-I_UNIT*getEpsilon(proton))
			  /(4.*PI*getBetasq(proton));
}

double FastParticle::getCTsigma(double z) const{
 return ((abs(z-getHitz()) > getLc())
			  ||(getHardScale() < getNkt_sq()))? 
			  1.:(abs(z-getHitz())/getLc() 
			  + getNkt_sq()/getHardScale()*
			  (1.-abs(z-getHitz()) / getLc())); 
  
}

double FastParticle::getCT_decay_sigma(double z, int level, MeanFieldNucleus *pnucleus) const{
  if((abs(z-getHitz()) > getLc())||(getHardScale() < getNkt_sq()))
    return getDecay_sigma(z,level,pnucleus);
  else 
    return abs(z-getHitz())/getLc() 
			  + getNkt_sq()/getHardScale()*(1.-abs(z-getHitz()) / getLc());

}
    
double FastParticle::getCT_decay_sigma(double z, bool proton) const{
  if((abs(z-getHitz()) > getLc())||(getHardScale() < getNkt_sq()))
    return getDecay_sigma(z,proton);
  else 
    return abs(z-getHitz())/getLc() 
			  + getNkt_sq()/getHardScale()*(1.-abs(z-getHitz()) / getLc());

}
double FastParticle::getDecay_sigma(double z, int level, MeanFieldNucleus *pnucleus) const{
    double a=exp(-getDecay_dil()*(z-getHitz())*INVHBARC);
    return a + getSigma_decay(level,pnucleus)/getSigma(level,pnucleus)*(1.-a);
  
}
    
double FastParticle::getDecay_sigma(double z, bool proton) const{
    double a=exp(-getDecay_dil()*(z-getHitz())*INVHBARC);
    return a + getSigma_decay(proton)/getSigma(proton)*(1.-a);
  
}
    



void FastParticle::printParticle() const{

  cout << "****************************************" << endl << "Type: " << getParticletype() << " = ";
 switch(getParticletype()){
   case(0):
     cout << "proton" << endl;
     break;
   case(1):
     cout << "neutron" << endl;
     break;
   case(2):
     cout << "pi+" << endl;
     break;
   case(3):
     cout << "pi-" << endl;
     break;
   case(4):
     cout << "rho0" << endl;
     break;
   default:
     cerr << "Particletypeerror!!!" << endl;
     exit(1);
 }
 if(incoming) cout << "This is the incoming particle" << endl;
 cout << "Momentum, norm[MeV]=" << getP() << ", theta[deg]=" << getTheta()*RADTODEGR << ", phi[deg]=" << getPhi()*RADTODEGR << endl;
 cout << "Scattering Parameters with proton:" << endl << "sigma[fm^2]: " << getSigmap() << 
  ", beta^2[fm^2]: " << getBeta2p() << ", epsilon: " << getEpsilonp() << endl; 
 cout << "Scattering Parameters with neutron:" << endl << "sigma[fm^2]: " << getSigman() << 
  ", beta^2[fm^2]: " << getBeta2n() << ", epsilon: " << getEpsilonn() << endl;
 cout << "Scattering hard scale [GeV^2]: " << hardscale << endl;
 cout << "****************************************" << endl << endl;
}




/////////////////////////////////////////////
//calculate parameters for profile function//
/////////////////////////////////////////////
void FastParticle::setGlauberParameters(double mom, double &sigmapp, double &beta2pp, double &epspp, 
					double &sigmapn, double &beta2pn, double &epspn){

  double p_gev = mom/1000.;
  double sigmappel, sigmapnel;

  // sigmapp

  if (p_gev > 1.4) sigmapp = 48 + 0.522*pow(log(p_gev),2) - 4.51*log(p_gev);
  else {
    if (p_gev > 0.9) sigmapp = 23 + 47.1 * (p_gev - 0.9);
    else {
      if (p_gev > 0.5 ) sigmapp = 25 - 5 * (p_gev - 0.5);
      else  sigmapp = 8.667/p_gev/p_gev - 23.114/p_gev + 36.559;
    }
  }	
  
  sigmapp /= 10;

  // sigmapp elastic

  if (p_gev > 1.7) 
    {
      sigmappel = 11.9 + 26.9 * pow(p_gev,-1.21) + 0.169 * pow(log(p_gev),2) - 
	1.85 * log(p_gev);
    }
  else{
    if (p_gev > 0.9) sigmappel = 23 + 2.65 * (p_gev - 0.9);
    else {
      if (p_gev > 0.5) sigmappel = 25 - 5 * (p_gev - 0.5);
      else  sigmappel = 8.667/p_gev/p_gev - 23.114/p_gev + 36.559;
    }
  }
  
  sigmappel /= 10;

  // sigmapn

  if (p_gev > 2) sigmapn = 40;
  else{
    if (p_gev > 1) sigmapn = 33.7 + 6.3 * (p_gev - 1);
    else{
      if (p_gev > 0.4) sigmapn = 19.7/p_gev/p_gev - 34.083/p_gev + 48.083;
      else  sigmapn = 26.426/p_gev/p_gev - 48.595/p_gev + 42.325;
    }
  }
  sigmapn /= 10;

  // sigmapnel
 
  if (p_gev > 1) sigmapnel = 31 - 10.85 * log(p_gev);
  else{
    if (p_gev > 0.4) sigmapnel = 19.7/p_gev/p_gev - 34.083/p_gev + 48.083;
    else  sigmapnel = 26.426/p_gev/p_gev - 48.595/p_gev + 42.325;
  }
  sigmapnel /= 10;

  // epspp

  if (p_gev > 0.32) epspp = -0.18/p_gev/p_gev + 1.45/p_gev - 0.77;
  else epspp = 0.8 + 6.49 * (p_gev - 0.135);
  
  // epspn

  if (p_gev > 3) epspn = 6.61/p_gev/p_gev - 3.51/p_gev - 0.0655;
  else{
    if (p_gev > 0.50) epspn = 0.27/p_gev/p_gev + 0.39/p_gev - 0.66;
    else  epspn = 0.25 + 1.04 * (log(p_gev) + 1.61);
  }
  // betapp

  beta2pp = (pow(sigmapp,2)*(1+pow(epspp,2)))/(16*PI*sigmappel);
  beta2pn = (pow(sigmapn,2)*(1+pow(epspn,2)))/(16*PI*sigmapnel);
}


void FastParticle::setPionGlauberData(double mom, double &sigmatotp, double &beta2p, double &epsp, 
		     double &sigmatotn, double &beta2n, double &epsn, const string dir){
  double logp=log10(mom/1000.);
  double sigmaelp, sigmaeln;
  double betap2, betan2;
  //read in all the data
  const string files[6]={dir+"/input/pimp_elastic.dat",dir+"/input/pimp_elastic.reim",dir+"/input/pimp_total.dat",
			       dir+"/input/pipp_elastic.dat",dir+"/input/pipp_elastic.reim",dir+"/input/pipp_total.dat"};
  int *lines=new int[6];
  double ***data = new double**[6];  //data is 3 dim array.  dim1 = 6 files, dim2 = 3 (momentum, y-value, sigma), dim3 = grid
  double dummy,dummy2; //used to read in redundant info etc
  for(int i=0;i<6;i++){
    data[i] = new double*[3];
    ifstream file(files[i].c_str(),ios::in);
    
    if(file.is_open()){

      //ignore first 10 lines, no usefull information
      for(int j=0;j<10;j++) file.ignore(500,'\n');

      //read in number of data points
      file >> lines[i];
      for(int j=0;j<3;j++) data[i][j]=new double[lines[i]];
      
      //read in data points
      for(int k=0;k<lines[i];k++){
	file >> dummy;
	file >> dummy;
	data[i][0][k]=log10(dummy); //momentum
	file >> dummy; file >> dummy;
	file >> dummy2;
	if(i==1||i==4)	data[i][1][k]=dummy2; //value(total,elastic or epsilon)
	else data[i][1][k]=log10(dummy2);
	file >> dummy;
	data[i][2][k]=1.;//dummy/log(10); //sigma for value
	file.ignore(500,'\n');
      }
      file.close();
    }
    
    else{
      cerr << "File " << files[i] << " could not be opened" << endl;
      exit(1);
    }
  }
  

  //new data for epsilon! extra input
  const string files2[2]={dir+"/input/epspim.inp",dir+"/input/epspip.inp"};
  
  double ***epsdata = new double**[2];
  int lines2[2]={164,134};
  for(int i=0;i<2;i++){
    epsdata[i] = new double*[3];
    
    ifstream file(files2[i].c_str(),ios::in);    
    if(file.is_open()){
      for(int j=0;j<3;j++) epsdata[i][j]=new double[lines2[i]];
      for(int k=0;k<lines2[i];k++){      
	file >> dummy;
	epsdata[i][0][k] = log10(dummy/1000.);
	file >> epsdata[i][1][k];
	epsdata[i][2][k] = 1.;
      }
      file.close();
    }
    
    else{
      cerr << "File " << files2[i] << " could not be opened" << endl;
      exit(1);
    }
  }


  //new data for beta, extra input
  const string files3[2]={dir+"/input/betapim.inp",dir+"/input/betapip.inp"};

  double ***betadata = new double**[2];
  int lines3[2]={130,96};
  for(int i=0;i<2;i++){
    betadata[i] = new double*[3];
    
    ifstream file(files3[i].c_str(),ios::in);    
    if(file.is_open()){
      for(int j=0;j<3;j++) betadata[i][j]=new double[lines3[i]];
      for(int k=0;k<lines3[i];k++){      
	file >> dummy;
	betadata[i][0][k] = log10(dummy/1000.);
	file >> betadata[i][1][k];
	betadata[i][2][k] = 1.;
	file.ignore(500,'\n');
      }
      file.close();
    }
    
    else{
      cerr << "File " << files3[i] << " could not be opened" << endl;
      exit(1);
    }
  }

  double *a;
  double chi2;

  //pi- scattering
    
  //sigmatot in fm^2
  if(logp<=-0.19){
    a=fitdata(data[2][0],data[2][1],data[2][2],0,195,9,chi2);
    sigmatotn=0.;
    for(int i=0;i<9;i++) sigmatotn+=a[i]*pow(logp,i);
    sigmatotn=pow(10.,sigmatotn)/10;
    delete [] a;
  }
  else if(logp<=0.03){
    a=fitdata(data[2][0],data[2][1],data[2][2],183,330,8,chi2);
    sigmatotn=0.;
    for(int i=0;i<8;i++) sigmatotn+=a[i]*pow(logp,i);
    sigmatotn=pow(10.,sigmatotn)/10;
      delete [] a;
    }
  else if(logp<=0.47){
    a=fitdata(data[2][0],data[2][1],data[2][2],315,450,10,chi2);
    sigmatotn=0.;
    for(int i=0;i<10;i++) sigmatotn+=a[i]*pow(logp,i);
    sigmatotn=pow(10.,sigmatotn)/10;
    delete [] a;
  }
  else {
    a=fitdata(data[2][0],data[2][1],data[2][2],440,lines[2]-1,5,chi2);
    sigmatotn=0.;
    for(int i=0;i<5;i++) sigmatotn+=a[i]*pow(logp,i);
    sigmatotn=pow(10.,sigmatotn)/10;
    delete [] a;
    }

   //sigmael
  if(logp<=-0.19){
    a=fitdata(data[0][0],data[0][1],data[0][2],0,55,9,chi2);
    sigmaeln=0.;
    for(int i=0;i<9;i++) sigmaeln+=a[i]*pow(logp,i);
    sigmaeln=pow(10.,sigmaeln)/10;
    delete [] a;
  }
  else if(logp<=0.1){
    a=fitdata(data[0][0],data[0][1],data[0][2],45,120,10,chi2);
    sigmaeln=0.;
    for(int i=0;i<10;i++) sigmaeln+=a[i]*pow(logp,i);
    sigmaeln=pow(10.,sigmaeln)/10;
    delete [] a;
  }
  else {
    a=fitdata(data[0][0],data[0][1],data[0][2],110,lines[0]-1,4,chi2);
    sigmaeln=0.;
    for(int i=0;i<4;i++) sigmaeln+=a[i]*pow(logp,i);
    sigmaeln=pow(10.,sigmaeln)/10;
    delete [] a;
  }
    


  //eps
  /*a=fitdata(data[1][0],data[1][1],data[1][2],0,lines[1]-1,4,chi2);
  epsn=0.;
  for(int i=0;i<4;i++) epsn+=a[i]*pow(logp,i);
  delete [] a;*/

  //eps2, fit to new data
  if(logp<=0.04){
    double *tempdata=new double[41];
    for(int i=0;i<41;i++) tempdata[i]=pow(10.,epsdata[0][0][i])*1000.;
    a=fitdata(tempdata,epsdata[0][1],epsdata[0][2],0,40,20,chi2);
    epsn=0.;
    for(int i=0;i<20;i++) epsn+=a[i]*pow(mom,i);
    delete [] a;
    delete [] tempdata;
  }

  else if(logp<=0.57){
    a=fitdata(epsdata[0][0],epsdata[0][1],epsdata[0][2],35,130,12,chi2);
    epsn=0.;
    for(int i=0;i<12;i++) epsn+=a[i]*pow(logp,i);
    delete [] a;
  }
  else{
    a=fitdata(epsdata[0][0],epsdata[0][1],epsdata[0][2],120,lines2[0]-1,4,chi2);
    epsn=0.;
    for(int i=0;i<4;i++) epsn+=a[i]*pow(logp,i);
    delete [] a;
  }

  //beta, 2 possibilites:
  // 1. perform fit
  // 2. calc from other parameters
  if(logp<=0.09){
    a=fitdata(betadata[0][0],betadata[0][1],betadata[0][2],0,47,11,chi2);
    beta2n=0.;
    for(int i=0;i<11;i++) beta2n+=a[i]*pow(logp,i);
    beta2n = beta2n*pow(HBARC/1000.,2.);
    delete [] a;
  }
  else if(logp<=0.6){
    a=fitdata(betadata[0][0],betadata[0][1],betadata[0][2],40,105,11,chi2);
    beta2n=0.;
    for(int i=0;i<11;i++) beta2n+=a[i]*pow(logp,i);
    beta2n = beta2n*pow(HBARC/1000.,2.);
    delete [] a;
  }
  else{
    a=fitdata(betadata[0][0],betadata[0][1],betadata[0][2],102,lines3[0]-1,5,chi2);
    beta2n=0.;
    for(int i=0;i<5;i++) beta2n+=a[i]*pow(logp,i);
    beta2n = beta2n*pow(HBARC/1000.,2.);
    delete [] a;
  }    

  betan2=sigmatotn*sigmatotn*(1+epsn*epsn)/(16*PI*sigmaeln);

  //epsn=sqrt((16*PI*sigmaeln*betan*betan)/sigmatotn/sigmatotn-1);
  
  //pi+ scattering
  //sigmatot in fm^2
  if(logp<=-0.22){
    a=fitdata(data[5][0],data[5][1],data[5][2],0,135,10,chi2);
    sigmatotp=0.;
    for(int i=0;i<10;i++) sigmatotp+=a[i]*pow(logp,i);
    sigmatotp=pow(10.,sigmatotp)/10;
    delete [] a;
  }
  else if(logp<=0.32){
    a=fitdata(data[5][0],data[5][1],data[5][2],120,283,11,chi2);
    sigmatotp=0.;
    for(int i=0;i<11;i++) sigmatotp+=a[i]*pow(logp,i);
    sigmatotp=pow(10.,sigmatotp)/10;
    delete [] a;
  }
  else {
    a=fitdata(data[5][0],data[5][1],data[5][2],277,lines[5]-1,10,chi2);
    sigmatotp=0.;
    for(int i=0;i<10;i++) sigmatotp+=a[i]*pow(logp,i);
    sigmatotp=pow(10.,sigmatotp)/10;
    delete [] a;
  }

  //sigmael
  if(logp<=-0.18){
    /*a=fitdata(data[3][0],data[3][1],data[3][2],0,40,8,chi2);
      sigmael=0.;
      for(int i=0;i<8;i++) sigmael+=a[i]*pow(logp,i);
      sigmael=pow(10.,sigmael);
      delete [] a;*/
    sigmaelp=sigmatotp; //data coincide, so no need to do second fit
    }
  else if(logp<=0.26){
    a=fitdata(data[3][0],data[3][1],data[3][2],30,85,10,chi2);
    sigmaelp=0.;
    for(int i=0;i<10;i++) sigmaelp+=a[i]*pow(logp,i);
    sigmaelp=pow(10.,sigmaelp)/10;
    delete [] a;
  }
  else {
    a=fitdata(data[3][0],data[3][1],data[3][2],77,lines[3]-1,7,chi2);
    sigmaelp=0.;
    for(int i=0;i<7;i++) sigmaelp+=a[i]*pow(logp,i);
    sigmaelp=pow(10.,sigmaelp)/10;
    delete [] a;
  }

  
  //eps
  /*a=fitdata(data[4][0],data[4][1],data[4][2],0,lines[4]-1,2,chi2);
  epsp=0.;
  for(int i=0;i<2;i++) epsp+=a[i]*pow(logp,i);
  delete [] a;*/

  //eps2, fit to new data
  if(logp<=0.15){
    a=fitdata(epsdata[1][0],epsdata[1][1],epsdata[1][2],0,57,10,chi2);
    epsp=0.;
    for(int i=0;i<10;i++) epsp+=a[i]*pow(logp,i);
    delete [] a;
  }

  else if(logp<=0.50){
    a=fitdata(epsdata[1][0],epsdata[1][1],epsdata[1][2],50,123,10,chi2);
    epsp=0.;
    for(int i=0;i<10;i++) epsp+=a[i]*pow(logp,i);
    delete [] a;
  }
  else{
    a=fitdata(epsdata[1][0],epsdata[1][1],epsdata[1][2],105,lines2[1]-1,3,chi2);
    epsp=0.;
    for(int i=0;i<3;i++) epsp+=a[i]*pow(logp,i);
    delete [] a;
  }

  
  //beta
  //again 2 possibilities
  if(logp<=0.255){
    a=fitdata(betadata[1][0],betadata[1][1],betadata[1][2],0,52,8,chi2);
    beta2p=0.;
    for(int i=0;i<8;i++) beta2p+=a[i]*pow(logp,i);
    beta2p=beta2p*pow(HBARC/1000.,2.);
    delete [] a;
  }
  else{
    a=fitdata(betadata[1][0],betadata[1][1],betadata[1][2],46,lines3[1]-1,8,chi2);
    beta2p=0.;
    for(int i=0;i<8;i++) beta2p+=a[i]*pow(logp,i);
    beta2p=beta2p*pow(HBARC/1000.,2.);
    delete [] a;
  }

  betap2=sigmatotp*sigmatotp*(1+epsp*epsp)/(16*PI*sigmaelp);
  //epsp=sqrt((16*PI*sigmaelp*betap*betap)/sigmatotp/sigmatotp-1);
      
  //maintenance
  delete []lines; 
  for(int i=0;i<6;i++){
    for(int j=0;j<3;j++){
      delete [] data[i][j];
    }
  }
  for(int i=0;i<6;i++) delete [] data[i];
  delete [] data;
  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
      delete [] epsdata[i][j];
      delete [] betadata[i][j];
    }
  }
  for(int i=0;i<2;i++) {
    delete [] epsdata[i];
    delete [] betadata[i];
  }
  delete [] epsdata; delete [] betadata;

  
}

//used to calc (b-b')^2 for the glauberphase in the ref frame of this particle
double FastParticle::getBdist(double r, double costheta, double sintheta, double cosphi, double sinphi, double zmom){
  return pow(getHitbvec()[0]-r*sintheta*cosphi+zmom*getEx(),2.)
	 +pow(getHitbvec()[1]-r*sintheta*sinphi+zmom*getEy(),2.)
	 +pow(getHitbvec()[2]-r*costheta+zmom*getEz(),2.);
}

void FastParticle::interpPionGlauberData(int particletype, double mom, double &sigmap, double &beta2p, double &epsilonp, 
		     double &sigman, double &beta2n, double &epsilonn){
  switch(particletype){
  case(2):
    //cout <<"pi+" << endl;
    sigmap = interpolate(sigmap_array,log10(mom),0.01,151,250.);
    beta2p = interpolate(beta2p_array,log10(mom),0.01,151,250.);
    epsilonp = interpolate(epsp_array,log10(mom),0.01,151,250.);
    sigman = interpolate(sigman_array,log10(mom),0.01,151,250.);
    beta2n = interpolate(beta2n_array,log10(mom),0.01,151,250.);
    epsilonn = interpolate(epsn_array,log10(mom),0.01,151,250.);
    break;
  case(3):
    //cout <<"pi-" << endl;
    sigmap = interpolate(sigmap_array,log10(mom),0.01,151,250.);
    beta2p = interpolate(beta2p_array,log10(mom),0.01,151,250.);
    epsilonp = interpolate(epsp_array,log10(mom),0.01,151,250.);
    sigman = interpolate(sigman_array,log10(mom),0.01,151,250.);
    beta2n = interpolate(beta2n_array,log10(mom),0.01,151,250.);
    epsilonn = interpolate(epsn_array,log10(mom),0.01,151,250.);
    break;
  case(4):
    //cout <<"rho0" << endl;
    sigman=sigmap = (interpolate(sigman_array,log10(mom),0.01,151,250.)+interpolate(sigmap_array,log10(mom),0.01,151,250.))*0.5;
    beta2n=beta2p = (interpolate(beta2n_array,log10(mom),0.01,151,250.)+interpolate(beta2p_array,log10(mom),0.01,151,250.))*0.5;
    epsilonn = epsilonp = (interpolate(epsn_array,log10(mom),0.01,151,250.)+interpolate(epsp_array,log10(mom),0.01,151,250.))*0.5;
    break;
  default:
    cerr << "Particle type is not yet supported!!: " << particletype << endl;
    exit(1);
  }
  return;
}


void FastParticle::setScatter(double sigma, double beta, double eps, double p_dil){
  setSigma(sigma);
  beta2n=beta2p=beta*pow(HBARC/1000.,2.);
  epsilonn=epsilonp=eps;
  decay_dil = Gamma*sqrt(1.-p_dil*p_dil/(mass*mass+p_dil*p_dil));
  userset=1;
}


const double FastParticle::sigmap_array[] = { 17.2781, 16.306, 15.2518, 14.1561, 13.0553, 11.9802, 10.9544, 9.99498, 9.11227, 8.31128,
				      7.59246, 6.95286, 6.3871, 5.8882, 5.44828, 5.05905, 4.71225, 4.39995, 4.1149, 3.85078,
				      3.60257, 3.36681, 3.14198, 2.92868, 2.72983, 2.55078, 2.3995, 2.28718, 2.21259, 2.07051,
				      1.93132, 1.80425, 1.69566, 1.60908, 1.54576, 1.50547, 1.48705, 1.48886, 1.50895, 1.54515,
				      1.59509, 1.65623, 1.72583, 1.80112, 1.87943, 1.95834, 2.03596, 2.11099, 2.18295, 2.2521,
				      2.31941, 2.38646, 2.4552, 2.52776, 2.60623, 2.69242, 2.7877, 2.89274, 3.00733, 3.1302,
				      3.25887, 3.38956, 3.51728, 3.63598, 3.73898, 3.81961, 3.87195, 3.89164, 3.8766, 3.82745,
				      3.74765, 3.64312, 3.52165, 3.39199, 3.26299, 3.14286, 3.03859, 2.95563, 2.89769, 2.86662,
				      2.86223, 2.88189, 2.91401, 2.94705, 2.97388, 2.99496, 3.01079, 3.02185, 3.02861, 3.03154,
				      3.03106, 3.02761, 3.02157, 3.0133, 3.00315, 2.99142, 2.9784, 2.96433, 2.94946, 2.93398,
				      2.91807, 2.90191, 2.88562, 2.86933, 2.85315, 2.83717, 2.82146, 2.80608, 2.7911, 2.77654,
				      2.76244, 2.74883, 2.73572, 2.72312, 2.71103, 2.69947, 2.68841, 2.67785, 2.66779, 2.6582,
				      2.64907, 2.64038, 2.63212, 2.62425, 2.61676, 2.60963, 2.60282, 2.59633, 2.59013, 2.58418,
				      2.57848, 2.573, 2.56771, 2.5626, 2.55764, 2.55282, 2.54812, 2.54353, 2.53902, 2.53457,
				      2.53019, 2.52585, 2.52154, 2.51725, 2.51297, 2.50869, 2.50442, 2.50013, 2.49582, 2.4915,
				      2.48716} ;

const double FastParticle::beta2p_array[] =	{ -0.0485909, 0.143607, 0.286071, 0.387333, 0.454878, 0.495238, 0.51407, 0.516242, 0.505898, 0.486539,
  0.46108, 0.431917, 0.400986, 0.369813, 0.339567, 0.311108, 0.285027, 0.261692, 0.241277, 0.223805,
  0.209168, 0.197164, 0.187519, 0.179907, 0.173971, 0.169342, 0.165654, 0.162554, 0.159713, 0.15684,
  0.153682, 0.150031, 0.145731, 0.140676, 0.134811, 0.128134, 0.120689, 0.11257, 0.10391, 0.0948794,
  0.0856826, 0.076548, 0.0677235, 0.0594694, 0.0520508, 0.0457307, 0.0407623, 0.037382, 0.0358023, 0.036205,
  0.0387351, 0.0434948, 0.0505389, 0.0598702, 0.0714361, 0.0851269, 0.100774, 0.118149, 0.136968, 0.156892,
  0.177532, 0.198456, 0.219196, 0.239257, 0.258134, 0.275319, 0.290324, 0.302695, 0.312037, 0.318036,
  0.320487, 0.319323, 0.314648, 0.306772, 0.296252, 0.28393, 0.259021, 0.251092, 0.245262, 0.241262,
  0.23884, 0.237765, 0.237824, 0.238821, 0.240579, 0.242935, 0.245743, 0.248871, 0.252202, 0.255632,
  0.259069, 0.262434, 0.26566, 0.268689, 0.271475, 0.27398, 0.276176, 0.278042, 0.279565, 0.28074,
  0.281567, 0.282052, 0.282208, 0.28205, 0.281599, 0.28088, 0.27992, 0.278749, 0.2774, 0.275906,
  0.274303, 0.272626, 0.270914, 0.269202, 0.267527, 0.265925, 0.264431, 0.263077, 0.261896, 0.260918,
  0.260171, 0.25968, 0.259467, 0.259552, 0.259954, 0.260685, 0.261755, 0.263173, 0.264941, 0.267061,
  0.269528, 0.272335, 0.275473, 0.278927, 0.282679, 0.28671, 0.290994, 0.295505, 0.300213, 0.305084,
  0.310085, 0.315177, 0.320321, 0.325476, 0.3306, 0.33565, 0.340583, 0.345354, 0.349921, 0.354241,
  0.358274} ;

const double FastParticle::epsp_array[] = { 40.2787, 27.7732, 18.5131, 11.7765, 6.97527, 3.63519, 1.37815, -0.0934038, -1.01015, -1.54774,
  -1.83704, -1.973, -2.02214, -2.0289, -2.02098, -2.01371, -2.01368, -2.02155, -2.0344, -2.04739,
  -2.05512, -2.05251, -2.0354, -2.00089, -1.94752, -1.87523, -1.78524, -1.67986, -1.5622, -1.43595,
  -1.30504, -1.17342, -1.0448, -0.922449, -0.809057, -0.706621, -0.616385, -0.538837, -0.473743, -0.420227,
  -0.376878, -0.341891, -0.313216, -0.288726, -0.266379, -0.244371, -0.221272, -0.196136, -0.168578, -0.138811,
  -0.107641, -0.0764157, -0.0469264, -0.0212679, -0.00165783, 0.00977428, 0.0112181, 0.00141802, -0.0201046, -0.0528668,
  -0.0953037, -0.144745, -0.197527, -0.249278, -0.295439, -0.332052, -0.350241, -0.371676, -0.38652, -0.393895,
  -0.39396, -0.387572, -0.376013, -0.360774, -0.343385, -0.325295, -0.307782, -0.291901, -0.278453, -0.267978,
  -0.260759, -0.256849, -0.256101, -0.258199, -0.262705, -0.269098, -0.276811, -0.285273, -0.29394, -0.302321,
  -0.310004, -0.316671, -0.322104, -0.32619, -0.328917, -0.330365, -0.33069, -0.330105, -0.328861, -0.327221,
  -0.325437, -0.327752, -0.325751, -0.323747, -0.32174, -0.319731, -0.31772, -0.315706, -0.313689, -0.31167,
  -0.309648, -0.307624, -0.305597, -0.303568, -0.301536, -0.299502, -0.297465, -0.295426, -0.293384, -0.29134,
  -0.289293, -0.287243, -0.285191, -0.283137, -0.28108, -0.27902, -0.276958, -0.274893, -0.272826, -0.270756,
  -0.268684, -0.266609, -0.264532, -0.262452, -0.26037, -0.258285, -0.256197, -0.254107, -0.252015, -0.24992,
  -0.247822, -0.245722, -0.24362, -0.241515, -0.239407, -0.237297, -0.235184, -0.233069, -0.230951, -0.228831,
  -0.226708} ;

const double FastParticle::sigman_array[] = { 6.01242, 5.70712, 5.38087, 5.04634, 4.71471, 4.39525, 4.09517, 3.81963, 3.57196, 3.35396,
  3.16623, 3.00846, 2.87971, 2.77865, 2.70364, 2.65289, 2.62443, 2.61612, 2.62557, 2.65014,
  2.68683, 2.73236, 2.78323, 2.83602, 2.88772, 2.93638, 2.98186, 3.02666, 3.0771, 3.14463,
  3.24794, 3.51793, 3.75167, 4.06071, 4.35431, 4.55696, 4.62355, 4.54971, 4.36772, 4.13054,
  3.8935, 3.70243, 3.58992, 3.57727, 3.67775, 3.89777, 4.23356, 4.66272, 5.13272, 5.55444,
  5.81446, 5.81624, 5.53835, 5.02546, 4.57165, 4.21632, 3.96029, 3.78954, 3.68587, 3.63125,
  3.60932, 3.60609, 3.61021, 3.61323, 3.60977, 3.5973, 3.57588, 3.54748, 3.5153, 3.48311,
  3.45458, 3.43288, 3.42036, 3.41836, 3.42719, 3.44609, 3.47332, 3.50628, 3.54167, 3.57582,
  3.60495, 3.62559, 3.63497, 3.63136, 3.61435, 3.58487, 3.54516, 3.49847, 3.44867, 3.39982,
  3.35566, 3.31922, 3.29238, 3.27558, 3.26757, 3.26535, 3.26429, 3.25874, 3.23093, 3.21291,
  3.19527, 3.17801, 3.16112, 3.14459, 3.12842, 3.11259, 3.0971, 3.08194, 3.0671, 3.05258,
  3.03836, 3.02445, 3.01083, 2.99749, 2.98444, 2.97166, 2.95915, 2.9469, 2.93491, 2.92316,
  2.91167, 2.90041, 2.88939, 2.87859, 2.86802, 2.85767, 2.84753, 2.8376, 2.82788, 2.81836,
  2.80903, 2.7999, 2.79095, 2.78219, 2.77361, 2.7652, 2.75696, 2.74889, 2.74099, 2.73325,
  2.72566, 2.71823, 2.71095, 2.70382, 2.69683, 2.68998, 2.68328, 2.6767, 2.67026, 2.66395,
  2.65777} ;

const double FastParticle::beta2n_array[] = { -1.87567, -0.90376, -0.318568, 0.0128662, 0.184313, 0.259972, 0.282194, 0.277594, 0.26181, 0.243137,
  0.225239, 0.209102, 0.194393, 0.180349, 0.166306, 0.151965, 0.137469, 0.12336, 0.110459, 0.0997147,
  0.0920482, 0.0882169, 0.0887069, 0.0936658, 0.102877, 0.115773, 0.131486, 0.148926, 0.166879, 0.184116,
  0.199505, 0.212112, 0.221296, 0.226768, 0.228634, 0.227401, 0.223952, 0.219488, 0.21544, 0.213362,
  0.214793, 0.22112, 0.233436, 0.252405, 0.278151, 0.310178, 0.347331, 0.387808, 0.429233, 0.468787,
  0.503414, 0.530072, 0.546046, 0.54929, 0.538782, 0.514848, 0.479421, 0.436178, 0.390469, 0.348979,
  0.330735, 0.30365, 0.292228, 0.293364, 0.302434, 0.314748, 0.326398, 0.334669, 0.338128, 0.336513,
  0.330494, 0.321373, 0.310778, 0.300387, 0.291688, 0.285819, 0.283458, 0.28479, 0.289526, 0.296978,
  0.306158, 0.31591, 0.325046, 0.33248, 0.337336, 0.339045, 0.337397, 0.332561, 0.325063, 0.315732,
  0.305605, 0.29582, 0.287483, 0.281535, 0.278633, 0.27905, 0.282615, 0.288696, 0.296247, 0.303914,
  0.3102, 0.313691, 0.313314, 0.30861, 0.299981, 0.288849, 0.277669, 0.269695, 0.268391, 0.27635,
  0.293187, 0.292679, 0.292596, 0.292911, 0.293596, 0.294623, 0.295967, 0.297603, 0.299505, 0.301651,
  0.304016, 0.306579, 0.309318, 0.312212, 0.315241, 0.318385, 0.321626, 0.324946, 0.328327, 0.331753,
  0.335209, 0.338679, 0.342149, 0.345605, 0.349035, 0.352428, 0.35577, 0.359053, 0.362265, 0.365398,
  0.368444, 0.371395, 0.374245, 0.376986, 0.379614, 0.382123, 0.384511, 0.386773, 0.388908, 0.390914,
  0.392789} ;

const double FastParticle::epsn_array[] = { 0.220782, 0.39322, 0.489897, 0.528353, 0.524508, 0.492505, 0.44458, 0.390976, 0.339892, 0.297477,
  0.267874, 0.253323, 0.254303, 0.269736, 0.297233, 0.333387, 0.374099, 0.414926, 0.45144, 0.479586,
  0.496009, 0.498344, 0.485439, 0.457499, 0.416129, 0.364263, 0.305979, 0.246191, 0.190234, 0.143366,
  0.110213, 0.0942141, 0.097116, 0.118584, 0.155994, 0.204469, 0.257203, 0.306086, 0.342616, 0.359026,
  0.349499, 0.311294, 0.245565, 0.157647, 0.0566115, -0.0459978, -0.138182, -0.209941, -0.255599, -0.275228,
  -0.274377, -0.261664, -0.244639, -0.22576, -0.174621, -0.179313, -0.159157, -0.135739, -0.118666, -0.11035,
  -0.109356, -0.112648, -0.116993, -0.119757, -0.119255, -0.11482, -0.106675, -0.0957192, -0.0832652, -0.0707917,
  -0.0597236, -0.0512651, -0.0462866, -0.0452688, -0.0482955, -0.0550886, -0.0650731, -0.0774616, -0.0913482, -0.105801,
  -0.119946, -0.133034, -0.144492, -0.153949, -0.161242, -0.166403, -0.169629, -0.171237, -0.171621, -0.171194,
  -0.170344, -0.169394, -0.168575, -0.168015, -0.167737, -0.16768, -0.167728, -0.167745, -0.167615, -0.167279,
  -0.166756, -0.166148, -0.165627, -0.165385, -0.165576, -0.166225, -0.167154, -0.161663, -0.161403, -0.161117,
  -0.160805, -0.160468, -0.160105, -0.159718, -0.159306, -0.158869, -0.158409, -0.157924, -0.157416, -0.156884,
  -0.156328, -0.15575, -0.155149, -0.154526, -0.15388, -0.153212, -0.152522, -0.151811, -0.151078, -0.150325,
  -0.14955, -0.148755, -0.14794, -0.147104, -0.146248, -0.145373, -0.144479, -0.143565, -0.142632, -0.141681,
  -0.140711, -0.139723, -0.138717, -0.137694, -0.136653, -0.135594, -0.134519, -0.133427, -0.132319, -0.131194,
  -0.130053} ;
  