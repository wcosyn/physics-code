#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std; 

#include "FastParticle.hpp"
#include <constants.hpp>
#include <Utilfunctions.hpp>

FastParticle::FastParticle(const int type, const int inc, const double momentum,
			   const double ptheta, const double pphi, const double hard_scale, const double Gamma, const string dir)
:particletype(type),
incoming(inc), 
p(momentum),
theta(ptheta),
phi(pphi),
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
    default:
      cerr << "Particle type is not yet supported!!: " << particletype << endl;
      exit(1);
  }
  E = sqrt(mass*mass+p*p);
  decay_dil = Gamma*sqrt(1.-p*p/(E*E));
  if(hard_scale<nkt_sq) cout << "Warning: hard scale is too low for color transparency effects for this particle, just so you know!" << endl;
  
}



FastParticle::FastParticle(const int type, const int inc, const TVector3 &pvec, 
			   const double hard_scale, const double Gamma, const string dir)
:particletype(type),
incoming(inc), 
p(pvec.Mag()),
theta(pvec.Theta()),
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
    default:
      cerr << "Particle type is not yet supported!!: " << particletype << endl;
      exit(1);
  }
  E = sqrt(mass*mass+p*p);
  decay_dil = Gamma/sqrt(1.-p*p/E*E);
  if(hard_scale<nkt_sq) cout << "Warning: hard scale is too low for color transparency effects for this particle, just so you know!" << endl;
  
}

FastParticle::FastParticle(const FastParticle &Copy):
particletype(Copy.getParticletype()),
incoming(Copy.getIncoming()),
p(Copy.getP()),
theta(Copy.getTheta()),
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
     (*this)=FastParticle(rhs);
     
     
     
  }

  return *this;
}



FastParticle::~FastParticle(){
  //cout << "Deleting Particle object..." << endl;
}

//getters
int FastParticle::getParticletype() const{
  return particletype;
}

bool FastParticle::getIncoming() const{
  return incoming;
}

double FastParticle::getP() const{
  return p;
}

double FastParticle::getTheta() const{
  return theta;
}

double FastParticle::getPhi() const{
  return phi;
}

double FastParticle::getEx() const{
  return ex;
}

double FastParticle::getEy() const{
  return ey;
}

double FastParticle::getEz() const{
  return ez;
}

double FastParticle::getHitz() const{
  return hitz;
}

double FastParticle::getHitbnorm() const{
  return hitbnorm;
}

const double* FastParticle::getHitbvec() const{
  return hitb;
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

double FastParticle::getSigmap() const{
  return sigmap;
}

double FastParticle::getBeta2p() const{
  return beta2p;
}

double FastParticle::getEpsilonp() const{
  return epsilonp;
}

double FastParticle::getSigman() const{
  return sigman;
}

double FastParticle::getBeta2n() const{
  return beta2n;
}

double FastParticle::getEpsilonn() const{
  return epsilonn;
}

double FastParticle::getSigma(bool proton) const{
  return proton? sigmap:sigman;
}

double FastParticle::getBetasq(bool proton) const{
  return proton? beta2p:beta2n;
}

double FastParticle::getEpsilon(bool proton) const{
  return proton? epsilonp:epsilonn;
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
			  *(1.-I*getEpsilon(level,pnucleus))
			  /(4.*PI*getBetasq(level,pnucleus));
}

complex<double> FastParticle::getScatterfront(bool proton) const{
  return getSigma(proton)
			  *(1.-I*getEpsilon(proton))
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
    
double FastParticle::getHardScale() const{
  return hardscale;
}

double FastParticle::getLc() const{
  return lc;
}
  
double FastParticle::getNkt_sq() const{
  return nkt_sq;
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


void FastParticle::setScatter(double sigma, double beta, double eps){
  setSigma(sigma);
  beta2n=beta2p=beta*pow(HBARC/1000.,2.);
  epsilonn=epsilonp=eps;
  userset=1;
}
  