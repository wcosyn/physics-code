#include "PionTCross.hpp"
#include <AuxFunction.hpp>
#include <TRotation.h>
#include <Utilfunctions.hpp>

using namespace std;

PionTCross::PionTCross(const int nucleus, const double p_max, const string dir,
  const bool user_set, const double user_sigma, const double precision, const int integr, const int max_Eval
):
homedir(dir),
pmax(p_max),
userset(user_set),
usersigma(user_sigma),
nucleusthick(nucleus,dir),
pdistgrid(NULL),
pfsigrid(NULL),
prec(precision),
integrator(integr),
abserror(1.E-012),
maxEval(max_Eval){

  pfsigrid = new GlauberGridThick*[nucleusthick.getTotalLevels()];
//   pdistgrid = new DistMomDistrGrid*[nucleusthick.getTotalLevels()];
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
    pfsigrid[i] = new GlauberGridThick(60,20,5,&nucleusthick,prec,2,homedir);
//     pdistgrid[i] = new DistMomDistrGrid(i, pmax, 30,20,5,pfsigrid[i],1.E-03,2,2E04,0.,homedir);
  }
  nrofcross=pfsigrid[0]->getNumber_of_grids()+1;
  cout << "Number of cross sections: " << nrofcross << endl;
}

PionTCross::~PionTCross(){
  
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
//     delete pdistgrid[i];
    delete pfsigrid[i];
  }
//   delete [] pdistgrid;
  delete [] pfsigrid;
}

void PionTCross::getCross(double *results, const double Ebeam, const double Eout, const double theta_e){

    
  
  int res=90;
  unsigned count=0;
    numint::vector_d ret(nrofcross,0.);
    numint::array<double,4> lower = {{Eout-200,0.,-1.,0.}};
    numint::array<double,4> upper = {{Eout+200,pmax,1.,2.*PI}};

    PionTCross::Ftor_pion F;
    F.cross = this;
    F.Ebeam=Ebeam;
    F.theta_e=theta_e;
    numint::mdfunction<numint::vector_d,4> mdf;
    mdf.func = &Ftor_pion::exec;
    mdf.param = &F;
    F.f=PionTCross::dist_momdistr_integral;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,abserror,1.E-03,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,abserror,1.E-03,1E03,maxEval,ret,count,0);
    for(int i=0;i<nrofcross;++i) results[i]=ret[i];
  for(int i=0;i<nrofcross;i++) results[i]*= 1.;
}



//all energies in MeV please
void PionTCross::getMomdistr(double *results, double ppi, double thetapi, double Q2, int shell, 
			    double pm, double pmcostheta, double pmphi){
  //if userset there's only one possible glaubergrid, so we don't have to check every time!
  if(!(userset&&pfsigrid[shell]->getFilledallgrid())){
    FastParticle pion(2, 0, ppi,0.,0.,Q2,0.,1.,1.,homedir);
    if(userset) pion.setScatter(usersigma,6.,-0.2,0.);
    pfsigrid[shell]->clearParticles();
    pfsigrid[shell]->addParticle(pion);
    pfsigrid[shell]->updateGrids();
    pfsigrid[shell]->clearKnockout();
  }
  //update the grid if necessary
  TRotation rot;
  rot.Rotate(-thetapi,TVector3(0.,1.,0.));
  double key=double(shell)/getNucleusthick().getTotalLevels()+int(round(thetapi*10.));
  map<double,DistMomDistrGrid>::iterator it=distgridmap.find(key);
//   cout << key << " " << thetarho << " " << shell  << endl;
  if(it==distgridmap.end()){
    distgridmap[key]=DistMomDistrGrid(shell, pmax, 30,20,5,pfsigrid[shell],1.E-03,2,5E04,0.,homedir);
//     distgridmap.insert(pair<double,DistMomDistrGrid>(key,DistMomDistrGrid(shell, pmax, 30,20,5,pfsigrid[shell],1.E-03,2,2E04,0.,homedir)));
    distgridmap[key].updateGrids(pfsigrid[shell],shell,rot);
//     distgridmap[key].printRho_grid(0);
//     cout << endl << endl;
    it=distgridmap.find(key);
  }
  else distgridmap[key].updateGrids(pfsigrid[shell],shell,rot);

  results[nrofcross-1] = distgridmap[key].getRhopwGridFull_interp(pm);
  //fsi
  for(int i=0;i<nrofcross-1;i++) results[i] = distgridmap[key].getRhoGridFull_interp3(i, pm, pmcostheta, pmphi);
//    if(pm>200.&&results[0]>50.) {distgridmap[key].printRho_grid(0); cout  << endl << endl << endl;}
  
  
//   pdistgrid[shell]->updateGrids(pfsigrid[shell],shell,rot);
//   //plane-wave
//   results[nrofcross-1] = pdistgrid[shell]->getRhopwGridFull_interp(pm);
//   //fsi
//   for(int i=0;i<nrofcross-1;i++) results[i]= pdistgrid[shell]->getRhoGridFull_interp3(i, pm, pmcostheta, pmphi);
}
  
  

void PionTCross::dist_momdistr_integral(numint::vector_d & results, double Eout, double pm, double costheta, double phi, 
                      PionTCross & cross, double Ebeam, double theta_e){
  results=numint::vector_d(cross.getNrofcross(),0.);
  if(pm==0.){
    return;
  }
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  double sintheta = sqrt(1.-costheta*costheta);
  double Ep = sqrt(pm*pm+MASSP*MASSP);
  
  double omega = Ebeam-Eout;
  double Q2 = 2.*Ebeam*Eout*(1.-cos(theta_e));
  double qvec = sqrt(omega*omega+Q2);
  cout << "Q2 " << Q2 << " omega " << omega << " qvec " << qvec << endl;
  double fScatterAngle = costheta;;
  double epsilon = 1./(1.+2.*(1.+omega*omega/Q2)
	 	     *((1.-fScatterAngle)/(1.+fScatterAngle)));

  double massA = cross.getNucleusthick().getMassA();
  double massX=massA+100;

  double CC=pow(omega+massA,2);
  double DD=CC+massX*massX+qvec*qvec-MASSPI*MASSPI;
  double a=4*(CC-qvec*qvec);
  double b=4*qvec*(DD-2*CC);
  double c=4*CC*(massX*massX+qvec*qvec)-DD*DD;
  double discr=b*b-4*a*c;
  if(discr<0) { cout << "discr <0 " << discr << endl; return;}
  double p_pi = (-b+sqrt(discr))*0.5/a;
  double E_pi = sqrt(p_pi*p_pi+MASSPI*MASSPI);
  double t = (Q2 - MASSPI*MASSPI+2*omega*E_pi-2*qvec*p_pi)/1e06; //in GeV^2
  double W = sqrt((-Q2+MASSP*MASSP+2*omega*MASSP)/1e06); // in GeV
  double sfree = -Q2+MASSP*MASSP+2*omega*sqrt(MASSP*MASSP+pm*pm)-2*qvec*pm*costheta; //mandelstam variable for the free nucleon process

  for(int i=0;i<cross.getNucleusthick().getPLevels();i++){
    double intresults[cross.getNrofcross()];
    cross.getMomdistr(intresults,p_pi,0.,Q2,i,pm,costheta,phi);	  
    double front=ALPHA/PI*Eout/Ebeam/Q2/(1-epsilon)*crosselectronfree(-t, Q2*1.E-06, W, 0., 0., epsilon)  //free cross section
                  *MASSP/(Ep-pm*costheta) //flux correction for moving proton
                  *p_pi/(8.*PI*PI)*massX/(MASSP*MASSN/pow(sfree-MASSP*MASSP,2)/4./PI);  //other factors
    //correct for not completely full shells
    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
    for(int dd=0;dd<cross.getNrofcross();dd++) { //cout << i << " " << dd<< " " << front << " " << intresults[dd] << endl;
                                                  results[dd]+=front*intresults[dd];}
  }
  for(int i=0;i<cross.getNrofcross();i++){
    results[i]*=pm*pm;
  }
  cout << pm << " " << costheta << " " << phi << " " << Eout << " "  << p_pi << " " << results[0] << " " << results[cross.getNrofcross()-1] << endl;
}

double PionTCross::crosselectronfree(double t, double Q2, double W, double phipq, double thetapq, double epsilon){

  double sigmaT, sigmaL, sigmaTT, sigmaLT;
  double A=350, B=16, C=7.5, D=4.5, E=2, F=5, G=0.79, H=3.4, K=1.1, L=3.6;

  sigmaT=D/Q2+E/Q2/Q2;
  sigmaL=A*Q2/pow(1+1.77*Q2+0.05*Q2*Q2,2)*exp((B-C*log(Q2))*abs(t));
  sigmaTT=-F*abs(t)*sin(thetapq)*sin(thetapq)/Q2/Q2/pow(abs(t)+MASSPI*MASSPI/1e06,2);
  sigmaLT=sin(thetapq)*(exp(G-H*t/sqrt(Q2))+K-L/Q2/Q2);

  double sigma= sigmaT+epsilon*(sigmaL+sigmaTT*cos(phipq))+sqrt(epsilon*(epsilon+1))*sigmaLT*cos(2*phipq);
  sigma*=8.539*0.5/PI/pow(W-MASSP*MASSP/1e06,2);

  //extra corrections
  if(abs(Q2-1.1)<=0.4) {
    return sigma*((-47.5984)+(43.4145)*W+(-9.64264)*W*W)
      *((1.32289)+(-0.698424)*Q2+(0.35561)*Q2*Q2)
      *((1.17152)+(-7.03367)*t+(52.053)*t*t)
      *((1.0612)+(0.147858)*cos(phipq)+(-0.0430268)*cos(2.0*phipq));
  }
  if(abs(Q2-2.15)<=0.4){
    return sigma*((-23.1723)+(20.6505)*W+(-4.37408)*W*W)
      *((2.29646)+(-1.11745)*Q2+(0.229736)*Q2*Q2)
      *((0.704879)+(1.61954)*t+(0.0859429)*t*t)
      *((0.979176)+(0.044882)*cos(phipq)+(-0.0743073)*cos(2.0*phipq)) ;
  }
  if(abs(Q2-3)<=0.4){
    return sigma*((-6.14191)+(5.64149)*W+(-1.0843)*W*W)
      *((2.43486)+(-0.888779)*Q2+(0.136267)*Q2*Q2)
      *((0.745356)+(1.22215)*t+(-1.24105)*t*t)
      *((0.962609)+(-0.0608404)*cos(phipq)+(-0.0084712)*cos(2.0*phipq));
  }
  if(abs(Q2-4)<=0.4){
    return sigma*((-7.8696)+(6.48878)*W+(-1.16624)*W*W)
      *((-0.703888)+(0.814839)*Q2+(-0.0957087)*Q2*Q2)
      *((0.723372)+(0.140101)*t+(0.809151)*t*t)
      *((1.00054)+(-0.100002)*cos(phipq)+(0.00780768)*cos(2.0*phipq));
  }
  if(abs(Q2-4.8)<=0.4){
    return sigma*((-11.1202)+(9.4995)*W+(-1.86788)*W*W)
      *((2.27339)+(-0.469771)*Q2+(0.0421723)*Q2*Q2)
      *((1.08961)+(-1.06851)*t+(1.36125)*t*t)
      *((0.89789)+(-0.118188)*cos(phipq)+(-0.0350948)*cos(2.0*phipq));
  }
  //12 GeV, no parametrization, just return 1.
  if((Q2>=4.9)&&(Q2<=10)) return 1.;
  else cout << "other corr needed" << endl;
  return 0.;
}
