#include "Model.hpp"
#include <GlauberGridThick.hpp>
#include <OneGlauberGrid.hpp>
#include <TSpinor.h>
#include <TMFSpinor.hpp>
#include <Utilfunctions.hpp>





Model::Model(MeanFieldNucleusThick *pnucleus, int setSRC, int setthick, string dir)
:SRC(setSRC), thick(setthick), pnucl(pnucleus), homedir(dir){
}


Model::~Model(){
  
}


complex<double> Model::getFreeMatrixEl(TKinematics2to2 &tk, int spinin, int spinout, int photonpol){
  double costheta=tk.GetCosthYlab();
  double sintheta=sqrt(1.-costheta*costheta);
  static FourVector<complex<double> > polVectorPlus(0.,
						    -1./sqrt(2.),
						    complex<double>(0.,-1./sqrt(2.)),
						 0.);
  static FourVector<complex<double> > polVectorMin(0.,
						   1./sqrt(2.),
						   complex<double>(0.,-1./sqrt(2.)),
						   0.);
  static FourVector<complex<double> > polVector0(1.,0.,0.,0.);
  double pout=tk.GetPYlab();
  double qvec = tk.GetKlab();
  double Eout = sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pout*pout);
  double pin = sqrt(qvec*qvec+pout*pout-2.*pout*qvec*costheta);
  double Ein = sqrt(pin*pin+tk.GetHyperonMass()*tk.GetHyperonMass());
/*  double wbar = Eout-Ein;
  double Q2bar = qvec*qvec-wbar*wbar;*/
  FourVector<double> q(tk.GetWlab(),0.,0.,qvec);
  FourVector<double> pf(Eout,pout*sintheta,0.,pout*costheta);
  FourVector<double> Pin(Ein,pf[1],0.,pf[3]-q[3]); 
  J= new NucleonEMOperator(tk.GetQsquared(),1,0);
  GammaStructure Jcontr;
  if(photonpol==0)  Jcontr = J->getCC2(q)*polVector0;
  else if(photonpol==-1) Jcontr= J->getCC2(q)*polVectorMin;
  else if(photonpol==1) Jcontr=J->getCC2(q)*polVectorPlus;
/*  if(photonpol==0)  Jcontr = J.getCC3(q,Pin,pf)*polVector0;
  else if(photonpol==-1) Jcontr= J.getCC3(q,Pin,pf)*polVectorMin;
  else if(photonpol==1) Jcontr=J.getCC3(q,Pin,pf)*polVectorPlus;*/
/*  if(photonpol==0) Jcontr = J.getCC1(Pin,pf)*polVector0;
  else if(photonpol==-1) Jcontr= J.getCC1(Pin,pf)*polVectorMin;
  else if(photonpol==1) Jcontr=J.getCC1(Pin,pf)*polVectorPlus;*/
   
  delete J;
  if(spinout==-1){
    if(spinin==-1) return TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr*
		      TSpinor(Pin,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity);
    else if(spinin==1) return TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr*
		      TSpinor(Pin,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity);
    else {cerr << "invalid spinin " << spinin << endl; exit(1);}
  }
  else if(spinout==1){
    if(spinin==-1) return TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr*
		      TSpinor(Pin,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity);
    else if(spinin==1) return TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr*
		      TSpinor(Pin,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity);
    else {cerr << "invalid spinin " << spinin << endl; exit(1);}
  }
  else {cerr << "invalid spinout " << spinout << endl; exit(1);}
  
  
}

complex<double> Model::getMatrixEl(TKinematics2to2 &tk, int spinout, int photonpol, int shellindex, int m, int CT, int pw){

  double costheta=tk.GetCosthYlab();
  double sintheta=sqrt(1.-costheta*costheta);
  static FourVector<complex<double> > polVectorPlus(0.,
						    -1./sqrt(2.)*costheta,
						    complex<double>(0.,-1./sqrt(2.)),
						 -1/sqrt(2.)*sintheta);
  static FourVector<complex<double> > polVectorMin(0.,
						   1./sqrt(2.)*costheta,
						   complex<double>(0.,-1./sqrt(2.)),
						   1./sqrt(2.)*sintheta);
  static FourVector<complex<double> > polVector0(1.,0.,0.,0.);
  static FourVector<complex<double> > polVectorX(0.,costheta,0.,sintheta);
  static FourVector<complex<double> > polVectorY(0.,0.,1.,0.);
  static FourVector<complex<double> > polVectorZ(0.,-sintheta,0.,costheta);
  
  shell=shellindex;
  mm=m;
  J= new NucleonEMOperator(tk.GetQsquared(),1,0);
  FastParticle proton(0, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()/1.e06,0.,homedir);
  if(!pw){
    if(SRC||thick) grid = new GlauberGridThick(60,18,5,pnucl,homedir);
    else grid = new OneGlauberGrid(60,18,pnucl,homedir);
    grid->addParticle(proton);
    grid->fillGrids();
    grid->clearKnockout();
    grid->addKnockout(shellindex,m);
    
  }
 
  GammaStructure Jcontr;
  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  /* FourVector<double> pi=pf-q;
  pi[0]=sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pi[1]*pi[1]+pi[2]*pi[2]+pi[3]*pi[3]);*/
  pm=TVector3(-q[1],0.,pf[3]-q[3]);
  if(photonpol==0)  Jcontr = J->getCC2(q)*polVector0;
  /* else if(photonpol==1)  Jcontr = J.getCC2(q)*polVectorX;
  else if(photonpol==2)  Jcontr = J.getCC2(q)*polVectorY;
  else if(photonpol==3)  Jcontr = J.getCC2(q)*polVectorZ;*/
  else if(photonpol==-1) Jcontr= J->getCC2(q)*polVectorMin;
  else if(photonpol==1) Jcontr=J->getCC2(q)*polVectorPlus;
  else if(photonpol==2) Jcontr=J->getCC2(q)*polVectorX;
  else{ cerr << "invalid photon pol" << endl;  exit(1); }
  if(spinout==-1) barcontract = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr;
  else if(spinout==1) barcontract = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr;
  else{ cerr << "invalid nucl pol" << endl;  exit(1); }
  
  complex<double> results[2];
  double restimate=0.,thestimate=0.,phiestimate=0.;
  rombergerN(this,&Model::intJR,0.,pnucl->getRange(),2,results,PREC,3,7,&restimate,pw,&thestimate, &phiestimate);
  if(!pw) delete grid;
  complex<double> result;
  delete J;
  if(CT) return results[1];
  else return results[0];
}


void Model::getMatrixEl(TKinematics2to2 &tk, Matrix<2,3> & matrixel, int shellindex, int m, int CT, int pw){

  double costheta=tk.GetCosthYlab();
  if(costheta>1.) costheta=1.;
  if(costheta<-1.) costheta=-1.;
  double sintheta=sqrt(1.-costheta*costheta);
  static FourVector<complex<double> > polVectorPlus(0.,
						    -1./sqrt(2.)*costheta,
						    complex<double>(0.,-1./sqrt(2.)),
						 -1/sqrt(2.)*sintheta);
  static FourVector<complex<double> > polVectorMin(0.,
						   1./sqrt(2.)*costheta,
						   complex<double>(0.,-1./sqrt(2.)),
						   1./sqrt(2.)*sintheta);
  static FourVector<complex<double> > polVector0(1.,0.,0.,0.);
 
  shell=shellindex;
  mm=m;
  J= new NucleonEMOperator(tk.GetQsquared(),1,0);
  FastParticle proton(0, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()/1.e06,0.,homedir);
  if(!pw){
    if(SRC||thick) grid = new GlauberGridThick(60,18,5,pnucl,homedir);
    else grid = new OneGlauberGrid(60,18,pnucl,homedir);
    grid->addParticle(proton);
    grid->fillGrids();
    grid->clearKnockout();
    grid->addKnockout(shellindex,m);
//     grid->printFsi_grid();
//     if(SRC) dynamic_cast<AbstractFsiGridThick *>(grid)->printFsi_src_grid();
    
  }
  GammaStructure Jcontr0, Jcontrmin, Jcontrplus;
  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  /* FourVector<double> pi=pf-q;
  pi[0]=sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pi[1]*pi[1]+pi[2]*pi[2]+pi[3]*pi[3]);*/
  pm=TVector3(-q[1],0.,pf[3]-q[3]);
  Jcontr0 = J->getCC2(q)*polVector0;
  Jcontrmin= J->getCC2(q)*polVectorMin;
  Jcontrplus=J->getCC2(q)*polVectorPlus;
  barcontract0down = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr0;
  barcontract0up = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr0;
  barcontractmindown = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontrmin;
  barcontractminup = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontrmin;
  barcontractplusdown = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontrplus;
  barcontractplusup = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontrplus;
  
  
  complex<double> results[12];
  double restimate=0.,thestimate=0.,phiestimate=0.;
  rombergerN(this,&Model::intJR12,0.,pnucl->getRange(),12,results,PREC,3,8,&restimate,pw,&thestimate, &phiestimate);
  if(!pw) delete grid;
  complex<double> result;
  if(CT){
    for(int j=0;j<2;++j){
      for(int i=0;i<3;i++) matrixel(j,i)=results[j*6+2*i+1];
    }
  }
  else{
    for(int j=0;j<2;++j){
      for(int i=0;i<3;i++) matrixel(j,i)=results[j*6+2*i];
    }
  }
  delete J;
}

void Model::intJR(const double r, complex<double> *results, va_list ap){
  
  int pw = va_arg(ap,int);
  double *pthetaestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);

  if(!pw) grid->setRinterp(r);
  rombergerN(this,&Model::intJCosTheta,-1.,1.,2,results,PREC,3,7,pthetaestimate, pw, r, pphiestimate);
  results[0]*=r;
  results[1]*=r;
  return;
}
void Model::intJR12(const double r, complex<double> *results, va_list ap){
  
  int pw = va_arg(ap,int);
  double *pthetaestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);

  if(!pw) grid->setRinterp(r);
  rombergerN(this,&Model::intJCosTheta12,-1.,1.,12,results,PREC,3,6,pthetaestimate, pw, r, pphiestimate);
  for(int i=0;i<12;i++) results[i]*=r;
  return;
}

void Model::intJCosTheta(const double costheta, complex<double> *results, va_list ap){
  
  int pw = va_arg(ap,int);
  double r = va_arg(ap,double);
  double *pphiestimate = va_arg(ap,double*);
  double sintheta=sqrt(1.-costheta*costheta);
  if(!pw) grid->setCthinterp(costheta);
  rombergerN(this,&Model::intJPhi,0.,2.*PI,2,results,PREC,3,6,pphiestimate, pw,r, costheta, sintheta);
  return;
}
void Model::intJCosTheta12(const double costheta, complex<double> *results, va_list ap){
  
  int pw = va_arg(ap,int);
  double r = va_arg(ap,double);
  double *pphiestimate = va_arg(ap,double*);
  double sintheta=sqrt(1.-costheta*costheta);
  if(!pw) grid->setCthinterp(costheta);
  rombergerN(this,&Model::intJPhi12,0.,2.*PI,12,results,PREC,3,6,pphiestimate, pw,r, costheta, sintheta);
  return;
}

void Model::intJPhi(const double phi, complex<double> *results, va_list ap){

  int pw = va_arg(ap,int);
  double r = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  
  results[1]=results[0]= barcontract*TMFSpinor(*pnucl,shell,mm,r,costheta,phi)*exp(-INVHBARC*pm*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*I);
   if(!pw){
    if(SRC){
      results[0]*=dynamic_cast<AbstractFsiGridThick *>(grid)->getFsiSrcGridFull_interp1(phi);
      results[1]*=dynamic_cast<AbstractFsiCTGridThick *>(grid)->getFsiSrcCtGridFull_interp1(phi);
    }
    else{
      results[0]*=grid->getFsiGridFull_interp1(phi);
      results[1]*=grid->getFsiCtGridFull_interp1(phi);
    }
  }
  return;
}

void Model::intJPhi12(const double phi, complex<double> *results, va_list ap){

  int pw = va_arg(ap,int);
  double r = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*pnucl,shell,mm,r,costheta,phi);
  complex<double> exp_pr=exp(-INVHBARC*pm*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*I);
  results[1]=results[0]= barcontract0down*wave*exp_pr;
  results[2]=results[3]= barcontractmindown*wave*exp_pr;
  results[4]=results[5]= barcontractplusdown*wave*exp_pr;
  results[6]=results[7]= barcontract0up*wave*exp_pr;
  results[8]=results[9]= barcontractminup*wave*exp_pr;
  results[10]=results[11]= barcontractplusup*wave*exp_pr;
  if(!pw){
    if(SRC){
      complex<double> src=dynamic_cast<AbstractFsiGridThick *>(grid)->getFsiSrcGridFull_interp1(phi);
      complex<double> srcct=dynamic_cast<AbstractFsiCTGridThick *>(grid)->getFsiSrcCtGridFull_interp1(phi);
      for(int i=0;i<6;i++){
	results[2*i]*=src;
	results[2*i+1]*=srcct;
      }      
    }
    else{
      complex<double> fsi=grid->getFsiGridFull_interp1(phi);
      complex<double> fsict=grid->getFsiCtGridFull_interp1(phi);
      for(int i=0;i<6;i++){
	results[2*i]*=fsi;
	results[2*i+1]*=fsict;
      }      
    }
  }
  return;
}




