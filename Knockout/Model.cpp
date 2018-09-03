#include "Model.hpp"
#include <TMFSpinor.hpp>
#include <AuxFunction.hpp>

using namespace std;



Model::Model(MeanFieldNucleusThick *pnucleus, double precision, int integr, string dir,
	      int max_Eval, bool user_sigma, double sigma_screening)
:pnucl(pnucleus), prec(precision), integrator(integr), homedir(dir), maxEval(max_Eval),
usersigma(user_sigma), sigmascreening(sigma_screening),
gridthick(GlauberGridThick(120,36,5,pnucleus,precision,2,dir)),
onegrid(OneGlauberGrid(120,36,pnucleus,precision,2,dir)){
}


Model::~Model(){
  
}

complex<double> Model::getFreeMatrixEl(TKinematics2to2 &tk, int current, int spinin, int spinout, int photonpol){
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=tk.GetCosthYlab();
  double sintheta=sqrt(1.-costheta*costheta);
  double pout=tk.GetPYlab(); //momentum final nucleon
  double qvec = tk.GetKlab(); //virtual photon momentum
  double Eout = sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pout*pout); //energy final nucleon
  double pin = sqrt(qvec*qvec+pout*pout-2.*pout*qvec*costheta); //momentum initial nucleon
  double Ein = sqrt(pin*pin+tk.GetHyperonMass()*tk.GetHyperonMass()); //energy initial nucleon
/*  double wbar = Eout-Ein;
  double Q2bar = qvec*qvec-wbar*wbar;*/
// polarization four vectors of the virtual photon
  static FourVector<complex<double> > polVectorPlus(0.,
                                                    -1./sqrt(2.),
                                                    complex<double>(0.,-1./sqrt(2.)),
                                                    0.);
  static FourVector<complex<double> > polVectorMin(0.,
                                                   1./sqrt(2.),
                                                   complex<double>(0.,-1./sqrt(2.)),
                                                   0.);
  static FourVector<complex<double> > polVector0(1.,0.,0.,0.);
  //Frame has z-axis along q!!!
  FourVector<double> q(tk.GetWlab(),0.,0.,qvec); //photon fourvector
  FourVector<double> pf(Eout,pout*sintheta,0.,pout*costheta); //final nucleon fourvector
  FourVector<double> Pin(Ein,pf[1],0.,pf[3]-q[3]);  //initial nucleon fourvector
  J= new NucleonEMOperator(tk.GetQsquared(),1,0);
  GammaStructure Jcontr; //nucleon-photon vertex contracted with the polarization fourvector of the photon
  if(photonpol==0)  Jcontr = J->getCC(current, q, Pin, pf, 0.,0,*pnucl)*polVector0;
  else if(photonpol==-1) Jcontr= J->getCC(current, q, Pin, pf, 0.,0,*pnucl)*polVectorMin;
  else if(photonpol==1) Jcontr=J->getCC(current, q, Pin, pf, 0.,0,*pnucl)*polVectorPlus;

  /*  if(photonpol==0)  Jcontr = J.getCC3(q,Pin,pf)*polVector0;
  else if(photonpol==-1) Jcontr= J.getCC3(q,Pin,pf)*polVectorMin;
  else if(photonpol==1) Jcontr=J.getCC3(q,Pin,pf)*polVectorPlus;*/
/*  if(photonpol==0) Jcontr = J.getCC1(Pin,pf)*polVector0;
  else if(photonpol==-1) Jcontr= J.getCC1(Pin,pf)*polVectorMin;
  else if(photonpol==1) Jcontr=J.getCC1(Pin,pf)*polVectorPlus;*/
   
  delete J;
  //contract everything
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


// complex<double> Model::getMatrixEl(TKinematics2to2 &tk, int spinout, int photonpol, int shellindex, 
// 				   int m, int CT, int pw, int current, int SRC, int thick){
// 
//   //to translate the 2to2kinematics language to our particles:
//   //p is A
//   //hyperon Y is fast final nucleon
//   //kaon is residual A-1 nucleus
//   double costheta=tk.GetCosthYlab();
//   double sintheta=sqrt(1.-costheta*costheta);
//   //careful!! we compute in a frame with p_f along z-axis, so have to rotate the photon polarization 4-vectors!
// Model::polVectorPlus(0.,
// 						    -1./sqrt(2.)*costheta,
// 						    complex<double>(0.,-1./sqrt(2.)),
// 						 -1/sqrt(2.)*sintheta);
// Model::polVectorMin(0.,
// 						   1./sqrt(2.)*costheta,
// 						   complex<double>(0.,-1./sqrt(2.)),
// 						   1./sqrt(2.)*sintheta);
// Model::polVector0(1.,0.,0.,0.);
// Model::polVectorX(0.,costheta,0.,sintheta);
// Model::polVectorY(0.,0.,1.,0.);
// Model::polVectorZ(0.,-sintheta,0.,costheta);
//   
//   shell=shellindex;
//   mm=m;
//    int prot = (shellindex < pnucl->getPLevels()?0:1);
//   J= new NucleonEMOperator(tk.GetQsquared(),!prot,0);
//   FastParticle proton(prot, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()/1.e06,0.,homedir);
//   if(getUsersigma()) proton.setScreening(getSigmascreening());
//   if(!pw){
//     if(SRC||thick) grid = new GlauberGridThick(60,18,5,pnucl,prec,integrator,homedir);
//     else grid = new OneGlauberGrid(60,18,pnucl,prec,integrator,homedir);
//     grid->clearParticles();
//     grid->addParticle(proton);
//     grid->updateGrids();
//     grid->clearKnockout();
//     grid->addKnockout(shellindex,m);
//     
//   }
//  
//   GammaStructure Jcontr;
//   FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
//   FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
//   FourVector<double> pi=pf-q;
//   pi[0]=sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pi[1]*pi[1]+pi[2]*pi[2]+pi[3]*pi[3]);
//   pm=TVector3(-q[1],0.,pf[3]-q[3]);
//   if(photonpol==0)  Jcontr = J->getCC(current, q, pi, pf)*polVector0;
//   /* else if(photonpol==1)  Jcontr = J.getCC2(q)*polVectorX;
//   else if(photonpol==2)  Jcontr = J.getCC2(q)*polVectorY;
//   else if(photonpol==3)  Jcontr = J.getCC2(q)*polVectorZ;*/
//   else if(photonpol==-1) Jcontr= J->getCC(current, q, pi, pf)*polVectorMin;
//   else if(photonpol==1) Jcontr=J->getCC(current, q, pi, pf)*polVectorPlus;
// //   else if(photonpol==2) Jcontr=J->getCC(current, q, pi, pf)*polVectorX;
//   else{ cerr << "invalid photon pol" << endl;  exit(1); }
//   if(spinout==-1) barcontract = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr;
//   else if(spinout==1) barcontract = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr;
//   else{ cerr << "invalid nucl pol" << endl;  exit(1); }
//   
//   complex<double> results[2];
//   double restimate=0.,thestimate=0.,phiestimate=0.;
//   rombergerN(this,&Model::intJR,0.,pnucl->getRange(),2,results,getPrec(),3,8,&restimate,pw,SRC,&thestimate, &phiestimate);
//   if(!pw) delete grid;
//   complex<double> result;
//   delete J;
//   if(CT) return results[1];
//   else return results[0];
// }


void Model::getMatrixEl(TKinematics2to2 &tk, Matrix<2,3> & matrixel, int shellindex, 
			int m, int CT, int pw, int current, int SRC, int thick){

  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=tk.GetCosthYlab();
  if(costheta>1.) costheta=1.;
  if(costheta<-1.) costheta=-1.;
  double sintheta=sqrt(1.-costheta*costheta);
  //careful!! we compute in a frame with p_f along z-axis, so have to rotate the photon polarization 4-vectors!
  FourVector<std::complex<double> > polVectorPlus(0.,
						    -1./sqrt(2.)*costheta,
						    complex<double>(0.,-1./sqrt(2.)),
						 -1/sqrt(2.)*sintheta);
  FourVector<std::complex<double> > polVectorMin(0.,
						   1./sqrt(2.)*costheta,
						   complex<double>(0.,-1./sqrt(2.)),
						   1./sqrt(2.)*sintheta);
 FourVector<std::complex<double> >  polVector0(1.,0.,0.,0.);
 
  shell=shellindex;
  mm=m;
  int prot = (shellindex < pnucl->getPLevels()?0:1);
  J= new NucleonEMOperator(tk.GetQsquared(),!prot,0);
  FastParticle proton(prot, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()/1.e06,0.,homedir);
  if(getUsersigma()) proton.setScreening(getSigmascreening());
  if(!pw){
    if(SRC||thick) grid = &gridthick;
    else grid = &onegrid;
    grid->clearParticles();
    grid->addParticle(proton);
    grid->updateGrids();
    grid->clearKnockout();
    grid->addKnockout(shellindex,m);    
  }
  GammaStructure Jcontr0, Jcontrmin, Jcontrplus;
  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  FourVector<double> pi=pf-q;
  pi[0]=sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pi[1]*pi[1]+pi[2]*pi[2]+pi[3]*pi[3]);
  pm=TVector3(-q[1],0.,pf[3]-q[3]);
  Jcontr0 = J->getCC(current, q, pi, pf, 0.,0,*pnucl)*polVector0;
  Jcontrmin= J->getCC(current, q, pi, pf, 0.,0,*pnucl)*polVectorMin;
  Jcontrplus=J->getCC(current, q, pi, pf, 0.,0,*pnucl)*polVectorPlus;
  barcontract0down = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr0;
  barcontract0up = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr0;
  barcontractmindown = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontrmin;
  barcontractminup = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontrmin;
  barcontractplusdown = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontrplus;
  barcontractplusup = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontrplus;
  
  int res=90;
  unsigned count=0;
  if(integrator==0){
//     complex<double> results[12];
//     double restimate=0.,thestimate=0.,phiestimate=0.;
//     rombergerN(this,&Model::intJR12,0.,pnucl->getRange(),12,results,getPrec(),3,8,&restimate,pw,SRC,&thestimate, &phiestimate);
//     //if(!pw) delete grid;
//     if(CT){
//       for(int j=0;j<2;++j){
// 	for(int i=0;i<3;i++) matrixel(j,i)=results[j*6+2*i+1];
//       }
//     }
//     else{
//       for(int j=0;j<2;++j){
// 	for(int i=0;i<3;i++) matrixel(j,i)=results[j*6+2*i];
//       }
//     }
  }
  else if(integrator==1||integrator==2){

    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{pnucl->getRange(),1.,2.*PI}};
    Model::Ftor_one F;
    F.model = this;
    F.SRC = SRC;
    F.pw = pw;
    numint::mdfunction<numint::vector_z,3> mdf;
    mdf.func = &Ftor_one::exec;
    mdf.param = &F;
    numint::vector_z ret(12,0.);
    F.f=Model::klaas_one_amp;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,1E03,maxEval,ret,count,0);
    if(CT){
      for(int j=0;j<2;++j){
	for(int i=0;i<3;i++) matrixel(j,i)=ret[j*6+2*i+1];
      }
    }
    else{
      for(int j=0;j<2;++j){
	for(int i=0;i<3;i++) matrixel(j,i)=ret[j*6+2*i];
      }
    }
  }      
  else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
      
//   cout << shellindex << " " << m << " ";
//   for(int i=0;i<6;i++) cout << matrixel(i/3,i%3) << " ";
//   cout << res << " " << count << endl << endl;
  
  delete J;
}

void Model::getAllMatrixElMult(TKinematics2to2 &tk, Matrix<2,3> *matrixel, int shellindex, int m, 
			       int current, int thick){
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=tk.GetCosthYlab();
  if(costheta>1.) costheta=1.;
  if(costheta<-1.) costheta=-1.;
  double sintheta=sqrt(1.-costheta*costheta);
  //careful!! we compute in a frame with p_f along z-axis, so have to rotate the photon polarization 4-vectors!
  shell=shellindex;
  mm=m;
  int prot = (shellindex < pnucl->getPLevels()?0:1);
  J= new NucleonEMOperator(tk.GetQsquared(),!prot,0);
  FastParticle proton(prot, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()/1.e06,0.,homedir);
  if(getUsersigma()) proton.setScreening(getSigmascreening());
  FourVector<std::complex<double> > polVectorPlus = FourVector< complex<double> >(0.,
						    -1./sqrt(2.)*costheta,
						    complex<double>(0.,-1./sqrt(2.)),
						 -1/sqrt(2.)*sintheta);
  FourVector<std::complex<double> > polVectorMin = FourVector< complex<double> >(0.,
						   1./sqrt(2.)*costheta,
						   complex<double>(0.,-1./sqrt(2.)),
						   1./sqrt(2.)*sintheta);
  FourVector<std::complex<double> > polVector0 = FourVector< complex<double> >(1.,0.,0.,0.);

  if(thick) grid = &gridthick;
  else grid= &onegrid;
  grid->clearParticles();
  grid->addParticle(proton);
  grid->updateGrids();
  grid->clearKnockout();
  grid->addKnockout(shellindex,m);    

  GammaStructure Jcontr;
  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  FourVector<double> pi=pf-q;
  pi[0]=sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pi[1]*pi[1]+pi[2]*pi[2]+pi[3]*pi[3]);
  pm=TVector3(-q[1],0.,pf[3]-q[3]);
  for(int spinout=0;spinout<2;++spinout){
    for(int photopol=0;photopol<3;++photopol){
      if(photopol==0) Jcontr = J->getCC(current, q, pi, pf, 0.,0,*pnucl)*polVector0;
      if(photopol==1) Jcontr= J->getCC(current, q, pi, pf, 0.,0,*pnucl)*polVectorMin;
      if(photopol==2) Jcontr=J->getCC(current, q, pi, pf, 0.,0,*pnucl)*polVectorPlus;
      //spinors with quantization axis along z-axis correspond to the helicity spinors here!!!
      if(spinout==0) barcontract = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr;
      else if(spinout==1) barcontract = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr;
      int res=90;
      unsigned count=0;
      int total=thick?5:3; //in thickness we need 5 diff FSI results, otherwise 3
      if(integrator==0){
// 	complex<double> results[total];
// 	double restimate=0.,thestimate=0.,phiestimate=0.;
// 	rombergerN(this,&Model::intJR,0.,pnucl->getRange(),total,results,getPrec(),3,8,
// 		   &restimate,total,&thestimate, &phiestimate);
// 	for(int k=0;k<total;++k) matrixel[k](spinout,photopol) = results[k];
// 	
      }
      else if(integrator==1||integrator==2){

	numint::array<double,3> lower = {{0.,-1.,0.}};
	numint::array<double,3> upper = {{pnucl->getRange(),1.,2.*PI}};
	Model::Ftor_one F;
	F.model = this;
	F.SRC = thick;
	F.pw = total;
	numint::mdfunction<numint::vector_z,3> mdf;
	mdf.func = &Ftor_one::exec;
	mdf.param = &F;
	numint::vector_z ret(total,0.);
	F.f=Model::klaas_mult_amp;
	if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
	else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,1E03,maxEval,ret,count,0);
	for(int k=0;k<total;++k) matrixel[k](spinout,photopol) = ret[k];
      }      
      else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
//       cout << shellindex << " " << m << " " << spinout << " " << photopol << " " << matrixel[0](spinout,photopol) << 
// 	" " << matrixel[1](spinout,photopol) << " " << matrixel[total-1](spinout,photopol) << " " << res << " " << count << endl;
    }
  }
        
  delete J;
  
}


void Model::getAllMatrixEl(TKinematics2to2 &tk, Matrix<2,3> *matrixel, int shellindex, int m, int current, int thick, int medium){
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=tk.GetCosthYlab();
  if(costheta>1.) costheta=1.;
  if(costheta<-1.) costheta=-1.;
  double sintheta=sqrt(1.-costheta*costheta);
  //careful!! we compute in a frame with p_f along z-axis, so have to rotate the photon polarization 4-vectors!
 
  shell=shellindex;
  mm=m;
  int prot = (shellindex < pnucl->getPLevels()?0:1);
  J= new NucleonEMOperator(tk.GetQsquared(),!prot,0);
  FastParticle proton(prot, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()*1.E-06,0.,homedir);
  if(getUsersigma()) proton.setScreening(getSigmascreening());
  
  if(thick) grid = &gridthick;
  else grid= &onegrid;
  grid->clearParticles();
  grid->addParticle(proton);
  grid->updateGrids();
  grid->clearKnockout();
  grid->addKnockout(shellindex,m);    

  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  FourVector<double> pi=pf-q;
  pi[0]=sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pi[1]*pi[1]+pi[2]*pi[2]+pi[3]*pi[3]);
  pm=TVector3(-q[1],0.,pf[3]-q[3]);
  FourVector<std::complex<double> > polVectorPlus = FourVector< complex<double> >(0.,
						    -1./sqrt(2.)*costheta,
						    complex<double>(0.,-1./sqrt(2.)),
						 -1/sqrt(2.)*sintheta);
  FourVector<std::complex<double> > polVectorMin = FourVector< complex<double> >(0.,
						   1./sqrt(2.)*costheta,
						   complex<double>(0.,-1./sqrt(2.)),
						   1./sqrt(2.)*sintheta);
  FourVector<std::complex<double> > polVector0 = FourVector< complex<double> >(1.,0.,0.,0.);
  if(!medium){
    GammaStructure Jcontr0, Jcontrmin, Jcontrplus;
    Jcontr0 = J->getCC(current, q, pi, pf, 0.,0,*pnucl)*polVector0;
    Jcontrmin= J->getCC(current, q, pi, pf, 0.,0,*pnucl)*polVectorMin;
    Jcontrplus=J->getCC(current, q, pi, pf, 0.,0,*pnucl)*polVectorPlus;
    Matrix<1,4> spindown=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
    Matrix<1,4> spinup=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity));
    barcontract0down = spindown*Jcontr0;
    barcontract0up = spinup*Jcontr0;
    barcontractmindown = spindown*Jcontrmin;
    barcontractminup = spinup*Jcontrmin;
    barcontractplusdown = spindown*Jcontrplus;
    barcontractplusup = spinup*Jcontrplus;
    int res=90;
    unsigned count=0;
  //   if(integrator==0){
  //     complex<double> results[12];
  //     double restimate=0.,thestimate=0.,phiestimate=0.;
  //     rombergerN(this,&Model::intJR12,0.,pnucl->getRange(),12,results,getPrec(),3,8,&restimate,pw,SRC,&thestimate, &phiestimate);
  //     //if(!pw) delete grid;
  //     if(CT){
  //       for(int j=0;j<2;++j){
  // 	for(int i=0;i<3;i++) matrixel(j,i)=results[j*6+2*i+1];
  //       }
  //     }
  //     else{
  //       for(int j=0;j<2;++j){
  // 	for(int i=0;i<3;i++) matrixel(j,i)=results[j*6+2*i];
  //       }
  //     }
  //   }
    int total=thick?5:3; //in thickness we need 5 diff FSI results, otherwise 3
    
    if(integrator==1||integrator==2){

      numint::array<double,3> lower = {{0.,-1.,0.}};
      numint::array<double,3> upper = {{pnucl->getRange(),1.,2.*PI}};
      Model::Ftor_one F;
      F.model = this;
      F.SRC = thick;
      F.pw = total;
      numint::mdfunction<numint::vector_z,3> mdf;
      mdf.func = &Ftor_one::exec;
      mdf.param = &F;
      numint::vector_z ret(total*6,0.);
      F.f=Model::klaas_all_amp;

  //     numint::array<double,1> lower = {{0.,}};
  //     numint::array<double,1> upper = {{pnucl->getRange()}};
  //     Model::Ftor_r F;
  //     F.model = this;
  //     F.thick = thick;
  //     F.total = total;
  //     numint::mdfunction<numint::vector_z,1> mdf;
  //     mdf.func = &Ftor_r::exec;
  //     mdf.param = &F;
  //     numint::vector_z ret(total*6,0.);
  //     F.f=Model::klaas_all_radial;

      if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
      else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,1E03,maxEval,ret,count,0);
      for(int k=0;k<total;++k){
	for(int j=0;j<2;++j){
	  for(int i=0;i<3;i++) matrixel[k](j,i)=ret[j*total*3+total*i+k];
	}
      }
    }      
    else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
	
  //   cout << shellindex << " " << m << " ";
  //   for(int i=0;i<6;i++) cout << matrixel[0](i/3,i%3) << " ";
  //   cout << res << " " << count << endl;
  //   for(int i=0;i<6;i++) cout << matrixel[1](i/3,i%3) << " ";
  //   cout << res << " " << count << endl ;
  //    for(int i=0;i<6;i++) cout << matrixel[total-1](i/3,i%3) << " ";
  //   cout << res << " " << count << endl << endl;
  }
  else{
    
   int res=90;
    unsigned count=0;
    int total=thick?5:3; //in thickness we need 5 diff FSI results, otherwise 3    
    if(integrator==1||integrator==2){

      numint::array<double,3> lower = {{0.,-1.,0.}};
      numint::array<double,3> upper = {{pnucl->getRange(),1.,2.*PI}};
      Model::Ftor_medium F;
      F.model = this;
      F.thick = thick;
      F.total = total;
      F.spindown=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),
				      TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
      F.spinup=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),
				    TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity));
      F.medium=medium;
      F.current=current;
      F.q=q;
      F.pi=pi;
      F.pf=pf;
      F.polmin=polVectorMin;
      F.polplus=polVectorPlus;
      F.pol0=polVector0;
      numint::mdfunction<numint::vector_z,3> mdf;
      mdf.func = &Ftor_medium::exec;
      mdf.param = &F;
      numint::vector_z ret(total*6,0.);
      F.f=Model::klaas_all_amp_medium;

      if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
      else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,1E03,maxEval,ret,count,0);
      for(int k=0;k<total;++k){
	for(int j=0;j<2;++j){
	  for(int i=0;i<3;i++) matrixel[k](j,i)=ret[j*total*3+total*i+k];
	}
      }
    }      
    else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
    
    
    
  }
  
  
  delete J;
  
}


// void Model::intJR(const double r, complex<double> *results, va_list ap){
//   
//   int total = va_arg(ap,int);
//   double *pthetaestimate = va_arg(ap,double*);
//   double *pphiestimate = va_arg(ap,double*);
// 
//   gridthick.setRinterp(r);
//   rombergerN(this,&Model::intJCosTheta,-1.,1.,total,results,getPrec(),3,8,pthetaestimate, r, total, pphiestimate);
//   //MFSpinor already contains a *r!!!
//   for(int i=0;i<5;++i) results[i]*=r;
//   return;
// }
// void Model::intJR12(const double r, complex<double> *results, va_list ap){
//   
//   int pw = va_arg(ap,int);
//   int SRC = va_arg(ap,int);
//   double *pthetaestimate = va_arg(ap,double*);
//   double *pphiestimate = va_arg(ap,double*);
// 
//   if(!pw) grid->setRinterp(r);
//   rombergerN(this,&Model::intJCosTheta12,-1.,1.,12,results,getPrec(),3,8,pthetaestimate, pw, SRC, r, pphiestimate);
//   //MFSpinor already contains a *r!!!
//   for(int i=0;i<12;i++) results[i]*=r;
//   return;
// }
// 
// void Model::intJCosTheta(const double costheta, complex<double> *results, va_list ap){
//   
//   double r = va_arg(ap,double);
//   int total = va_arg(ap,int);
//   double *pphiestimate = va_arg(ap,double*);
//   double sintheta=sqrt(1.-costheta*costheta);
//   gridthick.setCthinterp(costheta);
//   rombergerN(this,&Model::intJPhi,0.,2.*PI,total,results,getPrec(),3,8,pphiestimate, r, costheta, sintheta, total);
//   return;
// }
// void Model::intJCosTheta12(const double costheta, complex<double> *results, va_list ap){
//   
//   int pw = va_arg(ap,int);
//   int SRC = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double *pphiestimate = va_arg(ap,double*);
//   double sintheta=sqrt(1.-costheta*costheta);
//   if(!pw) grid->setCthinterp(costheta);
//   rombergerN(this,&Model::intJPhi12,0.,2.*PI,12,results,getPrec(),3,8,pphiestimate, pw, SRC, r, costheta, sintheta);
//   return;
// }
// 
// void Model::intJPhi(const double phi, complex<double> *results, va_list ap){
// 
//   double r = va_arg(ap,double);
//   double costheta = va_arg(ap,double);
//   double sintheta = va_arg(ap,double);
//   int total = va_arg(ap,int);
//   double cosphi,sinphi;
//   sincos(phi,&sinphi,&cosphi);
//   
//   results[total-1]= barcontract*TMFSpinor(*pnucl,shell,mm,r,costheta,phi)*
//     exp(-INVHBARC*pm*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*I_UNIT);
//   for(int i=0;i<total-1;++i) results[i]=results[total-1]*gridthick.getFsiGridN_interp1(i,phi);
//   
//   return;
// }
// 
// void Model::intJPhi12(const double phi, complex<double> *results, va_list ap){
// 
//   int pw = va_arg(ap,int);
//   int SRC = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double costheta = va_arg(ap,double);
//   double sintheta = va_arg(ap,double);
//   double cosphi,sinphi;
//   sincos(phi,&sinphi,&cosphi);
//   TMFSpinor wave(*pnucl,shell,mm,r,costheta,phi);
//   complex<double> exp_pr=exp(-INVHBARC*pm*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*I_UNIT);
//   results[1]=results[0]= barcontract0down*wave*exp_pr;
//   results[2]=results[3]= barcontractmindown*wave*exp_pr;
//   results[4]=results[5]= barcontractplusdown*wave*exp_pr;
//   results[6]=results[7]= barcontract0up*wave*exp_pr;
//   results[8]=results[9]= barcontractminup*wave*exp_pr;
//   results[10]=results[11]= barcontractplusup*wave*exp_pr;
//   if(!pw){
//     if(SRC){
//       complex<double> src=dynamic_cast<AbstractFsiGridThick *>(grid)->getFsiSrcGridFull_interp1(phi);
//       complex<double> srcct=dynamic_cast<AbstractFsiCTGridThick *>(grid)->getFsiSrcCtGridFull_interp1(phi);
//       for(int i=0;i<6;i++){
// 	results[2*i]*=src;
// 	results[2*i+1]*=srcct;
//       }      
//     }
//     else{
//       complex<double> fsi=grid->getFsiGridFull_interp1(phi);
//       complex<double> fsict=grid->getFsiCtGridFull_interp1(phi);
//       for(int i=0;i<6;i++){
// 	results[2*i]*=fsi;
// 	results[2*i+1]*=fsict;
//       }      
//     }
//   }
//   return;
// }
// 


void Model::klaas_one_amp(numint::vector_z & results, double r, double costheta, double phi, 
			  Model & model, int SRC, int pw){
  results=numint::vector_z(12,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*(model.getPnucleus()),model.getShell(),model.getM(),r,costheta,phi);
  complex<double> exp_pr=exp(-I_UNIT*INVHBARC*(model.getPm()*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  results[1]=results[0]= exp_pr*(model.getBarcontract0down()*wave);
  results[2]=results[3]= exp_pr*(model.getBarcontractmindown()*wave);
  results[4]=results[5]= exp_pr*(model.getBarcontractplusdown()*wave);
  results[6]=results[7]= exp_pr*(model.getBarcontract0up()*wave);
  results[8]=results[9]= exp_pr*(model.getBarcontractminup()*wave);
  results[10]=results[11]= exp_pr*(model.getBarcontractplusup()*wave);
  if(!pw){
    if(SRC){
      complex<double> src=dynamic_cast<AbstractFsiGridThick *>(model.getGrid())
			      ->getFsiSrcGridFull_interp3(r,costheta,phi);
      complex<double> srcct=dynamic_cast<AbstractFsiCTGridThick *>(model.getGrid())
			      ->getFsiSrcCtGridFull_interp3(r,costheta,phi);
      for(int i=0;i<6;i++){
	results[2*i]*=src;
	results[2*i+1]*=srcct;
      }      
    }
    else{
      complex<double> fsi=model.getGrid()->getFsiGridFull_interp3(r,costheta,phi);
      complex<double> fsict=model.getGrid()->getFsiCtGridFull_interp3(r,costheta,phi);
      for(int i=0;i<6;i++){
	results[2*i]*=fsi;
	results[2*i+1]*=fsict;
      }      
    }
  }
  for(int i=0;i<12;i++) results[i]*=r;
  return;

  
}
  

void Model::klaas_all_amp(numint::vector_z & results, double r, double costheta, double phi, 
			  Model & model, int thick, int total_grid){
  results=numint::vector_z(total_grid*6,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*(model.getPnucleus()),model.getShell(),model.getM(),r,costheta,phi);
  complex<double> exp_pr=r*exp(-I_UNIT*INVHBARC*(model.getPm()*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  results[0*total_grid]= exp_pr*(model.getBarcontract0down()*wave);
  results[1*total_grid]= exp_pr*(model.getBarcontractmindown()*wave);
  results[2*total_grid]= exp_pr*(model.getBarcontractplusdown()*wave);
  results[3*total_grid]= exp_pr*(model.getBarcontract0up()*wave);
  results[4*total_grid]= exp_pr*(model.getBarcontractminup()*wave);
  results[5*total_grid]= exp_pr*(model.getBarcontractplusup()*wave);
  for(int k=0;k<6;++k) for(int i=1;i<total_grid; ++i) results[k*total_grid+i] = results[k*total_grid];
//   complex<double> glauberphase[total_grid-1];
//   vector<complex<double> >phases(4,0.);
//   dynamic_cast<GlauberGridThick *>(model.getGrid())->getFsiphaseAll(phases,r,costheta,phi);
  for(int i=0;i<(total_grid-1);++i){
    complex<double> glauberphase=model.getGrid()->getFsiGridN_interp3(i,r,costheta,phi);
    for(int k=0;k<6;++k) results[k*total_grid+i]*=glauberphase;
//     for(int k=0;k<6;++k) results[k*total_grid+i]*=phases[i];
  }
  return; 
}

void Model::klaas_all_amp_medium(numint::vector_z & results, double r, double costheta, double phi, Model & model,
				   int thick, int total,Matrix<1,4> &spinup, Matrix<1,4> &spindown, const int medium,
				   const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				   const FourVector<std::complex<double> > &polmin, const FourVector<std::complex<double> > &polplus, 
				   const FourVector<std::complex<double> > &pol0,
				   int current){  
  results=numint::vector_z(total*6,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  GammaStructure Jcontr0, Jcontrmin, Jcontrplus;
  Jcontr0 = model.getJ().getCC(current, q, pi, pf, r,medium,*model.getPnucleus())*pol0;
  Jcontrmin= model.getJ().getCC(current, q, pi, pf, r,medium,*model.getPnucleus())*polmin;
  Jcontrplus=model.getJ().getCC(current, q, pi, pf, r,medium,*model.getPnucleus())*polplus;
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*(model.getPnucleus()),model.getShell(),model.getM(),r,costheta,phi);
  //includes phase space r
  complex<double> exp_pr=r*exp(-I_UNIT*INVHBARC*(model.getPm()*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  results[0*total]= exp_pr*(spindown*Jcontr0*wave);
  results[1*total]= exp_pr*(spindown*Jcontrmin*wave);
  results[2*total]= exp_pr*(spindown*Jcontrplus*wave);
  results[3*total]= exp_pr*(spinup*Jcontr0*wave);
  results[4*total]= exp_pr*(spinup*Jcontrmin*wave);
  results[5*total]= exp_pr*(spinup*Jcontrplus*wave);
  for(int k=0;k<6;++k) for(int i=1;i<total; ++i) results[k*total+i] = results[k*total];
//   complex<double> glauberphase[total-1];
//   vector<complex<double> >phases(4,0.);
//   dynamic_cast<GlauberGridThick *>(model.getGrid())->getFsiphaseAll(phases,r,costheta,phi);
  for(int i=0;i<(total-1);++i){
    complex<double> glauberphase=model.getGrid()->getFsiGridN_interp3(i,r,costheta,phi);
    for(int k=0;k<6;++k) results[k*total+i]*=glauberphase;
//     for(int k=0;k<6;++k) results[k*total+i]*=phases[i];
  }
  return; 
  
}


void Model::klaas_all_radial(numint::vector_z & results, double r,
			  Model & model, int thick, int total){
  results=numint::vector_z(total*6,0.);
  numint::array<double,2> lower = {{-1.,0.}};
  numint::array<double,2> upper = {{1.,2.*PI}};
  Model::Ftor_angles F;
  F.model = &model;
  F.thick = thick;
  F.total = total;
  F.r=r;
  numint::mdfunction<numint::vector_z,2> mdf;
  mdf.func = &Ftor_angles::exec;
  mdf.param = &F;
//   numint::vector_z ret(total*6,0.);
  F.f=Model::klaas_all_angular;
  unsigned count=0;
  if(model.getIntegrator()==1) numint::cube_romb(mdf,lower,upper,1.E-08,model.getPrec(),results,count,0);
  else numint::cube_adaptive(mdf,lower,upper,1.E-08,model.getPrec(),1E03,2E06,results,count,0);
  for(int i=0;i<total*6;++i) results[i]*=r;  
  return;
}
  
  
  

void Model::klaas_all_angular(numint::vector_z & results, double costheta, double phi, 
			  Model & model, int thick, int total, double r){
  results=numint::vector_z(total*6,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*(model.getPnucleus()),model.getShell(),model.getM(),r,costheta,phi);
  complex<double> exp_pr=exp(-I_UNIT*INVHBARC*(model.getPm()*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  results[0*total]= exp_pr*(model.getBarcontract0down()*wave);
  results[1*total]= exp_pr*(model.getBarcontractmindown()*wave);
  results[2*total]= exp_pr*(model.getBarcontractplusdown()*wave);
  results[3*total]= exp_pr*(model.getBarcontract0up()*wave);
  results[4*total]= exp_pr*(model.getBarcontractminup()*wave);
  results[5*total]= exp_pr*(model.getBarcontractplusup()*wave);
  for(int k=0;k<6;++k) for(int i=1;i<total; ++i) results[k*total+i] = results[k*total];
//   complex<double> glauberphase[total-1];
  for(int i=0;i<(total-1);++i){
    complex<double> glauberphase=model.getGrid()->getFsiGridN_interp3(i,r,costheta,phi);
    for(int k=0;k<6;++k) results[k*total+i]*=glauberphase;
  }
  return;

  
}


void Model::klaas_mult_amp(numint::vector_z & results, double r, double costheta, double phi, 
			  Model & model, int thick, int total){
  results=numint::vector_z(total,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  results[total-1]= r*exp(-I_UNIT*INVHBARC*(model.getPm()*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)))
	      *(model.getBarcontract()*TMFSpinor(*(model.getPnucleus()),model.getShell(),model.getM(),r,costheta,phi));
  for(int i=0;i<total-1;++i) results[i]=results[total-1]*model.getThickGrid()->getFsiGridN_interp3(i,r,costheta,phi);
  return;

  
}
  
