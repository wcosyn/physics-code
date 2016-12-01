#include "WeakQEHadronCurrent.hpp"
#include <TSpinor.h>
#include <TMFSpinor.hpp>
#include <AuxFunction.hpp>

using namespace std;



WeakQEHadronCurrent::WeakQEHadronCurrent(MeanFieldNucleusThick *pnucleus, double precision, int integr, string dir,
	      int max_Eval, bool cc, double M_A_in, bool user_sigma, double gA_s_in, double r_s2_in, 
	      double mu_s_in, double sigma_screening)
:pnucl(pnucleus), prec(precision), integrator(integr), homedir(dir), maxEval(max_Eval),
charged(cc), M_A(M_A_in), usersigma(user_sigma), gA_s(gA_s_in), r_s2(r_s2_in), mu_s(mu_s_in), sigmascreening(sigma_screening),
gridthick(GlauberGridThick(120,36,5,pnucleus,precision,2,dir)),
onegrid(OneGlauberGrid(120,36,pnucleus,precision,2,dir)){
}


WeakQEHadronCurrent::~WeakQEHadronCurrent(){
  
}



complex<double> WeakQEHadronCurrent::getFreeMatrixEl(double Q2, bool proton, 
						     int current, int spinin, int spinout, int photonpol){
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=1.;
  double sintheta=0.;
  double omega=0.5*(MASSN*MASSN-MASSP*MASSP+Q2)/MASSN;
  double qvec = sqrt(omega*omega+Q2);
  double pout = qvec;
  double Eout = sqrt(MASSN*MASSN+pout*pout);
  double pin = 0.;  // sqrt(qvec*qvec+pout*pout-2.*pout*qvec*costheta);
  double Ein = MASSP;
/*double wbar = Eout-Ein;
  double Q2bar = qvec*qvec-wbar*wbar;*/
  static FourVector<complex<double> > polVector0(1.,0.,0.,0.);
  static FourVector<complex<double> > polVector1(0.,1.,0.,0.);
  static FourVector<complex<double> > polVector2(0.,0.,1.,0.);
  static FourVector<complex<double> > polVector3(0.,0.,0.,1.);
  //Frame has z-axis along q!!!
  FourVector<double> q(omega,0.,0.,qvec);
  FourVector<double> pf(Eout,pout*sintheta,0.,pout*costheta);
  FourVector<double> Pin(Ein,pf[1],0.,pf[3]-q[3]); 
//cout << "q:" << q << endl;
//cout << "pf:" << pf << endl;
//cout << "Pin:" << Pin << endl;
  J= new NucleonWeakOperator(Q2,proton,0,charged,M_A,r_s2,mu_s,gA_s);
  GammaStructure Jcontr;
  if(photonpol==0)      Jcontr= J->getCC_weak(current, q, Pin, pf, 0.,0,*pnucl)*polVector0;
  else if(photonpol==1) Jcontr=-J->getCC_weak(current, q, Pin, pf, 0.,0,*pnucl)*polVector1;
  else if(photonpol==2) Jcontr=-J->getCC_weak(current, q, Pin, pf, 0.,0,*pnucl)*polVector2;
  else if(photonpol==3) Jcontr=-J->getCC_weak(current, q, Pin, pf, 0.,0,*pnucl)*polVector3;//minus sign to get J^z!!!
  
//   cout << "out " << Jcontr << endl;
  /*cout << "spinorout " << TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity) << endl;
  cout << "current " << Jcontr << endl;
  cout << "spinorin " << TSpinor(Pin,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity) << endl;
  */
 /* cout << "pin " << Pin << endl;
  cout << "pout " << pf << endl;
  cout << "kin " << pout << " " << qvec << " " << pin << " " << tk.GetHyperonMass() << endl;
 */ 
  delete J;
  //contract everything
  if(spinout==-1){
    if(spinin==-1) return TSpinor::Bar(TSpinor(pf,MASSN,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr*
		      TSpinor(Pin,MASSN,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity);
    else if(spinin==1) return TSpinor::Bar(TSpinor(pf,MASSN,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr*
		      TSpinor(Pin,MASSN,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity);
    else {cerr << "invalid spinin " << spinin << endl; exit(1);}
  }
  else if(spinout==1){
    if(spinin==-1) return TSpinor::Bar(TSpinor(pf,MASSN,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr*
		      TSpinor(Pin,MASSN,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity);
    else if(spinin==1) return TSpinor::Bar(TSpinor(pf,MASSN,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr*
		      TSpinor(Pin,MASSN,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity);
    else {cerr << "invalid spinin " << spinin << endl; exit(1);}
  }
  else {cerr << "invalid spinout " << spinout << endl; exit(1);}
  
  
}



void WeakQEHadronCurrent::getMatrixEl(TKinematics2to2 &tk, Matrix<2,4> & matrixel, int shellindex, 
			int m, int CT, int pw, int current, int SRC, int thick, bool vector_or_axial){

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
  FourVector<std::complex<double> >  polVectorz(0.,-sintheta,0.,costheta);
  
  
  bool proton = shellindex<pnucl->getPLevels()? 1:0; //determine if we're dealing with initial p or n
  
  shell=shellindex;
  mm=m;
  J= new NucleonWeakOperator(tk.GetQsquared(),proton,0,charged,M_A,r_s2,mu_s,gA_s);
  //in fastparticle outgoing nucleon (0=proton,1=neutron)!! changes for CC
  if(!pw){
    FastParticle nucleon(charged?proton:!proton, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()/1.e06,0.,homedir);
    if(getUsersigma()) nucleon.setScreening(getSigmascreening());
    if(SRC||thick) grid = &gridthick;
    else grid = &onegrid;
    grid->clearParticles();
    grid->addParticle(nucleon);
    grid->updateGrids();
    grid->clearKnockout();
    grid->addKnockout(shellindex,m);    
  }
  GammaStructure Jcontr0, Jcontrmin, Jcontrplus, Jcontrz;
  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  FourVector<double> pi=pf-q;
  pi[0]=sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pi[1]*pi[1]+pi[2]*pi[2]+pi[3]*pi[3]);
  pm=TVector3(-q[1],0.,pf[3]-q[3]);
  if(vector_or_axial){
    Jcontr0 = J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVector0;
    Jcontrmin= J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVectorMin;
    Jcontrplus=J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVectorPlus;
    Jcontrz=-J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVectorz; //minus sign to get J^z!!!
  }
  else{
    Jcontr0 = J->getAxial(q)*polVector0;
    Jcontrmin = J->getAxial(q)*polVectorMin;
    Jcontrplus =J->getAxial(q)*polVectorPlus;
    Jcontrz =-J->getAxial(q)*polVectorz; //minus sign to get J^z!!!
  }
  barcontract0down = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr0;
  barcontract0up = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr0;
  barcontractmindown = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontrmin;
  barcontractminup = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontrmin;
  barcontractplusdown = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontrplus;
  barcontractplusup = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontrplus;
  barcontractzdown = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontrz;
  barcontractzup = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontrz;
  
  int res=90;
  unsigned count=0;
  if(integrator==1||integrator==2){

    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{pnucl->getRange(),1.,2.*PI}};
    WeakQEHadronCurrent::Ftor_one F;
    F.model = this;
    F.SRC = SRC;
    F.pw = pw;
    numint::mdfunction<numint::vector_z,3> mdf;
    mdf.func = &Ftor_one::exec;
    mdf.param = &F;
    numint::vector_z ret(16,0.);
    F.f=WeakQEHadronCurrent::klaas_one_amp;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,1E03,maxEval,ret,count,0);
    if(CT){
      for(int j=0;j<2;++j){
	for(int i=0;i<4;i++) matrixel(j,i)=ret[j*8+2*i+1];
      }
    }
    else{
      for(int j=0;j<2;++j){
	for(int i=0;i<4;i++) matrixel(j,i)=ret[j*8+2*i];
      }
    }
  }      
  else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
      
//   cout << shellindex << " " << m << " ";
//   for(int i=0;i<6;i++) cout << matrixel(i/3,i%3) << " ";
//   cout << res << " " << count << endl << endl;
  
  delete J;
}

void WeakQEHadronCurrent::getAllMatrixElMult(TKinematics2to2 &tk, Matrix<2,4> *matrixel, int shellindex, int m, 
			       int current, int thick){
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=tk.GetCosthYlab();
  if(costheta>1.) costheta=1.;
  if(costheta<-1.) costheta=-1.;
  double sintheta=sqrt(1.-costheta*costheta);

  bool proton = shellindex<pnucl->getPLevels()? 1:0; //determine if we're dealing with initial p or n
  
  shell=shellindex;
  mm=m;
  J= new NucleonWeakOperator(tk.GetQsquared(),proton,0,charged,M_A,r_s2,mu_s,gA_s);
  FastParticle nucleon(charged?proton:!proton, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()/1.e06,0.,homedir);
  if(getUsersigma()) nucleon.setScreening(getSigmascreening());
  //careful!! we compute in a frame with p_f along z-axis, so have to rotate the photon polarization 4-vectors!
  FourVector<std::complex<double> > polVectorPlus = FourVector< complex<double> >(0.,
						    -1./sqrt(2.)*costheta,
						    complex<double>(0.,-1./sqrt(2.)),
						 -1/sqrt(2.)*sintheta);
  FourVector<std::complex<double> > polVectorMin = FourVector< complex<double> >(0.,
						   1./sqrt(2.)*costheta,
						   complex<double>(0.,-1./sqrt(2.)),
						   1./sqrt(2.)*sintheta);
  FourVector<std::complex<double> > polVector0 = FourVector< complex<double> >(1.,0.,0.,0.);
  FourVector<std::complex<double> >  polVectorz(0.,-sintheta,0.,costheta);

  if(thick) grid = &gridthick;
  else grid= &onegrid;
  grid->clearParticles();
  grid->addParticle(nucleon);
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
      if(photopol==0) Jcontr = J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVector0;
      if(photopol==1) Jcontr= J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVectorMin;
      if(photopol==2) Jcontr=J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVectorPlus;
      if(photopol==3) Jcontr=J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVectorz;
      //spinors with quantization axis along z-axis correspond to the helicity spinors here!!!
      if(spinout==0) barcontract = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity))*Jcontr;
      else if(spinout==1) barcontract = TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity))*Jcontr;
      int res=90;
      unsigned count=0;
      int total=thick?5:3; //in thickness we need 5 diff FSI results, otherwise 3
      if(integrator==1||integrator==2){

	numint::array<double,3> lower = {{0.,-1.,0.}};
	numint::array<double,3> upper = {{pnucl->getRange(),1.,2.*PI}};
	WeakQEHadronCurrent::Ftor_one F;
	F.model = this;
	F.SRC = thick;
	F.pw = total;
	numint::mdfunction<numint::vector_z,3> mdf;
	mdf.func = &Ftor_one::exec;
	mdf.param = &F;
	numint::vector_z ret(total,0.);
	F.f=WeakQEHadronCurrent::klaas_mult_amp;
	if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
	else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,1E03,2E06,ret,count,0);
	for(int k=0;k<total;++k) matrixel[k](spinout,photopol) = ret[k];
      }      
      else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
//       cout << shellindex << " " << m << " " << spinout << " " << photopol << " " << matrixel[0](spinout,photopol) << 
// 	" " << matrixel[1](spinout,photopol) << " " << matrixel[total-1](spinout,photopol) << " " << res << " " << count << endl;
    }
  }
        
  delete J;
  
}


void WeakQEHadronCurrent::getAllMatrixEl(TKinematics2to2 &tk, Matrix<2,4> *matrixel, 
					 int shellindex, int m, int current, int thick, int medium){
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=tk.GetCosthYlab();
  if(costheta>1.) costheta=1.;
  if(costheta<-1.) costheta=-1.;
  double sintheta=sqrt(1.-costheta*costheta);
 
  bool proton = shellindex<pnucl->getPLevels()? 1:0; //determine if we're dealing with initial p or n
  
  shell=shellindex;
  mm=m;
  J= new NucleonWeakOperator(tk.GetQsquared(),proton,0,charged,M_A,r_s2,mu_s,gA_s);
  FastParticle nucleon(charged?proton:!proton, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()/1.e06,0.,homedir);
  if(getUsersigma()) nucleon.setScreening(getSigmascreening());
  
  if(thick) grid = &gridthick;
  else grid= &onegrid;
  grid->clearParticles();
  grid->addParticle(nucleon);
  grid->updateGrids();
  grid->clearKnockout();
  grid->addKnockout(shellindex,m);    

  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  FourVector<double> pi=pf-q;
  pi[0]=sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+pi[1]*pi[1]+pi[2]*pi[2]+pi[3]*pi[3]);
  pm=TVector3(-q[1],0.,pf[3]-q[3]);
  //careful!! we compute in a frame with p_f along z-axis, so have to rotate the photon polarization 4-vectors!
  FourVector<std::complex<double> > polVectorPlus = FourVector< complex<double> >(0.,
						    -1./sqrt(2.)*costheta,
						    complex<double>(0.,-1./sqrt(2.)),
						 -1/sqrt(2.)*sintheta);
  FourVector<std::complex<double> > polVectorMin = FourVector< complex<double> >(0.,
						   1./sqrt(2.)*costheta,
						   complex<double>(0.,-1./sqrt(2.)),
						   1./sqrt(2.)*sintheta);
  FourVector<std::complex<double> > polVector0 = FourVector< complex<double> >(1.,0.,0.,0.);
  FourVector<std::complex<double> >  polVectorz(0.,-sintheta,0.,costheta);

  if(!medium){
    GammaStructure Jcontr0, Jcontrmin, Jcontrplus, Jcontrz;
    Jcontr0 = J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVector0;
    Jcontrmin= J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVectorMin;
    Jcontrplus=J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVectorPlus;
    Jcontrz=J->getCC_weak(current, q, pi, pf, 0.,0,*pnucl)*polVectorz;
    Matrix<1,4> spindown=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
    Matrix<1,4> spinup=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity));
    barcontract0down = spindown*Jcontr0;
    barcontract0up = spinup*Jcontr0;
    barcontractmindown = spindown*Jcontrmin;
    barcontractminup = spinup*Jcontrmin;
    barcontractplusdown = spindown*Jcontrplus;
    barcontractplusup = spinup*Jcontrplus;
    barcontractzdown = spindown*Jcontrz;
    barcontractzup = spinup*Jcontrz;
    int res=90;
    unsigned count=0;
    int total=thick?5:3; //in thickness we need 5 diff FSI results, otherwise 3
    
    if(integrator==1||integrator==2){

      numint::array<double,3> lower = {{0.,-1.,0.}};
      numint::array<double,3> upper = {{pnucl->getRange(),1.,2.*PI}};
      WeakQEHadronCurrent::Ftor_one F;
      F.model = this;
      F.SRC = thick;
      F.pw = total;
      numint::mdfunction<numint::vector_z,3> mdf;
      mdf.func = &Ftor_one::exec;
      mdf.param = &F;
      numint::vector_z ret(total*8,0.);
      F.f=WeakQEHadronCurrent::klaas_all_amp;


      if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
      else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,1E03,maxEval,ret,count,0);
      for(int k=0;k<total;++k){
	for(int j=0;j<2;++j){
	  for(int i=0;i<4;i++) matrixel[k](j,i)=ret[j*total*4+total*i+k];
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
      WeakQEHadronCurrent::Ftor_medium F;
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
      F.polz=polVectorz;
      numint::mdfunction<numint::vector_z,3> mdf;
      mdf.func = &Ftor_medium::exec;
      mdf.param = &F;
      numint::vector_z ret(total*8,0.);
      F.f=WeakQEHadronCurrent::klaas_all_amp_medium;

      if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
      else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,1E03,maxEval,ret,count,0);
      for(int k=0;k<total;++k){
	for(int j=0;j<2;++j){
	  for(int i=0;i<4;i++) matrixel[k](j,i)=ret[j*total*4+total*i+k];
	}
      }
    }      
    else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
    
    
    
  }
  
  
  delete J;
  
}




void WeakQEHadronCurrent::klaas_one_amp(numint::vector_z & results, double r, double costheta, double phi, 
			  WeakQEHadronCurrent & model, int SRC, int pw){
  results=numint::vector_z(16,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*(model.getPnucleus()),model.getShell(),model.getM(),r,costheta,phi);
  //no r squared because wave function already is multiplied by r
  complex<double> exp_pr=r*exp(-I_UNIT*INVHBARC*(model.getPm()*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  results[1]=results[0]= exp_pr*(model.getBarcontract0down()*wave);
  results[2]=results[3]= exp_pr*(model.getBarcontractmindown()*wave);
  results[4]=results[5]= exp_pr*(model.getBarcontractplusdown()*wave);
  results[6]=results[7]= exp_pr*(model.getBarcontractzdown()*wave);
  results[8]=results[9]= exp_pr*(model.getBarcontract0up()*wave);
  results[10]=results[11]= exp_pr*(model.getBarcontractminup()*wave);
  results[12]=results[13]= exp_pr*(model.getBarcontractplusup()*wave);
  results[14]=results[15]= exp_pr*(model.getBarcontractzup()*wave);
  if(!pw){
    if(SRC){
      complex<double> src=dynamic_cast<AbstractFsiGridThick *>(model.getGrid())
			      ->getFsiSrcGridFull_interp3(r,costheta,phi);
      complex<double> srcct=dynamic_cast<AbstractFsiCTGridThick *>(model.getGrid())
			      ->getFsiSrcCtGridFull_interp3(r,costheta,phi);
      for(int i=0;i<8;i++){
	results[2*i]*=src;
	results[2*i+1]*=srcct;
      }      
    }
    else{
      complex<double> fsi=model.getGrid()->getFsiGridFull_interp3(r,costheta,phi);
      complex<double> fsict=model.getGrid()->getFsiCtGridFull_interp3(r,costheta,phi);
      for(int i=0;i<8;i++){
	results[2*i]*=fsi;
	results[2*i+1]*=fsict;
      }      
    }
  }
  return;

  
}
  

void WeakQEHadronCurrent::klaas_all_amp(numint::vector_z & results, double r, double costheta, double phi, 
			  WeakQEHadronCurrent & model, int thick, int total_grid){
  results=numint::vector_z(total_grid*8,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*(model.getPnucleus()),model.getShell(),model.getM(),r,costheta,phi);
  complex<double> exp_pr=r*exp(-I_UNIT*INVHBARC*(model.getPm()*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  results[0*total_grid]= exp_pr*(model.getBarcontract0down()*wave);
  results[1*total_grid]= exp_pr*(model.getBarcontractmindown()*wave);
  results[2*total_grid]= exp_pr*(model.getBarcontractplusdown()*wave);
  results[3*total_grid]= exp_pr*(model.getBarcontractzdown()*wave);
  results[4*total_grid]= exp_pr*(model.getBarcontract0up()*wave);
  results[5*total_grid]= exp_pr*(model.getBarcontractminup()*wave);
  results[6*total_grid]= exp_pr*(model.getBarcontractplusup()*wave);
  results[7*total_grid]= exp_pr*(model.getBarcontractzup()*wave);
  for(int k=0;k<8;++k) for(int i=1;i<total_grid; ++i) results[k*total_grid+i] = results[k*total_grid];
//   complex<double> glauberphase[total_grid-1];
//   vector<complex<double> >phases(4,0.);
//   dynamic_cast<GlauberGridThick *>(model.getGrid())->getFsiphaseAll(phases,r,costheta,phi);
  for(int i=0;i<(total_grid-1);++i){
    complex<double> glauberphase=model.getGrid()->getFsiGridN_interp3(i,r,costheta,phi);
    for(int k=0;k<8;++k) results[k*total_grid+i]*=glauberphase;
//     for(int k=0;k<6;++k) results[k*total_grid+i]*=phases[i];
  }
  return; 
}

void WeakQEHadronCurrent::klaas_all_amp_medium(numint::vector_z & results, double r, double costheta, double phi, WeakQEHadronCurrent & model,
				   int thick, int total,Matrix<1,4> &spinup, Matrix<1,4> &spindown, const int medium,
				   const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				   const FourVector<std::complex<double> > &polmin, const FourVector<std::complex<double> > &polplus, 
				   const FourVector<std::complex<double> > &pol0, const FourVector<std::complex<double> > &polz,
				   int current){  
  results=numint::vector_z(total*8,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  GammaStructure Jcontr0, Jcontrmin, Jcontrplus, Jcontrz;
  Jcontr0 = model.getJ().getCC_weak(current, q, pi, pf, r,medium,*model.getPnucleus())*pol0;
  Jcontrmin= model.getJ().getCC_weak(current, q, pi, pf, r,medium,*model.getPnucleus())*polmin;
  Jcontrplus=model.getJ().getCC_weak(current, q, pi, pf, r,medium,*model.getPnucleus())*polplus;
  Jcontrplus=model.getJ().getCC_weak(current, q, pi, pf, r,medium,*model.getPnucleus())*polz;
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*(model.getPnucleus()),model.getShell(),model.getM(),r,costheta,phi);
  //includes phase space r
  complex<double> exp_pr=r*exp(-I_UNIT*INVHBARC*(model.getPm()*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  results[0*total]= exp_pr*(spindown*Jcontr0*wave);
  results[1*total]= exp_pr*(spindown*Jcontrmin*wave);
  results[2*total]= exp_pr*(spindown*Jcontrplus*wave);
  results[3*total]= exp_pr*(spindown*Jcontrz*wave);
  results[4*total]= exp_pr*(spinup*Jcontr0*wave);
  results[5*total]= exp_pr*(spinup*Jcontrmin*wave);
  results[6*total]= exp_pr*(spinup*Jcontrplus*wave);
  results[7*total]= exp_pr*(spinup*Jcontrz*wave);
  for(int k=0;k<8;++k) for(int i=1;i<total; ++i) results[k*total+i] = results[k*total];
//   complex<double> glauberphase[total-1];
//   vector<complex<double> >phases(4,0.);
//   dynamic_cast<GlauberGridThick *>(model.getGrid())->getFsiphaseAll(phases,r,costheta,phi);
  for(int i=0;i<(total-1);++i){
    complex<double> glauberphase=model.getGrid()->getFsiGridN_interp3(i,r,costheta,phi);
    for(int k=0;k<8;++k) results[k*total+i]*=glauberphase;
//     for(int k=0;k<6;++k) results[k*total+i]*=phases[i];
  }
  return; 
  
}




void WeakQEHadronCurrent::klaas_mult_amp(numint::vector_z & results, double r, double costheta, double phi, 
			  WeakQEHadronCurrent & model, int thick, int total){
  results=numint::vector_z(total,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  results[total-1]= r*exp(-I_UNIT*INVHBARC*(model.getPm()*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)))
	      *(model.getBarcontract()*TMFSpinor(*(model.getPnucleus()),model.getShell(),model.getM(),r,costheta,phi));
  for(int i=0;i<total-1;++i) results[i]=results[total-1]*model.getThickGrid()->getFsiGridN_interp3(i,r,costheta,phi);
  return;

  
}
  
