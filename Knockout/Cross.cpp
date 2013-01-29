#include "Cross.hpp"

#include <cmath>
#include <TSpinor.h>
#include <OneGlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include <TMFSpinor.hpp>

using namespace std;

Cross::Cross(TElectronKinematics &elec, MeanFieldNucleusThick *pnucleus, double precision, int integr,
	     string dir, bool user_sigma, double sigma_screening)
:homedir(dir),electron(elec),pnucl(pnucleus),reacmodel(NULL), prec(precision),integrator(integr),
usersigma(user_sigma),sigmascreening(sigma_screening){
  
}
Cross::~Cross(){
  
}

double Cross::getElCross(TKinematics2to2 &kin, int current, double phi){
  double Q2=kin.GetQsquared();
  double qvec=kin.GetKlab();
  double Q2overkk=Q2/qvec/qvec;
  double tan2=electron.GetTan2HalfAngle(kin);
  mott=(ALPHA*ALPHA*(electron.GetCosScatterAngle(kin)+1.)*pow(electron.GetBeamEnergy(kin)-kin.GetWlab(),2.))/Q2/Q2*2.;
  reacmodel=new Model(pnucl,prec,integrator,homedir,getUsersigma(),getSigmascreening());

  //electron kinematic factors
  kinfactors[0]=pow(Q2overkk,2.);
  kinfactors[1]=tan2+Q2overkk/2.;
  kinfactors[2]=-Q2overkk/2.;
  kinfactors[3]=-Q2overkk*sqrt((tan2+Q2overkk)/2.);

  //compute response[0] functions
  for(int i=0;i<6;i++) response[0][i]=0.;
  for(int spinin=-1;spinin<=1;spinin+=2){
    for(int spinout=-1;spinout<=1;spinout+=2){
      complex<double> jmin=reacmodel->getFreeMatrixEl(kin,current, spinin,spinout,-1);
      complex<double> j0=reacmodel->getFreeMatrixEl(kin,current, spinin,spinout,0);
      complex<double> jplus=reacmodel->getFreeMatrixEl(kin,current, spinin,spinout,1);
//       cout << jmin << " " << j0 << " " << jplus << endl;
      response[0][0]+=norm(j0);
      response[0][1]+=norm(jmin)+norm(jplus);
      response[0][2]+=2.*real(conj(jplus)*jmin);
      response[0][3]+=2.*real(conj(j0)*(jplus-jmin));
    }
  }
  
  double Einon=sqrt(kin.GetHyperonMass()*kin.GetHyperonMass()+kin.GetKlab()*kin.GetKlab()+kin.GetPYlab()*kin.GetPYlab()-2.*kin.GetPYlab()*kin.GetKlab()*kin.GetCosthYlab());
  response[0][0]/=2.;
  response[0][1]/=2.;
  response[0][2]/=2.;
  response[0][3]/=2.;
  delete reacmodel;
  //combine everything
  return kin.GetPYlab()*mott*kin.GetHyperonMass()*kin.GetHyperonMass()/(Einon)*(kinfactors[0]*response[0][0]
	  +kinfactors[1]*response[0][1]+kinfactors[2]*response[0][2]*cos(2.*phi)+kinfactors[3]*response[0][3]*cos(phi));
  
}


double Cross::getDiffCross(TKinematics2to2 &kin, int current, int thick, int SRC, int CT, int pw, int shellindex, double phi){
  
  //electron kinematics
  double Q2=kin.GetQsquared();
  double qvec=kin.GetKlab();
  double Q2overkk=Q2/qvec/qvec;
  double tan2=pow(tan(acos(electron.GetCosScatterAngle(kin))/2.),2.);
  //kinematical front factor
  frontfactor=kin.GetHyperonMass()*kin.GetMesonMass()*kin.GetPYlab()/(8*pow(PI,3.))
      /abs(sqrt(kin.GetMesonMass()*kin.GetMesonMass()+kin.GetPklab()+kin.GetPklab())+sqrt(kin.GetHyperonMass()*kin.GetHyperonMass()+kin.GetPYlab()*kin.GetPYlab())
	    *(1-qvec*kin.GetCosthYlab()/kin.GetPYlab()));
  mott=(ALPHA*ALPHA*(electron.GetCosScatterAngle(kin)+1.)*pow(electron.GetBeamEnergy(kin)-kin.GetWlab(),2.))/Q2/Q2*2.;
  reacmodel=new Model(pnucl,prec,integrator,homedir,getUsersigma(),getSigmascreening());
  
  kinfactors[0]=pow(Q2overkk,2.);
  kinfactors[1]=tan2+Q2overkk/2.;
  kinfactors[2]=-Q2overkk/2.;
  kinfactors[3]=-Q2overkk*sqrt((tan2+Q2overkk)/2.);
  kinfactors[4]=sqrt(tan2*(tan2+Q2overkk));
  kinfactors[5]=-1./sqrt(2)*Q2overkk*sqrt(tan2);
  
  //compute response functions
  for(int i=0;i<6;i++) response[0][i]=0.;
  for(int m=-pnucl->getJ_array()[shellindex];m<=pnucl->getJ_array()[shellindex];m+=2){
      Matrix<2,3> J;
      reacmodel->getMatrixEl(kin,J,shellindex,m,CT,pw, current, SRC, thick);
      for(int i=0;i<2;i++){
	response[0][0]+=norm(J(i,0));
	response[0][1]+=norm(J(i,1))+norm(J(i,2));
	response[0][2]+=2.*real(conj(J(i,2))*J(i,1));
	response[0][3]+=2.*real(conj(J(i,0))*(J(i,2)-J(i,1)));
	response[0][4]+=norm(J(i,1))-norm(J(i,2));
	response[0][5]+=2.*imag(J(i,0)*conj(J(i,2)-J(i,1)));
      }
  }
  double result=0.;
  //combine everything
  result=kinfactors[0]*response[0][0]+kinfactors[1]*response[0][1]+kinfactors[2]*response[0][2]*cos(2.*phi)+kinfactors[3]*response[0][3]*cos(phi);
  delete reacmodel;
  return mott*frontfactor*result/HBARC;
  
}

void  Cross::getAllDiffCross(std::vector<double> &cross, TKinematics2to2 &kin, int current, 
			     int shellindex, int thick, double phi, int maxEval, bool lab){
  
  cross=vector<double>(5,0.);
  //electron kinematics
  double Q2=kin.GetQsquared();
  double qvec=kin.GetKlab();
  double Q2overkk=Q2/qvec/qvec;
  double tan2=pow(tan(acos(electron.GetCosScatterAngle(kin))/2.),2.);
  //kinematical front factor
  frontfactor=lab? (kin.GetHyperonMass()*kin.GetMesonMass()*kin.GetPYlab()/(8*pow(PI,3.))
      /abs(sqrt(kin.GetMesonMass()*kin.GetMesonMass()+kin.GetPklab()+kin.GetPklab())+sqrt(kin.GetHyperonMass()*kin.GetHyperonMass()+kin.GetPYlab()*kin.GetPYlab())
	    *(1-qvec*kin.GetCosthYlab()/kin.GetPYlab()))) :
	    kin.GetHyperonMass()*kin.GetMesonMass()*kin.GetPkcm()/(8*pow(PI,3.)*kin.GetW());
  mott=(ALPHA*ALPHA*(electron.GetCosScatterAngle(kin)+1.)*pow(electron.GetBeamEnergy(kin)-kin.GetWlab(),2.))/Q2/Q2*2.;
  reacmodel=new Model(pnucl,prec,integrator,homedir, maxEval, getUsersigma(),getSigmascreening());
  
  kinfactors[0]=pow(Q2overkk,2.);
  kinfactors[1]=tan2+Q2overkk/2.;
  kinfactors[2]=-Q2overkk/2.;
  kinfactors[3]=-Q2overkk*sqrt((tan2+Q2overkk)/2.);
  kinfactors[4]=-1./sqrt(2)*Q2overkk*sqrt(tan2);
  kinfactors[5]=sqrt(tan2*(tan2+Q2overkk));
  
  int total=thick?5:3;
  //compute response functions
  for(int i=0;i<6;i++) for(int j=0;j<total;j++) response[j][i]=0.;
  //only half of the m values due to symmetry, careful!!! this symmetry is only valid for unpolarized cross sections!!!!!
  //R_TT does not have this symmetry!!!
  for(int m=-pnucl->getJ_array()[shellindex];m<=0/*pnucl->getJ_array()[shellindex]*/;m+=2){
    Matrix<2,3> J[total];
    reacmodel->getAllMatrixEl(kin,J,shellindex,m,current,thick,0);
//     cout << m << " " << J[total-1] << endl;
    for(int i=0;i<2;i++){
      for(int j=0;j<total;j++){
	response[j][0]+=norm(J[j](i,0));//W_L
	response[j][1]+=norm(J[j](i,1))+norm(J[j](i,2)); //W_T
	response[j][2]+=2.*real(conj(J[j](i,2))*J[j](i,1)); //W_TT
	response[j][3]+=2.*real(conj(J[j](i,0))*(J[j](i,2)-J[j](i,1))); //W_LT
	response[j][4]+=2.*imag(J[j](i,0)*conj(J[j](i,2)-J[j](i,1))); //W_LT'
	response[j][5]+=norm(J[j](i,2))-norm(J[j](i,1)); //W_T'
      }
    }
  }
  //combine everything, factor 2 because of symmetry!
  for(int j=0;j<total;j++) cross[j]=2.*(kinfactors[0]*response[j][0]+kinfactors[1]*response[j][1]
		    +kinfactors[2]*response[j][2]*cos(2.*phi)+kinfactors[3]*response[j][3]*cos(phi))*mott*frontfactor/HBARC;
  delete reacmodel;
  
}


void  Cross::getAllObs(std::vector<double> &obs, TKinematics2to2 &kin, int current, 
			     int shellindex, int thick, int medium, double phi, int maxEval, bool lab){
  
  //electron kinematics
  double Q2=kin.GetQsquared();
  double qvec=kin.GetKlab();
  double Q2overkk=Q2/qvec/qvec;
  double tan2=pow(tan(acos(electron.GetCosScatterAngle(kin))/2.),2.);
  //kinematical front factor
  frontfactor=lab? (kin.GetHyperonMass()*kin.GetMesonMass()*kin.GetPYlab()/(8*pow(PI,3.))
      /abs(sqrt(kin.GetMesonMass()*kin.GetMesonMass()+kin.GetPklab()+kin.GetPklab())+sqrt(kin.GetHyperonMass()*kin.GetHyperonMass()+kin.GetPYlab()*kin.GetPYlab())
	    *(1-qvec*kin.GetCosthYlab()/kin.GetPYlab()))) :
	    kin.GetHyperonMass()*kin.GetMesonMass()*kin.GetPkcm()/(8*pow(PI,3.)*kin.GetW());
  mott=(ALPHA*ALPHA*(electron.GetCosScatterAngle(kin)+1.)*pow(electron.GetBeamEnergy(kin)-kin.GetWlab(),2.))/Q2/Q2*2.;
  reacmodel=new Model(pnucl,prec,integrator,homedir, maxEval, getUsersigma(),getSigmascreening());
  
  kinfactors[0]=pow(Q2overkk,2.); //v_LL
  kinfactors[1]=tan2+Q2overkk/2.; //v_T
  kinfactors[2]=-Q2overkk/2.; //v_TT
  kinfactors[3]=-Q2overkk*sqrt((tan2+Q2overkk)/2.); //v_LT
  kinfactors[4]=-1./sqrt(2)*Q2overkk*sqrt(tan2); //v_LT'
  kinfactors[5]=sqrt(tan2*(tan2+Q2overkk)); //v_T'
  int total=thick?5:3;
  obs=vector<double>(total*8,0.);
  //compute response functions
  //for(int i=0;i<6;i++) for(int j=0;j<total;j++) response[j][i]=0.;
  Matrix<2,2> responsematrix[total][6];
  for(int i=0;i<total;++i) for(int j=0;j<6;++j) responsematrix[i][j]=Matrix<2,2>();
 
  //only half of the m values due to symmetry
  for(int m=-pnucl->getJ_array()[shellindex];m<=pnucl->getJ_array()[shellindex];m+=2){
      Matrix<2,3> J[total];
      reacmodel->getAllMatrixEl(kin,J,shellindex,m,current,thick,medium);
      for(int j=0;j<total;j++){
	  responsematrix[j][0]+=Matrix<2,2>(norm(J[j](1,0)),J[j](1,0)*conj(J[j](0,0)),
					    J[j](0,0)*conj(J[j](1,0)),norm(J[j](1,0))); //W_L
	  responsematrix[j][1]+=Matrix<2,2>(norm(J[j](1,1))+norm(J[j](1,2)),J[j](1,1)*conj(J[j](0,1))+J[j](1,2)*conj(J[j](0,2)),
						J[j](0,1)*conj(J[j](1,1))+J[j](0,2)*conj(J[j](1,2)),norm(J[j](0,1))+norm(J[j](0,2)));//W_T
	  responsematrix[j][2]+=Matrix<2,2>(J[j](1,2)*conj(J[j](1,1)),J[j](1,2)*conj(J[j](0,1)),
					    J[j](0,2)*conj(J[j](1,1)),J[j](0,2)*conj(J[j](0,1))); //W_TT
	  responsematrix[j][3]+=Matrix<2,2>(J[j](1,0)*conj(J[j](1,1)),J[j](1,0)*conj(J[j](0,1)),
					    J[j](0,0)*conj(J[j](1,1)),J[j](0,0)*conj(J[j](0,1))); //part of W_LT & W_LT' (w_0,-1)
	  responsematrix[j][4]+=Matrix<2,2>(J[j](1,0)*conj(J[j](1,2)),J[j](1,0)*conj(J[j](0,2)),
					    J[j](0,0)*conj(J[j](1,2)),J[j](0,0)*conj(J[j](0,2))); //part of W_LT & W_LT' (w_0,1)
	  responsematrix[j][5]+=Matrix<2,2>(norm(J[j](1,2))-norm(J[j](1,1)),J[j](1,2)*conj(J[j](0,2))-J[j](1,1)*conj(J[j](0,1)),
						J[j](0,2)*conj(J[j](1,2))-J[j](0,1)*conj(J[j](1,1)),norm(J[j](0,2))-norm(J[j](0,1)));//W_T'
	
      }
  }
  complex<double> responsespin[total][6][4];
  for(int i=0;i<total;++i) 
    for(int j=0;j<6;++j) 
      for(int k=0;k<4;++k) responsespin[i][j][k] = Trace(responsematrix[i][j]*TSpinor::kSigmaPauli[k]);
      
  for(int i=0;i<total;i++){
    obs[8*i] = real(kinfactors[0]*responsespin[i][0][0]+kinfactors[1]*responsespin[i][1][0])+
		   2.*kinfactors[2]*real(responsespin[i][2][0])*cos(2.*phi)
		   +2.*kinfactors[3]*real(responsespin[i][4][0]-responsespin[i][3][0])*cos(phi); //total cross section
    obs[8*i+1] = (-2.*kinfactors[4]*imag(responsespin[i][4][0]-responsespin[i][3][0])*sin(phi))/obs[8*i]; //A
    obs[8*i+2] = (2.*kinfactors[2]*imag(responsespin[i][2][1])*sin(2.*phi)
		   -2.*kinfactors[3]*imag(responsespin[i][4][1]+responsespin[i][3][1])*sin(phi))/obs[8*i]; //Px
    obs[8*i+3] = (real(kinfactors[0]*responsespin[i][0][2]+kinfactors[1]*responsespin[i][1][2])+
		  2.*kinfactors[2]*real(responsespin[i][2][2])*cos(2.*phi)
		   +2.*kinfactors[3]*real(responsespin[i][4][2]-responsespin[i][3][2])*cos(phi))/obs[8*i]; //Py
    obs[8*i+4] = (2.*kinfactors[2]*imag(responsespin[i][2][3])*sin(2.*phi)
		   -2.*kinfactors[3]*imag(responsespin[i][4][3]+responsespin[i][3][3])*sin(phi))/obs[8*i]; //Pz
    obs[8*i+5] = (2.*kinfactors[4]*real(responsespin[i][4][1]+responsespin[i][3][1])*cos(phi)
		  +kinfactors[5]*real(responsespin[i][5][1]))/obs[8*i]; //P'x
    obs[8*i+6] = (-2.*kinfactors[4]*imag(responsespin[i][4][2]-responsespin[i][3][2])*sin(phi))/obs[8*i]; //P'y
    obs[8*i+7] = (2.*kinfactors[4]*real(responsespin[i][4][3]+responsespin[i][3][3])*cos(phi)
		  +kinfactors[5]*real(responsespin[i][5][3]))/obs[8*i]; //P'z
    obs[8*i]*=mott*frontfactor/HBARC;
  }
  //combine everything, factor 2 because of symmetry!
//   for(int j=0;j<total;j++) cross[j]=2.*(kinfactors[0]*response[j][0]+kinfactors[1]*response[j][1]
// 		    +kinfactors[2]*response[j][2]*cos(2.*phi)+kinfactors[3]*response[j][3]*cos(phi))*mott*frontfactor/HBARC;
  delete reacmodel;
}

			       
void Cross::getDensr(std::vector<double> &densr, const TKinematics2to2 &tk, const int shellindex, 
		const int thick, const double r, const int maxEval){
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=tk.GetCosthYlab();
  if(costheta>1.) costheta=1.;
  if(costheta<-1.) costheta=-1.;
  double sintheta=sqrt(1.-costheta*costheta);
  //careful!! we compute in a frame with p_f along z-axis
  
  FastParticle proton(0, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()*1.E-06,0.,homedir);
  AbstractFsiGrid *grid; 
  if(thick) grid = new GlauberGridThick(120,36,5,pnucl,prec,2,homedir);
  else grid= new OneGlauberGrid(120,36,pnucl,prec,2,homedir);
  grid->clearParticles();
  grid->addParticle(proton);
  grid->updateGrids();

  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  TVector3 pm(-q[1],0.,pf[3]-q[3]);
  Matrix<1,4> spinoroutup, spinoroutdown;
  spinoroutdown=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
  spinoroutup=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity)); 

  int res=90;
  unsigned count=0;
  int total=thick?5:3; //in thickness we need 5 diff FSI results, otherwise 3
  complex<double> ** phid = new complex<double>* [2*abs(getPnucl()->getKappas()[shellindex])];
  for(int i=0;i<2*abs(getPnucl()->getKappas()[shellindex]);i++) phid[i]=new complex<double>[total];
  for(int i=0;i<2*abs(getPnucl()->getKappas()[shellindex]);i++){
    vector<complex<double> >phidd(total,0.);
    getPhid(phidd,tk,shellindex,2*(i/2)+1,2*(i%2)-1,thick,maxEval);
    for(int j=0;j<total;j++) {phid[i][j]=phidd[j];}
  }
  if(integrator==1||integrator==2){

    numint::array<double,2> lower = {{-1.,0.}};
    numint::array<double,2> upper = {{1.,2.*PI}};
    Cross::Ftor_densr F;
    F.cross = this;
    F.total = total;
    F.grid = grid;
    F.pm = pm;
    F.spinorup = spinoroutup;
    F.spinordown = spinoroutdown;
    F.r = r;
    F.phid=phid;
    F.shell = shellindex;
    numint::mdfunction<numint::vector_d,2> mdf;
    mdf.func = &Ftor_densr::exec;
    mdf.param = &F;
    densr = vector<double>(total,0.);
    F.f=Cross::klaas_densr;

    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,densr,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,maxEval,densr,count,0);
  }      
  else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
  
  delete grid;
  for(int i=0;i<2*abs(getPnucl()->getKappas()[shellindex]);i++) delete [] phid[i];
  delete [] phid;
  //factor of 2 because of symmetry in m_j
  for(int i=0;i<total;i++) densr[i]*=2./pow(2.*PI,3.);
  
}

void Cross::klaas_densr(numint::vector_d & results, double costheta, double phi, Cross & cross, int total,
	      AbstractFsiGrid *grid, TVector3 &pm, Matrix<1,4> &spinordown, Matrix<1,4> &spinorup,
	      std::complex<double>  **phid, double r, int shellindex){
  results=numint::vector_d(total,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  //includes phase space!
  complex<double> exp_pr=r*exp(-I_UNIT*INVHBARC*(pm*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  for(int m=1;m<=cross.getPnucl()->getJ_array()[shellindex];m+=2){
    TMFSpinor wave(*(cross.getPnucl()),shellindex,m,r,costheta,phi);
    grid->clearKnockout();
    grid->addKnockout(shellindex,m);    
    for(int i=0;i<total;i++){
      complex<double> gl = ( i<(total-1)? grid->getFsiGridN_interp3(i,r,costheta,phi) : 1.);
      complex<double> temp1 = exp_pr*(spinordown*wave)*gl;
      complex<double> temp2 = exp_pr*(spinorup*wave)*gl;
      results[i]+= real(temp1*conj(phid[m/2][i])+conj(temp1)*phid[m/2][i]
      +temp2*conj(phid[m/2+1][i])+conj(temp2)*phid[m/2+1][i])*0.5;
    }
  }
}

void Cross::getDensr_ctheta(std::vector<double> &densr, const TKinematics2to2 &tk, const int shellindex, 
		const int thick, const double r, const double denscostheta, const int maxEval){
  
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=tk.GetCosthYlab();
  if(costheta>1.) costheta=1.;
  if(costheta<-1.) costheta=-1.;
  double sintheta=sqrt(1.-costheta*costheta);
  //careful!! we compute in a frame with p_f along z-axis
  
  FastParticle proton(0, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()*1.E-06,0.,homedir);
  AbstractFsiGrid *grid; 
  if(thick) grid = new GlauberGridThick(120,36,5,pnucl,prec,2,homedir);
  else grid= new OneGlauberGrid(120,36,pnucl,prec,2,homedir);
  grid->clearParticles();
  grid->addParticle(proton);
  grid->updateGrids();

  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  TVector3 pm(-q[1],0.,pf[3]-q[3]);
  Matrix<1,4> spinoroutup, spinoroutdown;
  spinoroutdown=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
  spinoroutup=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity)); 

  int res=90;
  unsigned count=0;
  int total=thick?5:3; //in thickness we need 5 diff FSI results, otherwise 3
  complex<double> ** phid = new complex<double>* [2*abs(getPnucl()->getKappas()[shellindex])];
  for(int i=0;i<2*abs(getPnucl()->getKappas()[shellindex]);i++) phid[i]=new complex<double>[total];
  for(int i=0;i<2*abs(getPnucl()->getKappas()[shellindex]);i++){
    vector<complex<double> >phidd(total,0.);
    getPhid(phidd,tk,shellindex,i/2,2*(i%2)-1,thick,maxEval);
    for(int j=0;j<total;j++) phid[i][j]=phidd[j];
  }
  if(integrator==1||integrator==2){

    numint::array<double,1> lower = {{0.}};
    numint::array<double,1> upper = {{2.*PI}};
    Cross::Ftor_densr_ctheta F;
    F.cross = this;
    F.total = total;
    F.grid = grid;
    F.pm = pm;
    F.spinorup = spinoroutup;
    F.spinordown = spinoroutdown;
    F.r = r;
    F.ctheta=denscostheta;
    F.phid=phid;
    F.shell = shellindex;
    numint::mdfunction<numint::vector_d,1> mdf;
    mdf.func = &Ftor_densr_ctheta::exec;
    mdf.param = &F;
    densr = vector<double>(total,0.);
    F.f=Cross::klaas_densr_ctheta;

    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,densr,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,maxEval,densr,count,0);
  }      
  else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
  
  delete grid;
  for(int i=0;i<2*abs(getPnucl()->getKappas()[shellindex]);i++) delete [] phid[i];
  delete [] phid;
  //factor of 2 due to symmetry in m_j
  for(int i=0;i<total;i++) densr[i]*=2./pow(2.*PI,3.);
  
}

void Cross::klaas_densr_ctheta(numint::vector_d & results, double phi, Cross & cross, int total,
	      AbstractFsiGrid *grid, TVector3 &pm, Matrix<1,4> &spinordown, Matrix<1,4> &spinorup,
	      std::complex<double>  **phid, double r, double costheta, int shellindex){
  results=numint::vector_d(total,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  //includes phase space!
  complex<double> exp_pr=r*sintheta*exp(-I_UNIT*INVHBARC*(pm*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  for(int m=1;m<=cross.getPnucl()->getJ_array()[shellindex];m+=2){
    TMFSpinor wave(*(cross.getPnucl()),shellindex,m,r,costheta,phi);
    grid->clearKnockout();
    grid->addKnockout(shellindex,m);    
    for(int i=0;i<total;i++){
      complex<double> gl = ( i<(total-1)? grid->getFsiGridN_interp3(i,r,costheta,phi) : 1.);
      complex<double> temp1 = exp_pr*(spinordown*wave)*gl;
      complex<double> temp2 = exp_pr*(spinorup*wave)*gl;
      results[i]+= real(temp1*conj(phid[m/2][i])+conj(temp1)*phid[m/2][i]
      +temp2*conj(phid[m/2+1][i])+conj(temp2)*phid[m/2+1][i])*0.5;
    }
  }
}

void Cross::getPhid(std::vector<complex<double> > &phid, const TKinematics2to2 &tk, const int shellindex, const int m,
		const int ms, const int thick, const int maxEval){
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  double costheta=tk.GetCosthYlab();
  if(costheta>1.) costheta=1.;
  if(costheta<-1.) costheta=-1.;
  double sintheta=sqrt(1.-costheta*costheta);
  //careful!! we compute in a frame with p_f along z-axis
  
  FastParticle proton(0, 0, tk.GetPYlab(),0.,0.,tk.GetQsquared()*1.E-06,0.,homedir);
  AbstractFsiGrid *grid; 
  if(thick) grid = new GlauberGridThick(120,36,5,pnucl,prec,2,homedir);
  else grid= new OneGlauberGrid(120,36,pnucl,prec,2,homedir);
  grid->clearParticles();
  grid->addParticle(proton);
  grid->updateGrids();
  grid->clearKnockout();
  grid->addKnockout(shellindex,m);    

  FourVector<double> q(tk.GetWlab(),-tk.GetKlab()*sintheta,0.,tk.GetKlab()*costheta);
  FourVector<double> pf(sqrt(tk.GetHyperonMass()*tk.GetHyperonMass()+tk.GetPYlab()*tk.GetPYlab()),0.,0.,tk.GetPYlab());
  TVector3 pm(-q[1],0.,pf[3]-q[3]);
  Matrix<1,4> spinorout;
  if(ms==-1) spinorout=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
  else spinorout=TSpinor::Bar(TSpinor(pf,tk.GetHyperonMass(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity)); 

  int res=90;
  unsigned count=0;
  int total=thick?5:3; //in thickness we need 5 diff FSI results, otherwise 3
  if(integrator==1||integrator==2){

    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{pnucl->getRange(),1.,2.*PI}};
    Cross::Ftor_phid F;
    F.cross = this;
    F.total = total;
    F.grid = grid;
    F.pm = pm;
    F.spinor = spinorout;
    F.shell = shellindex;
    F.m = m;
    numint::mdfunction<numint::vector_z,3> mdf;
    mdf.func = &Ftor_phid::exec;
    mdf.param = &F;
    phid = vector<complex<double> >(total,0.);
    F.f=Cross::klaas_phid;

    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,phid,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,maxEval,phid,count,0);
  }      
  else {cerr  << "integrator type not implemented " << integrator << endl; exit(1);}
  
  delete grid;

}

void Cross::klaas_phid(numint::vector_z & results, double r, double costheta, double phi, Cross & cross, int total,
	      AbstractFsiGrid *grid, TVector3 &pm, Matrix<1,4> &spinor, int shell, int m){
  results=numint::vector_z(total,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*(cross.getPnucl()),shell,m,r,costheta,phi);
  //includes phase space!
  complex<double> exp_pr=r*exp(-I_UNIT*INVHBARC*(pm*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  results[0]= exp_pr*(spinor*wave);
  for(int i=1;i<total; ++i) results[i] = results[0];
  for(int i=0;i<(total-1);++i) results[i]*=grid->getFsiGridN_interp3(i,r,costheta,phi);
  

}


void Cross::printDensity_profile(const TKinematics2to2 &kin, const int shellindex, 
		const int thick, const int maxEval){
  vector<double > densr;
  double total=0.,totalpw=0.,avg_dens=0.,avg_denspw=0., avg_r=0., avg_rpw=0.;
  
  for(int k=0;k<100;k++){
    double r=pnucl->getRange()/100.*k;
    double dens=pnucl->getTotalDensity(r)*pnucl->getA();
//     double dens =  (pnucl->getF()[shellindex][k*12]*pnucl->getF()[shellindex][k*12]
//       +pnucl->getG()[shellindex][k*12]*pnucl->getG()[shellindex][k*12])
//       *abs(pnucl->getKappas()[shellindex])/(2.*PI);
    getDensr(densr,kin,shellindex,thick,r,maxEval);
    total+=densr[1];
    totalpw+=densr[4];
    avg_r+=densr[1]*r;
    avg_rpw+=densr[4]*r;
    if(k!=0.) avg_dens+=densr[1]*dens/r/r;
    if(k!=0.) avg_denspw+=densr[4]*dens/r/r;
    
    cout << kin.GetPklab() << " " << r << " " << densr[1] << " " << densr[4] << " " << dens
   << endl;
  }
//   double phi_pm=0.,phi_pm_pw=0.;
//   for(int m=1;m<=pnucl->getJ_array()[shellindex];m+=2){
//     for(int ms=-1;ms<=1;ms+=2){
//       vector<complex<double> > phid;
//       getPhid(phid,kin,shellindex,m,ms,thick,maxEval);
//       phi_pm+=norm(phid[1]);
//       phi_pm_pw+=norm(phid[4]);
//     }
//   }
//   //factor of 2 because of symmetry in m_j
//   phi_pm*=2./pow(2.*PI,3.);
//   phi_pm_pw*=2./pow(2.*PI,3.);
  
  cout << endl << endl;
  cout << kin.GetPklab() << " " << avg_dens/total << " " << avg_denspw/totalpw << " " 
    << avg_r/total << " " << avg_rpw/totalpw << " " 
    << total*pnucl->getRange()/100. << " " << totalpw*pnucl->getRange()/100. << endl << endl;
  cout << endl << endl;

}
