#include "Cross.hpp"

#include <cmath>

Cross::Cross(TElectronKinematics &elec, MeanFieldNucleusThick *pnucleus, double precision, string dir,
	     bool user_sigma, double sigma_screening)
:homedir(dir),electron(elec),pnucl(pnucleus),reacmodel(NULL), prec(precision),
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
  reacmodel=new Model(pnucl,0,0,getPrec(),homedir,getUsersigma(),getSigmascreening());

  //electron kinematic factors
  kinfactors[0]=pow(Q2overkk,2.);
  kinfactors[1]=tan2+Q2overkk/2.;
  kinfactors[2]=-Q2overkk/2.;
  kinfactors[3]=-Q2overkk*sqrt((tan2+Q2overkk)/2.);

  //compute response functions
  for(int i=0;i<6;i++) response[i]=0.;
  for(int spinin=-1;spinin<=1;spinin+=2){
    for(int spinout=-1;spinout<=1;spinout+=2){
      complex<double> jmin=reacmodel->getFreeMatrixEl(kin,current, spinin,spinout,-1);
      complex<double> j0=reacmodel->getFreeMatrixEl(kin,current, spinin,spinout,0);
      complex<double> jplus=reacmodel->getFreeMatrixEl(kin,current, spinin,spinout,1);
//       cout << jmin << " " << j0 << " " << jplus << endl;
      response[0]+=norm(j0);
      response[1]+=norm(jmin)+norm(jplus);
      response[2]+=2.*real(conj(jplus)*jmin);
      response[3]+=2.*real(conj(j0)*(jplus-jmin));
    }
  }
  
  double Einon=sqrt(kin.GetHyperonMass()*kin.GetHyperonMass()+kin.GetKlab()*kin.GetKlab()+kin.GetPYlab()*kin.GetPYlab()-2.*kin.GetPYlab()*kin.GetKlab()*kin.GetCosthYlab());
  response[0]/=2.;
  response[1]/=2.;
  response[2]/=2.;
  response[3]/=2.;
  delete reacmodel;
  //combine everything
  return kin.GetPYlab()*mott*kin.GetHyperonMass()*kin.GetHyperonMass()/(Einon)*(kinfactors[0]*response[0]
	  +kinfactors[1]*response[1]+kinfactors[2]*response[2]*cos(2.*phi)+kinfactors[3]*response[3]*cos(phi));
  
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
  reacmodel=new Model(pnucl,SRC,thick,prec,homedir,getUsersigma(),getSigmascreening());
  
  kinfactors[0]=pow(Q2overkk,2.);
  kinfactors[1]=tan2+Q2overkk/2.;
  kinfactors[2]=-Q2overkk/2.;
  kinfactors[3]=-Q2overkk*sqrt((tan2+Q2overkk)/2.);
  kinfactors[4]=sqrt(tan2*(tan2+Q2overkk));
  kinfactors[5]=-1./sqrt(2)*Q2overkk*sqrt(tan2);
  
  //compute response functions
  for(int i=0;i<6;i++) response[i]=0.;
  for(int m=-pnucl->getJ_array()[shellindex];m<=pnucl->getJ_array()[shellindex];m+=2){
      Matrix<2,3> J;
      reacmodel->getMatrixEl(kin,J,shellindex,m,CT,pw, current);
      for(int i=0;i<2;i++){
	response[0]+=norm(J(i,0));
	response[1]+=norm(J(i,1))+norm(J(i,2));
	response[2]+=2.*real(conj(J(i,2))*J(i,1));
	response[3]+=2.*real(conj(J(i,0))*(J(i,2)-J(i,1)));
	response[4]+=norm(J(i,2))-norm(J(i,1));
	response[5]+=2.*imag(J(i,0)*conj(J(i,2)-J(i,1)));
      }
  }
  double result=0.;
  //combine everything
  result=kinfactors[0]*response[0]+kinfactors[1]*response[1]+kinfactors[2]*response[2]*cos(2.*phi)+kinfactors[3]*response[3]*cos(phi);
  delete reacmodel;
  return mott*frontfactor*result/HBARC;
  
}

