#include "Cross.hpp"

#include <cmath>
#include <TSpinor.h>

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
    reacmodel->getAllMatrixEl(kin,J,shellindex,m,current,thick);
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
			     int shellindex, int thick, double phi, int maxEval, bool lab){
  
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
      reacmodel->getAllMatrixEl(kin,J,shellindex,m,current,thick);
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
