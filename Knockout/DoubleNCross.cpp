#include "DoubleNCross.hpp"

#include <cmath>

 DoubleNCross::DoubleNCross(TElectronKinematics &elec, MeanFieldNucleusThick *pnucleus, 
			    double precision, int integr, string dir)
:homedir(dir),electron(elec),pnucl(pnucleus),reacmodel(NULL),prec(precision),integrator(integr){
  
}
 DoubleNCross::~DoubleNCross(){
  
}



double DoubleNCross::getDiffCross(const TKinematics2to3 &kin, bool SRC, bool CT, bool pw, 
				  bool corr, int shellindex1, int shellindex2,double phi){
  int particletype1=(shellindex1 < pnucl->getPLevels()? 0:1);
  int particletype2=(shellindex2 < pnucl->getPLevels()? 0:1);
  double Q2=kin.GetQsquared();
  double qvec=kin.GetKlab();
  double Q2overkk=Q2/qvec/qvec;
  double tan2=pow(tan(acos(electron.GetCosScatterAngle(kin))/2.),2.);
//   frontfactor=kin.GetHyperonMass()*kin.GetMesonMass()*kin.GetPYlab()/(8*pow(PI,3.))
//       /abs(sqrt(kin.GetMesonMass()*kin.GetMesonMass()+kin.GetPklab()+kin.GetPklab())+sqrt(kin.GetHyperonMass()*kin.GetHyperonMass()+kin.GetPYlab()*kin.GetPYlab())
// 	    *(1-qvec*kin.GetCosthYlab()/kin.GetPYlab()));
  mott=(ALPHA*ALPHA*(electron.GetCosScatterAngle(kin)+1.)*pow(electron.GetBeamEnergy(kin)-kin.GetWlab(),2.))/Q2/Q2*2.; 
  reacmodel=new DoubleNModel(pnucl,pw,SRC,CT, corr, particletype1, particletype2, prec,integrator,homedir);
  kinfactors[0]=pow(Q2overkk,2.);
  kinfactors[1]=tan2+Q2overkk/2.;
  kinfactors[2]=-Q2overkk/2.;
  kinfactors[3]=-Q2overkk*sqrt((tan2+Q2overkk)/2.);
  kinfactors[4]=sqrt(tan2*(tan2+Q2overkk));
  kinfactors[5]=-1./sqrt(2)*Q2overkk*sqrt(tan2);
  for(int i=0;i<9;i++) response[i]=0.;
  double test=0.;
  for(int m1=-pnucl->getJ_array()[shellindex1];m1<=pnucl->getJ_array()[shellindex1];m1+=2){
    int maxm2 = (shellindex2==shellindex1? m1-2 : pnucl->getJ_array()[shellindex2]);
    for(int m2=-pnucl->getJ_array()[shellindex2];m2<=maxm2;m2+=2){
      for(int spin1=-1;spin1<=1;spin1+=2){
	for(int spin2=-1;spin2<=1;spin2+=2){
	  complex<double> jmin=reacmodel->getMatrixEl(kin,spin1, spin2,-1,shellindex1,shellindex2,m1,m2);
	  complex<double> j0=reacmodel->getMatrixEl(kin,spin1, spin2, 0,shellindex1,shellindex2,m1,m2);
	  complex<double> jplus=reacmodel->getMatrixEl(kin,spin1, spin2, 1, shellindex1,shellindex2,m1,m2);
	  cout << m1 << " " << m2 << " " << spin1 << " " << spin2 << j0 << " " << jmin << " " << jplus << endl;
	  response[0]+=norm(j0);
	  response[1]+=norm(jmin)+norm(jplus);
	  response[2]+=2.*real(conj(jplus)*jmin);
	  response[3]+=2.*imag(conj(jplus)*jmin);
	  response[4]+=2.*real(j0*(jplus-jmin));
	  response[5]+=2.*imag(j0*(jplus-jmin));
	  response[6]+=norm(jplus)-norm(jmin);
	  response[7]+=2.*imag(j0*conj(jplus-jmin));
	  response[8]+=2.*real(j0*conj(jplus-jmin));
	  
	}
      }
    }
  }
  double result=0.;
  result=kinfactors[0]*response[0]+kinfactors[1]*response[1]+kinfactors[2]*(response[2]*cos(2.*phi)+response[3]*sin(2.*phi))
	  +kinfactors[3]*(response[4]*cos(phi)+response[5]*cos(phi));
  delete reacmodel;
  return mott*result/HBARC;//mott*frontfactor*result/HBARC; 
  
}

