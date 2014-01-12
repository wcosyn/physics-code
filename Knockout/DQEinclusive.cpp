#include "DQEinclusive.hpp"
#include <Utilfunctions.hpp>
#include <gsl/gsl_poly.h>
#include <TSpinor.h>
#include <FourVector.h>
#include <TVector3.h>


using namespace std;

DQEinclusive::DQEinclusive(bool proton_in, int ff_param, string wavename, TElectronKinematics &elec, 
			       int offshell, double betaoffin, double lambdain)
:proton(proton_in),
betaoff(betaoffin*1.E-06),
lambda(lambdain*1.E+06),
massi(proton? MASSP:MASSN),
massr(proton? MASSN:MASSP),
offshellset(offshell),
electron(elec),
ffactors(NULL),
ffparam(ff_param){

  wf = TDeuteron::Wavefunction::CreateWavefunction(wavename);
  
}


DQEinclusive::~DQEinclusive(){
 delete wf; 
}


void DQEinclusive::calc_F2Dinc(double &contrib1, double &contrib2, double Q2,double x, int current){
  
  double result;
  ffactors=new NucleonEMOperator(Q2,proton,ffparam);
  
  numint::array<double,2> lower = {{0.,0.}};
  numint::array<double,2> upper = {{1.E03,2.*PI}};
  DQEinclusive::Ftor_planewave F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.current = current;
  numint::mdfunction<numint::vector_d,2> mdf;
  mdf.func = &Ftor_planewave::exec;
  mdf.param = &F;
  numint::vector_d ret(2,0.);
  F.f=DQEinclusive::planewave_int;
  int res=90;
  unsigned count=0;
//     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
  res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,2E04,ret,count,0);
//     cout << res << " " << count << endl;
  delete ffactors;
  contrib1= 2.*PI/3./MASSD*ret[0];
  contrib2= 2.*PI/3./MASSD*ret[1];    
}


void DQEinclusive::planewave_int(numint::vector_d & result, double pperp, double phi,
				   DQEinclusive& cross, double Q2, double x, int current){
  
  result=numint::vector_d(2,0.);
  if(pperp<1.E-03) return;
  double cosphi, sinphi;
  sincos(phi,&sinphi,&cosphi);
  vector<double> prz;
  double qvec=sqrt(Q2+pow(Q2/(2.*cross.massi*x),2.));
  double nu=Q2/(2.*cross.massi*x);
  static FourVector<complex<double> > polVectorPlus(0.,
                                                    -1./sqrt(2.),
                                                    complex<double>(0.,-1./sqrt(2.)),
                                                 0.);
  static FourVector<complex<double> > polVectorMin(0.,
                                                   1./sqrt(2.),
                                                   complex<double>(0.,-1./sqrt(2.)),
                                                   0.);
  static FourVector<complex<double> > polVector0(1.,0.,0.,0.);
  



  if(cross.get_prz(prz,pperp,Q2,nu,qvec)){ //prz evaluation is performed!!

    for(size_t it=0;it<prz.size();it++){      
      double prnorm=sqrt(pperp*pperp+prz[it]*prz[it]);
      double costheta=prz[it]/prnorm;
      double Ernorm=sqrt(prnorm*prnorm+cross.getMassr()*cross.getMassr());
      FourVector<double> q(nu,0.,0.,qvec);
      
      //direct term
      FourVector<double> pi(MASSD-Ernorm,-pperp*cosphi,-pperp*sinphi,prz[it]);
      FourVector<double> pn=pi+q;  //pi+q
      GammaStructure J0 = cross.getFFactors()->getCC(current, q, pi, pn)*polVector0;
      GammaStructure Jplus = cross.getFFactors()->getCC(current, q, pi, pn)*polVectorPlus;
      GammaStructure Jmin = cross.getFFactors()->getCC(current, q, pi, pn)*polVectorMin;
      
      
    }
    

  }
  

  
//   result[0]=pperp*(pow(cross.wf->GetUp(prnorm),2.)+pow(cross.wf->GetWp(prnorm),2.))/(4.*PI)*(MASSD/(2.*(MASSD-Er)))*structfactor;
  return;
  

}



bool DQEinclusive::get_prz(vector<double> &sol, double pt, double Q2, double nu, double qvec){
  sol.clear();
  double A=(massi*massi-getMassr()*getMassr()-MASSD*MASSD+Q2-2.*MASSD*nu)/(-2.*(MASSD+nu));
  double B= qvec/(MASSD+nu);
  double aa=B*B-1.;
  double bb=2.*A*B;
  double cc=A*A-getMassr()*getMassr()-pt*pt;
  double discr=bb*bb-4.*aa*cc;
  if(discr<0.) {return 0;}
  if(abs(discr)<1.E-09) {sol.push_back(-bb/(2.*aa)); return 1;}
  if(discr>0.){
    double p1 = (-bb+sqrt(discr))/(2.*aa);
    double Er=sqrt(getMassr()*getMassr()+pt*pt+p1*p1);
    if(SIGN(A+B*p1)==1) sol.push_back(p1);
    double p2 = (-bb-sqrt(discr))/(2.*aa);
    double Er2=sqrt(getMassr()*getMassr()+pt*pt+p2*p2);
    if(SIGN(A+B*p2)==1) sol.push_back(p2);
    cout << sol.size() << " " << A+B*p1 << " " << A+B*p2 << " " << Er << " " << Er2 << " " << p1 << " " << p2 << endl;
    return 1;
  }
}
