#include "DQEinclusive.hpp"
#include <Utilfunctions.hpp>
#include <gsl/gsl_poly.h>
#include <TSpinor.h>
#include <FourVector.h>
#include <TVector3.h>
#include <FastParticle.hpp>
#include <cassert>
#include <gsl/gsl_integration.h>

using namespace std;

const FourVector<complex<double> >DQEinclusive::polVectorPlus=FourVector<complex<double> >(0., -1./sqrt(2.),  														complex<double>(0.,-1./sqrt(2.)),0.);
const FourVector<complex<double> > DQEinclusive::polVectorMin=FourVector<complex<double> >(0.,
						  1./sqrt(2.),
						  complex<double>(0.,-1./sqrt(2.)),
						  0.);
const FourVector<complex<double> > DQEinclusive::polVector0 = FourVector<complex<double> >(1.,0.,0.,0.);


DQEinclusive::DQEinclusive(bool proton_in, int ff_param, string wavename, TElectronKinematics &elec, 
			       int offshell, double betaoffin, double lambdain)
:proton(proton_in),
betaoff(betaoffin*1.E-06),
lambda(lambdain*1.E+06),
massi(proton? MASSP:MASSN),
massr(proton? MASSN:MASSP),
offshellset(offshell),
electron(elec),
ffactorseq(NULL),
ffactorsdiff(NULL),
ffparam(ff_param){

  wf = TDeuteron::Wavefunction::CreateWavefunction(wavename);
  minpcm=1.E03;
}


DQEinclusive::~DQEinclusive(){
 delete wf; 
}


void DQEinclusive::calc_F2Dinc(double &contrib1, double &contrib2, double Q2,double x, int current, int integrator){
  
  ffactorseq=new NucleonEMOperator(Q2,proton,ffparam);
  ffactorsdiff=new NucleonEMOperator(Q2,!proton,ffparam);
  
  numint::array<double,1> lower = {{0.}};
  numint::array<double,1> upper = {{0.7E03}};
  DQEinclusive::Ftor_planewave F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.current = current;
  numint::mdfunction<numint::vector_d,1> mdf;
  mdf.func = &Ftor_planewave::exec;
  mdf.param = &F;
  numint::vector_d ret(2,0.);
  numint::vector_d error(2,0.);
  numint::vector_d prob(2,0.);
  F.f=DQEinclusive::planewave_int;
  int res=90;
  int count=0;
  unsigned int ccount =0;
  int fail     =0;
//     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
  switch(integrator){
    case(0):{  //numint adaptive cubature from MIT lib
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,2E04,ret,ccount,0);
      break;
    }
    case(1):{ //cuhre cubature from cuba lib
      int nvec=1;
      int flags=0x00;
      int key=11;
      char statefile[] = {""};
      int nregions =0;
      numint::cuhre(mdf,lower,upper,nvec,1.E-08,PREC,flags,1E03,2E04,key,statefile,nregions,count,fail,ret,error,prob);
      break;
    }
//     case(2):{ //suave MC from cuba lib
//       int nvec=1;
//       int flags = 0x00;
//       int seed = rand();
//       int nnew = 5000;
//       double flatness = 0.25;
//       char statefile[] = {""};
//       
//       // **** OUTPUT VARIABLES ****** //
//       int nregions    = 0;
//       numint::suave(mdf,lower,upper,nvec,1.E-08,PREC,flags,seed,1E03,2E04,nnew,flatness,
// 		    statefile,nregions,count,fail,ret,error,prob);
//       break;
//     }
//     case(3):{ //vegas MC from cuba lib
//       int nvec = 1;// unless you know what SIMD is leave this 1
//       int flags = 0x00;
//       int seed = rand();
//       int nstart = 100;
//       int nincrease = 1000;
//       int nbatch = 5000;
//       int gridno = 0;
//       char statefile[] = {""};
//       // **** OUTPUT VARIABLES ****** //
//       numint::vegas(mdf,lower,upper,nvec,1.E-08,PREC,flags,seed,
// 				1E03,2E04,nstart,nincrease,nbatch,gridno,
// 				statefile,count,fail,ret,error,prob);
//       break;
//     }
//     case(4):{ //divonne MC from cuba lib
// 	// **** INPUT VARIABLES ****** //
// 	int nvec = 1; // unless you know what SIMD is leave this 1
// 	int flags = 0x00;
// 	int seed  = rand();
// 	int key1 = 1;
// 	int key2 = 1;
// 	int key3 = 1;
// 	int maxpass = 50 ;
// 	double border = 0.;
// 	double maxchisq = 5;
// 	double mindeviation = 0.1;
// 	int ngiven = 0;
// 	int ldxgiven = 0;
// 	double* xgiven = NULL;
// 	int nextra = 0;
// 	peakfinder_t peakfinder = 0;
// 	char statefile[] = {""};
// 	
// 	// **** OUTPUT VARIABLES ****** //
// 	int nregions=0;
// 	numint::divonne(mdf,lower,upper,nvec,
// 		1.E-05,PREC,flags,seed,
// 		1E03,2E04,key1,key2,key3,
// 		maxpass,border,maxchisq,mindeviation,
// 		ngiven,ldxgiven,xgiven,nextra,peakfinder,
// 		statefile,nregions,count,fail,
// 		ret,error,prob);
// 	break;
//     }
    default:{
      cerr << "integrator not supported" << endl;
      assert(1==0);
    }
  }
//      cout << "pw " << res << " " << count << endl;
  delete ffactorseq;
  delete ffactorsdiff;
  contrib1= 2.*PI/3./MASSD*ret[0];
  contrib2= 2.*PI/3./MASSD*ret[1];    
}


void DQEinclusive::planewave_int(numint::vector_d & result, double pperp,
				   DQEinclusive& cross, double Q2, double x, int current){
  double phi=0.;
  result=numint::vector_d(2,0.);
  double res0=0.,res1=0.;
//   pperp=200.;
  if(pperp<1.E-03) return;
  double cosphi, sinphi;
  sincos(phi,&sinphi,&cosphi);
  vector<double> prz;
  double qvec=sqrt(Q2+pow(Q2/(2.*cross.massi*x),2.));
  double nu=Q2/(2.*cross.massi*x);



  if(cross.get_prz(prz,pperp,Q2,nu,qvec)){ //prz evaluation is performed!!
    for(size_t it=0;it<prz.size();it++){      
      double prnorm=sqrt(pperp*pperp+prz[it]*prz[it]);
//       double costheta=prz[it]/prnorm;
      double Ernorm=sqrt(pperp*pperp+prz[it]*prz[it]+cross.getMassr()*cross.getMassr());
      if(Ernorm<MASSD-200.){
	FourVector<double> q(nu,0.,0.,qvec);
	
	//direct term
	FourVector<double> pi(MASSD-Ernorm,-pperp*cosphi,-pperp*sinphi,-prz[it]); //pi=pD-ps
	FourVector<double> pn=pi+q;  //pn=pi+q
	GammaStructure J0 = cross.getFFactorseq()->getCC(current, q, pi, pn)*polVector0;
	GammaStructure Jplus = cross.getFFactorseq()->getCC(current, q, pi, pn)*polVectorPlus;
	GammaStructure Jmin = cross.getFFactorseq()->getCC(current, q, pi, pn)*polVectorMin;
	
	Matrix<1,4> un_down=TSpinor::Bar(TSpinor(pn,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
	Matrix<1,4> un_up=TSpinor::Bar(TSpinor(pn,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
	TSpinor ui_down(pi,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
	TSpinor ui_up(pi,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);

	//also needed for cross term
	FourVector<double> ps(Ernorm,pperp*cosphi,pperp*sinphi,prz[it]);
	FourVector<double> piprime=ps-q;  //pi'=ps-q
	double Erprimenorm=sqrt(pperp*pperp+(qvec-prz[it])*(qvec-prz[it])+cross.getMassi()*cross.getMassi());
	GammaStructure J0prime = cross.getFFactorsdiff()->getCC(current, q, piprime, ps)*polVector0;
	GammaStructure Jplusprime = cross.getFFactorsdiff()->getCC(current, q, piprime, ps)*polVectorPlus;
	GammaStructure Jminprime = cross.getFFactorsdiff()->getCC(current, q, piprime, ps)*polVectorMin;
	
	Matrix<1,4> us_down=TSpinor::Bar(TSpinor(ps,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
	Matrix<1,4> us_up=TSpinor::Bar(TSpinor(ps,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
	TSpinor uiprime_down(piprime,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
	TSpinor uiprime_up(piprime,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);

	res1=res0=0.;
	//summation over all the spin indices
	for(int spinn=-1;spinn<=1;spinn+=2){
	  for(int spini_in=-1;spini_in<=1;spini_in+=2){
	    complex<double> currentin0=(spinn==-1?un_down:un_up)*J0*(spini_in==-1?ui_down:ui_up);
	    complex<double> currentinplus=(spinn==-1?un_down:un_up)*Jplus*(spini_in==-1?ui_down:ui_up);
	    complex<double> currentinmin=(spinn==-1?un_down:un_up)*Jmin*(spini_in==-1?ui_down:ui_up);
	    for(int spini_out=-1;spini_out<=1;spini_out+=2){
	      //direct term
	      complex<double> currentout0=conj((spinn==-1?un_down:un_up)*J0*(spini_out==-1?ui_down:ui_up));
	      complex<double> currentoutplus=conj((spinn==-1?un_down:un_up)*Jplus*(spini_out==-1?ui_down:ui_up));
	      complex<double> currentoutmin=conj((spinn==-1?un_down:un_up)*Jmin*(spini_out==-1?ui_down:ui_up));
	      for(int spinr=-1;spinr<=1;spinr+=2){
		//cross term
		complex<double> currentout0prime=conj((spinr==-1?us_down:us_up)*J0prime*(spini_out==-1?uiprime_down:uiprime_up));
		complex<double> currentoutplusprime=conj((spinr==-1?us_down:us_up)*Jplusprime*(spini_out==-1?uiprime_down:uiprime_up));
		complex<double> currentoutminprime=conj((spinr==-1?us_down:us_up)*Jminprime*(spini_out==-1?uiprime_down:uiprime_up));
		for(int M=-2;M<=2;M+=2){
		  //direct
		  res0+=real(nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0
				      +Q2/(2.*qvec*qvec)*(currentinmin*currentoutmin+currentinplus*currentoutplus))*
				  cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp*cosphi,pperp*sinphi,prz[it]))
				  *conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinr, TVector3(pperp*cosphi,pperp*sinphi,prz[it]))));
		  //cross
		  if(Erprimenorm<MASSD-200.) res1+=real(nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0prime
				      +Q2/(2.*qvec*qvec)*(currentinmin*currentoutminprime+currentinplus*currentoutplusprime))*
				  cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp*cosphi,pperp*sinphi,prz[it]))
				  *conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinn, TVector3(-pperp*cosphi,-pperp*sinphi,qvec-prz[it]))));
									  
		}
	      }
	    }
	  }
	}
	res0*=pperp*MASSD/(2.*(MASSD-Ernorm))/2./abs(qvec-prz[it]/Ernorm*(MASSD+nu));
	//minus sign because of wave function!!
	if(Erprimenorm<MASSD-200.) res1*=-pperp*MASSD/(2.*sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm)))/2./abs(qvec-prz[it]/Ernorm*(MASSD+nu))
	  *sqrt(Erprimenorm/Ernorm);
	result[0]+=res0;
	result[1]+=res1;
      }

    }
  }
  

  return;
  

}


void DQEinclusive::calc_F2DincFSI(double &fsi1, double &fsi2, double &fsi1_off, double &fsi2_off, double Q2,double x, 
				  int current, int integrator){
  
  ffactorseq=new NucleonEMOperator(Q2,proton,ffparam);
  ffactorsdiff=new NucleonEMOperator(Q2,!proton,ffparam);

  
  fsi1=fsi2=0.;
  minpcm=1.E03;
  numint::array<double,3> lower = {{0.,0.,0.}};
  numint::array<double,3> upper = {{0.7E03,0.7E03,2.*PI}};
  DQEinclusive::Ftor_FSI F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.current = current;
  numint::mdfunction<numint::vector_d,3> mdf;
  mdf.func = &Ftor_FSI::exec;
  mdf.param = &F;
  numint::vector_d ret(4,0.);
  numint::vector_d error(4,0.);
  numint::vector_d prob(4,0.);
  F.f=DQEinclusive::FSI_int;
  int res=90;
  unsigned ccount=0;
  int count=0;
  int fail=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
  switch(integrator){
    case(0):{  //numint adaptive cubature from MIT lib
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,2E04,ret,ccount,0);
      break;
    }
    case(1):{ //cuhre cubature from cuba lib
      int nvec=1;
      int flags=0x00;
      int key=11;
      char statefile[] = {""};
      int nregions =0;
      numint::cuhre(mdf,lower,upper,nvec,1.E-08,PREC,flags,1E03,2E04,key,statefile,nregions,count,fail,ret,error,prob);
      break;
    }
 /*   case(2):{ //suave MC from cuba lib
      int nvec=1;
      int flags = 0x00;
      int seed = rand();
      int nnew = 5000;
      double flatness = 0.25;
      char statefile[] = {""};
      
      // **** OUTPUT VARIABLES ****** //
      int nregions    = 0;
      numint::suave(mdf,lower,upper,nvec,1.E-08,PREC,flags,seed,1E03,2E04,nnew,flatness,
		    statefile,nregions,count,fail,ret,error,prob);
      break;
    }
    case(3):{ //vegas MC from cuba lib
      int nvec = 1;// unless you know what SIMD is leave this 1
      int flags = 0x00;
      int seed = rand();
      int nstart = 100;
      int nincrease = 1000;
      int nbatch = 5000;
      int gridno = 0;
      char statefile[] = {""};
      // **** OUTPUT VARIABLES ****** //
      numint::vegas(mdf,lower,upper,nvec,1.E-08,PREC,flags,seed,
				1E03,2E04,nstart,nincrease,nbatch,gridno,
				statefile,count,fail,ret,error,prob);
      break;
    }
    case(4):{ //divonne MC from cuba lib
	// **** INPUT VARIABLES ****** //
	int nvec = 1; // unless you know what SIMD is leave this 1
	int flags = 0x00;
	int seed  = rand();
	int key1 = 1;
	int key2 = 1;
	int key3 = 1;
	int maxpass = 50 ;
	double border = 0.;
	double maxchisq = 5;
	double mindeviation = 0.1;
	int ngiven = 0;
	int ldxgiven = 0;
	double* xgiven = NULL;
	int nextra = 0;
	peakfinder_t peakfinder = 0;
	char statefile[] = {""};
	
	// **** OUTPUT VARIABLES ****** //
	int nregions=0;
	numint::divonne(mdf,lower,upper,nvec,
		1.E-05,PREC,flags,seed,
		1E03,2E04,key1,key2,key3,
		maxpass,border,maxchisq,mindeviation,
		ngiven,ldxgiven,xgiven,nextra,peakfinder,
		statefile,nregions,count,fail,
		ret,error,prob);
	break;
    }
 */   default:{
      cerr << "integrator not supported" << endl;
      assert(1==0);
    }
  }
//     cout << "fsi " << res << " " << count << endl;
  fsi1= 2.*PI/3./MASSD/(64.*PI*PI)*ret[0];
  fsi2= 2.*PI/3./MASSD/(64.*PI*PI)*ret[1];
  fsi1_off= 2.*PI/3./MASSD/(64.*PI*PI)*ret[2];
  fsi2_off= 2.*PI/3./MASSD/(64.*PI*PI)*ret[3];

  
  delete ffactorseq;
  delete ffactorsdiff;
  
  return;
}

void DQEinclusive::FSI_int(numint::vector_d & result, double pperp1, double qt, 
		      double qphi, DQEinclusive& cross, double Q2, double x, int current){
			
  double phi1=0.;
  result=numint::vector_d(4,0.);
  complex<double> res0=0.,res1=0.,res2=0.,res3=0.;
  if(pperp1<1.E-03||qt<1.E-03) { result[0]=result[1]=result[2]=result[3]=0.; return;}

  double cosphi1, sinphi1;
  sincos(phi1,&sinphi1,&cosphi1);
  double cosqphi, sinqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  vector<double> prz1, prz2;
  double qvec=sqrt(Q2+pow(Q2/(2.*cross.massi*x),2.));
  double nu=Q2/(2.*cross.massi*x);
  double s=(MASSD+nu)*(MASSD+nu)-qvec*qvec;
  double pcm=sqrt(s/4.-MASSP*MASSP);
  if(pcm<cross.minpcm) cross.minpcm=pcm;
  double chi=sqrt(s)*pcm*2.;

  FastParticle rescatter(0,0,pcm,0.,0.,Q2,0.,"");
  
  static FourVector<complex<double> > polVectorPlus(0.,
                                                    -1./sqrt(2.),
                                                    complex<double>(0.,-1./sqrt(2.)),
                                                 0.);
  static FourVector<complex<double> > polVectorMin(0.,
                                                   1./sqrt(2.),
                                                   complex<double>(0.,-1./sqrt(2.)),
                                                   0.);
  static FourVector<complex<double> > polVector0(1.,0.,0.,0.);
  
  double pperp2=sqrt(pow(pperp1+qt*cosqphi,2.)+pow(qt*sinqphi,2.));



  if(cross.get_prz(prz1,pperp1,Q2,nu,qvec)){ //prz evaluation is performed!!
    if(cross.get_prz(prz2,pperp2,Q2,nu,qvec)){ //prz evaluation is performed!!
      for(size_t it1=0;it1<prz1.size();it1++){      
	double prnorm=sqrt(pperp1*pperp1+prz1[it1]*prz1[it1]);
  //       double costheta=prz[it]/prnorm;
	double Ernorm=sqrt(pperp1*pperp1+prz1[it1]*prz1[it1]+cross.getMassr()*cross.getMassr());
	if(Ernorm<MASSD-200.){
	  FourVector<double> q(nu,0.,0.,qvec);
	  
	  //direct term
	  FourVector<double> pi1(MASSD-Ernorm,-pperp1*cosphi1,-pperp1*sinphi1,-prz1[it1]); //pi=pD-ps
	  FourVector<double> pn1=pi1+q;  //pn1=pi+q
	  GammaStructure J0 = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVector0;
	  GammaStructure Jplus = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVectorPlus;
	  GammaStructure Jmin = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVectorMin;
	  
	  Matrix<1,4> un1_down=TSpinor::Bar(TSpinor(pn1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
	  Matrix<1,4> un1_up=TSpinor::Bar(TSpinor(pn1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
	  TSpinor ui1_down(pi1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
	  TSpinor ui1_up(pi1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
	
	  for(size_t it2=0;it2<prz2.size();it2++){      
	    double Erprimenorm=sqrt(pperp2*pperp2+prz2[it2]*prz2[it2]+cross.getMassr()*cross.getMassr());
	    double Erprimenorm_cross=sqrt(pperp2*pperp2+prz2[it2]*prz2[it2]+cross.getMassi()*cross.getMassi());
	    if(Erprimenorm<MASSD-200.){

	      FourVector<double> pi2(MASSD-Erprimenorm,-pperp1*cosphi1-qt*cosqphi,-qt*sinqphi,-prz2[it2]); //pi2=pD-ps1+qt
	      FourVector<double> pn2=pi2+q; //pn2=pi2+q
	      double t=pow(Ernorm-Erprimenorm,2.)-pow(prz1[it1]-prz2[it2],2.)-qt*qt;
	      double t_cross=pow(pn2[0]-Ernorm,2.)-pow(pn2[3]-prz1[it1],2.)-qt*qt;
	      //double t_cross2=pow(pn1[0]-Erprimenorm,2.)-pow(pn1[3]-prz2[it2],2.)-qt*qt;
	      complex<double> directrescatt = rescatter.scatter(t,0);
	      complex<double> crossrescatt = rescatter.scatter(t_cross,0);
	      //direct term
	      GammaStructure J0prime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVector0;
	      GammaStructure Jplusprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorPlus;
	      GammaStructure Jminprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorMin;	     
	      
	      Matrix<1,4> un2_down=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
	      Matrix<1,4> un2_up=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
	      TSpinor ui2_down(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
	      TSpinor ui2_up(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
	      
	      //crossed term
	      FourVector<double> pi2_cross(MASSD-Erprimenorm_cross,pperp1*cosphi1+qt*cosqphi,qt*sinqphi,-prz2[it2]); //pi2=pD-ps1+qt
	      FourVector<double> pn2_cross=pi2_cross+q; //pn2=pi2+q
	      GammaStructure J0prime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVector0;
	      GammaStructure Jplusprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorPlus;
	      GammaStructure Jminprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorMin;	     
	      
	      Matrix<1,4> un2_down_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
	      Matrix<1,4> un2_up_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
	      TSpinor ui2_down_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
	      TSpinor ui2_up_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
	      
	      res0=res1=res2=res3=0.;
	      //summation over all the spin indices
	      for(int spinn=-1;spinn<=1;spinn+=2){
		for(int spini_in=-1;spini_in<=1;spini_in+=2){
		  complex<double> currentin0=(spinn==-1?un1_down:un1_up)*J0*(spini_in==-1?ui1_down:ui1_up);
		  complex<double> currentinplus=(spinn==-1?un1_down:un1_up)*Jplus*(spini_in==-1?ui1_down:ui1_up);
		  complex<double> currentinmin=(spinn==-1?un1_down:un1_up)*Jmin*(spini_in==-1?ui1_down:ui1_up);
		  for(int spini_out=-1;spini_out<=1;spini_out+=2){
		    //direct term
		    complex<double> currentout0=conj((spinn==-1?un2_down:un2_up)*J0prime*(spini_out==-1?ui2_down:ui2_up));
		    complex<double> currentoutplus=conj((spinn==-1?un2_down:un2_up)*Jplusprime*(spini_out==-1?ui2_down:ui2_up));
		    complex<double> currentoutmin=conj((spinn==-1?un2_down:un2_up)*Jminprime*(spini_out==-1?ui2_down:ui2_up));
		    //cross term
		    complex<double> currentout0_cross=conj((spinn==-1?un2_down_cross:un2_up_cross)*J0prime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
		    complex<double> currentoutplus_cross=conj((spinn==-1?un2_down_cross:un2_up_cross)*Jplusprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
		    complex<double> currentoutmin_cross=conj((spinn==-1?un2_down_cross:un2_up_cross)*Jminprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
		    for(int M=-2;M<=2;M+=2){
		      for(int spinr=-1;spinr<=1;spinr+=2){
			complex<double> wf=cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
					*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinr, 
						TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2])));
			complex<double> wf2=cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
					*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinn, 
						TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,prz2[it2])));
			complex<double> wf_off=cross.getDeutwf()->DeuteronPStateOff(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
					*conj(cross.getDeutwf()->DeuteronPStateOff(M, spini_out, spinr, 
						TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2])));
			complex<double> wf2_off=cross.getDeutwf()->DeuteronPStateOff(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
					*conj(cross.getDeutwf()->DeuteronPStateOff(M, spini_out, spinn, 
						TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,prz2[it2])));
			//direct onshell
			res0+=nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0+Q2/(2.*qvec*qvec)
			  *(currentinmin*currentoutmin+currentinplus*currentoutplus))*wf*directrescatt;
			//cross onshell
			res1+=nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0_cross+Q2/(2.*qvec*qvec)*
			  (currentinmin*currentoutmin_cross+currentinplus*currentoutplus_cross))*wf2*crossrescatt;
			//direct offshell
			res2+=nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0+Q2/(2.*qvec*qvec)
			  *(currentinmin*currentoutmin+currentinplus*currentoutplus))*wf_off*directrescatt;
			//cross offshell
			res3+=nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0_cross+Q2/(2.*qvec*qvec)*
			  (currentinmin*currentoutmin_cross+currentinplus*currentoutplus_cross))*wf2_off*crossrescatt;
		      }
		    }
		  }
		}
	      }
	      res0*=-pperp1*qt/sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))
	      /abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
	      res2*=pperp1*qt/sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))
	      /abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
	      //minus sign because of wave function!!
	      res1*=pperp1*qt*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm_cross))/sqrt(Ernorm*Erprimenorm_cross)
	      /abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))/abs(qvec-prz2[it2]/Erprimenorm_cross*(MASSD+nu))*chi;
	      res3*=-pperp1*qt*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm_cross))/sqrt(Ernorm*Erprimenorm_cross)
	      /abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))/abs(qvec-prz2[it2]/Erprimenorm_cross*(MASSD+nu))*chi;
	      result[0]+=imag(res0);
	      result[1]+=imag(res1);
	      result[2]+=imag(res2);
	      result[3]+=imag(res3);
//  	      cout << "BLAAAAAAAAA " << it1<< " " << it2 << " " << result[0] << " " << res0 << endl;

	    }
	  }
	}
      }

    }
  }
  
  
  return;
  

}


void DQEinclusive::calc_F2DincFSI_PVoff(double &fsi1_off, double &fsi2_off, double Q2,double x, 
				  int current, int integrator){
  
  ffactorseq=new NucleonEMOperator(Q2,proton,ffparam);
  ffactorsdiff=new NucleonEMOperator(Q2,!proton,ffparam);

  
  minpcm=1.E03;
  numint::array<double,3> lower = {{0.,0.,0.}};
  numint::array<double,3> upper = {{0.7E03,0.7E03,2.*PI}};
  DQEinclusive::Ftor_FSI F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.current = current;
  numint::mdfunction<numint::vector_d,3> mdf;
  mdf.func = &Ftor_FSI::exec;
  mdf.param = &F;
  numint::vector_d ret(4,0.);
  numint::vector_d error(4,0.);
  numint::vector_d prob(4,0.);
  F.f=DQEinclusive::FSI_PV;
  int res=90;
  unsigned ccount=0;
  int count=0;
  int fail=0;
  //     res = numint::cube_romb(mdf,lower,upper,1.E-08,PREC,ret,count,0);
  switch(integrator){
    case(0):{  //numint adaptive cubature from MIT lib
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,2E04,ret,ccount,0);
      break;
    }
    case(1):{ //cuhre cubature from cuba lib
      int nvec=1;
      int flags=0x00;
      int key=11;
      char statefile[] = {""};
      int nregions =0;
      numint::cuhre(mdf,lower,upper,nvec,1.E-08,PREC,flags,1E03,2E04,key,statefile,nregions,count,fail,ret,error,prob);
      break;
    }
 /*   case(2):{ //suave MC from cuba lib
      int nvec=1;
      int flags = 0x00;
      int seed = rand();
      int nnew = 5000;
      double flatness = 0.25;
      char statefile[] = {""};
      
      // **** OUTPUT VARIABLES ****** //
      int nregions    = 0;
      numint::suave(mdf,lower,upper,nvec,1.E-08,PREC,flags,seed,1E03,2E04,nnew,flatness,
		    statefile,nregions,count,fail,ret,error,prob);
      break;
    }
    case(3):{ //vegas MC from cuba lib
      int nvec = 1;// unless you know what SIMD is leave this 1
      int flags = 0x00;
      int seed = rand();
      int nstart = 100;
      int nincrease = 1000;
      int nbatch = 5000;
      int gridno = 0;
      char statefile[] = {""};
      // **** OUTPUT VARIABLES ****** //
      numint::vegas(mdf,lower,upper,nvec,1.E-08,PREC,flags,seed,
				1E03,2E04,nstart,nincrease,nbatch,gridno,
				statefile,count,fail,ret,error,prob);
      break;
    }
    case(4):{ //divonne MC from cuba lib
	// **** INPUT VARIABLES ****** //
	int nvec = 1; // unless you know what SIMD is leave this 1
	int flags = 0x00;
	int seed  = rand();
	int key1 = 1;
	int key2 = 1;
	int key3 = 1;
	int maxpass = 50 ;
	double border = 0.;
	double maxchisq = 5;
	double mindeviation = 0.1;
	int ngiven = 0;
	int ldxgiven = 0;
	double* xgiven = NULL;
	int nextra = 0;
	peakfinder_t peakfinder = 0;
	char statefile[] = {""};
	
	// **** OUTPUT VARIABLES ****** //
	int nregions=0;
	numint::divonne(mdf,lower,upper,nvec,
		1.E-05,PREC,flags,seed,
		1E03,2E04,key1,key2,key3,
		maxpass,border,maxchisq,mindeviation,
		ngiven,ldxgiven,xgiven,nextra,peakfinder,
		statefile,nregions,count,fail,
		ret,error,prob);
	break;
    }
 */   default:{
      cerr << "integrator not supported" << endl;
      assert(1==0);
    }
  }
//     cout << "fsi " << res << " " << count << endl;
  fsi1_off= 2.*PI/3./MASSD/(64.*PI*PI)*ret[0];
  fsi2_off= 2.*PI/3./MASSD/(64.*PI*PI)*ret[1];

  
  delete ffactorseq;
  delete ffactorsdiff;
  
  return;
}

void DQEinclusive::FSI_PV(numint::vector_d & result, double pperp1, double qt, 
		      double qphi, DQEinclusive& cross, double Q2, double x, int current){
  
  double phi1=0.;
  result=numint::vector_d(2,0.);
  if(pperp1<1.E-03||qt<1.E-03) { result[0]=result[1]=0.; return;}
  double cosphi1, sinphi1;
  sincos(phi1,&sinphi1,&cosphi1);
  double cosqphi, sinqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  vector<double> prz1, prz2;
  double qvec=sqrt(Q2+pow(Q2/(2.*cross.massi*x),2.));
  double nu=Q2/(2.*cross.massi*x);
  double s=(MASSD+nu)*(MASSD+nu)-qvec*qvec;
  double pcm=sqrt(s/4.-MASSP*MASSP);
  if(pcm<cross.minpcm) cross.minpcm=pcm;
  double chi=sqrt(s)*pcm*2.;
  FastParticle rescatter(0,0,pcm,0.,0.,Q2,0.,"");

  double pperp2=sqrt(pow(pperp1+qt*cosqphi,2.)+pow(qt*sinqphi,2.));

  Ftor_PV pvint;
  pvint.cross=&cross;
  pvint.Q2 = Q2;
  pvint.x = x;
  pvint.current = current;
  pvint.cosqphi=cosqphi;
  pvint.sinqphi=sinqphi;
  pvint.cosphi1=cosphi1;
  pvint.sinphi1=sinphi1;
  pvint.pperp1=pperp1;
  pvint.pperp2=pperp2;
  pvint.qt=qt;
  pvint.qvec=qvec;
  pvint.nu=nu;
  
  
  if(cross.get_prz(prz1,pperp1,Q2,nu,qvec)){ //prz evaluation is performed!!
    if(cross.get_prz(prz2,pperp2,Q2,nu,qvec)){ //prz evaluation is performed!!
      pvint.prz2poles = prz2;
      for(size_t it1=0;it1<prz1.size();it1++){      
	for(bool crossed=0;crossed<=1;crossed++){
	  pvint.crossed=crossed;
	  gsl_integration_workspace * w 
	  = gsl_integration_workspace_alloc (1E6);
	  double res, error;
	  double expected = 1.E-03;
	  gsl_function F1;
	  F1.function = &PV_int1;
	  F1.params = &pvint;
	  gsl_integration_qawc (&F1, -700., 700., prz1[it1], 1.e-08, PREC, 1E3,
				w, &res, &error); 
	  if(!crossed) result[0]+=res;
	  else result[1]+=res;
	  gsl_integration_workspace_free (w);

	}
      }
    }
  }
  result[0]*=chi*qt*pperp1;
  result[1]*=chi*qt*pperp1;
}

double DQEinclusive::PV_int1(double pz1, void * params){
  
  Ftor_PV pvint = *(Ftor_PV *) params;
  pvint.var1=pz1;
  double res=0.;
  for(size_t it2=0;it2<pvint.prz2poles.size();it2++){ 
    gsl_integration_workspace * w2 
      = gsl_integration_workspace_alloc (1E6);

    double result, error;
    double expected = -4.0;
    gsl_function F2;
    F2.function = &Ftor_PV::exec;
    F2.params = &pvint;
    gsl_integration_qawc (&F2, -700., 700.,pvint.prz2poles[it2] , 1.e-08, PREC, 1E3,
			  w2, &result, &error); 
    res+=result;
    gsl_integration_workspace_free (w2);
  }
  return res;
  
}



double DQEinclusive::FSI_intPV(double pz1, double pz2, double pperp1, double pperp2, double qt, double cosphi1, 
			       double sinphi1, double cosqphi, double sinqphi, vector<double> &prz2poles,  
			       DQEinclusive& cross, double Q2, double x, double qvec, double nu, int current, bool crossed){
			
//   complex<double> res0=0.,res1=0.,res2=0.,res3=0.;
// 
// 
//   
//   
// 
// 
// 
// 	double prnorm=sqrt(pperp1*pperp1+prz1[it1]*prz1[it1]);
//   //       double costheta=prz[it]/prnorm;
// 	double Ernorm=sqrt(pperp1*pperp1+prz1[it1]*prz1[it1]+cross.getMassr()*cross.getMassr());
// 	if(Ernorm<MASSD-200.){
// 	  FourVector<double> q(nu,0.,0.,qvec);
// 	  
// 	  //direct term
// 	  FourVector<double> pi1(MASSD-Ernorm,-pperp1*cosphi1,-pperp1*sinphi1,-prz1[it1]); //pi=pD-ps
// 	  FourVector<double> pn1=pi1+q;  //pn1=pi+q
// 	  GammaStructure J0 = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVector0;
// 	  GammaStructure Jplus = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVectorPlus;
// 	  GammaStructure Jmin = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVectorMin;
// 	  
// 	  Matrix<1,4> un1_down=TSpinor::Bar(TSpinor(pn1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
// 	  Matrix<1,4> un1_up=TSpinor::Bar(TSpinor(pn1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
// 	  TSpinor ui1_down(pi1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
// 	  TSpinor ui1_up(pi1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
// 	
// 	  for(size_t it2=0;it2<prz2.size();it2++){      
// 	    double Erprimenorm=sqrt(pperp2*pperp2+prz2[it2]*prz2[it2]+cross.getMassr()*cross.getMassr());
// 	    double Erprimenorm_cross=sqrt(pperp2*pperp2+prz2[it2]*prz2[it2]+cross.getMassi()*cross.getMassi());
// 	    if(Erprimenorm<MASSD-200.){
// 
// 	      FourVector<double> pi2(MASSD-Erprimenorm,-pperp1*cosphi1-qt*cosqphi,-qt*sinqphi,-prz2[it2]); //pi2=pD-ps1+qt
// 	      FourVector<double> pn2=pi2+q; //pn2=pi2+q
// 	      double t=pow(Ernorm-Erprimenorm,2.)-pow(prz1[it1]-prz2[it2],2.)-qt*qt;
// 	      double t_cross=pow(pn2[0]-Ernorm,2.)-pow(pn2[3]-prz1[it1],2.)-qt*qt;
// 	      //double t_cross2=pow(pn1[0]-Erprimenorm,2.)-pow(pn1[3]-prz2[it2],2.)-qt*qt;
// 	      complex<double> directrescatt = rescatter.scatter(t,0);
// 	      complex<double> crossrescatt = rescatter.scatter(t_cross,0);
// 	      //direct term
// 	      GammaStructure J0prime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVector0;
// 	      GammaStructure Jplusprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorPlus;
// 	      GammaStructure Jminprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorMin;	     
// 	      
// 	      Matrix<1,4> un2_down=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
// 	      Matrix<1,4> un2_up=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
// 	      TSpinor ui2_down(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
// 	      TSpinor ui2_up(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
// 	      
// 	      //crossed term
// 	      FourVector<double> pi2_cross(MASSD-Erprimenorm_cross,pperp1*cosphi1+qt*cosqphi,qt*sinqphi,-prz2[it2]); //pi2=pD-ps1+qt
// 	      FourVector<double> pn2_cross=pi2_cross+q; //pn2=pi2+q
// 	      GammaStructure J0prime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVector0;
// 	      GammaStructure Jplusprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorPlus;
// 	      GammaStructure Jminprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorMin;	     
// 	      
// 	      Matrix<1,4> un2_down_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
// 	      Matrix<1,4> un2_up_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
// 	      TSpinor ui2_down_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
// 	      TSpinor ui2_up_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
// 	      
// 	      res0=res1=res2=res3=0.;
// 	      //summation over all the spin indices
// 	      for(int spinn=-1;spinn<=1;spinn+=2){
// 		for(int spini_in=-1;spini_in<=1;spini_in+=2){
// 		  complex<double> currentin0=(spinn==-1?un1_down:un1_up)*J0*(spini_in==-1?ui1_down:ui1_up);
// 		  complex<double> currentinplus=(spinn==-1?un1_down:un1_up)*Jplus*(spini_in==-1?ui1_down:ui1_up);
// 		  complex<double> currentinmin=(spinn==-1?un1_down:un1_up)*Jmin*(spini_in==-1?ui1_down:ui1_up);
// 		  for(int spini_out=-1;spini_out<=1;spini_out+=2){
// 		    //direct term
// 		    complex<double> currentout0=conj((spinn==-1?un2_down:un2_up)*J0prime*(spini_out==-1?ui2_down:ui2_up));
// 		    complex<double> currentoutplus=conj((spinn==-1?un2_down:un2_up)*Jplusprime*(spini_out==-1?ui2_down:ui2_up));
// 		    complex<double> currentoutmin=conj((spinn==-1?un2_down:un2_up)*Jminprime*(spini_out==-1?ui2_down:ui2_up));
// 		    //cross term
// 		    complex<double> currentout0_cross=conj((spinn==-1?un2_down_cross:un2_up_cross)*J0prime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
// 		    complex<double> currentoutplus_cross=conj((spinn==-1?un2_down_cross:un2_up_cross)*Jplusprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
// 		    complex<double> currentoutmin_cross=conj((spinn==-1?un2_down_cross:un2_up_cross)*Jminprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
// 		    for(int M=-2;M<=2;M+=2){
// 		      for(int spinr=-1;spinr<=1;spinr+=2){
// 			complex<double> wf=cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
// 					*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinr, 
// 						TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2])));
// 			complex<double> wf2=cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
// 					*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinn, 
// 						TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,prz2[it2])));
// 			complex<double> wf_off=cross.getDeutwf()->DeuteronPStateOff(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
// 					*conj(cross.getDeutwf()->DeuteronPStateOff(M, spini_out, spinr, 
// 						TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2])));
// 			complex<double> wf2_off=cross.getDeutwf()->DeuteronPStateOff(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
// 					*conj(cross.getDeutwf()->DeuteronPStateOff(M, spini_out, spinn, 
// 						TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,prz2[it2])));
// 			//direct onshell
// 			res0+=nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0+Q2/(2.*qvec*qvec)
// 			  *(currentinmin*currentoutmin+currentinplus*currentoutplus))*wf*directrescatt;
// 			//cross onshell
// 			res1+=nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0_cross+Q2/(2.*qvec*qvec)*
// 			  (currentinmin*currentoutmin_cross+currentinplus*currentoutplus_cross))*wf2*crossrescatt;
// 			//direct offshell
// 			res2+=nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0+Q2/(2.*qvec*qvec)
// 			  *(currentinmin*currentoutmin+currentinplus*currentoutplus))*wf_off*directrescatt;
// 			//cross offshell
// 			res3+=nu*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0_cross+Q2/(2.*qvec*qvec)*
// 			  (currentinmin*currentoutmin_cross+currentinplus*currentoutplus_cross))*wf2_off*crossrescatt;
// 		      }
// 		    }
// 		  }
// 		}
// 	      }
// 	      res0*=-pperp1*qt/sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))
// 	      /abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
// 	      res2*=pperp1*qt/sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))
// 	      /abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
// 	      //minus sign because of wave function!!
// 	      res1*=pperp1*qt*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm_cross))/sqrt(Ernorm*Erprimenorm_cross)
// 	      /abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))/abs(qvec-prz2[it2]/Erprimenorm_cross*(MASSD+nu))*chi;
// 	      res3*=-pperp1*qt*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm_cross))/sqrt(Ernorm*Erprimenorm_cross)
// 	      /abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))/abs(qvec-prz2[it2]/Erprimenorm_cross*(MASSD+nu))*chi;
// 	      result[0]+=imag(res0);
// 	      result[1]+=imag(res1);
// 	      result[2]+=imag(res2);
// 	      result[3]+=imag(res3);
// //  	      cout << "BLAAAAAAAAA " << it1<< " " << it2 << " " << result[0] << " " << res0 << endl;
// 
// 	    }
// 	  }
// 	}
//       }
// 
//     }
//   }
//   
  
  return 0.;
  

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
//     cout << sol.size() << " " << A+B*p1 << " " << A+B*p2 << " " << Er << " " << Er2 << " " << p1 << " " << p2 << endl;
    return 1;
  }
}
