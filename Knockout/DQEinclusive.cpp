#include "DQEinclusive.hpp"
#include <Utilfunctions.hpp>
#include <AuxFunction.hpp>
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
massi(proton_in? MASSP:MASSN),
massr(proton_in? MASSN:MASSP),
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


void DQEinclusive::calc_Crossinc(vector<double> &pw, double Q2,double x, int current, int integrator, double thetapol){
  
  double nu=Q2/(x*2.*massi);
  ffactorseq=new NucleonEMOperator(Q2,proton,ffparam);
  ffactorsdiff=new NucleonEMOperator(Q2,!proton,ffparam);
  
  numint::array<double,1> lower = {{0.}};
  numint::array<double,1> upper = {{0.7E03}};
  DQEinclusive::Ftor_planewave F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.current = current;
  F.tanhalfth2=electron.GetTan2HalfAngle(Q2,nu);
  F.thetapol=thetapol;
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
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,2E04,ret,ccount,0);
      break;
    }
//     case(1):{ //cuhre cubature from cuba lib
//       int nvec=1;
//       int flags=0x00;
//       int key=11;
//       char statefile[] = {""};
//       int nregions =0;
//       numint::cuhre(mdf,lower,upper,nvec,1.E-08,PREC,flags,1E03,2E04,key,statefile,nregions,count,fail,ret,error,prob);
//       break;
//     }
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
      cerr << "integrator not supported " << integrator << endl;
      assert(1==0);
    }
  }
//      cout << "pw " << res << " " << count << endl;
  delete ffactorseq;
  delete ffactorsdiff;
  double sigmamott=ALPHA*ALPHA*(1+electron.GetCosScatterAngle(Q2,nu))/
    (2.*pow(electron.GetBeamEnergy(Q2,nu)*(1-electron.GetCosScatterAngle(Q2,nu)),2.))*HBARC*HBARC*1.E07;
//     cout << sigmamott << " " << acos(electron.GetCosScatterAngle(Q2,nu))/PI*180. << " " << 2.*PI/3./MASSD*1.E03*ret[0] << endl;
  for(int i=0;i<4;i++) pw[i]= 2.*PI/3./MASSD*1.E03*ret[i]*sigmamott;
}


void DQEinclusive::planewave_int(numint::vector_d & result, double pperp,
				   DQEinclusive& cross, double Q2, double x, int current, double tanhalfth2, double thetapol){
  double phi=0.;
  result=numint::vector_d(4,0.);
  double res0=0.,res1=0.,res2=0.,res3=0.;
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
//       cout << "pw " << it << " " << pperp << " " <<  2.*qvec*prz[it]-2.*(MASSD+nu)*Ernorm << " " << -MASSD*MASSD+Q2-2.*MASSD*nu -cross.getMassr()*cross.getMassr()+ cross.getMassi()*cross.getMassi() << endl;
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

	res3=res2=res1=res0=0.;
	//summation over all the spin indices
	for(int spinn=-1;spinn<=0;spinn+=2){  //exploit parity symmetry here!
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

		complex<double> wavemin=cross.getDeutwf()->DeuteronPState(-2, spini_in, spinr, TVector3(pperp*cosphi,pperp*sinphi,prz[it]));
		complex<double> wavezero=cross.getDeutwf()->DeuteronPState(0, spini_in, spinr, TVector3(pperp*cosphi,pperp*sinphi,prz[it]));
		complex<double> waveplus=cross.getDeutwf()->DeuteronPState(2, spini_in, spinr, TVector3(pperp*cosphi,pperp*sinphi,prz[it]));
		complex<double> wavemin_out=cross.getDeutwf()->DeuteronPState(-2, spini_out, spinr, TVector3(pperp*cosphi,pperp*sinphi,prz[it]));
		complex<double> wavezero_out=cross.getDeutwf()->DeuteronPState(0, spini_out, spinr, TVector3(pperp*cosphi,pperp*sinphi,prz[it]));
		complex<double> waveplus_out=cross.getDeutwf()->DeuteronPState(2, spini_out, spinr, TVector3(pperp*cosphi,pperp*sinphi,prz[it]));
		complex<double> wavemin_out_crossed=cross.getDeutwf()->DeuteronPState(-2, spini_out, spinn, TVector3(-pperp*cosphi,-pperp*sinphi,-prz[it]+qvec));
		complex<double> wavezero_out_crossed=cross.getDeutwf()->DeuteronPState(0, spini_out, spinn, TVector3(-pperp*cosphi,-pperp*sinphi,-prz[it]+qvec));
		complex<double> waveplus_out_crossed=cross.getDeutwf()->DeuteronPState(2, spini_out, spinn, TVector3(-pperp*cosphi,-pperp*sinphi,-prz[it]+qvec));
		  
		complex<double> dens0=wavemin*conj(wavemin_out)+waveplus*conj(waveplus_out)+wavezero*conj(wavezero_out);
		complex<double> dens1=wavemin*conj(wavemin_out)+waveplus*conj(waveplus_out)-2.*wavezero*conj(wavezero_out);
		complex<double> dens2=wavemin*conj(wavezero_out)+wavezero*conj(wavemin_out)-wavezero*conj(waveplus_out)-waveplus*conj(wavezero_out);
		complex<double> dens3=wavemin*conj(waveplus_out)+waveplus*conj(wavemin_out);
		complex<double> dens0_crossed=wavemin*conj(wavemin_out_crossed)+waveplus*conj(waveplus_out_crossed)+wavezero*conj(wavezero_out_crossed);
		complex<double> dens1_crossed=wavemin*conj(wavemin_out_crossed)+waveplus*conj(waveplus_out_crossed)-2.*wavezero*conj(wavezero_out_crossed);
		complex<double> dens2_crossed=wavemin*conj(wavezero_out_crossed)+wavezero*conj(wavemin_out_crossed)-wavezero*conj(waveplus_out_crossed)-waveplus*conj(wavezero_out_crossed);
		complex<double> dens3_crossed=wavemin*conj(waveplus_out_crossed)+waveplus*conj(wavemin_out_crossed);
		
		double Q2overkk=Q2/qvec/qvec;
		complex<double> TandL=  pow(Q2overkk,2.)*currentin0*currentout0
				    +(Q2overkk/2.+tanhalfth2)*(currentinmin*currentoutmin+currentinplus*currentoutplus);
		complex<double> TandL_crossed=  pow(Q2overkk,2.)*currentin0*currentout0prime
				    +(Q2overkk/2.+tanhalfth2)*(currentinmin*currentoutminprime+currentinplus*currentoutplusprime);
		complex<double> LT=-Q2overkk*sqrt((tanhalfth2+Q2overkk)/2.)*(currentin0*conj(currentoutplus-currentoutmin)+
					(currentinplus-currentinmin)*conj(currentout0));
		complex<double> LT_crossed=-Q2overkk*sqrt((tanhalfth2+Q2overkk)/2.)*(currentin0*conj(currentoutplusprime-currentoutminprime)+
					(currentinplus-currentinmin)*conj(currentout0prime));
		complex<double> TT=-Q2overkk/2.*(currentinmin*conj(currentoutplus)+currentinplus*conj(currentoutmin));
		complex<double> TT_crossed=-Q2overkk/2.*(currentinmin*conj(currentoutplusprime)+currentinplus*conj(currentoutminprime));
		//direct
		res0+=real(TandL*dens0);
		//cross
		if(Erprimenorm<MASSD-200.) res1+=real(TandL_crossed*dens0_crossed);
		res2+=real(TandL*dens1*(1.+3.*cos(2.*thetapol))/4.
		      +(3./(4.*sqrt(2.))*sin(2.*thetapol)*dens2*LT)
		+3./8.*(1.-cos(2.*thetapol))*dens3*TT);
		if(Erprimenorm<MASSD-200.) res3+=real(TandL_crossed*dens1_crossed*(1.+3.*cos(2.*thetapol))/4.
		      +(3./(4.*sqrt(2.))*sin(2.*thetapol)*dens2_crossed*LT_crossed)
		+3./8.*(1.-cos(2.*thetapol))*dens3_crossed*TT_crossed);
									  
	      }
	    }
	  }
	}
	//factors of 2 because of parity symmetry exploited
	res0*=2.*pperp*MASSD/(2.*(MASSD-Ernorm))/2./abs(qvec-prz[it]/Ernorm*(MASSD+nu));
	res2*=2.*pperp*MASSD/(2.*(MASSD-Ernorm))/2./abs(qvec-prz[it]/Ernorm*(MASSD+nu));
	//minus sign because other nucleon is on-shell now in the wave function for cross term
	if(Erprimenorm<MASSD-200.) {
	  res1*=-2.*pperp*MASSD/(2.*sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm)))/2./abs(qvec-prz[it]/Ernorm*(MASSD+nu))
	  *sqrt(Erprimenorm/Ernorm);
	  res3*=-2.*pperp*MASSD/(2.*sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm)))/2./abs(qvec-prz[it]/Ernorm*(MASSD+nu))
	  *sqrt(Erprimenorm/Ernorm);
	}
	result[0]+=res0;
	result[1]+=res1;
	result[2]+=res2;
	result[3]+=res3;
      }

    }
  }
  

  return;
  

}


void DQEinclusive::calc_CrossincFSI(vector<double> &fsi, double Q2,double x, 
				  int current, int integrator, int maxEval, bool nopt, double thetapol){
  
  double nu=Q2/(x*2.*MASSP);
  ffactorseq=new NucleonEMOperator(Q2,proton,ffparam);
  ffactorsdiff=new NucleonEMOperator(Q2,!proton,ffparam);

  
  for(int i=0;i<8;i++) fsi[i]=0.;
  minpcm=1.E03;
  numint::array<double,3> lower = {{0.,0.,0.}};
  numint::array<double,3> upper = {{0.7E03,0.7E03,2.*PI}};
  DQEinclusive::Ftor_FSI F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.current = current;
  F.tanhalfth2=electron.GetTan2HalfAngle(Q2,nu);
  F.nopt=nopt;
  F.thetapol=thetapol;
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
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E02,maxEval,ret,ccount,0);
      break;
    }
//     case(1):{ //cuhre cubature from cuba lib
//       int nvec=1;
//       int flags=0x00;
//       int key=11;
//       char statefile[] = {""};
//       int nregions =0;
//       numint::cuhre(mdf,lower,upper,nvec,1.E-08,PREC,flags,1E03,maxEval,key,statefile,nregions,count,fail,ret,error,prob);
//       break;
//     }
//     case(2):{ //suave MC from cuba lib
//       int nvec=1;
//       int flags = 0x00;
//       int seed = rand();
//       int nnew = 5000;
//       double flatness = 0.25;
//       char statefile[] = {""};
//       
// //       **** OUTPUT VARIABLES ****** //
//       int nregions    = 0;
//       numint::suave(mdf,lower,upper,nvec,1.E-08,PREC,flags,seed,1E03,maxEval,nnew,flatness,
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
// //       **** OUTPUT VARIABLES ****** //
//       numint::vegas(mdf,lower,upper,nvec,1.E-08,PREC,flags,seed,
// 				1E03,maxEval,nstart,nincrease,nbatch,gridno,
// 				statefile,count,fail,ret,error,prob);
//       break;
//     }
//     case(4):{ //divonne MC from cuba lib
// // 	**** INPUT VARIABLES ****** //
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
// // 	**** OUTPUT VARIABLES ****** //
// 	int nregions=0;
// 	numint::divonne(mdf,lower,upper,nvec,
// 		1.E-05,PREC,flags,seed,
// 		1E03,maxEval,key1,key2,key3,
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
//     cout << "fsi " << res << " " << count << endl;
  double sigmamott=ALPHA*ALPHA*(1+electron.GetCosScatterAngle(Q2,nu))/
    (2.*pow(electron.GetBeamEnergy(Q2,nu)*(1-electron.GetCosScatterAngle(Q2,nu)),2.))*HBARC*HBARC*1.E07;
  for(int i=0;i<8;i++) fsi[i]= 2.*PI/3./MASSD*1.E03/(64.*PI*PI)*ret[i]*sigmamott;

  
  delete ffactorseq;
  delete ffactorsdiff;
  
  return;
}

void DQEinclusive::FSI_int(numint::vector_d & result, double pperp1, double qt, 
		      double qphi, DQEinclusive& cross, double Q2, double x, int current, 
		      double tanhalfth2, bool nopt, double thetapol){
			
  double phi1=0.;
  result=numint::vector_d(8,0.);
  complex<double> res0=0.,res1=0.,res2=0.,res3=0., res4=0.,res5=0.,res6=0.,res7=0.;
  if(pperp1<1.E-03||qt<1.E-03) return;

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



  if(cross.get_prz(prz1,nopt?0.:pperp1,Q2,nu,qvec)){ //prz evaluation is performed!!
    if(cross.get_prz(prz2,nopt?0.:pperp2,Q2,nu,qvec)){ //prz evaluation is performed!!
      for(size_t it1=0;it1<prz1.size();it1++){      
	double prnorm=sqrt(pperp1*pperp1+prz1[it1]*prz1[it1]);
  //       double costheta=prz[it]/prnorm;
	double Ernorm=sqrt(pperp1*pperp1+prz1[it1]*prz1[it1]+cross.getMassr()*cross.getMassr());
	if(Ernorm<MASSD-200.){
// 	cout << "1 " << it1 << " " << pperp1 << " " <<  prz1[it1] << " " << Ernorm << " " << 2.*qvec*prz1[it1]-2.*(MASSD+nu)*Ernorm << " " << -MASSD*MASSD+Q2-2.*MASSD*nu -cross.getMassr()*cross.getMassr()+ cross.getMassi()*cross.getMassi()<< endl;
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
	    if(Erprimenorm<MASSD-200.){
// 	    cout << "2 " << it2 << " " << pperp2 << " " <<  prz2[it2] << " " << Erprimenorm << " " << 2.*qvec*prz2[it2]-2.*(MASSD+nu)*Erprimenorm << " " << -MASSD*MASSD+Q2-2.*MASSD*nu  -cross.getMassr()*cross.getMassr()+ cross.getMassi()*cross.getMassi()<< endl;

	      FourVector<double> pi2(MASSD-Erprimenorm,-pperp1*cosphi1-qt*cosqphi,-qt*sinqphi,-prz2[it2]); //pi2=pD-ps2, in perp dir ps2=ps1+qt
	      FourVector<double> pn2=pi2+q; //pn2=pi2+q
	      double t=pow(Ernorm-Erprimenorm,2.)-pow(prz1[it1]-prz2[it2],2.)-qt*qt;
	      //double t_cross2=pow(pn1[0]-Erprimenorm,2.)-pow(pn1[3]-prz2[it2],2.)-qt*qt;
	      complex<double> directrescatt = rescatter.scatter(t,0);

	      //direct term
	      GammaStructure J0prime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVector0;
	      GammaStructure Jplusprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorPlus;
	      GammaStructure Jminprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorMin;	     
	      
	      Matrix<1,4> un2_down=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
	      Matrix<1,4> un2_up=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
	      TSpinor ui2_down(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
	      TSpinor ui2_up(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
	      
	      //crossed term
	      //pr2 remains the first spectator after rescattering
	      //it's on shell through pr2^0 integration (before sending out the photon)
	      //prz pole puts the nucleon in the upper leg (that absorbed the photon) after the rescattering on-shell
// 	      FourVector<double> pr2_cross(0.,pn1-qt*cosqphi,pn1-qt*sinqhi,-prz
	      FourVector<double> pi2_cross(Erprimenorm-nu,pperp1*cosphi1+qt*cosqphi,qt*sinqphi,prz2[it2]-qvec); //pi2_cross = pr2 - qvec
	      FourVector<double> pn2_cross=pi2_cross+q; //pn2=pi2+q
// 	      double t_cross=pow(pn2_cross[0]-Ernorm,2.)-pow(pn2_cross[3]-prz1[it1],2.)-qt*qt;
// 	      complex<double> crossrescatt = rescatter.scatter(t_cross,0);

	      GammaStructure J0prime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVector0;
	      GammaStructure Jplusprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorPlus;
	      GammaStructure Jminprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorMin;	     
	      
	      Matrix<1,4> un2_down_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
	      Matrix<1,4> un2_up_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
	      TSpinor ui2_down_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
	      TSpinor ui2_up_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
	      
	      res0=res1=res2=res3=res4=res5=res6=res7=0.;
	      //summation over all the spin indices
	      for(int spinn=-1;spinn<=0;spinn+=2){ //exploiting parity symmetry here (see factor 2 below)
		for(int spini_in=-1;spini_in<=1;spini_in+=2){
		  complex<double> currentin0=(spinn==-1?un1_down:un1_up)*J0*(spini_in==-1?ui1_down:ui1_up);
		  complex<double> currentinplus=(spinn==-1?un1_down:un1_up)*Jplus*(spini_in==-1?ui1_down:ui1_up);
		  complex<double> currentinmin=(spinn==-1?un1_down:un1_up)*Jmin*(spini_in==-1?ui1_down:ui1_up);
		  for(int spini_out=-1;spini_out<=1;spini_out+=2){
		    //direct term
		    complex<double> currentout0=conj((spinn==-1?un2_down:un2_up)*J0prime*(spini_out==-1?ui2_down:ui2_up));
		    complex<double> currentoutplus=conj((spinn==-1?un2_down:un2_up)*Jplusprime*(spini_out==-1?ui2_down:ui2_up));
		    complex<double> currentoutmin=conj((spinn==-1?un2_down:un2_up)*Jminprime*(spini_out==-1?ui2_down:ui2_up));
		    for(int spinr=-1;spinr<=1;spinr+=2){
		  //cross term
		      complex<double> currentout0_cross=conj((spinr==-1?un2_down_cross:un2_up_cross)*J0prime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
		      complex<double> currentoutplus_cross=conj((spinr==-1?un2_down_cross:un2_up_cross)*Jplusprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
		      complex<double> currentoutmin_cross=conj((spinr==-1?un2_down_cross:un2_up_cross)*Jminprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));

		      complex<double> wavemin=cross.getDeutwf()->DeuteronPState(-2, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]));
		      complex<double> wavezero=cross.getDeutwf()->DeuteronPState(0, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]));
		      complex<double> waveplus=cross.getDeutwf()->DeuteronPState(2, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]));
		      
		      complex<double> wavemin_out=cross.getDeutwf()->DeuteronPState(-2, spini_out, spinr, TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2]));
		      complex<double> wavezero_out=cross.getDeutwf()->DeuteronPState(0, spini_out, spinr, TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2]));
		      complex<double> waveplus_out=cross.getDeutwf()->DeuteronPState(2, spini_out, spinr, TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2]));

		      complex<double> wavemin_out_crossed=cross.getDeutwf()->DeuteronPState(-2, spini_out, spinn, TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-prz2[it2]+qvec));
		      complex<double> wavezero_out_crossed=cross.getDeutwf()->DeuteronPState(0, spini_out, spinn, TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-prz2[it2]+qvec));
		      complex<double> waveplus_out_crossed=cross.getDeutwf()->DeuteronPState(2, spini_out, spinn, TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-prz2[it2]+qvec));

		      complex<double> wavemin_off=cross.getDeutwf()->DeuteronPStateOff(-2, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]));
		      complex<double> wavezero_off=cross.getDeutwf()->DeuteronPStateOff(0, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]));
		      complex<double> waveplus_off=cross.getDeutwf()->DeuteronPStateOff(2, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]));
		      
		      complex<double> wavemin_out_off=cross.getDeutwf()->DeuteronPStateOff(-2, spini_out, spinr, TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2]));
		      complex<double> wavezero_out_off=cross.getDeutwf()->DeuteronPStateOff(0, spini_out, spinr, TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2]));
		      complex<double> waveplus_out_off=cross.getDeutwf()->DeuteronPStateOff(2, spini_out, spinr, TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2]));

		      complex<double> wavemin_out_crossed_off=cross.getDeutwf()->DeuteronPStateOff(-2, spini_out, spinn, TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-prz2[it2]+qvec));
		      complex<double> wavezero_out_crossed_off=cross.getDeutwf()->DeuteronPStateOff(0, spini_out, spinn, TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-prz2[it2]+qvec));
		      complex<double> waveplus_out_crossed_off=cross.getDeutwf()->DeuteronPStateOff(2, spini_out, spinn, TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-prz2[it2]+qvec));

		      complex<double> dens0=wavemin*conj(wavemin_out)+waveplus*conj(waveplus_out)+wavezero*conj(wavezero_out);
		      complex<double> dens1=wavemin*conj(wavemin_out)+waveplus*conj(waveplus_out)-2.*wavezero*conj(wavezero_out);
		      complex<double> dens2=wavemin*conj(wavezero_out)+wavezero*conj(wavemin_out)-wavezero*conj(waveplus_out)-waveplus*conj(wavezero_out);
		      complex<double> dens3=wavemin*conj(waveplus_out)+waveplus*conj(wavemin_out);
		      complex<double> dens0_crossed=wavemin*conj(wavemin_out_crossed)+waveplus*conj(waveplus_out_crossed)+wavezero*conj(wavezero_out_crossed);
		      complex<double> dens1_crossed=wavemin*conj(wavemin_out_crossed)+waveplus*conj(waveplus_out_crossed)-2.*wavezero*conj(wavezero_out_crossed);
		      complex<double> dens2_crossed=wavemin*conj(wavezero_out_crossed)+wavezero*conj(wavemin_out_crossed)-wavezero*conj(waveplus_out_crossed)-waveplus*conj(wavezero_out_crossed);
		      complex<double> dens3_crossed=wavemin*conj(waveplus_out_crossed)+waveplus*conj(wavemin_out_crossed);
		      
		      complex<double> dens0_off=wavemin_off*conj(wavemin_out_off)+waveplus_off*conj(waveplus_out_off)+wavezero_off*conj(wavezero_out_off);
		      complex<double> dens1_off=wavemin_off*conj(wavemin_out_off)+waveplus_off*conj(waveplus_out_off)-2.*wavezero_off*conj(wavezero_out_off);
		      complex<double> dens2_off=wavemin_off*conj(wavezero_out_off)+wavezero_off*conj(wavemin_out_off)-wavezero_off*conj(waveplus_out_off)-waveplus_off*conj(wavezero_out_off);
		      complex<double> dens3_off=wavemin_off*conj(waveplus_out_off)+waveplus_off*conj(wavemin_out_off);
		      complex<double> dens0_crossed_off=wavemin_off*conj(wavemin_out_crossed_off)+waveplus_off*conj(waveplus_out_crossed_off)+wavezero_off*conj(wavezero_out_crossed_off);
		      complex<double> dens1_crossed_off=wavemin_off*conj(wavemin_out_crossed_off)+waveplus_off*conj(waveplus_out_crossed_off)-2.*wavezero_off*conj(wavezero_out_crossed_off);
		      complex<double> dens2_crossed_off=wavemin_off*conj(wavezero_out_crossed_off)+wavezero_off*conj(wavemin_out_crossed_off)-wavezero_off*conj(waveplus_out_crossed_off)-waveplus_off*conj(wavezero_out_crossed_off);
		      complex<double> dens3_crossed_off=wavemin_off*conj(waveplus_out_crossed_off)+waveplus_off*conj(wavemin_out_crossed_off);

		      double Q2overkk=Q2/qvec/qvec;
		      complex<double> TandL=  pow(Q2overkk,2.)*currentin0*currentout0
					  +(Q2overkk/2.+tanhalfth2)*(currentinmin*currentoutmin+currentinplus*currentoutplus);
		      complex<double> TandL_crossed=  pow(Q2overkk,2.)*currentin0*currentout0_cross
					  +(Q2overkk/2.+tanhalfth2)*(currentinmin*currentoutmin_cross+currentinplus*currentoutplus_cross);
		      complex<double> LT=-Q2overkk*sqrt((tanhalfth2+Q2overkk)/2.)*(currentin0*conj(currentoutplus-currentoutmin)+
					      (currentinplus-currentinmin)*conj(currentout0));
		      complex<double> LT_crossed=-Q2overkk*sqrt((tanhalfth2+Q2overkk)/2.)*(currentin0*conj(currentoutplus_cross-currentoutmin_cross)+
					      (currentinplus-currentinmin)*conj(currentout0_cross));
		      complex<double> TT=-Q2overkk/2.*(currentinmin*conj(currentoutplus)+currentinplus*conj(currentoutmin));
		      complex<double> TT_crossed=-Q2overkk/2.*(currentinmin*conj(currentoutplus_cross)+currentinplus*conj(currentoutmin_cross));
		      
		      res0+=TandL*dens0*directrescatt;
		      res1+=TandL_crossed*dens0_crossed*directrescatt;
		      
		      res2+=TandL*dens0_off*directrescatt;
		      res3+=TandL_crossed*dens0_crossed_off*directrescatt;
		      
		      res4+=(TandL*dens1*(1.+3.*cos(2.*thetapol))/4.
			+(3./(4.*sqrt(2.))*sin(2.*thetapol)*dens2*LT)
			  +3./8.*(1.-cos(2.*thetapol))*dens3*TT)*directrescatt;
		      res5+=(TandL_crossed*dens1_crossed*(1.+3.*cos(2.*thetapol))/4.
			    +(3./(4.*sqrt(2.))*sin(2.*thetapol)*dens2_crossed*LT_crossed)
			    +3./8.*(1.-cos(2.*thetapol))*dens3_crossed*TT_crossed)*directrescatt;
		      res6+=(TandL*dens1_off*(1.+3.*cos(2.*thetapol))/4.
			+(3./(4.*sqrt(2.))*sin(2.*thetapol)*dens2_off*LT)
			  +3./8.*(1.-cos(2.*thetapol))*dens3_off*TT)*directrescatt;
		      res7+=(TandL_crossed*dens1_crossed_off*(1.+3.*cos(2.*thetapol))/4.
			    +(3./(4.*sqrt(2.))*sin(2.*thetapol)*dens2_crossed_off*LT_crossed)
			    +3./8.*(1.-cos(2.*thetapol))*dens3_crossed_off*TT_crossed)*directrescatt;
		      
// 		      complex<double> ares0=0.,ares1=0.,ares2=0.,ares3=0.;		      
// 		      for(int M=-2;M<=2;M+=2){
// 			complex<double> wf=/*(M==0? -2.:1.)**/cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
// 					*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinr, 
// 						TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2])));
// 			complex<double> wf2=/*(M==0? -2.:1.)**/cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
// 					*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinn, 
// 						TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-prz2[it2]+qvec)));
// 			complex<double> wf_off=/*(M==0? -2.:1.)**/cross.getDeutwf()->DeuteronPStateOff(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
// 					*conj(cross.getDeutwf()->DeuteronPStateOff(M, spini_out, spinr, 
// 						TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,prz2[it2])));
// 			complex<double> wf2_off=/*(M==0? -2.:1.)**/cross.getDeutwf()->DeuteronPStateOff(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,prz1[it1]))
// 					*conj(cross.getDeutwf()->DeuteronPStateOff(M, spini_out, spinn, 
// 						TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-prz2[it2]+qvec)));
// 			//direct onshell
// 			ares0+=(Q2*Q2/pow(qvec,4.)*currentin0*currentout0+(Q2/(2.*qvec*qvec)+tanhalfth2)
// 			  *(currentinmin*currentoutmin+currentinplus*currentoutplus))*wf*directrescatt;
// 			//cross onshell
// 			ares1+=(Q2*Q2/pow(qvec,4.)*currentin0*currentout0_cross+(Q2/(2.*qvec*qvec)+tanhalfth2)*
// 			  (currentinmin*currentoutmin_cross+currentinplus*currentoutplus_cross))*wf2*directrescatt;
// 			//direct offshell
// 			ares2+=(Q2*Q2/pow(qvec,4.)*currentin0*currentout0+(Q2/(2.*qvec*qvec)+tanhalfth2)
// 			  *(currentinmin*currentoutmin+currentinplus*currentoutplus))*wf_off*directrescatt;
// 			//cross offshell
// 			ares3+=(Q2*Q2/pow(qvec,4.)*currentin0*currentout0_cross+(Q2/(2.*qvec*qvec)+tanhalfth2)*
// 			  (currentinmin*currentoutmin_cross+currentinplus*currentoutplus_cross))*wf2_off*directrescatt;
// 		      }
// 		      cout << res0 << " " << res1 << " " << res2 << " " << res3 << endl;
// 		      cout << ares0 << " " << ares1 << " " << ares2 << " " << ares3 << endl;

		    }
		  }
		}
	      }	  
	      //factor 2 because of exploiting parity symmetry
	      res0*=-2.*pperp1*qt/sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))
	      /abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
	      res2*=2.*pperp1*qt/sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))
	      /abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
	      //minus sign because of wave function!!
	      res1*=2.*pperp1*qt*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm))/sqrt(Ernorm*Erprimenorm)
	      /abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))/abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
	      res3*=-2.*pperp1*qt*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm))/sqrt(Ernorm*Erprimenorm)
	      /abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))/abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;

	      res4*=-2.*pperp1*qt/sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))
	      /abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
	      res5*=2.*pperp1*qt/sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))
	      /abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
	      //minus sign because of wave function!!
	      res6*=2.*pperp1*qt*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm))/sqrt(Ernorm*Erprimenorm)
	      /abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))/abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
	      res7*=-2.*pperp1*qt*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm))/sqrt(Ernorm*Erprimenorm)
	      /abs(qvec-prz1[it1]/Ernorm*(MASSD+nu))/abs(qvec-prz2[it2]/Erprimenorm*(MASSD+nu))*chi;
	      
	      result[0]+=imag(res0);
	      result[1]+=imag(res1);
	      result[2]+=imag(res2);
	      result[3]+=imag(res3);
	      result[4]+=imag(res4);
	      result[5]+=imag(res5);
	      result[6]+=imag(res6);
	      result[7]+=imag(res7);
//  	      cout << "BLAAAAAAAAA " << it1<< " " << it2 << " " << result[0] << " " << res0 << endl;

	    }
	  }
	}
      }

    }
  }
  
  
  return;
  

}


void DQEinclusive::calc_CrossincFSI_PVoff(double &fsi1_off, double &fsi2_off, double Q2,double x, 
				  int current, int integrator, int maxEval){
  
  double nu=Q2/(x*2.*MASSP);
  ffactorseq=new NucleonEMOperator(Q2,proton,ffparam);
  ffactorsdiff=new NucleonEMOperator(Q2,!proton,ffparam);
  //first we perform the regular 4d integration, inside we do the PV integrations.
  
  minpcm=1.E03;
  numint::array<double,3> lower = {{0.,0.,0.}};
  numint::array<double,3> upper = {{0.7E03,0.7E03,2.*PI}};
  DQEinclusive::Ftor_FSI F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.current = current;
  F.tanhalfth2=electron.GetTan2HalfAngle(Q2,nu);
  F.nopt=1;
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
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E01,maxEval,ret,ccount,0);
      break;
    }
//     case(1):{ //cuhre cubature from cuba lib
//       int nvec=1;
//       int flags=0x00;
//       int key=11;
//       char statefile[] = {""};
//       int nregions =0;
//       numint::cuhre(mdf,lower,upper,nvec,1.E-08,PREC,flags,1E03,2E04,key,statefile,nregions,count,fail,ret,error,prob);
//       break;
//     }
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
//     cout << "fsi " << res << " " << count << endl;
  double sigmamott=ALPHA*ALPHA*(1+electron.GetCosScatterAngle(Q2,nu))/
    (2.*pow(electron.GetBeamEnergy(Q2,nu)*(1-electron.GetCosScatterAngle(Q2,nu)),2.))*HBARC*HBARC*1.E07;
  fsi1_off= 2.*PI/3./MASSD*1.E03/(64.*PI*PI*PI)*ret[0]*sigmamott;
  fsi2_off= 2.*PI/3./MASSD*1.E03/(64.*PI*PI*PI)*ret[1]*sigmamott;

  
  delete ffactorseq;
  delete ffactorsdiff;
  
  return;
}

void DQEinclusive::FSI_PV(numint::vector_d & result, double pperp1, double qt, 
		      double qphi, DQEinclusive& cross, double Q2, double x, int current, double tanhalfth2, bool nopt, double thetapol){
  
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
  pvint.tanhalfth2 = tanhalfth2;
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
  pvint.rescatter=&rescatter;
  pvint.f=DQEinclusive::FSI_intPV;
  
  cross.get_prz(prz1,nopt?0.:pperp1,Q2,nu,qvec);
  cross.get_prz(prz2,nopt?0.:pperp2,Q2,nu,qvec);
  cout << prz1.size() << " " << prz2.size() << " ";
  for(size_t i=0;i<prz1.size();i++) cout << prz1[i] << " ";
  for(size_t i=0;i<prz2.size();i++) cout << prz2[i] << " ";
  cout << (-MASSD*MASSD+Q2-2.*MASSD*nu+2.*(MASSD+nu)*cross.getMassr())/2./qvec << endl;
  
  if(cross.get_prz(prz1,nopt?0.:pperp1,Q2,nu,qvec)){ //prz evaluation is performed!!
    if(cross.get_prz(prz2,nopt?0.:pperp2,Q2,nu,qvec)){ //prz evaluation is performed!!
      pvint.prz2poles = prz2;
      for(size_t it1=0;it1<prz1.size();it1++){      
	for(int crossed=1;crossed<=1;crossed++){
	  pvint.crossed=bool(crossed);
	  gsl_integration_workspace * w 
	  = gsl_integration_workspace_alloc (1E6);
	  double res, error;
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
  cout << pperp1 << " " << qt << " " << qphi << " " << result[0] << " " << result[1] << endl;
}

double DQEinclusive::PV_int1(double pz1, void * params){
  

  Ftor_PV pvint = *(Ftor_PV *) params;
  pvint.var1=pz1;
  double res=0.;
  for(size_t it2=0;it2<pvint.prz2poles.size();it2++){ 
    gsl_integration_workspace * w2 
      = gsl_integration_workspace_alloc (1E6);

    double result, error;
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
			       double sinphi1, double cosqphi, double sinqphi,
			       DQEinclusive& cross, double Q2, double x, double tanhalfth2, double qvec, double nu, int current, bool crossed,
			       FastParticle &rescatter
			      ){
  complex<double>result=0.;

  double prnorm=sqrt(pperp1*pperp1+pz1*pz1);
  double Ernorm=sqrt(pperp1*pperp1+pz1*pz1+cross.getMassr()*cross.getMassr());
  if(Ernorm<MASSD-200.){
    FourVector<double> q(nu,0.,0.,qvec);
    
    //direct term
    FourVector<double> pi1(MASSD-Ernorm,-pperp1*cosphi1,-pperp1*sinphi1,-pz1); //pi=pD-ps
    FourVector<double> pn1=pi1+q;  //pn1=pi+q
    GammaStructure J0 = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVector0;
    GammaStructure Jplus = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVectorPlus;
    GammaStructure Jmin = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVectorMin;
    
    Matrix<1,4> un1_down=TSpinor::Bar(TSpinor(pn1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
    Matrix<1,4> un1_up=TSpinor::Bar(TSpinor(pn1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
    TSpinor ui1_down(pi1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
    TSpinor ui1_up(pi1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
  
    double Erprimenorm=sqrt(pperp2*pperp2+pz2*pz2+cross.getMassr()*cross.getMassr());
    if(Erprimenorm<MASSD-200.){

      FourVector<double> pi2(MASSD-Erprimenorm,-pperp1*cosphi1-qt*cosqphi,-qt*sinqphi,-pz2); //pi2=pD-ps1+qt
      FourVector<double> pn2=pi2+q; //pn2=pi2+q
      if(crossed) {pi2[1]*=-1.;pi2[2]*=-1.;pn2[1]*=-1.;pn2[2]*=-1.;} //transverse components opposite for crossed term
      double t=pow(Ernorm-Erprimenorm,2.)-pow(pz1-pz2,2.)-qt*qt;

      complex<double> rescatt = rescatter.scatter(t,0);

      //direct term
      GammaStructure J0prime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVector0;
      GammaStructure Jplusprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorPlus;
      GammaStructure Jminprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorMin;	     
      
      Matrix<1,4> un2_down=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
      Matrix<1,4> un2_up=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
      TSpinor ui2_down(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
      TSpinor ui2_up(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);

      FourVector<double> pi2_cross(Erprimenorm-nu,pperp1*cosphi1+qt*cosqphi,qt*sinqphi,pz2-qvec); //pi2_cross = pr2 - qvec
      FourVector<double> pn2_cross=pi2_cross+q; //pn2=pi2+q (or pr2)
// 	      double t_cross=pow(pn2_cross[0]-Ernorm,2.)-pow(pn2_cross[3]-prz1[it1],2.)-qt*qt;
// 	      complex<double> crossrescatt = rescatter.scatter(t_cross,0);

      GammaStructure J0prime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVector0;
      GammaStructure Jplusprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorPlus;
      GammaStructure Jminprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorMin;	     
      
      Matrix<1,4> un2_down_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
      Matrix<1,4> un2_up_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
      TSpinor ui2_down_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
      TSpinor ui2_up_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
      
      //summation over all the spin indices
      //check and exploit parity symmetry!!!!
      for(int spinn=-1;spinn<=0;spinn+=2){
	for(int spini_in=-1;spini_in<=1;spini_in+=2){
	  complex<double> currentin0=(spinn==-1?un1_down:un1_up)*J0*(spini_in==-1?ui1_down:ui1_up);
	  complex<double> currentinplus=(spinn==-1?un1_down:un1_up)*Jplus*(spini_in==-1?ui1_down:ui1_up);
	  complex<double> currentinmin=(spinn==-1?un1_down:un1_up)*Jmin*(spini_in==-1?ui1_down:ui1_up);
	  for(int spini_out=-1;spini_out<=1;spini_out+=2){
	    for(int spinr=-1;spinr<=1;spinr+=2){
	      complex<double> currentout0=0.,currentoutplus=0.,currentoutmin=0.;
	      if(!crossed){
		currentout0=conj((spinn==-1?un2_down:un2_up)*J0prime*(spini_out==-1?ui2_down:ui2_up));
		currentoutplus=conj((spinn==-1?un2_down:un2_up)*Jplusprime*(spini_out==-1?ui2_down:ui2_up));
		currentoutmin=conj((spinn==-1?un2_down:un2_up)*Jminprime*(spini_out==-1?ui2_down:ui2_up));
	      }
	      else{
		currentout0=conj((spinr==-1?un2_down_cross:un2_up_cross)*J0prime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
		currentoutplus=conj((spinr==-1?un2_down_cross:un2_up_cross)*Jplusprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
		currentoutmin=conj((spinr==-1?un2_down_cross:un2_up_cross)*Jminprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
	      }
	      complex<double> wftemp=0.;
	      for(int M=-2;M<=2;M+=2){
		if(!crossed) wftemp+=cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,pz1))
				*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinr, 
					TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,pz2)));
		else wftemp+=cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,pz1))
				*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinn, 
					TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-pz2+qvec)));
	      }
	      result+=wftemp*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0+(Q2/(2.*qvec*qvec)+tanhalfth2)
		  *(currentinmin*currentoutmin+currentinplus*currentoutplus))*rescatt;
	      
	    }	    
	  }
	}
      }
      //exploited parity symmetry, hence extra factor 2
      if(!crossed) result*=2./sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-pz1/Ernorm*(MASSD+nu))
      /abs(qvec-pz2/Erprimenorm*(MASSD+nu));
      else result*=-2.*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm))/sqrt(Ernorm*Erprimenorm)
      /abs(qvec-pz1/Ernorm*(MASSD+nu))/abs(qvec-pz2/Erprimenorm*(MASSD+nu));
    }
  }

  
  
  return imag(result);
  

}

void DQEinclusive::calc_CrossincFSI_PVoff2(double &fsi1_off, double &fsi2_off, double Q2,double x, 
				  int current, int integrator, int maxEval){
  
  double nu=Q2/(x*2.*MASSP);
  ffactorseq=new NucleonEMOperator(Q2,proton,ffparam);
  ffactorsdiff=new NucleonEMOperator(Q2,!proton,ffparam);
  //first we perform the two PV integrations, inside we do the other (regular) integrations over the perp coordinates
  //no pt dependence in the pole values for the prz now of course! (neglected).
  
  minpcm=1.E03;
  DQEinclusive::Ftor_PVfirst F;
  F.cross=this;
  F.Q2 = Q2;
  F.x = x;
  F.current = current;
  F.tanhalfth2=electron.GetTan2HalfAngle(Q2,nu);
  F.maxEval=maxEval;
  F.integrator=integrator;
  
  F.f=DQEinclusive::FSI_PV_pz2;
  vector<double>prz1;
  
  double qvec=sqrt(Q2+pow(Q2/(2.*massi*x),2.));
  get_prz(prz1,0.,Q2,nu,qvec);
  cout << prz1.size() << " ";
  for(size_t i=0;i<prz1.size();i++) cout << prz1[i] << " ";
  cout << (-MASSD*MASSD+Q2-2.*MASSD*nu+2.*(MASSD+nu)*massr+massi*massi-massr*massr)/2./qvec << endl;

  if(prz1.size()==0){
    cerr<< "implement 0 case"<<endl;
    assert(1==0);
  }
  if(prz1.size()==1)   F.pzpole=prz1[0];
  if(prz1.size()==2) {F.pzpole=prz1[0]<prz1[1]?prz1[0]:prz1[1];}
  if(abs(F.pzpole)>700.){    cerr<< "pole outside boundaries, implement regular integral! "<< F.pzpole << endl;
    assert(1==0);
  }

  F.crossed=0;
  gsl_integration_workspace * w2 
    = gsl_integration_workspace_alloc (1E6);

  double result, error;
  gsl_function F2;
  F2.function = &FSI_PV_pz1;
  F2.params = &F;
  gsl_integration_qawc (&F2, -700., 700.,F.pzpole, 1.e-08, PREC, 1E3,
			w2, &result, &error); 
  double sigmamott=ALPHA*ALPHA*(1+electron.GetCosScatterAngle(Q2,nu))/
    (2.*pow(electron.GetBeamEnergy(Q2,nu)*(1-electron.GetCosScatterAngle(Q2,nu)),2.))*HBARC*HBARC*1.E07;
  fsi1_off= 2.*PI/3./MASSD*1.E03/(64.*PI*PI*PI)*result*sigmamott;
  F.crossed=1;
  gsl_integration_qawc (&F2, -700., 700.,F.pzpole, 1.e-08, PREC, 1E3,
			w2, &result, &error); 
  fsi2_off= 2.*PI/3./MASSD*1.E03/(64.*PI*PI*PI)*result*sigmamott;

  gsl_integration_workspace_free (w2);
  
  
  //     cout << "fsi " << res << " " << count << endl;

  
  delete ffactorseq;
  delete ffactorsdiff;
  
  return;
}


double DQEinclusive::FSI_PV_pz1(double pz1, void * params){
  Ftor_PVfirst pvint = *(Ftor_PVfirst *) params;
  pvint.var1=pz1;
  double res=0.;
  gsl_integration_workspace * w2 
    = gsl_integration_workspace_alloc (1E6);

  double result, error;
  gsl_function F2;
  F2.function = &Ftor_PVfirst::exec;
  F2.params = &pvint;
  gsl_integration_qawc (&F2, -700., 700.,pvint.pzpole , 1.e-08, PREC, 1E3,
			w2, &result, &error); 
  res=result;
  gsl_integration_workspace_free (w2);
  return res;

}

double DQEinclusive::FSI_PV_pz2(double pz1, double pz2, DQEinclusive& cross, 
	      double Q2, double x, double tanhalfth2, int current, bool crossed, int maxEval, int integrator){
  
  numint::array<double,3> lower = {{0.,0.,0.}};
  numint::array<double,3> upper = {{0.7E03,0.7E03,2.*PI}};
  DQEinclusive::Ftor_PVfirst_inner F;
  F.pz1=pz1;
  F.pz2=pz2;
  F.cross= &cross;
  F.Q2 = Q2;
  F.x = x;
  F.current = current;
  F.tanhalfth2=tanhalfth2;
  F.crossed=crossed;
  F.maxEval=maxEval;
  
  numint::mdfunction<numint::vector_d,3> mdf;
  mdf.func = &Ftor_PVfirst_inner::exec;
  mdf.param = &F;
  F.f=DQEinclusive::FSI_PVfirst_perp;
  int res=90;
  unsigned ccount=0;
  numint::vector_d ret = numint::vector_d(1,0.);
  switch(integrator){
    case(0):{  //numint adaptive cubature from MIT lib
      res = numint::cube_adaptive(mdf,lower,upper,1.E-08,PREC,1E01,maxEval,ret,ccount,0);
      break;
    }
    default:{
      cerr << "integrator not supported" << endl;
      assert(1==0);
    }
  }
//   cout << pz1 << " " << pz2 << " " << ret[0] << endl;
//   return ret[0];
}

  
void DQEinclusive::FSI_PVfirst_perp(numint::vector_d & res, double pz1, double pz2, double pperp1, double qt, double qphi, DQEinclusive& cross, 
	    double Q2, double x, int current, double tanhalfth2, bool crossed){
  res= numint::vector_d(1,0.);
  double phi1=0.;
  if(pperp1<1.E-03||qt<1.E-03) { res[0]=0.; return;}
  double cosphi1, sinphi1;
  sincos(phi1,&sinphi1,&cosphi1);
  double cosqphi, sinqphi;
  sincos(qphi,&sinqphi,&cosqphi);
  double qvec=sqrt(Q2+pow(Q2/(2.*cross.getMassi()*x),2.));
  double nu=Q2/(2.*cross.getMassi()*x);
  double s=(MASSD+nu)*(MASSD+nu)-qvec*qvec;
  double pcm=sqrt(s/4.-MASSP*MASSP);
  if(pcm<cross.getMinpcm()) cross.minpcm=pcm;
  double chi=sqrt(s)*pcm*2.;
  FastParticle rescatter(0,0,pcm,0.,0.,Q2,0.,"");

  double pperp2=sqrt(pow(pperp1+qt*cosqphi,2.)+pow(qt*sinqphi,2.));

  complex<double>result=0.;

  double prnorm=sqrt(pperp1*pperp1+pz1*pz1);
  double Ernorm=sqrt(pperp1*pperp1+pz1*pz1+cross.getMassr()*cross.getMassr());
  if(Ernorm<MASSD-200.){
    FourVector<double> q(nu,0.,0.,qvec);
    
    //direct term
    FourVector<double> pi1(MASSD-Ernorm,-pperp1*cosphi1,-pperp1*sinphi1,-pz1); //pi=pD-ps
    FourVector<double> pn1=pi1+q;  //pn1=pi+q
    GammaStructure J0 = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVector0;
    GammaStructure Jplus = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVectorPlus;
    GammaStructure Jmin = cross.getFFactorseq()->getCC(current, q, pi1, pn1)*polVectorMin;
    
    Matrix<1,4> un1_down=TSpinor::Bar(TSpinor(pn1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
    Matrix<1,4> un1_up=TSpinor::Bar(TSpinor(pn1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
    TSpinor ui1_down(pi1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
    TSpinor ui1_up(pi1,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
  
    double Erprimenorm=sqrt(pperp2*pperp2+pz2*pz2+cross.getMassr()*cross.getMassr());
    if(Erprimenorm<MASSD-200.){

      FourVector<double> pi2(MASSD-Erprimenorm,-pperp1*cosphi1-qt*cosqphi,-qt*sinqphi,-pz2); //pi2=pD-ps1+qt
      FourVector<double> pn2=pi2+q; //pn2=pi2+q
      double t=pow(Ernorm-Erprimenorm,2.)-pow(pz1-pz2,2.)-qt*qt;

      complex<double> rescatt = rescatter.scatter(t,0);

      //direct term
      GammaStructure J0prime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVector0;
      GammaStructure Jplusprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorPlus;
      GammaStructure Jminprime = cross.getFFactorseq()->getCC(current, q, pi2, pn2)*polVectorMin;	     
      
      Matrix<1,4> un2_down=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
      Matrix<1,4> un2_up=TSpinor::Bar(TSpinor(pn2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
      TSpinor ui2_down(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
      TSpinor ui2_up(pi2,cross.getMassi(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);

      FourVector<double> pi2_cross(Erprimenorm-nu,pperp1*cosphi1+qt*cosqphi,qt*sinqphi,pz2-qvec); //pi2_cross = pr2 - qvec
      FourVector<double> pn2_cross=pi2_cross+q; //pn2=pi2+q
// 	      double t_cross=pow(pn2_cross[0]-Ernorm,2.)-pow(pn2_cross[3]-prz1[it1],2.)-qt*qt;
// 	      complex<double> crossrescatt = rescatter.scatter(t_cross,0);

      GammaStructure J0prime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVector0;
      GammaStructure Jplusprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorPlus;
      GammaStructure Jminprime_cross = cross.getFFactorsdiff()->getCC(current, q, pi2_cross,pn2_cross)*polVectorMin;	     
      
      Matrix<1,4> un2_down_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass));
      Matrix<1,4> un2_up_cross=TSpinor::Bar(TSpinor(pn2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass));
      TSpinor ui2_down_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kDoubleMass);
      TSpinor ui2_up_cross(pi2_cross,cross.getMassr(),TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kDoubleMass);
      
      
      //summation over all the spin indices
      //parity symmetry exploited in first summation (checkeD)
      for(int spinn=-1;spinn<=0;spinn+=2){
	for(int spini_in=-1;spini_in<=1;spini_in+=2){
	  complex<double> currentin0=(spinn==-1?un1_down:un1_up)*J0*(spini_in==-1?ui1_down:ui1_up);
	  complex<double> currentinplus=(spinn==-1?un1_down:un1_up)*Jplus*(spini_in==-1?ui1_down:ui1_up);
	  complex<double> currentinmin=(spinn==-1?un1_down:un1_up)*Jmin*(spini_in==-1?ui1_down:ui1_up);
	  for(int spini_out=-1;spini_out<=1;spini_out+=2){
	    for(int spinr=-1;spinr<=1;spinr+=2){
	      complex<double> currentout0=0.,currentoutplus=0.,currentoutmin=0.;
	      if(!crossed){
		currentout0=conj((spinn==-1?un2_down:un2_up)*J0prime*(spini_out==-1?ui2_down:ui2_up));
		currentoutplus=conj((spinn==-1?un2_down:un2_up)*Jplusprime*(spini_out==-1?ui2_down:ui2_up));
		currentoutmin=conj((spinn==-1?un2_down:un2_up)*Jminprime*(spini_out==-1?ui2_down:ui2_up));
	      }
	      else{
		currentout0=conj((spinr==-1?un2_down_cross:un2_up_cross)*J0prime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
		currentoutplus=conj((spinr==-1?un2_down_cross:un2_up_cross)*Jplusprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
		currentoutmin=conj((spinr==-1?un2_down_cross:un2_up_cross)*Jminprime_cross*(spini_out==-1?ui2_down_cross:ui2_up_cross));
	      }
	      complex<double> wftemp=0.;
	      for(int M=-2;M<=2;M+=2){
		if(!crossed) wftemp+=cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,pz1))
				*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinr, 
					TVector3(pperp1*cosphi1+qt*cosqphi,pperp1*sinphi1+qt*sinqphi,pz2)));
		else wftemp+=cross.getDeutwf()->DeuteronPState(M, spini_in, spinr, TVector3(pperp1*cosphi1,pperp1*sinphi1,pz1))
				*conj(cross.getDeutwf()->DeuteronPState(M, spini_out, spinn, 
					TVector3(-pperp1*cosphi1-qt*cosqphi,-pperp1*sinphi1-qt*sinqphi,-pz2+qvec)));
	      }
	      result+=wftemp*(Q2*Q2/pow(qvec,4.)*currentin0*currentout0+(Q2/(2.*qvec*qvec)+tanhalfth2)
		*(currentinmin*currentoutmin+currentinplus*currentoutplus))*rescatt;
	      
	    }	    
	  }
	}
      }
      //extra factor of 2 because of the parity symmetry used.
      if(!crossed) result*=2./sqrt(Ernorm*Erprimenorm)*MASSD/(2.*(MASSD-Ernorm))/abs(qvec-pz1/Ernorm*(MASSD+nu))
      /abs(qvec-pz2/Erprimenorm*(MASSD+nu));
      else result*=-2.*MASSD/2./sqrt((MASSD-Ernorm)*(MASSD-Erprimenorm))/sqrt(Ernorm*Erprimenorm)
      /abs(qvec-pz1/Ernorm*(MASSD+nu))/abs(qvec-pz2/Erprimenorm*(MASSD+nu));
    }
  }

  res[0]=imag(result);
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
//      cout << std::setprecision(9)<< sol.size() << " " << A+B*p1 << " " << A+B*p2 << " " << Er << " " << Er2 << " " << p1 << " " << p2 << endl;
    return 1;
  }
}
