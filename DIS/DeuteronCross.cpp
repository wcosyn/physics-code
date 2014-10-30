#include "DeuteronCross.hpp"
#include <fstream>
#include <TLorentzRotation.h>
#include <FourVector.h>
#include <NuclStructure.hpp>
#include <TVector3.h>
#include <cassert>

using namespace std;

DeuteronCross::DeuteronCross(string name, bool proton, string struc_name,
			     double sigmain, double betain, double epsilonin, double betaoffin, double lambdain, int offshellset,
			     int looplimit
			    ):
massi(proton? MASSP:MASSN),
momdistr(name,massi,offshellset,sigmain,betain,epsilonin,betaoffin,lambdain,looplimit),
structure(proton,struc_name),
strucname(struc_name)
{
//   for(int i=0;i<200;i+=5){
//     TVector3 p(i,0.,0.);
//     double E=sqrt(p.Mag2()+MASSP*MASSP);
//     double alpha=2.*(E-p[2])/MASSD; //lightcone alpha_s=(E-p_z)/M_n
//     double pt2=p[0]*p[0]+p[1]*p[1];
//     double k=sqrt((M_NUCL*M_NUCL+pt2)/(alpha*(2.-alpha))-M_NUCL*M_NUCL); //lightcone momentum rescaling
//     cout << i << " " << alpha << " " << k << " " << momdistr.getMomDistrpw(p) << " " << momdistr.getLCMomDistrpw(p)/(2.-alpha) << endl;
//       cout << i << " " << alpha << " " << k << " " << pow(MASSD-E,2.)-MASSP*MASSP <<" " << momdistr.getMomDistrpw(p)/pow(momdistr.getDeuteronwf()->getResidu(),2.)*pow(i*i+0.2316*0.2316*HBARC*HBARC,2.) << " " 
// 	<< momdistr.getLCMomDistrpw(p)  << endl;
//   }
//   exit(1);
}


DeuteronCross::~DeuteronCross(){
  
}

DeuteronCross::DeuteronCross(const DeuteronCross& rhs){
  massi=rhs.massi;
  momdistr=rhs.momdistr;
  structure=rhs.structure;
  strucname=rhs.strucname;
}


DeuteronCross& DeuteronCross::operator=(const DeuteronCross& rhs){
  if(this!=&rhs) { // avoid self-assignment
    massi=rhs.massi;
    momdistr=rhs.momdistr;
    structure=rhs.structure;
    strucname=rhs.strucname;
  }
  return *this;

}


double DeuteronCross::getavgBonus(TKinematics2to2 &kin,TElectronKinematics &elec, bool lc){
  
  double dens=/*lc? momdistr.getLCMomDistrpw(kin) :*/ momdistr.getMomDistrpw(kin); 
  FourVector<double> k_in, k_out;
  elec.GetLeptonVectors(kin,k_in,k_out);
  double costhetap=kin.GetCosthklab();
  double sinthetap=sqrt(1.-costhetap*costhetap);
  TLorentzRotation transf(kin.GetPklab()*sinthetap/kin.GetEklab(),0.,kin.GetPklab()*costhetap/kin.GetEklab());
  FourVector<double> k_in_transf=transf*k_in;
  FourVector<double> k_out_transf=transf*k_out;
  double thetae_transf=acos((k_in_transf[1]*k_out_transf[1]+k_in_transf[2]*k_out_transf[2]
		     +k_in_transf[3]*k_out_transf[3])/(k_in_transf[0]*k_out_transf[0])); //Bonus xsection formulas in rest frame of struck nucleon
  double thetae=acos((k_in[1]*k_out[1]+k_in[2]*k_out[2]
		     +k_in[3]*k_out[3])/(k_in[0]*k_out[0])); //Bonus xsection formulas in rest frame of struck nucleon
  FourVector<double> q(kin.GetWlab(),0.,0.,kin.GetKlab());
  FourVector<double> q_transf=transf*q;
  double xtransf=kin.GetQsquared()/(2.*massi*q_transf[0]); //get modified x
  NuclStructure nucl(massi==MASSP?1:0,kin.GetQsquared(),xtransf,0,structure.getName());
  double F1,F2;
  nucl.getF(F1,F2);
  double res=2.*PI*ALPHA*ALPHA*pow(cos(thetae_transf/2.),2.)/(4.*pow(k_in_transf[0],2.)*pow(sin(thetae_transf/2.),4.))*dens
    *(F2+2.*q_transf[0]/massi*F1*pow(tan(thetae_transf/2.),2.))/q_transf[0];
  return F2==0.? 0.: res;
  
}
void DeuteronCross::getBonusextrapolate(double Q2, double W, double Ein, double pr, double costhetar, bool proton, bool lc, 
			 double xref,double norm, double Rdata, double error){
  double massi=proton? MASSP:MASSN;
  double massr=proton? MASSN:MASSP;
  
  double Er=sqrt(massr*massr+pr*pr);
  double Einoff=MASSD-Er;
  double massoff=sqrt(Einoff*Einoff-pr*pr);
  
  double xprime=Q2/(W*W-massoff*massoff+Q2);
  
  //calc nu
  double prz=pr*costhetar;
  double aaa=Einoff*Einoff-prz*prz;
  double bbb=-Einoff*Q2/xprime;
  double ccc=Q2*Q2/(4.*xprime*xprime)-Q2*prz*prz;
  
  double discr=sqrt(bbb*bbb-4.*aaa*ccc);
  double nu1=(-bbb+discr)/(2.*aaa);
  double nu2=(-bbb-discr)/(2.*aaa);
  //cout << nu1 << " " << nu2 << endl;
  double nu=nu2;
  if(costhetar<0.) nu=nu1;
  double qvec=sqrt(Q2+nu*nu);
/*  double xx=Q2/2./(Einoff*nu+prz*qvec);
  cout << xprime << " " << xx << endl;*/
  double x=Q2/(2.*massi*nu);
    
  TKinematics2to2 kin("","",MASSD,MASSP,W,"qsquared:wlab:pklab",Q2,nu,pr);
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
//   DeuteronCross test(*elec,"paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  double alphar=(Er-prz)/(MASSD/2.);

  double Eout=Ein-nu;
  //unphysical kinematics!!!
//   if(Eout<0.){ /*cout << Eout << endl;*/ MCresult=modelresultpw=modelresultfsi=0.;return;}
//   if(std::isnan(asin(sqrt(Q2/(4.*Ein*Eout))))){ /*cout << sqrt(Q2/(4.*Ein*Eout)) << endl; */MCresult=modelresultpw=modelresultfsi=0.;return;}
  
  double MCresult= getavgBonus(kin,*elec,lc);
  double cross_data_extr=Rdata*MCresult/norm; //MeV-6
  double error_extr=error*MCresult/norm; // MeV-6
  double y=kin.GetWlab()/elec->GetBeamEnergy(kin); 
  double front=(4.*PI*ALPHA*ALPHA)/(x*kin.GetQsquared()*kin.GetQsquared())
		*(1-y-(x*x*y*y*massi*massi)/kin.GetQsquared());
  double denspw=momdistr.getMomDistrpw(kin);
  double densfsi=momdistr.getMomDistrfsi(kin,0.);
  double Dstrucs=structure.getavgStructure(kin,*elec,Einoff);
  double strucprefactor=structure.getavgPrefactor(kin,*elec,Einoff);
//   cout << Dstrucs << endl;
//   front*dens*Dstrucs*kin.GetEklab()*1.E18;
//   modelresultpw= getavgCross(kin,*elec,1, Einoff)*2.*Ein*Eout*x/nu/Er;  //go to dEdOmegaed^3ps in MeV-6
//   modelresultfsi= pw? modelresultpw: getavgCross(kin,*elec,0, Einoff)/HBARC/HBARC/10.*2.*Ein*Eout*x/nu/Er;  //go to dEdOmegaed^3ps in MeV-6
// //   double residu=sqrt(2.)*c[0]*sqrt(HBARC)/PI; //[MeV^1/2]
  double corr_residu = pow((massoff*massoff-massi*massi)/momdistr.getDeuteronwf()->getResidu(),2.); // [MeV^3]
  cross_data_extr/=front*2.*Ein*Eout*x/nu*strucprefactor/corr_residu; //dimensionless
  error_extr/=front*2.*Ein*Eout*x/nu*strucprefactor/corr_residu; //dimensionless
  double F2pw=Dstrucs*denspw*corr_residu/strucprefactor; // dimensionless
  double F2fsi=Dstrucs*densfsi*corr_residu/strucprefactor; //dimensionless
  
  NuclStructure strfunction(proton,kin.GetQsquared(),xref,0,strucname);
  double F2ref;
  F2ref=strfunction.getF2();
//   cout << xref <<  " " << kin.GetQsquared() << " " << F2ref << endl;
  
  cout << xprime << " " << W << " " << pr << " " << costhetar << " " 
	<< alphar << " " << -(massoff*massoff-massi*massi)*1.E-06 << " " << F2pw << " " 
	<< F2fsi << " " << F2ref << " " << cross_data_extr << " " << error_extr << endl;
  
}

void DeuteronCross::getBonusMCresult(double &MCresult, double &modelresultpw, double &modelresultfsi, 
				     double Q2, double W, double Ein, double pr, double costhetar, bool proton,
				     bool pw, bool lc
				    ){

  double massi=proton? MASSP:MASSN;
  double massr=proton? MASSN:MASSP;
  
  double Er=sqrt(massr*massr+pr*pr);
  double Einoff=MASSD-Er;
  double massoff=sqrt(Einoff*Einoff-pr*pr);
  
  double xprime=Q2/(W*W-massoff*massoff+Q2);
  
  //calc nu
  double prz=pr*costhetar;
  double aaa=Einoff*Einoff-prz*prz;
  double bbb=-Einoff*Q2/xprime;
  double ccc=Q2*Q2/(4.*xprime*xprime)-Q2*prz*prz;
  
  double discr=sqrt(bbb*bbb-4.*aaa*ccc);
  double nu1=(-bbb+discr)/(2.*aaa);
  double nu2=(-bbb-discr)/(2.*aaa);
  //cout << nu1 << " " << nu2 << endl;
  double nu=nu2;
  if(costhetar<0.) nu=nu1;
  double qvec=sqrt(Q2+nu*nu);
/*  double xx=Q2/2./(Einoff*nu+prz*qvec);
  cout << xprime << " " << xx << endl;*/
  double x=Q2/(2.*massi*nu);
    
  TKinematics2to2 kin("","",MASSD,MASSP,W,"qsquared:wlab:pklab",Q2,nu,pr);
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
//   DeuteronCross test(*elec,"paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  
  double Eout=Ein-nu;
  //unphysical kinematics!!!
  if(Eout<0.){ /*cout << Eout << endl;*/ MCresult=modelresultpw=modelresultfsi=0.;return;}
  if(std::isnan(asin(sqrt(Q2/(4.*Ein*Eout))))){ /*cout << sqrt(Q2/(4.*Ein*Eout)) << endl; */MCresult=modelresultpw=modelresultfsi=0.;return;}
  
  MCresult= getavgBonus(kin,*elec,lc)*1.E18; //GeV^-6
  modelresultpw= getavgVNALabCross(kin,*elec,1, Einoff)/HBARC/HBARC/10.*2.*Ein*Eout*x/nu/Er;  //go to dEdOmegaed^3ps in GeV-6
  modelresultfsi= pw? modelresultpw: getavgVNALabCross(kin,*elec,0, Einoff)/HBARC/HBARC/10.*2.*Ein*Eout*x/nu/Er;  //go to dEdOmegaed^3ps in GeV-6
}

double DeuteronCross::getavgVNALabCross(TKinematics2to2 &kin,TElectronKinematics &elec, bool pw, double Einoff){
  
  double y=kin.GetWlab()/elec.GetBeamEnergy(kin); 
  double x=kin.GetQsquared()/(2.*massi*kin.GetWlab());
  double front=(4.*PI*ALPHA*ALPHA)/(x*kin.GetQsquared()*kin.GetQsquared())
		*(1-y-(x*x*y*y*massi*massi)/kin.GetQsquared());
  double dens=pw?momdistr.getMomDistrpw(kin):momdistr.getMomDistrfsi(kin,0.);
  double Dstrucs=structure.getavgStructure(kin,elec,Einoff);
//   cout << Dstrucs << endl;
//   cout << "VNA " << front*kin.GetQsquared()*kin.GetWlab()/kin.GetKlab()/kin.GetKlab() << " " << dens*kin.GetEklab() << " " << Dstrucs/(kin.GetQsquared()*kin.GetWlab()/kin.GetKlab()/kin.GetKlab()) << endl;
  return front*dens*Dstrucs*kin.GetEklab()*HBARC*HBARC*1.E19;
  
}

double DeuteronCross::getLCCross(LightConeKin2to2 &kin, bool pw){
//   if(!pw){
//     cerr << "DeuteronCross::getavgLCCross with FSI still in debugging mode. Not useable for now." << endl;
//     assert(1==0);
//   }
  double front=2.*pow(kin.getYA()*ALPHA/kin.getQ2(),2.);
  double dens=pw?momdistr.getMomDistrpwLC(kin):momdistr.getMomDistrfsiLC(kin);
  double Dstrucs=structure.getStructureLC(kin);
  return front*dens*Dstrucs*HBARC*HBARC*1.E19*kin.getAlpha_s()/kin.getAlpha_i();
  
}

double DeuteronCross::getavgLCCross(LightConeKin2to2 &kin, bool pw){
  if(!kin.getIsCollinear()){
    cerr << "DeuteronCross::getavgLCCross only useable in collinear kinematics" << endl;
    assert(1==0);
  }
//   if(!pw){
//     cerr << "DeuteronCross::getavgLCCross with FSI still in debugging mode. Not useable for now." << endl;
//     assert(1==0);
//   }
  double Q2=kin.getQ2();
  double front=2.*PI*ALPHA*ALPHA*kin.getYA()*kin.getYA()/(Q2*Q2*(1.-kin.getEpsilon()));
  double dens=pw?momdistr.getMomDistrpwLC(kin):momdistr.getMomDistrfsiLC(kin);
  double Dstrucs=structure.getavgStructureLC(kin);
//   cout << "LC " << front*kin.getEpsilon()*MASSD << " " << dens*kin.getAlpha_s()/kin.getAlpha_i() << " " << Dstrucs/(kin.getEpsilon()*MASSD) << endl;
  return front*dens*Dstrucs*HBARC*HBARC*1.E19*kin.getAlpha_s()/kin.getAlpha_i();
  
}

double DeuteronCross::getVNACross(LightConeKin2to2 &kin, bool pw){
  double front=2.*pow(kin.getYA()*ALPHA/kin.getQ2(),2.);
  double nu_lab=(kin.getS()+kin.getQ2()-kin.getMassA()*kin.getMassA())/(2.*kin.getMassA());
  double pr=sqrt(pow(kin.getPs_mu()[1]-kin.getA_mu()[1]/2.,2.)+
		 pow(kin.getPs_mu()[2]-kin.getA_mu()[2]/2.,2.)+
		 pow(kin.getPs_mu()[3]-kin.getA_mu()[3]/2.,2.));
  TKinematics2to2 VNAkin("","",kin.getMassA(),kin.getMassN(),kin.getMassX(),"qsquared:wlab:pklab",kin.getQ2(),nu_lab,pr);
  double dens=pw?momdistr.getMomDistrpw(VNAkin):momdistr.getMomDistrfsi(VNAkin,0.);
  double Dstrucs=structure.getStructureLC(kin);
  return front*dens*Dstrucs*HBARC*HBARC*1.E19*sqrt(pr*pr+kin.getMassN()*kin.getMassN());
  
}

double DeuteronCross::getavgVNACross(LightConeKin2to2 &kin, bool pw){
  if(!kin.getIsCollinear()){
    cout << "DeuteronCross::getavgVNACross only useable in collinear kinematics" << endl;
    assert(1==0);
  }
  double Q2=kin.getQ2();
  double front=2.*PI*ALPHA*ALPHA*kin.getYA()*kin.getYA()/(Q2*Q2*(1.-kin.getEpsilon()));
  double nu_lab=(kin.getS()+Q2-kin.getMassA()*kin.getMassA())/(2.*kin.getMassA());
  double pr=sqrt(pow(kin.getPs_mu()[1]-kin.getA_mu()[1]/2.,2.)+
		 pow(kin.getPs_mu()[2]-kin.getA_mu()[2]/2.,2.)+
		 pow(kin.getPs_mu()[3]-kin.getA_mu()[3]/2.,2.));
  TKinematics2to2 VNAkin("","",kin.getMassA(),kin.getMassN(),kin.getMassX(),"qsquared:wlab:pklab",Q2,nu_lab,pr);
  double dens=pw?momdistr.getMomDistrpw(VNAkin):momdistr.getMomDistrfsi(VNAkin,0.);
  double Dstrucs=structure.getavgStructureLC(kin);
  /*cout << "VNA " << front*kin.getEpsilon()*MASSD << " " << dens*kin.getAlpha_s()/kin.getAlpha_i() << " " << Dstrucs/(kin.getEpsilon()*MASSD) 
  << " " << kin.getXN() << " " << endl;
 */ return front*dens*Dstrucs*HBARC*HBARC*1.E19*sqrt(pr*pr+kin.getMassN()*kin.getMassN());
  
}



double DeuteronCross::getavgAzz(TKinematics2to2 &kin,TElectronKinematics &elec, bool azz, bool pw, double Einoff){
  
  double y=kin.GetWlab()/elec.GetBeamEnergy(kin); 
  double x=kin.GetQsquared()/(2.*massi*kin.GetWlab());
  double front=(4.*PI*ALPHA*ALPHA)/(x*kin.GetQsquared()*kin.GetQsquared())
		*(1-y-(x*x*y*y*massi*massi)/kin.GetQsquared());
  double dens=azz?(pw?momdistr.getAzzDistrpw(kin):momdistr.getAzzDistrfsi(kin,0.)):
    (pw?momdistr.getMomDistrpw(kin):momdistr.getMomDistrfsi(kin,0.));
  double Dstrucs=structure.getavgStructure(kin,elec,Einoff);
  return front*dens*Dstrucs*kin.GetEklab()*HBARC*HBARC*1.E19;
  
}

void DeuteronCross::getDeepsresult(double Q2, double W, double Ein, double pr, double costhetar, bool proton,
    double &planewave, double &fsi){

  double massi=proton? MASSP:MASSN;
  double massr=proton? MASSN:MASSP;
  
  double Er=sqrt(massr*massr+pr*pr);
  double Einoff=MASSD-Er;
  double massoff=sqrt(Einoff*Einoff-pr*pr);
  
  double xprime=Q2/(W*W-massoff*massoff+Q2);
  
  //calc nu
  double prz=pr*costhetar;
  double aaa=Einoff*Einoff-prz*prz;
  double bbb=-Einoff*Q2/xprime;
  double ccc=Q2*Q2/(4.*xprime*xprime)-Q2*prz*prz;
  
  double discr=sqrt(bbb*bbb-4.*aaa*ccc);
  double nu1=(-bbb+discr)/(2.*aaa);
  double nu2=(-bbb-discr)/(2.*aaa);
  //cout << nu1 << " " << nu2 << endl;
  double nu=nu2;
  if(costhetar<0.) nu=nu1;
  double qvec=sqrt(Q2+nu*nu);
/*  double xx=Q2/2./(Einoff*nu+prz*qvec);
  cout << xprime << " " << xx << endl;*/
  double x=Q2/(2.*massi*nu);
    
  TKinematics2to2 kin("","",MASSD,MASSP,W,"qsquared:wlab:pklab",Q2,nu,pr);
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
//   DeuteronCross test(*elec,"paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  
  double Eout=Ein-nu;
  //unphysical kinematics!!!
  if(Eout<0.){ planewave=fsi=0.0/0.; return;}
  if(std::isnan(asin(sqrt(Q2/(4.*Ein*Eout))))){ planewave=fsi=0.0/0.; return;}

  double thetain=acos((Ein*Ein+qvec*qvec-Eout*Eout)/(2.*Ein*qvec));
  double yprime=(Einoff*nu+prz*qvec)/(Ein*Einoff+Ein*cos(thetain)*pr*costhetar/*+Ein*sin(thetain)*pr*sqrt(1.-costhetar*costhetar)*/);  
  double R=0.18;
  double frontdeeps=(4.*PI*ALPHA*ALPHA)/(xprime*Q2*Q2)*
  (yprime*yprime/(2.*(1.+R))+(1.-yprime)+(pow(massoff*xprime*yprime,2.)*(1.-R))/(Q2*(1+R)));
  double dxprimedx=-2.*xprime*xprime*nu/(x*qvec)*((Einoff+prz)/(nu-qvec)+1./(2.*xprime));
  
  
  planewave = getavgVNALabCross(kin,*elec,1,Einoff)/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10;
  fsi= getavgVNALabCross(kin,*elec,0,Einoff)/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10;
  return;
    
}

void DeuteronCross::getDeepsresultLC(double Q2, double W, double Ein, double pr, double costhetar, bool proton,
    double &planewave, double &fsi){

  double massi=proton? MASSP:MASSN;
  double massr=proton? MASSN:MASSP;
  
  double Er=sqrt(massr*massr+pr*pr);
  double Einoff=MASSD-Er;
  double massoff=sqrt(Einoff*Einoff-pr*pr);
  
  double xprime=Q2/(W*W-massoff*massoff+Q2);
  
  //calc nu
  double prz=pr*costhetar;
  double aaa=Einoff*Einoff-prz*prz;
  double bbb=-Einoff*Q2/xprime;
  double ccc=Q2*Q2/(4.*xprime*xprime)-Q2*prz*prz;
  
  double discr=sqrt(bbb*bbb-4.*aaa*ccc);
  double nu1=(-bbb+discr)/(2.*aaa);
  double nu2=(-bbb-discr)/(2.*aaa);
  //cout << nu1 << " " << nu2 << endl;
  double nu=nu2;
  if(costhetar<0.) nu=nu1;
  double qvec=sqrt(Q2+nu*nu);
/*  double xx=Q2/2./(Einoff*nu+prz*qvec);
  cout << xprime << " " << xx << endl;*/
  double x=Q2/(2.*massi*nu);
    
  TKinematics2to2 kin("","",MASSD,MASSP,W,"qsquared:wlab:pklab",Q2,nu,pr);
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
//   DeuteronCross test(*elec,"paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  
  double Eout=Ein-nu;
  //unphysical kinematics!!!
  if(Eout<0.){ planewave=fsi=0.0/0.; return;}
  if(std::isnan(asin(sqrt(Q2/(4.*Ein*Eout))))){ planewave=fsi=0.0/0.; return;}

  double thetain=acos((Ein*Ein+qvec*qvec-Eout*Eout)/(2.*Ein*qvec));
  double yprime=(Einoff*nu+prz*qvec)/(Ein*Einoff+Ein*cos(thetain)*pr*costhetar/*+Ein*sin(thetain)*pr*sqrt(1.-costhetar*costhetar)*/);  
  double R=0.18;
  double frontdeeps=(4.*PI*ALPHA*ALPHA)/(xprime*Q2*Q2)*
  (yprime*yprime/(2.*(1.+R))+(1.-yprime)+(pow(massoff*xprime*yprime,2.)*(1.-R))/(Q2*(1+R)));
  double dxprimedx=-2.*xprime*xprime*nu/(x*qvec)*((Einoff+prz)/(nu-qvec)+1./(2.*xprime));
  
  
  planewave = getavgVNALabCross(kin,*elec,1,Einoff)/*/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10*/;
  fsi= getavgVNALabCross(kin,*elec,0,Einoff)/*/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10*/;

  TVector3 vecps(pr*sqrt(1-costhetar*costhetar),0.,prz);
  TVector3 veckin(Ein*sin(thetain),0.,Ein*cos(thetain));
  
  LightConeKin2to2 kinLC(MASSD,Q2,massr,0.,qvec,vecps,veckin);
  double planewaveLC = getavgLCCross(kinLC,1);
  double fsiLC = getavgLCCross(kinLC,0);
  cout << RADTODEGR*acos(costhetar) << " " << planewave << " " << planewaveLC << " " << fsi << " " << fsiLC <<  endl;
  return;
    
}

void DeuteronCross::getDeepsAzz(double Q2, double W, double Ein, double pr, double costhetar, bool proton,
    double &Azz, double &Azzfsi, double &planewave, double &fsi){

  double massi=proton? MASSP:MASSN;
  double massr=proton? MASSN:MASSP;
  
  double Er=sqrt(massr*massr+pr*pr);
  double Einoff=MASSD-Er;
  double massoff=sqrt(Einoff*Einoff-pr*pr);
  
  double xprime=Q2/(W*W-massoff*massoff+Q2);
  
  //calc nu
  double prz=pr*costhetar;
  double aaa=Einoff*Einoff-prz*prz;
  double bbb=-Einoff*Q2/xprime;
  double ccc=Q2*Q2/(4.*xprime*xprime)-Q2*prz*prz;
  
  double discr=sqrt(bbb*bbb-4.*aaa*ccc);
  double nu1=(-bbb+discr)/(2.*aaa);
  double nu2=(-bbb-discr)/(2.*aaa);
  //cout << nu1 << " " << nu2 << endl;
  double nu=nu2;
  if(costhetar<0.) nu=nu1;
  double qvec=sqrt(Q2+nu*nu);
/*  double xx=Q2/2./(Einoff*nu+prz*qvec);
  cout << xprime << " " << xx << endl;*/
  double x=Q2/(2.*massi*nu);
    
  TKinematics2to2 kin("","",MASSD,MASSP,W,"qsquared:wlab:pklab",Q2,nu,pr);
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
//   DeuteronCross test(*elec,"paris",proton,"SLAC",36.3274,1.97948,-0.5,8.,1.2,4);
  
  double Eout=Ein-nu;
  //unphysical kinematics!!!
  if(Eout<0.){ planewave=fsi=0.0/0.; return;}
  if(std::isnan(asin(sqrt(Q2/(4.*Ein*Eout))))){ planewave=fsi=0.0/0.; return;}

  double thetain=acos((Ein*Ein+qvec*qvec-Eout*Eout)/(2.*Ein*qvec));
  double yprime=(Einoff*nu+prz*qvec)/(Ein*Einoff+Ein*cos(thetain)*pr*costhetar/*+Ein*sin(thetain)*pr*sqrt(1.-costhetar*costhetar)*/);  
  double R=0.18;
  double frontdeeps=(4.*PI*ALPHA*ALPHA)/(xprime*Q2*Q2)*
  (yprime*yprime/(2.*(1.+R))+(1.-yprime)+(pow(massoff*xprime*yprime,2.)*(1.-R))/(Q2*(1+R)));
  double dxprimedx=-2.*xprime*xprime*nu/(x*qvec)*((Einoff+prz)/(nu-qvec)+1./(2.*xprime));
  
  
  Azz = getavgAzz(kin,*elec,1,1,Einoff)/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10;
  Azzfsi= getavgAzz(kin,*elec,1,0,Einoff)/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10;
  planewave = getavgAzz(kin,*elec,0,1,Einoff);
  fsi= getavgAzz(kin,*elec,0,0,Einoff);
  
  return;
    
}

void DeuteronCross::readin_deeps(double ******pdeepsarray, string dir){
  
  int Q2index=2;
  int Windex=5;
  int psindex=5;
  int thetaindex=34;
 
  (*pdeepsarray)=new double ****[Q2index];
  for(int i=0;i<Q2index;i++){
    (*pdeepsarray)[i]=new double ***[Windex];
    for(int j=0;j<Windex;j++){
      (*pdeepsarray)[i][j]=new double **[psindex];
      for(int k=0;k<psindex;k++){
	(*pdeepsarray)[i][j][k]=new double*[thetaindex];
	for(int l=0;l<thetaindex;l++) (*pdeepsarray)[i][j][k][l]=new double[4];
      }
    }
  }
 
  string qq[2]={"18","28"};
  string ww[5]={"125","15","173","2","24"};
  string pp[5]={"300","340","390","460","560"};
  for(int i=0;i<Q2index;i++){
    for(int j=0;j<Windex;j++){
      for(int k=0;k<psindex;k++){ 
	string filename=dir+"/deepsdata/exp.Q"+qq[i]+".W"+ww[j]+".p"+pp[k]+".dat";
	ifstream file(filename.c_str(),ios::in);
	if (file.is_open()){
	  for(int l=0;l<thetaindex;l++) file >> (*pdeepsarray)[i][j][k][l][0]>>(*pdeepsarray)[i][j][k][l][1]
	    >>(*pdeepsarray)[i][j][k][l][2]>>(*pdeepsarray)[i][j][k][l][3];
	}
	else{cout << filename << endl; cerr << "could not read in deeps files" << endl;}
      }
    }
  }
  
  return;
}


void DeuteronCross::maint_deepsarray(double *****deepsarray){
  
  
  int Q2index=2;
  int Windex=5;
  int psindex=5;
  int thetaindex=34;
  
  for(int i=0;i<Q2index;i++){
    for(int j=0;j<Windex;j++){
      for(int k=0;k<psindex;k++){ 
	for(int l=0;l<thetaindex;l++) delete [] deepsarray[i][j][k][l];
	delete [] deepsarray[i][j][k];
      }
      delete [] deepsarray[i][j];
    }
    delete [] deepsarray[i];
  }
  delete [] deepsarray;
  
}

void DeuteronCross::setScatter(double sigmain, double betain, double epsin){
  momdistr.setScatter(sigmain,betain,epsin);
}


double DeuteronCross::sigmaparam(double W_sq, double Q2){
  /*cout << sqrt(W_sq) << " " << 65/10.*1.8E06/Q2*INVHBARC*INVHBARC << " " << (25.3*1.E-06*2.3+53*(sqrt(W_sq>5.76E06?5.76E06:W_sq)-MASSP)*1.E-03)
	/(1.E-05*Q2)*INVHBARC*INVHBARC << " " << endl;
  */if(abs(W_sq-1.232*1.232E06)<2.5E5) return 65/10.*1.8E06/Q2*INVHBARC*INVHBARC;
  return (25.3*1.E-06*2.3+53*(sqrt(W_sq>5.76E06?5.76E06:W_sq)-MASSP)*1.E-03)
	/(1.E-05*Q2)*INVHBARC*INVHBARC;
}