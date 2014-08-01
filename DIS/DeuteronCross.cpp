#include "DeuteronCross.hpp"
#include <fstream>
#include <TLorentzRotation.h>
#include <FourVector.h>
#include <NuclStructure.hpp>
#include <TVector3.h>

using namespace std;

DeuteronCross::DeuteronCross(string name, bool proton, string strucname,
			     double sigmain, double betain, double epsilonin, double betaoffin, double lambdain, int offshellset,
			     int looplimit
			    ):
massi(proton? MASSP:MASSN),
momdistr(name,massi,offshellset,sigmain,betain,epsilonin,betaoffin,lambdain,looplimit),
structure(proton,strucname)
{
//   for(int i=0;i<200;i+=5){
//     TVector3 p(0.,0.,i);
//     double E=sqrt(p.Mag2()+MASSP*MASSP);
//     double alpha=2.*(E-p[2])/MASSD; //lightcone alpha_s=(E-p_z)/M_n
//     double pt2=p[0]*p[0]+p[1]*p[1];
//     double k=sqrt((M_NUCL*M_NUCL+pt2)/(alpha*(2.-alpha))-M_NUCL*M_NUCL); //lightcone momentum rescaling
// 
//       cout << i << " " << alpha << " " << k << " " << pow(MASSD-E,2.)-MASSP*MASSP <<" " << momdistr.getMomDistrpw(p)/pow(momdistr.getDeuteronwf()->getResidu(),2.)*pow(i*i+0.2316*0.2316*HBARC*HBARC,2.) << " " 
// 	<< momdistr.getLCMomDistrpw(p)  << endl;
//   }
//   exit(1);
}


DeuteronCross::~DeuteronCross(){
  
}

double DeuteronCross::getavgBonus(TKinematics2to2 &kin,TElectronKinematics &elec, bool lc){
  
  double dens=lc? momdistr.getLCMomDistrpw(kin) : momdistr.getMomDistrpw(kin); 
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
  
  double Eout=Ein-nu;
  //unphysical kinematics!!!
//   if(Eout<0.){ /*cout << Eout << endl;*/ MCresult=modelresultpw=modelresultfsi=0.;return;}
//   if(std::isnan(asin(sqrt(Q2/(4.*Ein*Eout))))){ /*cout << sqrt(Q2/(4.*Ein*Eout)) << endl; */MCresult=modelresultpw=modelresultfsi=0.;return;}
  
  double MCresult= getavgBonus(kin,*elec,lc)*1.E18;
  double cross_data_extr=Rdata*MCresult/norm;
  double y=kin.GetWlab()/elec->GetBeamEnergy(kin); 
  double front=(4.*PI*ALPHA*ALPHA)/(x*kin.GetQsquared()*kin.GetQsquared())
		*(1-y-(x*x*y*y*massi*massi)/kin.GetQsquared());
  double denspw=momdistr.getMomDistrpw(kin);
  double densfsi=momdistr.getMomDistrfsi(kin,0.);
  double Dstrucs=structure.getavgStructure(kin,*elec,Einoff);
//   cout << Dstrucs << endl;
//   front*dens*Dstrucs*kin.GetEklab()*HBARC*HBARC*1.E19;
//   modelresultpw= getavgCross(kin,*elec,1, Einoff)/HBARC/HBARC/10.*2.*Ein*Eout*x/nu/Er;  //go to dEdOmegaed^3ps in MeV-6
//   modelresultfsi= pw? modelresultpw: getavgCross(kin,*elec,0, Einoff)/HBARC/HBARC/10.*2.*Ein*Eout*x/nu/Er;  //go to dEdOmegaed^3ps in MeV-6
// //   double residu=sqrt(2.)*c[0]*sqrt(HBARC)/PI; //[MeV^1/2]
//   double corr = pow((massoff*massoff-massi*massi)/residu,2.) // [MeV^2]
  
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
  
  MCresult= getavgBonus(kin,*elec,lc)*1.E18;
  modelresultpw= getavgCross(kin,*elec,1, Einoff)/HBARC/HBARC/10.*2.*Ein*Eout*x/nu/Er;  //go to dEdOmegaed^3ps in MeV-6
  modelresultfsi= pw? modelresultpw: getavgCross(kin,*elec,0, Einoff)/HBARC/HBARC/10.*2.*Ein*Eout*x/nu/Er;  //go to dEdOmegaed^3ps in MeV-6
}

double DeuteronCross::getavgCross(TKinematics2to2 &kin,TElectronKinematics &elec, bool pw, double Einoff){
  
  double y=kin.GetWlab()/elec.GetBeamEnergy(kin); 
  double x=kin.GetQsquared()/(2.*massi*kin.GetWlab());
  double front=(4.*PI*ALPHA*ALPHA)/(x*kin.GetQsquared()*kin.GetQsquared())
		*(1-y-(x*x*y*y*massi*massi)/kin.GetQsquared());
  double dens=pw?momdistr.getMomDistrpw(kin):momdistr.getMomDistrfsi(kin,0.);
  double Dstrucs=structure.getavgStructure(kin,elec,Einoff);
//   cout << Dstrucs << endl;
  return front*dens*Dstrucs*kin.GetEklab()*HBARC*HBARC*1.E19;
  
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
  
  
  planewave = getavgCross(kin,*elec,1,Einoff)/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10;
  fsi= getavgCross(kin,*elec,0,Einoff)/frontdeeps/dxprimedx/Er/HBARC/HBARC/1.E10;
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
