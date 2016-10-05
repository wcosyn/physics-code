#include "NucleonStructure.hpp"
#include "constants.hpp"
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include "mstwpdf.h"

using namespace std;

//nucleon structure parametrization from Christy & Bosted (fortran)
extern"C"{
    void f1f2in09_(double *Z, double *A, double *QSQ, double *Wsq, double *F1, double *F2, double *rc);
}

//nucleon structure parametrization from Alekhin leading twist (fortran)
extern"C"{
  void alekhin_(double *xb,double *q2,double PDFS[22],int *IPAR,int *ICOL);
}

//data structures needed in Alekhin parametrization
extern "C"{
  extern struct{
    char homedir[128];
    char a09file1[128];
    char a09file2[128];
  } dir_;
}

//nucleon structure parametrization from CTEQ (based on PDFs)
extern "C"{
  double f2_ctq_(int *proton, double *x, double *Q2);
}

//g1N structure parametrization from GRSV ( hep-ph/0011215, grsv2000pdf_g1.f)
extern "C"{
  void parpol_(int *iset, double *x, double *q2, double *u, 
               double *d, double *ub, double *db, double *st, double *gl, double *g1p, double *g1n);
}


//fortran prototypes for structs needed in grsv g1 parametrization
extern "C"{
  extern struct{
    char fileslo[128];
  } grsvdir_;
}

extern "C"{
  extern struct{
    int iini=7;
  } intini_;
}

//HMRS F1,2 from LO pdfs
extern "C"{
 double f1p_hmrs_(double *x, double *Q2);
}

extern "C"{
 double f1n_hmrs_(double *x, double *Q2);
}

extern "C"{
 double f2n_hmrs_(double *x, double *Q2);
}

extern "C"{
 double f2p_hmrs_(double *x, double *Q2);
}


//r1990 SLAC parametrization
extern "C"{
  double r1990_(double *x, double *q2);
}

//r1990 SLAC parametrization
extern "C"{
  double r1998_(double *x, double *q2);
}


NucleonStructure::NucleonStructure():
name("SLAC"),dir(HOMEDIR){
  if(!name.compare("MSTW")){
    string file= string(HOMEDIR)+"/mstw_grids/mstw2008lo.00.dat";
    mstw = new c_mstwpdf(file,false,true);
  }
  else mstw=NULL;

  
  return;

}


NucleonStructure::NucleonStructure(const string nm):
name(nm),dir(HOMEDIR){
  if(!name.compare("MSTW")){
    string file= string(HOMEDIR)+"/mstw_grids/mstw2008lo.00.dat";
    mstw = new c_mstwpdf(file,false,true);
  }
  else mstw=NULL;

  
  return;

}


NucleonStructure::~NucleonStructure(){
  if(mstw) delete mstw;
}

NucleonStructure::NucleonStructure(const NucleonStructure& rhs){
  name=rhs.name;
  dir=rhs.dir;
  if(!name.compare("MSTW")){
    string file= string(HOMEDIR)+"/mstw_grids/mstw2008lo.00.dat";
    mstw = new c_mstwpdf(file,false,true);
  }
  else mstw=NULL;

  
  return;
}

NucleonStructure& NucleonStructure::operator=(const NucleonStructure& rhs){
  if(this!=&rhs) { // avoid self-assignment
    name=rhs.name;
    dir=rhs.dir;
    if(!name.compare("MSTW")){
      string file= string(HOMEDIR)+"/mstw_grids/mstw2008lo.00.dat";
      mstw = new c_mstwpdf(file,false,true);
    }
    else mstw=NULL;
    
  }
  return *this;
 
  
}



void NucleonStructure::getF_Alekhin(double &F1, double &F2, const bool proton, const double x, const double Q2)const{
  if(x>1){
    F1=0.;
    F2=0.;
    return;
  }
  strcpy(dir_.homedir,dir.c_str());
  strcpy(dir_.a09file1,dir.c_str());
  strcpy(dir_.a09file2,dir.c_str());
  strcat(dir_.a09file1,"/alekhin/a09.sfs_lNNC");
  strcat(dir_.a09file2,"/alekhin/a09.dsfs_lNNC");
  //cout <<dir_.a09file1 << endl;
  int ONE=1;
  int ZERO=0;
  double SF[22];
  double Qsq=Q2*1.E-06;
  double xx=x;
  alekhin_(&xx,&Qsq,SF,&ZERO,&ONE);
  if(proton){
    F1=SF[7]/(2.*x);
    F2=SF[1];
  }
  else{
    F1=SF[10]/(2.*x);
    F2=SF[4];
  }    
}

double NucleonStructure::getF1_Alekhin(const bool proton, const double x, const double Q2) const{
  double F1,F2;
  getF_Alekhin(F1,F2,proton,x,Q2);
  return F1;
}

double NucleonStructure::getF2_Alekhin(const bool proton, const double x, const double Q2) const{
  double F1,F2;
  getF_Alekhin(F1,F2,proton,x,Q2);
  return F2;  
}

void NucleonStructure::getF_SLAC(double &F1, double &F2, const bool proton, const double x, const double Q2) const{
  F2=getF2_SLAC(proton,x,Q2);
  double R=getr1998(x,Q2);
  double gamma2=pow(2*(proton? MASSP:MASSN)*x,2.)/Q2;
  F1=F2/2./x/(R+1)*(1+gamma2);
}


double NucleonStructure::getF2_SLAC(const bool proton, const double x, const double Q2) const{
  double mass = proton? MASSP:MASSN;
  double Wsq=Q2*(1./x-1.)+mass*mass;
  return proton? f2p_b(mass*1.E-03, x,Q2*1.E-06, sqrt(Wsq*1.E-06)):f2n_b(mass*1.E-03, x,Q2*1.E-06, sqrt(Wsq*1.E-06));
}

void NucleonStructure::getF_CTEQ(double &F1, double &F2, const bool proton, const double x, const double Q2) const{
  F2=getF2_CTEQ(proton,x,Q2);
  double R=getr1998(x,Q2);
  double gamma2=pow(2*(proton? MASSP:MASSN)*x,2.)/Q2;
  F1=F2/2./x/(R+1)*(1+gamma2);

  //F1=F2*2.*x/(1+R)*(pow(alphai/alphaq+1/(2.*x),2.);
}


double NucleonStructure::getF2_CTEQ(const bool proton, const double x, const double Q2) const{
  double q2=Q2*1.E-06;
  int pr = int(proton);
  double xx=x;
  return f2_ctq_(&pr,&xx,&q2);
}

void NucleonStructure::getF_CB(double &F1, double &F2, const bool proton, const double x, const double Q2) const{
  double mass=proton? MASSP:MASSN;
  double Wsq=Q2*(1./x-1.)+mass*mass;
  if(Wsq>25.E06){ 
    cerr << "Invariant mass too big for Christy & Bosted Parametrization!! " << Wsq << endl;
    exit(1);
  }
  double R_CB;
  double ONE=1.;
  double ZERO=0.;
  double QSQ=Q2*1.E-06;
  double WSQ=Wsq*1.E-06;
  if(proton) f1f2in09_(&ONE, &ZERO, &QSQ, &WSQ, &F1, &F2, &R_CB);
  else f1f2in09_(&ZERO, &ONE, &QSQ, &WSQ, &F1, &F2, &R_CB);
//   cout << F1 << " " << F2 << " " << QSQ << " " << WSQ << " " << 1./(1.+(WSQ-MASSN*MASSN*1.E-06)/QSQ) << endl;
  return;
}

double NucleonStructure::getF1_CB(const bool proton, const double x, const double Q2) const{
  double F1,F2;
  getF_CB(F1,F2,proton,x,Q2);
  return F1;
  
}


double NucleonStructure::getF2_CB(const bool proton, const double x, const double Q2) const{
  double F1,F2;
  getF_CB(F1,F2,proton,x,Q2);
  return F2;
  
}

void NucleonStructure::getF_HMRS(double &F1,double &F2, const bool proton, const double x, const double Q2) const{
  double q2gev=Q2*1.E-06;
  double xx=x;
  if(proton){ F1=f1p_hmrs_(&xx, &q2gev); F2=f2p_hmrs_(&xx, &q2gev);}
  else{F1=f1n_hmrs_(&xx, &q2gev); F2=f2n_hmrs_(&xx, &q2gev);}
  
}


void NucleonStructure::getF_MSTW(double &F1,double &F2,const bool proton, const double x, const double Q2) const{
  double u,d,s,c,b,ubar,dbar,sbar,cbar,bbar;
  double qgev=sqrt(Q2*1.E-06);
  mstw->update(x,qgev);

  d=mstw->parton(1,x,qgev);
  u=mstw->parton(2,x,qgev);
  s=mstw->parton(3,x,qgev);
  c=mstw->parton(4,x,qgev);
  b=mstw->parton(5,x,qgev);
  dbar=mstw->parton(-1,x,qgev);
  ubar=mstw->parton(-2,x,qgev);
  sbar=mstw->parton(-3,x,qgev);
  cbar=mstw->parton(-4,x,qgev);
  bbar=mstw->parton(-5,x,qgev);
  if(proton){ 
    F2=1./9.*(4.*(u+ubar)+(d+dbar)+(s+sbar)+4*(c+cbar)+(b+bbar));
    F1=1./(2.*x)/(1.+getr1998(x,Q2))*(1.+4.*MASSP*MASSP*x*x/Q2)*F2;
  }
  else{
    F2=1./9.*((u+ubar)+4.*(d+dbar)+(s+sbar)+4*(c+cbar)+(b+bbar));
    F1=1./(2.*x)/(1.+getr1998(x,Q2))*(1.+4.*MASSN*MASSN*x*x/Q2)*F2;
  }
//   cout << x << " " << q2gev << " " << F1 << " " << F2 << endl;
}



void NucleonStructure::getF_xQ(double &F1, double &F2,const bool proton, const double x, const double Q2) const{
  if(!name.compare("SLAC")) getF_SLAC(F1,F2,proton,x,Q2);
  else if(!name.compare("Alekhin")) getF_Alekhin(F1,F2,proton,x,Q2);
  else if(!name.compare("CB")) getF_CB(F1,F2,proton,x,Q2);
  else if(!name.compare("CTEQ")) getF_CTEQ(F1,F2,proton,x,Q2);
  else if(!name.compare("HMRS")) getF_HMRS(F1,F2,proton,x,Q2);
  else if(!name.compare("MSTW")) getF_MSTW(F1,F2,proton,x,Q2);
  else { cerr << "invalid name for structure function" << endl; exit(1);}
  return;
}

double NucleonStructure::getF1_xQ(const bool proton, const double x, const double Q2) const{
  double F1,F2;
  getF_xQ(F1,F2,proton,x,Q2);
  return F1;
}

double NucleonStructure::getF2_xQ(const bool proton, const double x, const double Q2) const{
  double F1,F2;
  getF_xQ(F1,F2,proton,x,Q2);
  return F2;  
}


void NucleonStructure::getF_WQ(double &F1, double &F2,const bool proton, const double Wsq, const double Q2) const{

  double mass=proton? MASSP:MASSN;
  double x=Q2/(Wsq-mass*mass+Q2);
  getF_xQ(F1,F2,proton,x,Q2);
}

double NucleonStructure::getF2_WQ(const bool proton, const double Wsq, const double Q2) const{
  double mass=proton? MASSP:MASSN;
  double x=Q2/(Wsq-mass*mass+Q2);
  return getF2_xQ(proton,x,Q2);
}
  
double NucleonStructure::getF1_WQ(const bool proton, const double Wsq, const double Q2) const{
  double mass=proton? MASSP:MASSN;
  double x=Q2/(Wsq-mass*mass+Q2);
  return getF1_xQ(proton,x,Q2);
}
  


double NucleonStructure::f2p_b(double massi, double xp, double q2, double fm) const{

  double f2p,yn,wwi,t,gw;
  //f2p =  0.0;
  yn  =  q2/(2*massi*xp);//nuoffshell
  //fm2 = -q2 +2.0*massi*yn + wstar2;
  if ((xp > 1.0) | (fm < massi)){return 0;}
  //	fm=pow(fm2,.5);
  //fm  = sqrt(fm2);
  //cout << fm;
  wwi = (2.*yn*massi+1.642)/(q2 +0.376);
  t = 1.0 - 1.0/wwi;
  gw=0.256*pow(t,3)+2.178*pow(t,4)+
        0.898*pow(t,5)-6.716*pow(t,6)+3.756*pow(t,7);
  f2p = bodek(fm,q2)*gw*wwi*xp;
  return f2p;
}

// Bodek's parameterization of neutron's F2 
double NucleonStructure::f2n_b(double massi, double xp, double q2, double fm) const{
  
  
  double f2n,yn,wwi,t,gw;
  //f2n =  0.0;
  yn  =  q2/(2*massi*xp);
  //fm2 = -q2 + 2.0*massi*yn + wstar2;
  if ((xp > 1.0) | (fm < massi)){return 0;}
  //  fm=pow(fm2,2);
  //fm = sqrt(fm2);
  //cout << fm;
  wwi = (2*yn*massi + 1.642)/(q2 + 0.376);
  t   = 1.0-1.0/wwi;
  gw  = 0.064*pow(t,3) + 0.225*pow(t,4)+
  4.106*pow(t,5) - 7.079*pow(t,6)+3.055*pow(t,7);
  f2n = bodek(fm,q2)*gw*wwi*xp;
  return f2n;
}


/************************************************/
/*		BODEK PARAMETRIZATION		*/
/************************************************/

double NucleonStructure::bodek(double wm, double qsq) const
{
 const double c[25]={0., 1.0741163,  0.75531124, 3.3506491   , 1.7447015  ,
                    3.5102405,  1.040004  , 1.2299128   , 0.10625394 ,
                    0.48132786, 1.5101467 , 0.081661975 , 0.65587179 ,
                    1.7176216 , 0.12551987, 0.7473379   , 1.953819   ,
                    0.19891522,-0.17498537, 0.0096701919,-0.035256748, 
		    3.5185207 ,-0.59993696, 4.7615828   , 0.41167589}; 
        const int lspin[5]={0,1,2,3,2};
	      int nres=4,nbkg=5,index,i;
	
	double pmsq=.880324,pm2=1.876512,pm=0.938256;
	double b,b1,b2,wsq,omega,x,xpx,piemsq;
	double eb1,eb2,bbkg,bres,ressum;
	double ram,rma,rwd,qstarn,qstaro,j,k;
	double term,termo,gamres,brwig,res;
	b=0;

	//  cout <<  "c2 " << c[2]  <<  endl;

	if (wm < .94){ return b;}
	wsq=wm*wm;
	omega=1+(wsq-pmsq)/qsq;
	x=1/omega;
	xpx=c[22]+c[23]*(x-c[24])*(x-c[24]);
	//  cout <<  "xpx " << xpx  <<  endl;

	piemsq=(c[1]-pm)*(c[1]-pm);
	b1=0;
	if (wm == c[1]){goto label11;}
	b1=max(0. , (wm - c[1]))/(wm-c[1])*c[2];
	//  cout <<  "b1 " << b1  <<  endl;
	
label11:
	eb1=c[3]*(wm-c[1]);
	if (eb1 > 25.){goto label1;}
	b1=b1*(1-exp(-eb1));
	b2=0;
	if (wm == c[4]){goto label12;}
label1:
	b2=max(0., (wm - c[4]))/(wm-c[4])*(1.-c[2]);
	//  cout <<  "b2 " << b2  <<  endl;

label12:
	eb2=c[5]*(wsq-c[4]*c[4]);
	if (eb2 > 25.){goto label2;}
	b2=b2*(1.-exp(-eb2));
	//  cout <<  "b2 " << b2  <<  endl;

label2:
	bbkg=b1+b2;
	bres=c[2]+b2;
	ressum=0.;
	for(i=1;i<=nres;i++)
	{
		index=(i-1)*3+1+nbkg;
		ram=c[index];
		if (i==1) {ram=c[index]+c[18]*qsq+c[19]*qsq*qsq;}
		rma=c[index+1];
		if (i==3) {rma=rma*(1+c[20]/(1+c[21]*qsq));}
		rwd=c[index+2];
		qstarn=sqrt(max(0.,pow((wsq+pmsq-piemsq)/(2*wm),2)-pmsq));
		qstaro=sqrt(max(0.,pow((rma*rma-pmsq+piemsq)/(2*rma),2)-piemsq));
		if (qstaro <= 1E-10){goto label40;}
		term=6.08974*qstarn;
		termo=6.08974*qstaro;
		j=2*lspin[i];
		k=j+1;
		gamres=rwd*pow((term/termo),k)*(1+pow(termo,j))/(1+pow(term,j));
		gamres=gamres/2;
		brwig=gamres/((pow(wm-rma,2)+pow(gamres,2))*3.1415926);
		res=ram*brwig/pm2;
		goto label30;
	label40:
		res=0;
	label30:
		ressum=ressum+res;
	}
	b=bbkg*(1.+(1.-bbkg)*xpx)+ressum*(1.-bres);
	return b;
}


double NucleonStructure::getG1_grsv2000(const bool proton, const double x, const double Q2){
  
  strcpy(grsvdir_.fileslo,HOMEDIR);
  strcat(grsvdir_.fileslo,"/grsv/std2000_lo_g1.grid");
  if(intini_.iini!=1) intini_.iini=0;
  int iset=3;
  double dummy,g1p,g1n;
  double q2_ = Q2/1.E06;
  double xx=x;
  parpol_(&iset, &xx, &q2_, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &g1p, &g1n);
  return (proton? g1p:g1n);
}



double NucleonStructure::getG1plusg2_grsv2000(const bool proton, const double x, const double Q2) {
    
  double step=(1.-x)/100.;
  double result=0.;
  for(int i=0;i<=100;i++){
    double y=x+step*i;
    result+=NucleonStructure::getG1_grsv2000(proton,y,Q2)/y;
  }
  return result*step;
}

double NucleonStructure::getr1990(const double x, const double Q2){
    double q2gev=Q2*1.E-06;
    double xx=x;
    return r1990_(&xx,&q2gev);
}

double NucleonStructure::getr1998(const double x, const double Q2){
    double q2gev=Q2*1.E-06;
    double xx=x;
    return r1998_(&xx,&q2gev);
}