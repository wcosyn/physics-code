#include "NuclStructure.hpp"
#include "constants.hpp"
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <cstring>

extern"C"{
    void f1f2in09_(double *Z, double *A, double *QSQ, double *Wsq, double *F1, double *F2, double *rc);
}

extern"C"{
  void alekhin_(double *xb,double *q2,double PDFS[22],int *IPAR,int *ICOL);
}


extern "C"{
  extern struct{
    char homedir[128];
    char a09file1[128];
    char a09file2[128];
  } dir_;
}


NuclStructure::NuclStructure(bool pr, double var1, double var2, int switchvar, string nm="SLAC"):
name(nm),proton(pr), mass(pr?MASSP:MASSN),dir(HOMEDIR){
  switch(switchvar){
    case 0:
      Q2=var1;
      x=var2;
      Wsq=Q2*(1./x-1.)+mass*mass;
      break;
    case 1:
      Q2=var1;
      Wsq=var2;
      x=Q2/(Wsq-mass*mass+Q2);
      break;
    case 2:
      x=var1;
      Wsq=var2;
      Q2=(Wsq-mass*mass)/(1./x-1.);
      break;
    default:
      cerr << "non supported format in constructor" << endl;
      exit(1);
  }
  return;

}
NuclStructure::NuclStructure(bool pr, double Q2in, double xin, double Wsqin, string nm="SLAC"):
name(nm),proton(pr), mass(pr?MASSP:MASSN),dir(HOMEDIR),x(xin),Q2(Q2in),Wsq(Wsqin){
  
}


void NuclStructure::getF_Alekhin(double &F1, double &F2){
  strcpy(dir_.homedir,dir.c_str());
  strcpy(dir_.a09file1,dir.c_str());
  strcpy(dir_.a09file2,dir.c_str());
  strcat(dir_.a09file1,"/grids/a09.sfs_lNNC");
  strcat(dir_.a09file2,"/grids/a09.dsfs_lNNC");
  //cout <<dir_.a09file1 << endl;
  int ONE=1;
  int ZERO=0;
  double SF[22];
  double Qsq=Q2*1.E-06;
  alekhin_(&x,&Qsq,SF,&ZERO,&ONE);
  if(proton){
    F1=SF[7]/(2.*x);
    F2=SF[1];
  }
  else{
    F1=SF[10]/(2.*x);
    F2=SF[4];
  }    
}

double NuclStructure::getF1_Alekhin(){
  double F1,F2;
  getF_Alekhin(F1,F2);
  return F1;
}

double NuclStructure::getF2_Alekhin(){
  double F1,F2;
  getF_Alekhin(F1,F2);
  return F2;  
}

void NuclStructure::getF_SLAC(double &F1, double &F2) const{
  F2=getF2_SLAC();
  double R=0.18;
  //F1=F2*2.*x/(1+R)*(pow(alphai/alphaq+1/(2.*x),2.);
}

double NuclStructure::getF1_SLAC() const{
  return 0.;
}

double NuclStructure::getF2_SLAC() const{
  return proton? f2p_b(mass*1.E-03, x,Q2*1.E-06, sqrt(Wsq*1.E-06)):f2n_b(mass*1.E-03, x,Q2*1.E-06, sqrt(Wsq*1.E-06));
}

void NuclStructure::getF_CB(double &F1, double &F2) const{
  if(Wsq>25.E06){ 
    cerr << "Invariant mass too big for Christy & Bosted Parametrization!!" << endl;
    exit(1);
  }
  double R_CB;
  double ONE=1.;
  double ZERO=0.;
  double QSQ=Q2*1.E-06;
  double WSQ=Wsq*1.E-06;
  if(proton) f1f2in09_(&ONE, &ZERO, &QSQ, &WSQ, &F1, &F2, &R_CB);
  else f1f2in09_(&ZERO, &ONE, &QSQ, &WSQ, &F1, &F2, &R_CB);
  return;
}

double NuclStructure::getF1_CB() const{
  double F1,F2;
  getF_CB(F1,F2);
  return F1;
  
}


double NuclStructure::getF2_CB() const{
  double F1,F2;
  getF_CB(F1,F2);
  return F2;
  
}


void NuclStructure::getF(double &F1, double &F2){
  if(!name.compare("SLAC")) getF_SLAC(F1,F2);
  else if(!name.compare("Alekhin")) getF_Alekhin(F1,F2);
  else if(!name.compare("CB")) getF_CB(F1,F2);
  else { cerr << "invalid name for structure function" << endl; exit(1);}
  return;
}

double NuclStructure::getF1(){
  double F1,F2;
  getF(F1,F2);
  return F1;
}

double NuclStructure::getF2(){
  double F1,F2;
  getF(F1,F2);
  return F2;  
}






double NuclStructure::f2p_b(double massi, double xp, double q2, double fm) const{

  double pmo=massi;
  double f2p,yn,wwi,t,gw;
  //f2p =  0.0;
  yn  =  q2/(2*massi*xp);//nuoffshell
  //fm2 = -q2 +2.0*massi*yn + wstar2;
  if ((xp > 1.0) | (fm < massi)){return 0;}
  //	fm=pow(fm2,.5);
  //fm  = sqrt(fm2);
  //cout << fm;
  wwi = (2.*yn*pmo+1.642)/(q2 +0.376);
  t = 1.0 - 1.0/wwi;
  gw=0.256*pow(t,3)+2.178*pow(t,4)+
        0.898*pow(t,5)-6.716*pow(t,6)+3.756*pow(t,7);
  f2p = bodek(fm,q2)*gw*wwi*xp;
  return f2p;
}

// Bodek's parameterization of neutron's F2 
double NuclStructure::f2n_b(double massi, double xp, double q2, double fm) const{
  
  
  double pmo=massi;
  double f2n,yn,wwi,t,gw;
  //f2n =  0.0;
  yn  =  q2/(2*massi*xp);
  //fm2 = -q2 + 2.0*massi*yn + wstar2;
  if ((xp > 1.0) | (fm < massi)){return 0;}
  //  fm=pow(fm2,2);
  //fm = sqrt(fm2);
  //cout << fm;
  wwi = (2*yn*pmo + 1.642)/(q2 + 0.376);
  t   = 1.0-1.0/wwi;
  gw  = 0.064*pow(t,3) + 0.225*pow(t,4)+
  4.106*pow(t,5) - 7.079*pow(t,6)+3.055*pow(t,7);
  f2n = bodek(fm,q2)*gw*wwi*xp;
  return f2n;
}


/************************************************/
/*		BODEK PARAMETRIZATION		*/
/************************************************/

double NuclStructure::bodek(double wm, double qsq) const
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

