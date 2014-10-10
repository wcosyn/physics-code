#include "He3wf.hpp"
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <cassert>

using namespace std;

#define INVHBARC 0.00506770453255 /*!< \def defines (hbar*c)^-1 [(MeV*fm)^-1] */


//data type needed for complex (fortran data type) return 
typedef struct{double dr; double di;} complex_;

//fortran imports
extern"C"{
    //reads inputfile
    //DATNAME: import file
    //mbjset: spin state of he3 nucleus
    //usecdep: include PV parts or not (use FALSE) 
    void readhe3wave(int *mbjset, int *usecdep);
}

extern"C"{
    //resets he3 spin state
    //mbjset: spin state of he3 nucleus
    //usecdep: include PV parts or not (use FALSE) 
    void he3reset(int *mbjset, int *usecdep);
}

extern"C"{
    //returns wf state
    //input momenta are fm^-1 !!!!!!!!!!
     complex_ gethe3wave(double *p1x,double *p1y,double *p1z,int *ms1,int *mt1,
		      double *p2x,double *p2y,double *p2z,int *ms2,int *mt2,int *ms3,int *mt3);
}

extern"C"{
    char he3filename[256];
}




He3wf::He3wf(string wf_name,string dir):mA_set(-1){
  int usecdep=0;
  if(abs(mA_set)!=1) { cerr << "invalid he3 spin value" << endl; assert(1==0);}
  char *copy = new char[256];
  if(!wf_name.compare("AV18")) strcpy(he3filename,(dir+"/he_3.av18.form").c_str());
  else if(!wf_name.compare("AV18-UIX")) strcpy(he3filename,(dir+"/he_3.av18+urb-ix.jmax=6.esuche.cdep_app.fad.c.form").c_str());
  else if(!wf_name.compare("AV18-TML")) strcpy(he3filename,(dir+"/he_3.av18+tml5156.jmax=6.esuche.cdep_app.fad.c.form").c_str());
  else if(!wf_name.compare("CDBonn")) strcpy(he3filename,(dir+"/he_3.cdbonn.jmax=6.esuche.cdep_app.fad.form").c_str());
  else if(!wf_name.compare("CDBonn-TML")) strcpy(he3filename,(dir+"/he_3.cdbonn+tml4784.jmax=6.esuche.cdep_app.fad.form").c_str());
  else { cerr << "invalid name for he3 parametrisation function" << endl; assert(1==0);}
  readhe3wave(&mA_set,&usecdep);
  he3reset(&mA_set,&usecdep);
  delete []copy;
  
  return;

}
He3wf::He3wf():mA_set(-1){
}

He3wf::He3wf(const He3wf& rhs){
  mA_set=-1;
  int ex=0; he3reset(&mA_set,&ex);
}

He3wf& He3wf::operator=(const He3wf& rhs){
  if(this!=&rhs) { // avoid self-assignment
    mA_set=-1;
    int ex=0; he3reset(&mA_set,&ex); 
  }
  return *this;
  
}



complex<double> He3wf::getWF(TVector3 &p1, TVector3 &p2, int ms1, int ms2, int mt1, int mt2, int mA){
  //sanity checks
  int ex=0; he3reset(&mA,&ex);
  if(abs(mA)!=1) { cerr << "invalid he3 spin value" << endl; assert(1==0);}
  if(mA!=mA_set) {mA_set=mA; int ex=0; he3reset(&mA,&ex);}
  int ms3=mA_set-ms1-ms2;
  int mt3=1-mt1-mt2;
  if(abs(ms1)!=1||abs(ms2)!=1||abs(mt1)!=1||abs(mt2)!=1||abs(mt3)!=1) { 
    cerr << "something wrong with the spins or isospins" << endl; 
    assert(1==0);
  }
  if(abs(ms3)!=1) return 0.;
  double p1x=p1.X()*INVHBARC;
  double p1y=p1.Y()*INVHBARC;
  double p1z=p1.Z()*INVHBARC;
  double p2x=p2.X()*INVHBARC;
  double p2y=p2.Y()*INVHBARC;
  double p2z=p2.Z()*INVHBARC;  
  double p3x=-p1x-p2x;
  double p3y=-p1y-p2y;
  double p3z=-p1z-p2z;
  complex_ wfvalue1=gethe3wave(&p1x,&p1y,&p1z,&ms1,&mt1,
			    &p2x,&p2y,&p2z,&ms2,&mt2,&ms3,&mt3);
  complex_ wfvalue2=gethe3wave(&p2x,&p2y,&p2z,&ms2,&mt2,
			    &p3x,&p3y,&p3z,&ms3,&mt3,&ms1,&mt1);
  complex_ wfvalue3=gethe3wave(&p3x,&p3y,&p3z,&ms3,&mt3,
			    &p1x,&p1y,&p1z,&ms1,&mt1,&ms2,&mt2);
  return complex<double>((wfvalue1.dr+wfvalue2.dr+wfvalue3.dr)*pow(INVHBARC,3.)
			  ,(wfvalue1.di+wfvalue2.di+wfvalue3.di)*pow(INVHBARC,3.));
}

complex<double> He3wf::getWF(double p1[3], double p2[3], int ms1, int ms2, int mt1, int mt2, int mA){
  //sanity checks
  if(abs(mA)!=1) { cerr << "invalid he3 spin value" << endl; assert(1==0);}
  if(mA!=mA_set) {mA_set=mA; int ex=0; he3reset(&mA,&ex);}
  int ms3=mA_set-ms1-ms2;
  int mt3=1-mt1-mt2;
  if(abs(ms1)!=1||abs(ms2)!=1||abs(ms3)!=1||abs(mt1)!=1||abs(mt2)!=1||abs(mt3)!=1) { 
    cerr << "something wrong with the spins or isospins" << endl; 
    assert(1==0);
    
  }
  double p1x=p1[0]*INVHBARC;
  double p1y=p1[1]*INVHBARC;
  double p1z=p1[2]*INVHBARC;
  double p2x=p2[0]*INVHBARC;
  double p2y=p2[1]*INVHBARC;
  double p2z=p2[2]*INVHBARC;  
  double p3x=-p1x-p2x;
  double p3y=-p1y-p2y;
  double p3z=-p1z-p2z;
  complex_ wfvalue1=gethe3wave(&p1x,&p1y,&p1z,&ms1,&mt1,
			    &p2x,&p2y,&p2z,&ms2,&mt2,&ms3,&mt3);
  complex_ wfvalue2=gethe3wave(&p2x,&p2y,&p2z,&ms2,&mt2,
			    &p3x,&p3y,&p3z,&ms3,&mt3,&ms1,&mt1);
  complex_ wfvalue3=gethe3wave(&p3x,&p3y,&p3z,&ms3,&mt3,
			    &p1x,&p1y,&p1z,&ms1,&mt1,&ms2,&mt2);
  return complex<double>((wfvalue1.dr+wfvalue2.dr+wfvalue3.dr)*pow(INVHBARC,3.)
			  ,(wfvalue1.di+wfvalue2.di+wfvalue3.di)*pow(INVHBARC,3.));
  
}
