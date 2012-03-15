#include <complex>
#include <iostream>
#include <cstdlib>

using namespace std;

#include "OneGlauberGrid.hpp"
#include <Utilfunctions.hpp>

//constructor, calls glaubergrid's constructor
OneGlauberGrid::OneGlauberGrid(const int r_grid, const int cth_grid, MeanFieldNucleus *pnucl, string dir):
AbstractFsiGrid(r_grid,cth_grid,0,pnucl,dir),
GlauberGrid(r_grid,cth_grid,0,pnucl,dir){
}

OneGlauberGrid::~OneGlauberGrid(){
}

//2d interpolation of array
complex<double> OneGlauberGrid::getInterp(complex<double> ***grid){
  return (getComp_t_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()][0]
	+ getS_interp()*grid[getRindex()+1][getCthindex()][0]) 
	+ getT_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()+1][0] 
	+ getS_interp()*grid[getRindex()+1][getCthindex()+1][0]));
}

void OneGlauberGrid::addParticle(FastParticle& newparticle){
 if(getParticles().size()!=0){
   cerr << "You cannot add more than one particle to the Glauber phase for one particle!!!" << endl;
   exit(1);
 }
 else{
   if(newparticle.getPhi()!=0.||newparticle.getTheta()!=0.){
     cerr << "phi and theta angles of the particle should be 0 in this special case!!!" << endl;
     exit(1);
   }
   else AbstractFsiGrid::addParticle(newparticle);
 }
}


//set filenames, don't forget to add ending "/" to dir if necessary!!!
void OneGlauberGrid::setFilenames(string dir){
 AbstractFsiCTGrid::setFilenames(dir+"OneGl."); 
}




//calc glauberphases for one gridpoint
void OneGlauberGrid::calcGlauberphasesBoth(const int i, const int j, const int k){
  double *results = new double[2];
  for(int level=0;level<getPnucleus()->getTotalLevels();level++){
    for(int mm=0; mm<((getPnucleus()->getJ_array()[level]+1)/2);mm++){
      double restimate=0.,thetaestimate=0.,phiestimate=0.;
      rombergerN(this, &OneGlauberGrid::intGlauberB,1.e-06,getPnucleus()->getRange()-PREC,2,results,PREC,3,8,&restimate,level,
		 mm,&thetaestimate, &phiestimate); 
      fsi_grid[level][mm][i][j][k]=1.-getParticles()[0].getSigma(level,getPnucleus())
			  *(1.-I*getParticles()[0].getEpsilon(level,getPnucleus()))
			  /(4.*PI*getParticles()[0].getBetasq(level,getPnucleus()))
			  *results[0];
      fsi_ct_grid[level][mm][i][j][k]=1.-getParticles()[0].getSigma(level,getPnucleus())
			  *(1.-I*getParticles()[0].getEpsilon(level,getPnucleus()))
			  /(4.*PI*getParticles()[0].getBetasq(level,getPnucleus()))
			  *results[1];
    }
  }   
  delete []results;
}

//calc glauberphases for one gridpoint,only CT grid
void OneGlauberGrid::calcGlauberphasesCt(const int i, const int j, const int k){
  double result=0.;
  for(int level=0;level<getPnucleus()->getTotalLevels();level++){
    for(int mm=0; mm<((getPnucleus()->getJ_array()[level]+1)/2);mm++){
      double restimate=0.,thetaestimate=0.,phiestimate=0.;
      //2*mm+1 = m quantum number in [1/2...j] (times two)
      rombergerN(this, &OneGlauberGrid::intGlauberBCT,1.e-06,getPnucleus()->getRange(),1,&result,PREC,3,8,&restimate,level,
		 mm,&thetaestimate, &phiestimate); 
      fsi_ct_grid[level][mm][i][j][k]=1.-getParticles()[0].getSigma(level,getPnucleus())
			  *(1.-I*getParticles()[0].getEpsilon(level,getPnucleus()))
			  /(4.*PI*getParticles()[0].getBetasq(level,getPnucleus()))
			  *result;
    }
  }     
}


void OneGlauberGrid::intGlauberB(const double b, double *results, va_list ap){
  int level = va_arg(ap,int);
  int m = va_arg(ap,int);
  double *pthetaestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  
  double ceiling=sqrt(getPnucleus()->getRange()*getPnucleus()->getRange() - b*b);
  double bottom;
  //because of Heaviside function, bottom limit is z
  if(getParticles()[0].getHitz()>=-ceiling) bottom = getParticles()[0].getHitz();
  //if limit because of data range is more limiting
  else bottom = -ceiling;
  if(ceiling<bottom)  { //due to rounding errors
    results[0]= results[1] = 0.;
    return;
  }
  
  double phiint=0;
  rombergerN(this,&OneGlauberGrid::intGlauberPhi,0.,2.*PI,1,&phiint,PREC,3,8,pphiestimate,b,level);
  rombergerN(this,&OneGlauberGrid::intGlauberZ,bottom+1.E-08,ceiling-1.E-08,2,results,PREC,3,8,pthetaestimate,b, level,m);
  phiint*=b*exp(-power(b-getParticles()[0].getHitbnorm(),2.)
			      /(2.*getParticles()[0].getBetasq(level,getPnucleus())));
  results[0]*=phiint;
  results[1]*=phiint;
}

void OneGlauberGrid::intGlauberZ(const double z, double *results, va_list ap){
  double b = va_arg(ap,double);
  int level = va_arg(ap,int);
  int m = va_arg(ap,int);
  
  double r=sqrt(b*b+z*z);
  double costheta=z/r;
  
  double temp=power(getPnucleus()->getWave_F(level,r)/r,2)*getPnucleus()->getYminkappacos(level,m,costheta)
	      +power(getPnucleus()->getWave_G(level,r)/r,2)*getPnucleus()->getYkappacos(level,m,costheta);
  results[0]=temp;
  results[1]=temp*(((abs(z-getParticles()[0].getHitz()) > getParticles()[0].getLc())
			  ||(getParticles()[0].getHardScale() < getParticles()[0].getNkt_sq()))? 
			  1.:(abs(z-getParticles()[0].getHitz())/getParticles()[0].getLc() 
			  + getParticles()[0].getNkt_sq()/getParticles()[0].getHardScale()*
			  (1.-abs(z-getParticles()[0].getHitz()) / getParticles()[0].getLc())));


}

void OneGlauberGrid::intGlauberPhi(const double phi, double *result, va_list ap){
  double b = va_arg(ap,double);
  int level = va_arg(ap,int);
  
  *result=exp(-b*getParticles()[0].getHitbnorm()/getParticles()[0].getBetasq(level,getPnucleus())*2.
	     *power(sin(phi/2.),2.));
  
  
}


void OneGlauberGrid::intGlauberBCT(const double b, double *result, va_list ap){
  int level = va_arg(ap,int);
  int m = va_arg(ap,int);
  double *pthetaestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  
  double ceiling=sqrt(getPnucleus()->getRange()*getPnucleus()->getRange() - b*b);
  double bottom;
  //because of Heaviside function, bottom limit is z
  if(getParticles()[0].getHitz()>=-ceiling) bottom = getParticles()[0].getHitz();
  //if limit because of data range is more limiting
  else bottom = -ceiling;
  if(ceiling<bottom)  { //due to rounding errors
    *result = 0.;
    return;
  }
  
  double phiint=0;
  rombergerN(this,&OneGlauberGrid::intGlauberPhi,0.,2.*PI,1,&phiint,PREC,3,8,pphiestimate,b,level);
  rombergerN(this,&OneGlauberGrid::intGlauberZ,bottom+1.E-08,ceiling-1.E-08,1,result,PREC,3,8,pthetaestimate,b, level,m);
  phiint*=b*exp(-power(b-getParticles()[0].getHitbnorm(),2.)
			      /(2.*getParticles()[0].getBetasq(level,getPnucleus())));
  *result*=phiint;
}

void OneGlauberGrid::intGlauberZCT(const double z, double *result, va_list ap){
  double b = va_arg(ap,double);
  int level = va_arg(ap,int);
  int m = va_arg(ap,int);
  
  double r=sqrt(b*b+z*z);
  double costheta=z/r;
  
  double temp=power(getPnucleus()->getWave_F(level,r)/r,2)*getPnucleus()->getYminkappacos(level,m,costheta)
	      +power(getPnucleus()->getWave_G(level,r)/r,2)*getPnucleus()->getYkappacos(level,m,costheta);
  *result=temp*(((abs(z-getParticles()[0].getHitz()) > getParticles()[0].getLc())
			  ||(getParticles()[0].getHardScale() < getParticles()[0].getNkt_sq()))? 
			  1.:(abs(z-getParticles()[0].getHitz())/getParticles()[0].getLc() 
			  + getParticles()[0].getNkt_sq()/getParticles()[0].getHardScale()*
			  (1.-abs(z-getParticles()[0].getHitz()) / getParticles()[0].getLc())));
  
}



