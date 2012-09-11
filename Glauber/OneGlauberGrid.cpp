#include <complex>
#include <iostream>
#include <cstdlib>

using namespace std;

#include "OneGlauberGrid.hpp"
#include <Utilfunctions.hpp>

//constructor, calls glaubergrid's constructor
OneGlauberGrid::OneGlauberGrid(const int r_grid, const int cth_grid, MeanFieldNucleus *pnucl, 
			       double prec, int integrator, string dir):
AbstractFsiGrid(r_grid,cth_grid,0,pnucl,prec,integrator,dir),
GlauberGrid(r_grid,cth_grid,0,pnucl,prec,integrator,dir){
}

OneGlauberGrid::~OneGlauberGrid(){
}

//2d interpolation of array
complex<double> OneGlauberGrid::getInterp(complex<double> ***grid){
  return (
    complex<double>(getComp_t_interp()*(getComp_s_interp()*real(grid[getRindex()][getCthindex()][0])
	+ getS_interp()*real(grid[getRindex()+1][getCthindex()][0]))
	+ getT_interp()*(getComp_s_interp()*real(grid[getRindex()][getCthindex()+1][0])
	+ getS_interp()*real(grid[getRindex()+1][getCthindex()+1][0])),
	getComp_t_interp()*(getComp_s_interp()*imag(grid[getRindex()][getCthindex()][0])
	+ getS_interp()*imag(grid[getRindex()+1][getCthindex()][0])) 
	+ getT_interp()*(getComp_s_interp()*imag(grid[getRindex()][getCthindex()+1][0])
	+ getS_interp()*imag(grid[getRindex()+1][getCthindex()+1][0]))
 		  ));
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
  double results[2]={0.,0.};
  for(int level=0;level<getPnucleus()->getTotalLevels();level++){
    for(int mm=0; mm<((getPnucleus()->getJ_array()[level]+1)/2);mm++){
      int res=90;
      unsigned count=0;
      double deeserror=0.;
      if(integrator==0){
	double restimate=0.,thetaestimate=0.,phiestimate=0.;
	rombergerN(this, &OneGlauberGrid::intGlauberB,1.e-06,getPnucleus()->getRange()-getPrec(),2,results,getPrec(),3,8,&restimate,level,
		  mm,&thetaestimate, &phiestimate); 
      }
      else if(integrator==1||integrator==2){
	numint::array<double,3> lower = {{1.E-06,getParticles()[0].getHitz(),0.}};
	numint::array<double,3> upper = {{getPnucleus()->getRange(),getPnucleus()->getRange(),0.5*PI}};
	Ftor_one F;
	F.grid = this;
	F.level = level;
	F.mm = mm;
	numint::mdfunction<numint::vector_d,3> mdf;
	mdf.func = &Ftor_one::exec;
	mdf.param = &F;
	vector<double> ret(2,0.);
	F.f=klaas_one_bound;
	if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
	else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,ret,count,0);
	results[0]=ret[0]; results[1]=ret[1];
	
      }      
      else {cerr  << "integrator type not implemented" << endl; exit(1);}
      
      fsi_grid[level][mm][i][j][k]=1.-getParticles()[0].getScatterfront(level,getPnucleus())*results[0];
      fsi_ct_grid[level][mm][i][j][k]=1.-getParticles()[0].getScatterfront(level,getPnucleus())*results[1];
//       cout << i << " " << j << " " << k << " " << level << " " << mm << fsi_grid[level][mm][i][j][k] << " " << fsi_ct_grid[level][mm][i][j][k] << 
//       " " << res << " " << count << endl;
    }
  }   
  
}

//calc glauberphases for one gridpoint,only CT grid
void OneGlauberGrid::calcGlauberphasesCt(const int i, const int j, const int k){
  double results=0.;
  for(int level=0;level<getPnucleus()->getTotalLevels();level++){
    for(int mm=0; mm<((getPnucleus()->getJ_array()[level]+1)/2);mm++){
      int res=90;
      unsigned count=0;
      double deeserror=0.;
      if(integrator==0){
	double restimate=0.,thetaestimate=0.,phiestimate=0.;
	rombergerN(this, &OneGlauberGrid::intGlauberBCT,1.e-06,getPnucleus()->getRange()-getPrec(),1,&results,getPrec(),3,8,&restimate,level,
		  mm,&thetaestimate, &phiestimate); 
      }
      else if(integrator==1||integrator==2){
	numint::array<double,3> lower = {{1.E-06,getParticles()[0].getHitz(),0.}};
	numint::array<double,3> upper = {{getPnucleus()->getRange(),getPnucleus()->getRange(),0.5*PI}};
	Ftor_one F;
	F.grid = this;
	F.level = level;
	F.mm = mm;
	numint::mdfunction<numint::vector_d,3> mdf;
	mdf.func = &Ftor_one::exec;
	mdf.param = &F;
	vector<double> ret(2,0.);
	F.f=klaas_one_bound_ct;
	if(integrator==1) res = numint::cube_romb(mdf,lower,upper,1.E-08,prec,ret,count,0);
	else res = numint::cube_adaptive(mdf,lower,upper,1.E-08,prec,ret,count,0);
	results=ret[0];
	
      }      
      else {cerr  << "integrator type not implemented" << endl; exit(1);}
      
      fsi_ct_grid[level][mm][i][j][k]=1.-getParticles()[0].getScatterfront(level,getPnucleus())*results;
//       cout << i << " " << j << " " << k << " " << level << " " << mm << fsi_grid[level][mm][i][j][k] << " " << fsi_ct_grid[level][mm][i][j][k] << 
//       " " << res << " " << count << endl;
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
  rombergerN(this,&OneGlauberGrid::intGlauberPhi,0.,PI/2.,1,&phiint,getPrec(),3,8,pphiestimate,b,level);
  rombergerN(this,&OneGlauberGrid::intGlauberZ,bottom+1.E-08,ceiling-1.E-08,2,results,getPrec(),3,8,pthetaestimate,b, level,m);
  phiint*=4.*b*exp(-power(b-getParticles()[0].getHitbnorm(),2.)
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
  results[1]=temp*getParticles()[0].getCTsigma(z);


}

void OneGlauberGrid::intGlauberPhi(const double phi, double *result, va_list ap){
  double b = va_arg(ap,double);
  int level = va_arg(ap,int);
  
  *result=exp(-b*getParticles()[0].getHitbnorm()/getParticles()[0].getBetasq(level,getPnucleus())*2.
	     *power(sin(phi),2.));
  
  
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
  rombergerN(this,&OneGlauberGrid::intGlauberPhi,0.,PI/2.,1,&phiint,getPrec(),3,8,pphiestimate,b,level);
  rombergerN(this,&OneGlauberGrid::intGlauberZCT,bottom+1.E-08,ceiling-1.E-08,1,result,getPrec(),3,8,pthetaestimate,b, level,m);
  phiint*=4.*b*exp(-power(b-getParticles()[0].getHitbnorm(),2.)
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
  *result=temp*getParticles()[0].getCTsigma(z);
  
}

void klaas_one_bound(numint::vector_d & results, double b, double z, double phi, OneGlauberGrid & grid, 
		   int level, int m){
  //double sintheta=sqrt(1.-costheta*costheta);
  results=numint::vector_d(2,0.);
  double r=sqrt(b*b+z*z);
  if(r>grid.getPnucleus()->getRange()) return;
  double costheta=z/r;
  double sintheta=b/r;
  double cosphi, sinphi;
  sincos(phi,&sinphi,&cosphi);
  
  results[0]=power(grid.getPnucleus()->getWave_F(level,r)/r,2)*grid.getPnucleus()->getYminkappacos(level,m,costheta)
	      +power(grid.getPnucleus()->getWave_G(level,r)/r,2)*grid.getPnucleus()->getYkappacos(level,m,costheta)
	      *4.*b*exp(-(power(b-grid.getParticles()[0].getHitbnorm(),2)
	      +b*grid.getParticles()[0].getHitbnorm()*4.*sinphi*sinphi)
	      /(2.*grid.getParticles()[0].getBetasq(level,grid.getPnucleus())));
  results[1]=results[0]*grid.getParticles()[0].getCTsigma(z);
  return;
  
}

void klaas_one_bound_ct(numint::vector_d & results, double b, double z, double phi, OneGlauberGrid & grid, 
		   int level, int m){
  //double sintheta=sqrt(1.-costheta*costheta);
  results=numint::vector_d(1,0.);
  double r=sqrt(b*b+z*z);
  if(r>grid.getPnucleus()->getRange()) return;
  double costheta=z/r;
  double sintheta=b/r;
  double cosphi, sinphi;
  sincos(phi,&sinphi,&cosphi);
  
  results[0]=power(grid.getPnucleus()->getWave_F(level,r)/r,2)*grid.getPnucleus()->getYminkappacos(level,m,costheta)
	      +power(grid.getPnucleus()->getWave_G(level,r)/r,2)*grid.getPnucleus()->getYkappacos(level,m,costheta)
	      *4.*b*exp(-(power(b-grid.getParticles()[0].getHitbnorm(),2)
	      +b*grid.getParticles()[0].getHitbnorm()*4.*sinphi*sinphi)
	      /(2.*grid.getParticles()[0].getBetasq(level,grid.getPnucleus())))*grid.getParticles()[0].getCTsigma(z);
  return;
  
}



