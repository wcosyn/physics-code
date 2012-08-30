#include <cstdlib>
#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

#include "DistMomDistrGrid.hpp"
#include "AbstractFsiCTGrid.hpp"
#include <constants.hpp>
#include <Utilfunctions.hpp>
#include <FourVector.h>
#include <Matrix.h>
#include <TSpinor.h>
#include "TMFSpinor.hpp"

DistMomDistrGrid::DistMomDistrGrid(int shell, const double p_max, const int p_grid, const int cth_grid, const int phi_grid,
				   AbstractFsiCTGrid *pfsi_grid, double precision, string homedir):
shellindex(shell),
pmax(p_max),
filledgrid(0),
filledctgrid(0),
filledallgrid(0),
dir(homedir+"/grids/"),
pgrid(p_grid),
cthgrid(cth_grid),
phigrid(phi_grid),
pfsigrid(pfsi_grid),
prec(precision){
  
  mass = getShellindex()<getPfsigrid()->getPnucleus()->getPLevels()? MASSP:MASSN;
  invpstep=pgrid/p_max;
  invcthstep=cthgrid/2.;
  invphistep=0.5*phigrid/PI;
  rhogrid=new double***[getPfsigrid()->getNumber_of_grids()/2];
  rhoctgrid=new double***[getPfsigrid()->getNumber_of_grids()/2];
  for(int i=0;i<getPfsigrid()->getNumber_of_grids()/2;i++){
    rhogrid[i]=new double**[getPgrid()+1];
    rhoctgrid[i]=new double**[getPgrid()+1];
    for(int j=0;j<(getPgrid()+1);j++){
      rhogrid[i][j]=new double*[getCthgrid()+1];
      rhoctgrid[i][j]=new double*[getCthgrid()+1];
      for(int k=0;k<(getCthgrid()+1);k++){
	rhogrid[i][j][k]=new double[getPhigrid()+1];
	rhoctgrid[i][j][k]=new double[getPhigrid()+1];
      }
    }
  }
  rhopwgrid = new double[getPgrid()+1];
  constructpwGrid();
  //fillGrids();
  
  
}  
  
DistMomDistrGrid::~DistMomDistrGrid(){
  //cout << "Deleting FSI object" << endl;
  if(rhogrid!=NULL){
    for(int i=0;i<getPfsigrid()->getNumber_of_grids()/2;i++){
      for(int j=0;j<(getPgrid()+1);j++){
	for(int k=0;k<(getCthgrid()+1);k++){
	  delete [] rhogrid[i][j][k];
	  delete [] rhoctgrid[i][j][k];
	}
	delete [] rhogrid[i][j];
	delete [] rhoctgrid[i][j];
      }
      delete [] rhogrid[i];
      delete [] rhoctgrid[i];
    }
    delete [] rhogrid;
    delete [] rhoctgrid;
  }
  delete [] rhopwgrid;
  
}

//set up inteprolation help variables for r
void DistMomDistrGrid::setPinterp(const double p){
  if(p>getPmax()){
    cerr << "p out of range: " << p << endl;
    exit(1);
  }
  pindex = int(floor(p*getInvPstep()));
  if(pindex==getPgrid()) pindex-=1;
  s_interp = (p*getInvPstep() - (pindex));
  comp_s_interp=1.-s_interp;
}

//set up interpolation help variables for theta
void DistMomDistrGrid::setCthinterp(const double costheta){
  cthindex = int(floor((-costheta+1.)*getInvCthstep()));
  if(cthindex==getCthgrid()) cthindex-=1;
  t_interp = ((-costheta+1.)*getInvCthstep() - (cthindex));  
  comp_t_interp=1.-t_interp;
}

//set up interpolation help variables for phi
void DistMomDistrGrid::setPhiinterp(const double phi){
  phiindex = int(floor(phi*getInvPhistep()));
  if(phiindex==getPhigrid()) phiindex-=1;
  u_interp = (phi*getInvPhistep() - (phiindex));  
  comp_u_interp=1.-u_interp;
}

//interpolation functions
double DistMomDistrGrid::getRhoGridFull_interpvec(int gridindex, const TVector3 &pvec){
  return getRhoGridFull_interp3(gridindex, pvec.Mag(),pvec.CosTheta(),pvec.Phi());
}

double DistMomDistrGrid::getRhoGridFull_interp3(int gridindex, const double p, const double costheta, const double phi){
  setPinterp(p);
  setCthinterp(costheta);
  setPhiinterp(phi);
  return getRhoGridFull_interp(gridindex);
}

double DistMomDistrGrid::getRhoGridFull_interp2(int gridindex, const double costheta, const double phi){
  setCthinterp(costheta);
  setPhiinterp(phi);
  return getRhoGridFull_interp(gridindex);
}

double DistMomDistrGrid::getRhoGridFull_interp1(int gridindex, const double phi){
  setPhiinterp(phi);
  return getRhoGridFull_interp(gridindex);
}

double DistMomDistrGrid::getRhoGridFull_interp(int gridindex){
  return Interp3d((gridindex < (getPfsigrid()->getNumber_of_grids()/2))? 
		    rhogrid[gridindex] : rhoctgrid[gridindex-getPfsigrid()->getNumber_of_grids()/2],
		  getS_interp(), getT_interp(), getU_interp(), 
		  getComp_s_interp(), getComp_t_interp(), getComp_u_interp(), 
		  getPindex(), getCthindex(), getPhiindex());  
}



double DistMomDistrGrid::getRhopwGridFull_interp(double p){
  if(p>getPmax()){
    cerr << "p out of range in : getRhopwGridFull" << p << endl;
    exit(1);
  }
  return interpolate(rhopwgrid,p,1./getInvPstep(),getPgrid()+1,0); 
  
}



//set filenames
void DistMomDistrGrid::setFilenames(string homedir){
  rho_filename=getPfsigrid()->getFsi_Filename();
  rhoct_filename=getPfsigrid()->getFsi_Ct_Filename();
  size_t found;

  found=rho_filename.rfind("/");
  if (found!=string::npos){
    string addendum = "Rhogrid.l"+to_string(getShellindex())+".p"+to_string(getPgrid())
		    +".cth"+to_string(getCthgrid())+".phi"+to_string(getPhigrid())+".pmax"+
		    to_string(getPmax())+".Prec"+to_string(getPrec()*1.E05)+".";
    rho_filename.insert(found+1,addendum);
    rhoct_filename.insert(found+1,addendum);
  }
//   cout << rho_filename << endl << rhoct_filename << endl;
}
  

void DistMomDistrGrid::printRhopw_grid(){
  cout << "Printing plane-wave momentum array for " << endl;
//   double tot = 0.;
  for(int i=0; i<=getPgrid(); i++){
    double p = float(i)/getInvPstep();
    cout << p << " " << getRhopwGridFull_interp(p) << endl;    
//     tot +=p*p*getRhopwGridFull_interp(p);
  }
//   tot*=4.*PI/getInvPstep()*pow(INVHBARC,3.);
//   cout << "tot " << tot << endl;
}
  
void DistMomDistrGrid::printRho_grid(int gridindex){
  cout << "Printing Distorted momentum array for " << getRho_Filename() << endl << "Gridindex " << gridindex << endl;
  for(int i=0; i<=getPgrid(); i++){
    double p = float(i)/getInvPstep();
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = k/getInvPhistep();	
	cout << p << " " << costheta << " " << phi*RADTODEGR << " " << getRhoGridFull_interp3(gridindex,p,costheta,phi) << endl;
      }
    }
  }

  
}
  
//readin fsi grid
void DistMomDistrGrid::readinRhoGrid(ifstream &infile){
  for(int i=0;i<getPfsigrid()->getNumber_of_grids()/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  infile.read(reinterpret_cast<char *>(&rhogrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}
  
void DistMomDistrGrid::readinRhoCtGrid(ifstream &infile){
  for(int i=0;i<getPfsigrid()->getNumber_of_grids()/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  infile.read(reinterpret_cast<char *>(&rhoctgrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}
  
  
//write out rho grid
void DistMomDistrGrid::writeoutRhoGrid(ofstream &outfile){
  for(int i=0;i<getPfsigrid()->getNumber_of_grids()/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  outfile.write(reinterpret_cast<char *>(&rhogrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}


//write out rho grid
void DistMomDistrGrid::writeoutRhoCtGrid(ofstream &outfile){
  for(int i=0;i<getPfsigrid()->getNumber_of_grids()/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  outfile.write(reinterpret_cast<char *>(&rhoctgrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}


//calc both fsi and fsi+ct grid
void DistMomDistrGrid::constructAllGrids(){
  //fill the arrays!
  for(int i=0; i<=getPgrid(); i++){
    p_hit = float(i)/getInvPstep();
    for(int j=0;j<=getCthgrid();j++){
      costheta_hit = -2.*j/getCthgrid()+1.;
      sintheta_hit=sqrt(1.-costheta_hit*costheta_hit);
      for(int k=0;k<=getPhigrid();k++){
	phi_hit = 2.*PI*k/getPhigrid();
	sincos(phi_hit,&sinphi_hit,&cosphi_hit);
	pvec_hit= TVector3(p_hit*sintheta_hit*cosphi_hit,p_hit*sintheta_hit*sinphi_hit,p_hit*costheta_hit);	
	for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++) rhogrid[l][i][j][k]=rhoctgrid[l][i][j][k]=0.;
	//symmetry for (-m,-ms)(m,ms), so only half needed, multiply everything by 2 in the end.
	for(int m=1;m<=getPfsigrid()->getPnucleus()->getJ_array()[shellindex];m+=2){
	  getPfsigrid()->clearKnockout();
	  getPfsigrid()->addKnockout(getShellindex(),m);
	  for(int ms=0;ms<=1;ms++){
	    FourVector<double> pf(sqrt(p_hit*p_hit+mass*mass),pvec_hit.X(),pvec_hit.Y(),pvec_hit.Z());
	    if(ms) Upm_bar=TSpinor::Bar(TSpinor(pf,mass,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
	    else Upm_bar=TSpinor::Bar(TSpinor(pf,mass,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity));
	    complex<double> results[getPfsigrid()->getNumber_of_grids()];
	    double restimate=0.,thestimate=0.,phiestimate=0.;
	    rombergerN(this,&DistMomDistrGrid::intRhoR,0.,getPfsigrid()->getPnucleus()->getRange(),getPfsigrid()->getNumber_of_grids(),
		       results,getPrec(),3,7,&restimate,m,&thestimate, &phiestimate);
	    for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++){
	      rhogrid[l][i][j][k]+=norm(results[l]);
	      rhoctgrid[l][i][j][k]+=norm(results[l+getPfsigrid()->getNumber_of_grids()/2]);
// 	      if(l==0) cout << i << " " << j << " " << k << " " << m << " " << ms << " " << l << " " 
// 	      << norm(results[l]) << " "<< norm(results[l+getPfsigrid()->getNumber_of_grids()/2]) << endl;

	    }
	  }
	}
	for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++){
	  rhogrid[l][i][j][k]/=pow(2.*PI,3.)/2.;
	  rhoctgrid[l][i][j][k]/=pow(2.*PI,3.)/2.;
	}
	//r=0 symmetry shortcut
	if(i==0){
	  for(j=0;j<=getCthgrid();j++){
	    for(k=0;k<=getPhigrid();k++){
	      for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++){
		rhogrid[l][i][j][k]=rhogrid[l][0][0][0];
		rhoctgrid[l][i][j][k]=rhoctgrid[l][0][0][0];
	      }
	    }
	  }
	}
	//theta=0 or Pi symmetry shortcut
	else if(j==0||j==getCthgrid()){
	  for(k=1;k<=getPhigrid();k++){
	    for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++){
	      rhogrid[l][i][j][k]=rhogrid[l][i][j][0];
	      rhoctgrid[l][i][j][k]=rhoctgrid[l][i][j][0];
	    }
	  }
	}
		    

	
	
      }
    }
  }
}

//calc both fsi and fsi+ct grid
void DistMomDistrGrid::constructpwGrid(){
  //cout << "Constructing pw momentum distribution" << endl;
  //fill the arrays!
  for(int i=0; i<=getPgrid(); i++){
    p_hit = float(i)/getInvPstep();
    double restimate=0.;
    double result;
    rombergerN(this,&DistMomDistrGrid::intRhoRpw,0.,getPfsigrid()->getPnucleus()->getRange(),1,
		       &result,1.E-07,3,10,&restimate);
    rhopwgrid[i] = result*result*2.*getMass()/(sqrt(getMass()*getMass()+p_hit*p_hit)+getMass())*
		    (getPfsigrid()->getPnucleus()->getJ_array()[getShellindex()]+1)/(4.*PI)*(2./PI);  //(2j+1)/(4pi)*(2/pi)
  }
}

//calc both fsi and fsi+ct grid
void DistMomDistrGrid::constructCtGrid(){
  //fill the arrays!
  for(int i=0; i<=getPgrid(); i++){
    p_hit = float(i)/getInvPstep();
    for(int j=0;j<=getCthgrid();j++){
      costheta_hit = -2.*j/getCthgrid()+1.;
      sintheta_hit=sqrt(1.-costheta_hit*costheta_hit);
      for(int k=0;k<=getPhigrid();k++){
	phi_hit = 2.*PI*k/getPhigrid();
	sincos(phi_hit,&sinphi_hit,&cosphi_hit);
	pvec_hit= TVector3(p_hit*sintheta_hit*cosphi_hit,p_hit*sintheta_hit*sinphi_hit,p_hit*costheta_hit);	
	for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++) rhoctgrid[l][i][j][k]=0.;
	for(int m=-getPfsigrid()->getPnucleus()->getJ_array()[shellindex];
	    m<=getPfsigrid()->getPnucleus()->getJ_array()[shellindex];m+=2){
	  getPfsigrid()->clearKnockout();
	  getPfsigrid()->addKnockout(getShellindex(),m);
	  for(int ms=0;ms<=1;ms++){
	    FourVector<double> pf(sqrt(p_hit*p_hit+mass*mass),pvec_hit.X(),pvec_hit.Y(),pvec_hit.Z());
	    if(ms) Upm_bar=TSpinor::Bar(TSpinor(pf,mass,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
	    else Upm_bar=TSpinor::Bar(TSpinor(pf,mass,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity));
	    complex<double> results[getPfsigrid()->getNumber_of_grids()/2];
	    double restimate=0.,thestimate=0.,phiestimate=0.;
	    rombergerN(this,&DistMomDistrGrid::intRhoRCT,0.,getPfsigrid()->getPnucleus()->getRange(),getPfsigrid()->getNumber_of_grids()/2,
		       results,getPrec(),3,7,&restimate,m,&thestimate, &phiestimate);
	    for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++){
	      rhoctgrid[l][i][j][k]+=norm(results[l]);
//  	      cout << i << " " << j << " " << k << " " << ms << " " << m << " " << l << " " << norm(results[l]) << endl;
	    }
	  }
	}
	for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++){
	  rhoctgrid[l][i][j][k]/=pow(2.*PI,3.);
	}
	//r=0 symmetry shortcut
	if(i==0){
	  for(j=0;j<=getCthgrid();j++){
	    for(k=0;k<=getPhigrid();k++){
	      for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++){
		rhoctgrid[l][i][j][k]=rhoctgrid[l][0][0][0];
	      }
	    }
	  }
	}
	//theta=0 or Pi symmetry shortcut
	else if(j==0||j==getCthgrid()){
	  for(k=1;k<=getPhigrid();k++){
	    for(int l=0;l<getPfsigrid()->getNumber_of_grids()/2;l++){
	      rhoctgrid[l][i][j][k]=rhoctgrid[l][i][j][0];
	    }
	  }
	}
		    

	
	
      }
    }
  }
}

void DistMomDistrGrid::intRhoR(const double r, complex<double> *results, va_list ap){
  
  int m = va_arg(ap,int);
  double *pthetaestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);

  getPfsigrid()->setRinterp(r);
  rombergerN(this,&DistMomDistrGrid::intRhoCosTheta,-1.,1.,getPfsigrid()->getNumber_of_grids(),
	     results,getPrec(),3,6,pthetaestimate, m, r, pphiestimate);
  for(int i=0;i<getPfsigrid()->getNumber_of_grids();i++) results[i]*=r;
  return;
}

void DistMomDistrGrid::intRhoRCT(const double r, complex<double> *results, va_list ap){
  
  int m = va_arg(ap,int);
  double *pthetaestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);

  getPfsigrid()->setRinterp(r);
  rombergerN(this,&DistMomDistrGrid::intRhoCosThetaCT,-1.,1.,getPfsigrid()->getNumber_of_grids()/2,
	     results,getPrec(),3,6,pthetaestimate, m, r, pphiestimate);
  for(int i=0;i<getPfsigrid()->getNumber_of_grids()/2;i++) results[i]*=r;
  return;
}

void DistMomDistrGrid::intRhoRpw(const double r, double *result, va_list ap){
  
  *result= r*getPfsigrid()->getPnucleus()->getWave_G(getShellindex(),r)
	*gsl_sf_bessel_jl(getPfsigrid()->getPnucleus()->getL_array()[getShellindex()],p_hit*r*INVHBARC);
}
void DistMomDistrGrid::intRhoCosTheta(const double costheta, complex<double> *results, va_list ap){
  
  int m = va_arg(ap,int);
  double r = va_arg(ap,double);
  double *pphiestimate = va_arg(ap,double*);
  double sintheta=sqrt(1.-costheta*costheta);
  getPfsigrid()->setCthinterp(costheta);
  rombergerN(this,&DistMomDistrGrid::intRhoPhi,0.,2.*PI,getPfsigrid()->getNumber_of_grids(),
	     results,getPrec(),3,6,pphiestimate, m, r, costheta, sintheta);
  return;
}

void DistMomDistrGrid::intRhoCosThetaCT(const double costheta, complex<double> *results, va_list ap){
  
  int m = va_arg(ap,int);
  double r = va_arg(ap,double);
  double *pphiestimate = va_arg(ap,double*);
  double sintheta=sqrt(1.-costheta*costheta);
  getPfsigrid()->setCthinterp(costheta);
  rombergerN(this,&DistMomDistrGrid::intRhoPhiCT,0.,2.*PI,getPfsigrid()->getNumber_of_grids()/2,
	     results,getPrec(),3,6,pphiestimate, m, r, costheta, sintheta);
  return;
}


void DistMomDistrGrid::intRhoPhi(const double phi, complex<double> *results, va_list ap){

  int m = va_arg(ap,int);
  double r = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*getPfsigrid()->getPnucleus(),getShellindex(),m,r,costheta,phi);
  complex<double> exp_pr=exp(-INVHBARC*pvec_hit*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*I)
			  *(Upm_bar*wave);  

  for(int i=0;i<getPfsigrid()->getNumber_of_grids();i++) results[i]=exp_pr*getPfsigrid()->getFsiGridN_interp1(i,phi);
  return;
}

void DistMomDistrGrid::intRhoPhiCT(const double phi, complex<double> *results, va_list ap){

  int m = va_arg(ap,int);
  double r = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*getPfsigrid()->getPnucleus(),getShellindex(),m,r,costheta,phi);
  complex<double> exp_pr=exp(-INVHBARC*pvec_hit*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*I)
			  *(Upm_bar*wave);  

  for(int i=0;i<getPfsigrid()->getNumber_of_grids()/2;i++) results[i]=exp_pr*getPfsigrid()->getFsiGridN_interp1(i+getPfsigrid()->getNumber_of_grids()/2,phi);
  return;
}


void DistMomDistrGrid::fillGrids(){
  setFilenames(dir);
  ifstream infile(rho_filename.c_str(),ios::in|ios::binary);
  //check if object has been created sometime earlier and read it in
  if(infile.is_open()){
    //cout << "Reading in FSI grid from memory: " << rho_filename << endl;
    readinRhoGrid(infile);
    filledgrid=1;
    infile.close();
  }
  else{
    cout << "Constructing all grids" << endl;
    constructAllGrids();
    filledgrid=filledallgrid=filledctgrid=1;
    ofstream outfile(rho_filename.c_str(),ios::out|ios::binary);
    if(outfile.is_open()){
      cout << "Writing out FSI grid: " << rho_filename << endl;
      writeoutRhoGrid(outfile);
      outfile.close();
    }
    else{
      cerr << "could not open file for writing fsi grid output: " << rho_filename << endl;
    }
    ofstream outfile2(rhoct_filename.c_str(),ios::out|ios::binary);
    if(outfile2.is_open()){
      cout << "Writing out FSI+CT grid: " << rhoct_filename << endl;
      writeoutRhoCtGrid(outfile2);
      outfile2.close();
      return;
    }
    else{
      cerr << "could not open file for writing FSICT grid output: " << rhoct_filename << endl;
    }
  } 
  ifstream infile2(rhoct_filename.c_str(),ios::in|ios::binary);
  //check if object has been created sometime earlier and read it in
  if(infile2.is_open()){
    //cout << "Reading in FSI+CT grid from memory: " << rhoct_filename << endl;
    readinRhoCtGrid(infile2);
    filledctgrid=1;
    infile2.close();
  }
  else{
    cout << "Constructing FSI+CT grid" << endl;
    constructCtGrid();
    filledctgrid=1;
    ofstream outfile(rhoct_filename.c_str(),ios::out|ios::binary);
    if(outfile.is_open()){
      cout << "Writing out FSI+CT grid: " << rhoct_filename << endl;
      writeoutRhoCtGrid(outfile);
      outfile.close();
      return;
    }
    else{
      cerr << "could not open file for writing corrgrid output: " << rhoct_filename << endl;
    }    
  }
  
}


void DistMomDistrGrid::updateGrids(AbstractFsiCTGrid *pfsi_grid, int shell){
  if(!filledgrid){
  fillGrids();
//     cout << "Grids still empty " << endl;
  }
  else{
    string old_rho_filename = rho_filename;
    string old_rhoct_filename = rhoct_filename;
    pfsigrid = pfsi_grid;
    shellindex=shell;
    setFilenames(dir);
  
    if(rho_filename.compare(old_rho_filename)){
//       cout << "fsi grid not equal" << rho_filename << endl << old_rho_filename << endl;
      fillGrids();      
    }
//     else cout << "fsi grid equal to the earlier one, doing nothing" << endl << rho_filename << endl << old_rho_filename << endl;
    if(old_rhoct_filename.compare(rhoct_filename)){
//     cout << "fsict grid not equal " << endl << rhoct_filename << endl << old_rhoct_filename << endl;
    
      ifstream infile2(rhoct_filename.c_str(),ios::in|ios::binary);
      //check if object has been created sometime earlier and read it in
      if(infile2.is_open()){
      // cout << "Reading in FSI+CT grid from memory: " << fsi_ct_filename << endl;
	readinRhoCtGrid(infile2);
	filledctgrid=1;
	infile2.close();
      }
      else{
	cout << "Constructing FSI+CT grid" << endl;
	constructCtGrid();
	filledctgrid=1;
	ofstream outfile(rhoct_filename.c_str(),ios::out|ios::binary);
	if(outfile.is_open()){
	  //cout << "Writing out FSI+CT grid: " << fsi_ct_filename << endl;
	  writeoutRhoCtGrid(outfile);
	  outfile.close();
	  return;
	}
	else{
	  cerr << "could not open file for writing corrgrid output: " << rhoct_filename << endl;
	}    
      }
    }
//   else cout << "fsict grid equal to the earlier one, doing nothing" << endl << rhoct_filename << endl << old_rhoct_filename << endl;

  }
}
 