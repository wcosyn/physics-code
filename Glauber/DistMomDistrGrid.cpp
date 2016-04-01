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
#include <TRotation.h>
#include "TMFSpinor.hpp"

DistMomDistrGrid::DistMomDistrGrid():
pgrid(0),
cthgrid(0),
phigrid(0),
pfsigrid(NULL),
dir(""),
rhogrid(NULL),
rhoctgrid(NULL),
rhopwgrid(NULL){
//   cerr << "This DistMomDistrGrid default constructor shouldn't be called!" << endl;
//   exit(1);
}

DistMomDistrGrid::DistMomDistrGrid(const int shell, const double p_max, const int p_grid, const int cth_grid, const int phi_grid,
				   AbstractFsiCTGrid *pfsi_grid, const double precision, const int integr, 
				   const int max_Eval, const double theta_rot, const std::string homedir):
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
prec(precision),
integrator(integr),
maxEval(max_Eval),
thetarot(theta_rot),
rhogrid(NULL),
rhoctgrid(NULL),
rhopwgrid(NULL),
number_of_Grids(pfsi_grid->getNumber_of_grids()){
  
  invpstep=pgrid/p_max;
  invcthstep=cthgrid/2.;
  invphistep=0.5*phigrid/PI;
  if(pfsigrid!=NULL){
    mass = getShellindex()<getPfsigrid()->getPnucleus()->getPLevels()? MASSP:MASSN;
    rhogrid=new double***[number_of_Grids/2];
    rhoctgrid=new double***[number_of_Grids/2];
    for(int i=0;i<number_of_Grids/2;i++){
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
  }
   //setFilenames(dir);
  //fillGrids();
  
  
}  
  
DistMomDistrGrid::~DistMomDistrGrid(){
  //cout << "Deleting FSI object" << endl;
  if(rhogrid!=NULL){
    for(int i=0;i<number_of_Grids/2;i++){
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
  
    delete [] rhopwgrid;
  }
}

DistMomDistrGrid::DistMomDistrGrid(const DistMomDistrGrid &Copy):
shellindex(Copy.getShellindex()),
rho_filename(Copy.getRho_Filename()),
rhoct_filename(Copy.getRhoCT_Filename()),
dir(Copy.getDir()),
mass(Copy.getMass()),
pmax(Copy.getPmax()),
pgrid(Copy.getPgrid()),
cthgrid(Copy.getCthgrid()),
phigrid(Copy.getPhigrid()),
invpstep(Copy.getInvPstep()),
invcthstep(Copy.getInvCthstep()),
invphistep(Copy.getInvPhistep()),
s_interp(Copy.getS_interp()),
t_interp(Copy.getT_interp()),
u_interp(Copy.getU_interp()),
comp_s_interp(Copy.getComp_s_interp()),
comp_t_interp(Copy.getComp_t_interp()),
comp_u_interp(Copy.getComp_u_interp()),
pindex(Copy.getPindex()),
cthindex(Copy.getCthindex()),
phiindex(Copy.getPhiindex()),
pfsigrid(Copy.getPfsigrid()),
filledgrid(Copy.getFilledgrid()),
filledctgrid(Copy.getFilledctgrid()),
filledallgrid(Copy.getFilledallgrid()),
pvec_hit(Copy.getPvec_hit()),
p_hit(Copy.getP_hit()),
costheta_hit(Copy.getCostheta_hit()),
sintheta_hit(Copy.getSintheta_hit()),
phi_hit(Copy.getPhi_hit()),
cosphi_hit(Copy.getCosphi_hit()),
sinphi_hit(Copy.getSinphi_hit()),
Upm_bar(Copy.getUpm_bar()),
prec(Copy.getPrec()),
integrator(Copy.getIntegrator()),
maxEval(Copy.getMaxEval()),
thetarot(Copy.getThetarot()),
rhogrid(NULL),
rhoctgrid(NULL),
rhopwgrid(NULL),
number_of_Grids(Copy.number_of_Grids){
  
  if(pfsigrid!=NULL){
    rhogrid=new double***[Copy.number_of_Grids/2];
    rhoctgrid=new double***[Copy.number_of_Grids/2];
    for(int i=0;i<Copy.number_of_Grids/2;i++){
      rhogrid[i]=new double**[Copy.getPgrid()+1];
      rhoctgrid[i]=new double**[Copy.getPgrid()+1];
      for(int j=0;j<(Copy.getPgrid()+1);j++){
	rhogrid[i][j]=new double*[Copy.getCthgrid()+1];
	rhoctgrid[i][j]=new double*[Copy.getCthgrid()+1];
	for(int k=0;k<(Copy.getCthgrid()+1);k++){
	  rhogrid[i][j][k]=new double[Copy.getPhigrid()+1];
	  rhoctgrid[i][j][k]=new double[Copy.getPhigrid()+1];
	  for(int l=0;l<(Copy.getPhigrid()+1);++l){
	    rhogrid[i][j][k][l] = Copy.getRhogrid()[i][j][k][l];
	    rhoctgrid[i][j][k][l] = Copy.getRhoctgrid()[i][j][k][l];
	  }
	}
      }
    }
    rhopwgrid = new double[Copy.getPgrid()+1];
    for(int i=0;i<(Copy.getPgrid()+1);++i) rhopwgrid[i] = Copy.getRhopwgrid()[i];
  }
  
}


DistMomDistrGrid& DistMomDistrGrid::operator=(const DistMomDistrGrid& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
    shellindex = rhs.getShellindex();
    rho_filename = rhs.getRho_Filename();
    rhoct_filename = rhs.getRhoCT_Filename();
    dir = rhs.getDir();
    mass = rhs.getMass();
    pmax = rhs.getPmax();
    pgrid = rhs.getPgrid();
    cthgrid = rhs.getCthgrid();
    phigrid = rhs.getPhigrid();
    invpstep = rhs.getInvPstep();
    invcthstep = rhs.getInvCthstep();
    invphistep = rhs.getInvPhistep();
    s_interp = rhs.getS_interp();
    t_interp = rhs.getT_interp();
    u_interp = rhs.getU_interp();
    comp_s_interp = rhs.getComp_s_interp();
    comp_t_interp = rhs.getComp_t_interp();
    comp_u_interp = rhs.getComp_u_interp();
    pindex = rhs.getPindex();
    cthindex = rhs.getCthindex();
    phiindex = rhs.getPhiindex();
    pfsigrid = rhs.getPfsigrid();
    filledgrid = rhs.getFilledgrid();
    filledctgrid = rhs.getFilledctgrid();
    filledallgrid = rhs.getFilledallgrid();
    pvec_hit = rhs.getPvec_hit();
    p_hit = rhs.getP_hit();
    costheta_hit = rhs.getCostheta_hit();
    sintheta_hit = rhs.getSintheta_hit();
    phi_hit = rhs.getPhi_hit();
    cosphi_hit = rhs.getCosphi_hit();
    sinphi_hit = rhs.getSinphi_hit();
    Upm_bar = rhs.getUpm_bar();
    prec = rhs.getPrec();
    integrator = rhs.getIntegrator();
    maxEval = rhs.getMaxEval();
    thetarot = rhs.getThetarot();
    rhogrid = NULL;
    rhoctgrid = NULL;
    rhopwgrid = NULL;
    number_of_Grids=rhs.number_of_Grids;
    
    if(pfsigrid!=NULL){
      rhogrid=new double***[rhs.number_of_Grids/2];
      rhoctgrid=new double***[rhs.number_of_Grids/2];
      for(int i=0;i<rhs.number_of_Grids/2;i++){
	rhogrid[i]=new double**[rhs.getPgrid()+1];
	rhoctgrid[i]=new double**[rhs.getPgrid()+1];
	for(int j=0;j<(rhs.getPgrid()+1);j++){
	  rhogrid[i][j]=new double*[rhs.getCthgrid()+1];
	  rhoctgrid[i][j]=new double*[rhs.getCthgrid()+1];
	  for(int k=0;k<(rhs.getCthgrid()+1);k++){
	    rhogrid[i][j][k]=new double[rhs.getPhigrid()+1];
	    rhoctgrid[i][j][k]=new double[rhs.getPhigrid()+1];
	    for(int l=0;l<(rhs.getPhigrid()+1);++l){
	      rhogrid[i][j][k][l] = rhs.getRhogrid()[i][j][k][l];
	      rhoctgrid[i][j][k][l] = rhs.getRhoctgrid()[i][j][k][l];
	    }
	  }
	}
      }
      rhopwgrid = new double[rhs.getPgrid()+1];
      for(int i=0;i<(rhs.getPgrid()+1);++i) rhopwgrid[i] = rhs.getRhopwgrid()[i];
    }
     
     
  }

  return *this;
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
  return Interp3d((gridindex < (number_of_Grids/2))? 
		    rhogrid[gridindex] : rhoctgrid[gridindex-number_of_Grids/2],
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
    string addendum = "Rhogrid.l"+to_string(getShellindex())+".throt"+to_string(int(round(thetarot*10.)))
		    +".p"+to_string(getPgrid())
		    +".cth"+to_string(getCthgrid())+".phi"+to_string(getPhigrid())+".pmax"+
		    to_string(getPmax())+".prec"+to_string(getPrec()*1.E05)+".integr"+to_string(integrator)+".";
    rho_filename.insert(found+1,addendum);
    rhoct_filename.insert(found+1,addendum);
  }
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
  for(int i=0;i<number_of_Grids/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  infile.read(reinterpret_cast<char *>(&rhogrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}
  
void DistMomDistrGrid::readinRhoCtGrid(ifstream &infile){
  for(int i=0;i<number_of_Grids/2;i++){
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
  for(int i=0;i<number_of_Grids/2;i++){
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
  for(int i=0;i<number_of_Grids/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  outfile.write(reinterpret_cast<char *>(&rhoctgrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}


//calc both fsi and fsi+ct grid
void DistMomDistrGrid::constructAllGrids(TRotation & rot){
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
	pvec_hit=pvec_hit.Transform(rot);	
	for(int l=0;l<number_of_Grids/2;l++) rhogrid[l][i][j][k]=rhoctgrid[l][i][j][k]=0.;
	//symmetry for (-m,-ms)(m,ms), so only half needed, multiply everything by 2 in the end.
	double absprec=1.E-12;
	for(int m=1;m<=getPfsigrid()->getPnucleus()->getJ_array()[shellindex];m+=2){
	  getPfsigrid()->clearKnockout();
	  getPfsigrid()->addKnockout(getShellindex(),m);
	  for(int ms=0;ms<=1;ms++){
	    FourVector<double> pf(sqrt(p_hit*p_hit+mass*mass),pvec_hit.X(),pvec_hit.Y(),pvec_hit.Z());
	    if(ms) Upm_bar=TSpinor::Bar(TSpinor(pf,mass,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
	    else Upm_bar=TSpinor::Bar(TSpinor(pf,mass,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity));
	    int res=90;
	    unsigned count=0;
	    double deeserror=0.;
	    if(integrator==0){
// 	      complex<double> results[number_of_Grids];
// 	      double restimate=0.,thestimate=0.,phiestimate=0.;
// // 	      rombergerN(this,&DistMomDistrGrid::intRhoR,0.,getPfsigrid()->getPnucleus()->getRange(),number_of_Grids,
// // 			results,getPrec(),3,7,&restimate,m,&thestimate, &phiestimate);
// 	      for(int l=0;l<number_of_Grids/2;l++){
// 		rhogrid[l][i][j][k]+=norm(results[l]);
// 		rhoctgrid[l][i][j][k]+=norm(results[l+number_of_Grids/2]);
// 	      }
	    }
	    else if(integrator==2||integrator==1){
	      numint::array<double,3> lower = {{0.,-1.,0.}};
	      numint::array<double,3> upper = {{getPfsigrid()->getPnucleus()->getRange(),1.,2.*PI}};
	      DistMomDistrGrid::Ftor_one F;
	      F.grid = this;
	      F.m = m;
	      numint::mdfunction<numint::vector_z,3> mdf;
	      mdf.func = &DistMomDistrGrid::Ftor_one::exec;
	      mdf.param = &F;
	      vector<complex<double> > ret(number_of_Grids,0.);
	      F.f=DistMomDistrGrid::klaas_distint;
	      if(integrator==1) res = numint::cube_romb(mdf,lower,upper,absprec,prec,ret,count,0);
	      else res = numint::cube_adaptive(mdf,lower,upper,absprec,prec,2E03,maxEval,ret,count,0);
	      if(abs(ret[i])*prec*0.1>absprec) absprec=abs(ret[i])*prec*0.1;
	      for(int l=0;l<number_of_Grids/2;l++){
		rhogrid[l][i][j][k]+=norm(ret[l]);
		rhoctgrid[l][i][j][k]+=norm(ret[l+number_of_Grids/2]);
	      }
	    }      
	    else {cerr  << "integrator type not implemented" << endl; exit(1);}
//   	    cout << i << " " << j << " " << k << " " << m << " " << ms << " " << rhogrid[0][i][j][k] << " " << rhogrid[1][i][j][k]<< " " << rhogrid[2][i][j][k]<< " " << rhoctgrid[0][i][j][k]<< " " << res << " " << count << endl;
	    
	  }
	}
	for(int l=0;l<number_of_Grids/2;l++){
	  //extra factor of 2 due to symmetry
	  rhogrid[l][i][j][k]/=pow(2.*PI,3.)/2.;
	  rhoctgrid[l][i][j][k]/=pow(2.*PI,3.)/2.;
	}
	//r=0 symmetry shortcut
// 	cout << p_hit << " " << costheta_hit << " " << phi_hit*RADTODEGR << " " << rhogrid[0][i][j][k] << " " << rhogrid[1][i][j][k]<< " " << rhogrid[2][i][j][k]<< " " << rhoctgrid[0][i][j][k]<< endl;
	if(i==0){
	  for(j=0;j<=getCthgrid();j++){
	    for(k=0;k<=getPhigrid();k++){
	      for(int l=0;l<number_of_Grids/2;l++){
		rhogrid[l][i][j][k]=rhogrid[l][0][0][0];
		rhoctgrid[l][i][j][k]=rhoctgrid[l][0][0][0];
	      }
	    }
	  }
	}
	//theta=0 or Pi symmetry shortcut
	else if(j==0||j==getCthgrid()){
	  for(k=1;k<=getPhigrid();k++){
	    for(int l=0;l<number_of_Grids/2;l++){
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
//     rombergerN(this,&DistMomDistrGrid::intRhoRpw,0.,getPfsigrid()->getPnucleus()->getRange(),1,
// 		       &result,1.E-07,3,10,&restimate);
    rhopwgrid[i] = result*result*2.*getMass()/(sqrt(getMass()*getMass()+p_hit*p_hit)+getMass())*
		    (getPfsigrid()->getPnucleus()->getJ_array()[getShellindex()]+1)/(4.*PI)*(2./PI);  //(2j+1)/(4pi)*(2/pi)
  }
}

//calc both fsi and fsi+ct grid
void DistMomDistrGrid::constructCtGrid(TRotation & rot){
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
	pvec_hit=pvec_hit.Transform(rot);
	for(int l=0;l<number_of_Grids/2;l++) rhoctgrid[l][i][j][k]=0.;
	double absprec=1.E-12;
	//symmetry in m_j,m_s final proton needs only half the values
	for(int m=1;m<=getPfsigrid()->getPnucleus()->getJ_array()[shellindex];m+=2){
	  getPfsigrid()->clearKnockout();
	  getPfsigrid()->addKnockout(getShellindex(),m);
	  for(int ms=0;ms<=1;ms++){
	    FourVector<double> pf(sqrt(p_hit*p_hit+mass*mass),pvec_hit.X(),pvec_hit.Y(),pvec_hit.Z());
	    if(ms) Upm_bar=TSpinor::Bar(TSpinor(pf,mass,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kDown),TSpinor::kUnity));
	    else Upm_bar=TSpinor::Bar(TSpinor(pf,mass,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity));
	    
    	    int res=90;
	    unsigned count=0;
	    double deeserror=0.;
	    if(integrator==0){
// 	      complex<double> results[number_of_Grids/2];
// 	      double restimate=0.,thestimate=0.,phiestimate=0.;
// 	      rombergerN(this,&DistMomDistrGrid::intRhoR,0.,getPfsigrid()->getPnucleus()->getRange(),number_of_Grids/2,
// 			results,getPrec(),3,7,&restimate,m,&thestimate, &phiestimate);
// 	      for(int l=0;l<number_of_Grids/2;l++){
// 		//cout << results[l] << endl;
// 		rhoctgrid[l][i][j][k]+=norm(results[l]);
// 	      if(l==0) cout << i << " " << j << " " << k << " " << m << " " << ms << " " << l << " " 
// 	      << norm(results[l]) << " "<< norm(results[l+number_of_Grids/2]) << endl;
// 	      }
	    }
	    else if(integrator==2||integrator==1){
	      numint::array<double,3> lower = {{0.,-1.,0.}};
	      numint::array<double,3> upper = {{getPfsigrid()->getPnucleus()->getRange(),1.,2.*PI}};
	      DistMomDistrGrid::Ftor_one F;
	      F.grid = this;
	      F.m = m;
	      numint::mdfunction<numint::vector_z,3> mdf;
	      mdf.func = &Ftor_one::exec;
	      mdf.param = &F;
	      vector<complex<double> > ret(number_of_Grids/2,0.);
	      F.f=DistMomDistrGrid::klaas_distint_ct;
	      if(integrator==1) res = numint::cube_romb(mdf,lower,upper,absprec,prec,ret,count,0);
	      else res = numint::cube_adaptive(mdf,lower,upper,absprec,prec,2E03,maxEval,ret,count,0);
	      if(abs(ret[i])*prec*0.1>absprec) absprec=abs(ret[i])*prec*0.1;
	      for(int l=0;l<number_of_Grids/2;l++){
		rhoctgrid[l][i][j][k]+=norm(ret[l]);
	      }
	    }      
	    else {cerr  << "integrator type not implemented" << endl; exit(1);}
// 	    cout << i << " " << j << " " << k << " " << m << " " << ms << " " << rhoctgrid[0][i][j][k] << " " << rhoctgrid[1][i][j][k]<< " " << rhoctgrid[2][i][j][k]<< " " << res << " " << count << endl;
	  }
	}
	for(int l=0;l<number_of_Grids/2;l++){
	  rhoctgrid[l][i][j][k]/=pow(2.*PI,3.)/2.;  //factor 2 because of symmetry
	}
	//r=0 symmetry shortcut
	if(i==0){
	  for(j=0;j<=getCthgrid();j++){
	    for(k=0;k<=getPhigrid();k++){
	      for(int l=0;l<number_of_Grids/2;l++){
		rhoctgrid[l][i][j][k]=rhoctgrid[l][0][0][0];
	      }
	    }
	  }
	}
	//theta=0 or Pi symmetry shortcut
	else if(j==0||j==getCthgrid()){
	  for(k=1;k<=getPhigrid();k++){
	    for(int l=0;l<number_of_Grids/2;l++){
	      rhoctgrid[l][i][j][k]=rhoctgrid[l][i][j][0];
	    }
	  }
	}
		    

	
	
      }
    }
  }
}

// void DistMomDistrGrid::intRhoR(const double r, complex<double> *results, va_list ap){
//   
//   int m = va_arg(ap,int);
//   double *pthetaestimate = va_arg(ap,double*);
//   double *pphiestimate = va_arg(ap,double*);
//   
//   
//   getPfsigrid()->setRinterp(r);
//   rombergerN(this,&DistMomDistrGrid::intRhoCosTheta,-1.,1.,number_of_Grids,
// 	     results,getPrec(),3,6,pthetaestimate, m,  r, pphiestimate);
//   for(int i=0;i<number_of_Grids;i++) results[i]*=r;
//   return;
// }
// 
// void DistMomDistrGrid::intRhoRCT(const double r, complex<double> *results, va_list ap){
//   
//   int m = va_arg(ap,int);
//   double *pthetaestimate = va_arg(ap,double*);
//   double *pphiestimate = va_arg(ap,double*);
// 
//   getPfsigrid()->setRinterp(r);
//   rombergerN(this,&DistMomDistrGrid::intRhoCosThetaCT,-1.,1.,number_of_Grids/2,
// 	     results,getPrec(),3,6,pthetaestimate, m,  r, pphiestimate);
//   for(int i=0;i<number_of_Grids/2;i++) results[i]*=r;
//   return;
// }
// 
// void DistMomDistrGrid::intRhoRpw(const double r, double *result, va_list ap){
//   
//   *result= r*getPfsigrid()->getPnucleus()->getWave_G(getShellindex(),r)
// 	*gsl_sf_bessel_jl(getPfsigrid()->getPnucleus()->getL_array()[getShellindex()],p_hit*r*INVHBARC);
// }
// void DistMomDistrGrid::intRhoCosTheta(const double costheta, complex<double> *results, va_list ap){
//   
//   int m = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double *pphiestimate = va_arg(ap,double*);
//   double sintheta=sqrt(1.-costheta*costheta);
//   getPfsigrid()->setCthinterp(costheta);
//   rombergerN(this,&DistMomDistrGrid::intRhoPhi,0.,2.*PI,number_of_Grids,
// 	     results,getPrec(),3,6,pphiestimate, m, r, costheta, sintheta);
//   return;
// }
// 
// void DistMomDistrGrid::intRhoCosThetaCT(const double costheta, complex<double> *results, va_list ap){
//   
//   int m = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double *pphiestimate = va_arg(ap,double*);
//   double sintheta=sqrt(1.-costheta*costheta);
//   getPfsigrid()->setCthinterp(costheta);
//   rombergerN(this,&DistMomDistrGrid::intRhoPhiCT,0.,2.*PI,number_of_Grids/2,
// 	     results,getPrec(),3,6,pphiestimate, m, r, costheta, sintheta);
//   return;
// }
// 
// 
// void DistMomDistrGrid::intRhoPhi(const double phi, complex<double> *results, va_list ap){
// 
//   int m = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double costheta = va_arg(ap,double);
//   double sintheta = va_arg(ap,double);
//   double cosphi,sinphi;
//   sincos(phi,&sinphi,&cosphi);
//   
//   TMFSpinor wave(*getPfsigrid()->getPnucleus(),getShellindex(),m,r,costheta,phi);
//   complex<double> exp_pr=exp(-INVHBARC*pvec_hit*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*I_UNIT)
// 			  *(Upm_bar*wave);  
// 
//   for(int i=0;i<number_of_Grids;i++) results[i]=exp_pr*getPfsigrid()->getFsiGridN_interp1(i,phi);
//   return;
// }
// 
// void DistMomDistrGrid::intRhoPhiCT(const double phi, complex<double> *results, va_list ap){
// 
//   int m = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double costheta = va_arg(ap,double);
//   double sintheta = va_arg(ap,double);
//   double cosphi,sinphi;
//   sincos(phi,&sinphi,&cosphi);
//   TMFSpinor wave(*getPfsigrid()->getPnucleus(),getShellindex(),m,r,costheta,phi);
//   
//   
//   complex<double> exp_pr=exp(-INVHBARC*pvec_hit*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*I_UNIT)
// 			  *(Upm_bar*wave);  
// 
//   for(int i=0;i<number_of_Grids/2;i++) results[i]=exp_pr*getPfsigrid()->getFsiGridN_interp1(i+number_of_Grids/2,phi);
//   return;
// }


void DistMomDistrGrid::fillGrids(TRotation & rot){
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
    constructAllGrids(rot);
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
    constructCtGrid(rot);
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


void DistMomDistrGrid::updateGrids(AbstractFsiCTGrid *pfsi_grid, int shell, TRotation & rot){
  TVector3 axis;
  rot.AngleAxis(thetarot,axis);
  thetarot*=-SIGN(axis[1]);
  if(!filledgrid){
  fillGrids(rot);
//      cout << "Grids still empty " << endl;
  }
  else{
    string old_rho_filename = rho_filename;
    string old_rhoct_filename = rhoct_filename;
    pfsigrid = pfsi_grid;
    shellindex=shell;
    setFilenames(dir);
  
    if(rho_filename.compare(old_rho_filename)){
//        cout << "fsi grid not equal" << rho_filename << endl << old_rho_filename << endl;
      fillGrids(rot);      
    }
//      else cout << "fsi grid equal to the earlier one, doing nothing" << endl << rho_filename << endl << old_rho_filename << endl;
    if(old_rhoct_filename.compare(rhoct_filename)){
//      cout << "fsict grid not equal " << endl << rhoct_filename << endl << old_rhoct_filename << endl;
    
      ifstream infile2(rhoct_filename.c_str(),ios::in|ios::binary);
      //check if object has been created sometime earlier and read it in
      if(infile2.is_open()){
      // cout << "Reading in FSI+CT grid from memory: " << fsi_ct_filename << endl;
	readinRhoCtGrid(infile2);
	filledctgrid=1;
	infile2.close();
      }
      else{
// 	cout << "Constructing FSI+CT grid" << endl;
	constructCtGrid(rot);
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
 
 
/*! integrandum function (clean ones)*/
void DistMomDistrGrid::klaas_distint(numint::vector_z &results, double r, double costheta, double phi, DistMomDistrGrid & grid, 
			  int m){
  results=numint::vector_z(grid.number_of_Grids,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*grid.getPfsigrid()->getPnucleus(),grid.getShellindex(),m,r,costheta,phi);
  complex<double> exp_pr=exp(-INVHBARC*((grid.getPvec_hit())*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta))*I_UNIT)
			  *((grid.getUpm_bar())*wave);  
  for(int i=0;i<grid.number_of_Grids;i++) results[i]=exp_pr*grid.getPfsigrid()->getFsiGridN_interp3(i,r,costheta,phi);
  for(int i=0;i<grid.number_of_Grids;i++) results[i]*=r;
  return;
}
/*! integrandum function (clean ones), only CT*/
void DistMomDistrGrid::klaas_distint_ct(numint::vector_z &results, double r, double costheta, double phi, DistMomDistrGrid & grid, 
			      int m){

  results=numint::vector_z(grid.number_of_Grids/2,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*grid.getPfsigrid()->getPnucleus(),grid.getShellindex(),m,r,costheta,phi);
  complex<double> exp_pr=exp(-INVHBARC*((grid.getPvec_hit())*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta))*I_UNIT)
			  *((grid.getUpm_bar())*wave);  

  for(int i=0;i<grid.number_of_Grids/2;i++) results[i]=exp_pr*
  grid.getPfsigrid()->getFsiGridN_interp3(i+grid.number_of_Grids/2,r,costheta,phi);
  for(int i=0;i<grid.number_of_Grids/2;i++) results[i]*=r;
  return;

}




