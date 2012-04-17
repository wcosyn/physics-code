#include <cstdlib>
#include <iostream>
#include <fstream>
#include <complex>

using namespace std;

#include "GlauberDecayGridThick.hpp"
#include <Utilfunctions.hpp>

//constructor, calls abstractfsigrid's constructor
GlauberDecayGridThick::GlauberDecayGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick *pnucl,
				   string dir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,dir),
AbstractFsiCTGrid(r_grid,cth_grid,phi_grid,pnucl,dir),
AbstractFsiGridThick(r_grid,cth_grid,phi_grid,pnucl,dir),
AbstractFsiCTDecayGridThick(r_grid,cth_grid,phi_grid,pnucl,dir),fsi_grid(NULL),fsi_ct_grid(NULL),treshold(NULL){
}

//destructor
GlauberDecayGridThick::~GlauberDecayGridThick(){
  if(fsi_grid!=NULL){
    for(int i=0;i<4;i++){
      for(int j=0;j<2;j++){
	for(int k=0;k<(getRgrid()+1);k++){
	  for(int l=0;l<(getCthgrid()+1);l++){
	    delete [] fsi_grid[i][j][k][l];
	    delete [] fsi_ct_grid[i][j][k][l];
	  }
	  delete [] fsi_grid[i][j][k];
	  delete [] fsi_ct_grid[i][j][k];
	}
	delete [] fsi_grid[i][j];
	delete [] fsi_ct_grid[i][j];
      }
      delete [] fsi_grid[i];
      delete [] fsi_ct_grid[i];
    }
    delete [] fsi_grid;
    delete [] fsi_ct_grid;
  }

}  

//interpolation of the grid after r,th,phi have been set
complex<double> GlauberDecayGridThick::getFsiGridFull_interp(){
  return pow(getInterp(fsi_grid[0][0]),getPnucleusthick()->getN()-getProtonKnockout())
	  *pow(getInterp(fsi_grid[0][1]),getPnucleusthick()->getZ()-getNeutronKnockout());
}
  
//interpolation of the grid after r,th,phi have been set
complex<double> GlauberDecayGridThick::getFsiCtGridFull_interp(){
  return pow(getInterp(fsi_ct_grid[0][0]),getPnucleusthick()->getN()-getProtonKnockout())
	  *pow(getInterp(fsi_ct_grid[0][1]),getPnucleusthick()->getZ()-getNeutronKnockout());
}

//interpolation of the grid after r,th,phi have been set
complex<double> GlauberDecayGridThick::getFsiSrcGridFull_interp(){
  return pow(getInterp(fsi_grid[1][0]),getPnucleusthick()->getN()-getProtonKnockout())
	  *pow(getInterp(fsi_grid[1][1]),getPnucleusthick()->getZ()-getNeutronKnockout());
}
  
//interpolation of the grid after r,th,phi have been set
complex<double> GlauberDecayGridThick::getFsiSrcCtGridFull_interp(){
  return pow(getInterp(fsi_ct_grid[1][0]),getPnucleusthick()->getN()-getProtonKnockout())
	  *pow(getInterp(fsi_ct_grid[1][1]),getPnucleusthick()->getZ()-getNeutronKnockout());
}


//interpolation of the grid after r,th,phi have been set
complex<double> GlauberDecayGridThick::getFsiDecayGridFull_interp(){
  return pow(getInterp(fsi_grid[2][0]),getPnucleusthick()->getN()-getProtonKnockout())
	  *pow(getInterp(fsi_grid[2][1]),getPnucleusthick()->getZ()-getNeutronKnockout());
}
  
//interpolation of the grid after r,th,phi have been set
complex<double> GlauberDecayGridThick::getFsiCtDecayGridFull_interp(){
  return pow(getInterp(fsi_ct_grid[2][0]),getPnucleusthick()->getN()-getProtonKnockout())
	  *pow(getInterp(fsi_ct_grid[2][1]),getPnucleusthick()->getZ()-getNeutronKnockout());
}

//interpolation of the grid after r,th,phi have been set
complex<double> GlauberDecayGridThick::getFsiSrcDecayGridFull_interp(){
  return pow(getInterp(fsi_grid[3][0]),getPnucleusthick()->getN()-getProtonKnockout())
	  *pow(getInterp(fsi_grid[3][1]),getPnucleusthick()->getZ()-getNeutronKnockout());
}
  
//interpolation of the grid after r,th,phi have been set
complex<double> GlauberDecayGridThick::getFsiSrcCtDecayGridFull_interp(){
  return pow(getInterp(fsi_ct_grid[3][0]),getPnucleusthick()->getN()-getProtonKnockout())
	  *pow(getInterp(fsi_ct_grid[3][1]),getPnucleusthick()->getZ()-getNeutronKnockout());
}

complex<double> GlauberDecayGridThick::getFsiGridN_interp(int grid){
  switch(grid){
    case(0):
      return getFsiGridFull_interp();
      break;
    case(1):
      return getFsiSrcGridFull_interp();
      break;
    case(2):
      return getFsiDecayGridFull_interp();
      break;
    case(3):
      return getFsiSrcDecayGridFull_interp();
      break;
    case(4):
      return getFsiCtGridFull_interp();
      break;
    case(5):
      return getFsiSrcCtGridFull_interp();
      break;
    case(6):
      return getFsiCtDecayGridFull_interp();
      break;
    case(7):
      return getFsiSrcCtDecayGridFull_interp();
      break;
    default:
      cerr << "grid number not supported" << endl;
      exit(1);
  }
}


void GlauberDecayGridThick::print_grid(int grid){
  switch(grid){
    case(0):
      return printFsi_grid();
      break;
    case(1):
      return printFsi_src_grid();
      break;
    case(2):
      return printFsi_decay_grid();
      break;
    case(3):
      return printFsi_src_decay_grid();
      break;
    case(4):
      return printFsi_ct_grid();
      break;
    case(5):
      return printFsi_src_ct_grid();
      break;
    case(6):
      return printFsi_ct_decay_grid();
      break;
    case(7):
      return printFsi_src_ct_decay_grid();
      break;
    default:
      cerr << "grid number not supported" << endl;
      exit(1);
  }
}


void GlauberDecayGridThick::printFsi_grid(){
  cout << "Printing Glauberarray for " << getFsi_Filename() << endl;
  printKnockout();
  for(int i=0; i<=getRgrid(); i++){
    double r = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi)) phi=0.;
	complex<double> value=AbstractFsiGrid::getFsiGridFull_interp3(r,costheta,phi);
	cout << r << " " << costheta << " " << phi*RADTODEGR << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }

  
}
void GlauberDecayGridThick::printFsi_ct_grid(){
  cout << "Printing CT Glauberarray for " << getFsi_Ct_Filename() << endl;
  printKnockout();
  for(int i=0; i<=getRgrid(); i++){
    double r = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi)) phi=0.;
	complex<double> value=AbstractFsiCTGrid::getFsiCtGridFull_interp3(r,costheta,phi);
	cout << r << " " << costheta << " " << phi*RADTODEGR << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }

  
}


void GlauberDecayGridThick::printFsi_src_grid(){
  cout << "Printing SRC Glauberarray for " << getFsi_Filename() << endl;
  printKnockout();
  for(int i=0; i<=getRgrid(); i++){
    double r = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi)) phi=0.;
	complex<double> value=AbstractFsiGridThick::getFsiSrcGridFull_interp3(r,costheta,phi);
	cout << r << " " << costheta << " " << phi*RADTODEGR << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }

  
}
void GlauberDecayGridThick::printFsi_src_ct_grid(){
  cout << "Printing SRC CT Glauberarray for " << getFsi_Ct_Filename() << endl;
  printKnockout();
  for(int i=0; i<=getRgrid(); i++){
    double r = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi)) phi=0.;
	complex<double> value=AbstractFsiCTDecayGridThick::getFsiSrcCtGridFull_interp3(r,costheta,phi);
	cout << r << " " << costheta << " " << phi*RADTODEGR << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }
}


void GlauberDecayGridThick::printFsi_decay_grid(){
  cout << "Printing Decay Glauberarray for " << getFsi_Filename() << endl;
  printKnockout();
  for(int i=0; i<=getRgrid(); i++){
    double r = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi)) phi=0.;
	complex<double> value=AbstractFsiCTDecayGrid::getFsiDecayGridFull_interp3(r,costheta,phi);
	cout << r << " " << costheta << " " << phi*RADTODEGR << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }

  
}
void GlauberDecayGridThick::printFsi_ct_decay_grid(){
  cout << "Printing CT Decay Glauberarray for " << getFsi_Ct_Filename() << endl;
  printKnockout();
  for(int i=0; i<=getRgrid(); i++){
    double r = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi)) phi=0.;
	complex<double> value=AbstractFsiCTDecayGrid::getFsiCtDecayGridFull_interp3(r,costheta,phi);
	cout << r << " " << costheta << " " << phi*RADTODEGR << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }

  
}


void GlauberDecayGridThick::printFsi_src_decay_grid(){
  cout << "Printing SRC Decay Glauberarray for " << getFsi_Filename() << endl;
  printKnockout();
  for(int i=0; i<=getRgrid(); i++){
    double r = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi)) phi=0.;
	complex<double> value=AbstractFsiCTDecayGridThick::getFsiSrcDecayGridFull_interp3(r,costheta,phi);
	cout << r << " " << costheta << " " << phi*RADTODEGR << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }

  
}
void GlauberDecayGridThick::printFsi_src_ct_decay_grid(){
  cout << "Printing SRC CT Decay Glauberarray for " << getFsi_Ct_Filename() << endl;
  printKnockout();
  for(int i=0; i<=getRgrid(); i++){
    double r = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi)) phi=0.;
	complex<double> value=AbstractFsiCTDecayGridThick::getFsiSrcCtDecayGridFull_interp3(r,costheta,phi);
	cout << r << " " << costheta << " " << phi*RADTODEGR << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }
}


//set filenames, don't forget to add ending "/" to dir if necessary!!!
void GlauberDecayGridThick::setFilenames(string dir){
 AbstractFsiCTDecayGridThick::setFilenames(dir+"Gl."); 
}


//calc both fsi and fsi+ct grid
void GlauberDecayGridThick::constructAllGrids(){
  //mem reservation for the grids, index [shelllevel][mindex(postive only due to symmetry][rgrid][thgrid][phgrid]
  //fsi grid
  if(fsi_grid==NULL){
    //cout << "Mem reservation for fsi_grid" << endl;
    fsi_grid=new complex<double>****[4];
    for(int i=0;i<4;i++){
      fsi_grid[i]=new complex<double>***[2];
      for(int j=0;j<2;j++){
	fsi_grid[i][j]=new complex<double>**[getRgrid()+1];
	for(int k=0;k<(getRgrid()+1);k++){
	  fsi_grid[i][j][k]=new complex<double>*[getCthgrid()+1];
	  for(int l=0;l<(getCthgrid()+1);l++){
	    fsi_grid[i][j][k][l]=new complex<double>[getPhigrid()+1];
	  }
	}
      }
    }
  }
  else{
    cerr << "fsi_grid should not be initialized more than once!!!!" << endl;
    exit(1);
  }
  //ct grid
  if(fsi_ct_grid==NULL){
    //cout << "Mem reservation for fsi_ct_grid" << endl;
    fsi_ct_grid=new complex<double>****[4];
    for(int i=0;i<4;i++){
      fsi_ct_grid[i]=new complex<double>***[2];
      for(int j=0;j<2;j++){
	fsi_ct_grid[i][j]=new complex<double>**[getRgrid()+1];
	for(int k=0;k<(getRgrid()+1);k++){
	  fsi_ct_grid[i][j][k]=new complex<double>*[getCthgrid()+1];
	  for(int l=0;l<(getCthgrid()+1);l++){
	    fsi_ct_grid[i][j][k][l]=new complex<double>[getPhigrid()+1];
	  }
	}
      }
    }
  }
  else{
    cerr << "fsi_ct_grid should not be initialized more than once!!!!" << endl;
    exit(1);
  }    
  //mem reservation for treshold array
  if(treshold==NULL){
    //cout << "Mem reservation for treshold array" << endl;
    treshold = new int*[getCthgrid()+1];
    for(int i=0;i<=getCthgrid();i++){
      treshold[i] = new int[getPhigrid()+1];
      for(int j=0;j<=getPhigrid();j++) treshold[i][j]=0;
    }
  }
  else{
    cerr << "Treshold should only be initialized once!" << endl;
    exit(1);
  }
  
  //fill the arrays!
  for(int i=0; i<=getRgrid(); i++){
    r_hit = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      costheta_hit = -2.*j/getCthgrid()+1.;
      //sincos(theta_hit,&sintheta_hit,&costheta_hit);
      sintheta_hit = sqrt(1.-costheta_hit*costheta_hit);
      for(int k=0;k<=getPhigrid();k++){
	phi_hit = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi_hit)) phi_hit=0.;
	sincos(phi_hit,&sinphi_hit,&cosphi_hit);
	if(treshold[j][k]==0){ //if treshold is true, we're close enough to one and we don't have to compute everything again and again 
	  for(size_t it=0; it<getParticles().size();it++){
	    getParticles()[it].setHitcoord(r_hit, costheta_hit,sintheta_hit, cosphi_hit, sinphi_hit);
	  }
	  calcGlauberphasesBoth(i,j,k);
		  	  
	  if(abs(pow(fsi_ct_grid[0][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_ct_grid[0][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_ct_grid[1][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_ct_grid[1][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_ct_grid[2][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_ct_grid[2][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_ct_grid[3][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_ct_grid[3][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_grid[0][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_grid[0][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_grid[1][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_grid[1][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_grid[2][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_grid[2][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_grid[3][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_grid[3][1][i][j][k],getPnucleusthick()->getZ()))>0.99){
	    treshold[j][k]=1;
	    if(j==0||j==getCthgrid()){
	      for(int ll=1;ll<=getPhigrid();ll++) treshold[j][ll]=1;
	    }
	  }
	  //r=0 symmetry shortcut
	  if(i==0){
	    for(j=0;j<=getCthgrid();j++){
	      for(k=0;k<=getPhigrid();k++){
		for(int src=0;src<4;src++){
		  for(int proton=0; proton<2;proton++){
		    fsi_grid[src][proton][i][j][k]=fsi_grid[src][proton][0][0][0];
		    fsi_ct_grid[src][proton][i][j][k]=fsi_ct_grid[src][proton][0][0][0];		    
		  }
		}
		
	      }
	    }
	  }
	  //theta=0 or Pi symmetry shortcut
	  else if(j==0||j==getCthgrid()){
	    for(k=1;k<=getPhigrid();k++){
	      for(int src=0;src<4;src++){
		for(int proton=0; proton<2;proton++){
		  fsi_grid[src][proton][i][j][k]=fsi_grid[src][proton][i][j][0];
		  fsi_ct_grid[src][proton][i][j][k]=fsi_ct_grid[src][proton][i][j][0];
		}
	      }
	    }
	  }
		    

	}
	//treshold reached, fill with 1s
	else{
	  for(int s=0; s<4; s++){
	    for(int mm=0; mm<2; mm++){
	      fsi_grid[s][mm][i][j][k] = 1.;
	      fsi_ct_grid[s][mm][i][j][k] = 1.;
	    }
	  }
	
	}
	
      }
    }
  }
  //mem cleanup
  for(int i=0;i<=getCthgrid();i++) delete [] treshold[i];
  delete [] treshold;
  treshold=NULL;
}


//calc only ct grid
void GlauberDecayGridThick::constructCtGrid(){
  //ct grid
  if(fsi_ct_grid==NULL){
    cout << "Mem reservation for fsi_ct_grid" << endl;
    fsi_ct_grid=new complex<double>****[4];
    for(int i=0;i<4;i++){
      fsi_ct_grid[i]=new complex<double>***[2];
      for(int j=0;j<2;j++){
	fsi_ct_grid[i][j]=new complex<double>**[getRgrid()+1];
	for(int k=0;k<(getRgrid()+1);k++){
	  fsi_ct_grid[i][j][k]=new complex<double>*[getCthgrid()+1];
	  for(int l=0;l<(getCthgrid()+1);l++){
	    fsi_ct_grid[i][j][k][l]=new complex<double>[getPhigrid()+1];
	  }
	}
      }
    }
  }
  else{
    cerr << "fsi_ct_grid should not be initialized more than once!!!!" << endl;
    exit(1);
  }    
  //mem reservation for treshold array
  if(treshold==NULL){
    cout << "Mem reservation for treshold array" << endl;
    treshold = new int*[getCthgrid()+1];
    for(int i=0;i<=getCthgrid();i++){
      treshold[i] = new int[getPhigrid()+1];
      for(int j=0;j<=getPhigrid();j++) treshold[i][j]=0;
    }
  }
  else{
    cerr << "Treshold should only be initialized once!" << endl;
    exit(1);
  }
  
  //fill the arrays!
  for(int i=0; i<=getRgrid(); i++){
    r_hit = float(i)*getPnucleusthick()->getRange()/getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for(int j=0;j<=getCthgrid();j++){
      costheta_hit = -2.*j/getCthgrid()+1.;
      //sincos(theta_hit,&sintheta_hit,&costheta_hit);
      sintheta_hit = sqrt(1.-costheta_hit*costheta_hit);      
      for(int k=0;k<=getPhigrid();k++){
	phi_hit = (getAllinplane()?1.:2.)*PI*k/getPhigrid();
	if(isnan(phi_hit)) phi_hit=0.;
	sincos(phi_hit,&sinphi_hit,&cosphi_hit);
	if(treshold[j][k]==0){ //if treshold is true, we're close enough to one and we don't have to compute everything again and again 
	  for(size_t it=0; it<getParticles().size();it++){
	    getParticles()[it].setHitcoord(r_hit, costheta_hit,sintheta_hit, cosphi_hit, sinphi_hit);
	  }
	  calcGlauberphasesCt(i,j,k);
	
	  if(abs(pow(fsi_ct_grid[0][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_ct_grid[0][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_ct_grid[1][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_ct_grid[1][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_ct_grid[2][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_ct_grid[2][1][i][j][k],getPnucleusthick()->getZ()))>0.99
		  &&abs(pow(fsi_ct_grid[3][0][i][j][k],getPnucleusthick()->getN())*pow(fsi_ct_grid[3][1][i][j][k],getPnucleusthick()->getZ()))>0.99){
	    treshold[j][k]=1;
	    if(j==0||j==getCthgrid()){
	      for(int ll=1;ll<=getPhigrid();ll++) treshold[j][ll]=1;
	    }
	  }
	  //r=0 symmetry shortcut
	  if(i==0){
	    for(j=0;j<=getCthgrid();j++){
	      for(k=0;k<=getPhigrid();k++){
		for(int src=0;src<4;src++){
		  for(int proton=0; proton<2;proton++){
		    fsi_ct_grid[src][proton][i][j][k]=fsi_ct_grid[src][proton][0][0][0];		    
		  }
		}
	      }
	    }
	  }
	  //theta=0 or Pi symmetry shortcut
	  else if(j==0||j==getCthgrid()){
	    for(k=1;k<=getPhigrid();k++){
	      for(int src=0;src<4;src++){
		for(int proton=0; proton<2;proton++){
		  fsi_ct_grid[src][proton][i][j][k]=fsi_ct_grid[src][proton][i][j][0];
		}
	      }
	    }
	  }
		    

	}
	//treshold reached, fill with 1s
	else{
	  for(int s=0; s<4; s++){
	    for(int mm=0; mm<2; mm++){
	      fsi_ct_grid[s][mm][i][j][k] = 1.;
	    }
	  }
	}
	
      }
    }
  }
  //mem cleanup
  for(int i=0;i<=getCthgrid();i++) delete [] treshold[i];
  delete [] treshold;
  treshold=NULL;
}

//calc glauberphases for one gridpoint
void GlauberDecayGridThick::calcGlauberphasesBoth(const int i, const int j, const int k){
  complex<double>* results= new complex<double>[8];
  for(int proton=0;proton<2;proton++){
    double src=getFsiCorrelator().getCorrGrid_interp(r_hit,costheta_hit,proton);
    double restimate=0.,thetaestimate=0.,phiestimate=0.;
    rombergerN(this, &GlauberDecayGridThick::intGlauberR,0.,getPnucleusthick()->getRange(),8,results,PREC,3,8,&restimate,proton,&thetaestimate, &phiestimate); 
    fsi_grid[0][proton][i][j][k]=results[0];
    fsi_grid[1][proton][i][j][k]=results[1]*src;
    fsi_grid[2][proton][i][j][k]=results[2];
    fsi_grid[3][proton][i][j][k]=results[3]*src;
    fsi_ct_grid[0][proton][i][j][k]=results[4];
    fsi_ct_grid[1][proton][i][j][k]=results[5]*src;
    fsi_ct_grid[2][proton][i][j][k]=results[6];
    fsi_ct_grid[3][proton][i][j][k]=results[7]*src;    
  }   
  delete[] results;
}

//calc glauberphases for one gridpoint,only CT grid
void GlauberDecayGridThick::calcGlauberphasesCt(const int i, const int j, const int k){
  complex<double>* results= new complex<double>[4];
  for(int proton=0;proton<2;proton++){
    double restimate=0.,thetaestimate=0.,phiestimate=0.;
    rombergerN(this, &GlauberDecayGridThick::intGlauberRCT,0.,getPnucleusthick()->getRange(),4,results,PREC,3,8,&restimate,proton,&thetaestimate, &phiestimate); 
    double src=getFsiCorrelator().getCorrGrid_interp(r_hit,costheta_hit,proton);
    fsi_ct_grid[0][proton][i][j][k]=results[0];
    fsi_ct_grid[1][proton][i][j][k]=results[1]*src;
    fsi_ct_grid[2][proton][i][j][k]=results[2];
    fsi_ct_grid[3][proton][i][j][k]=results[3]*src;
  }   
  delete[] results;
}



//readin fsi grid
void GlauberDecayGridThick::readinFsiGrid(ifstream &infile){
  if(fsi_grid==NULL){
    fsi_grid=new complex<double>****[4];
    for(int i=0;i<4;i++){
      fsi_grid[i]=new complex<double>***[2];
      for(int j=0;j<2;j++){
	fsi_grid[i][j]=new complex<double>**[getRgrid()+1];
	for(int k=0;k<(getRgrid()+1);k++){
	  fsi_grid[i][j][k]=new complex<double>*[getCthgrid()+1];
	  for(int l=0;l<(getCthgrid()+1);l++){
	    fsi_grid[i][j][k][l]=new complex<double>[getPhigrid()+1];
	    for(int mm=0;mm<(getPhigrid()+1);mm++){
	      infile.read(reinterpret_cast<char *>(&fsi_grid[i][j][k][l][mm]),sizeof(complex<double>));
	    }
	  }
	}
      }
    }
  }
  else{
    cerr << "fsi_grid should not be initialized more than once!!!!" << endl;
    exit(1);
  }
}

//readin fsi+ct grid
void GlauberDecayGridThick::readinFsiCtGrid(ifstream &infile){
  if(fsi_ct_grid==NULL){
    fsi_ct_grid=new complex<double>****[4];
    for(int i=0;i<4;i++){
      fsi_ct_grid[i]=new complex<double>***[2];
      for(int j=0;j<2;j++){
	fsi_ct_grid[i][j]=new complex<double>**[getRgrid()+1];
	for(int k=0;k<(getRgrid()+1);k++){
	  fsi_ct_grid[i][j][k]=new complex<double>*[getCthgrid()+1];
	  for(int l=0;l<(getCthgrid()+1);l++){
	    fsi_ct_grid[i][j][k][l]=new complex<double>[getPhigrid()+1];
	    for(int mm=0;mm<(getPhigrid()+1);mm++){
	      infile.read(reinterpret_cast<char *>(&fsi_ct_grid[i][j][k][l][mm]),sizeof(complex<double>));
	    }
	  }
	}
      }
    }  
  }
  else{
    cerr << "fsi_ct_grid should not be initialized more than once!!!!" << endl;
    exit(1);
  }    
}

//write fsi grid to file
void GlauberDecayGridThick::writeoutFsiGrid(ofstream &outfile){
  for(int i=0;i<4;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<(getRgrid()+1);k++){
	for(int l=0;l<(getCthgrid()+1);l++){
	  for(int mm=0;mm<(getPhigrid()+1);mm++){
	    outfile.write(reinterpret_cast<char *>(&fsi_grid[i][j][k][l][mm]),sizeof(complex<double>));
	  }
	}
      }
    }
  }  
}

//write fsi+ct grid to file
void GlauberDecayGridThick::writeoutFsiCtGrid(ofstream &outfile){
  for(int i=0;i<4;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<(getRgrid()+1);k++){
	for(int l=0;l<(getCthgrid()+1);l++){
	  for(int mm=0;mm<(getPhigrid()+1);mm++){
	    outfile.write(reinterpret_cast<char *>(&fsi_ct_grid[i][j][k][l][mm]),sizeof(complex<double>));
	  }
	}
      }
    }
  }    
}
  
void GlauberDecayGridThick::intGlauberR(const double r, complex<double> *results, va_list ap){
  int proton = va_arg(ap,int);
  double *pthetaestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  rombergerN(this,&GlauberDecayGridThick::intGlauberCosTheta,-1.,1.,8,results,PREC,3,8,pthetaestimate,r, proton,pphiestimate);
  double dens=getPnucleusthick()->getDensity(r,proton);
  getFsiCorrelator().setRinterp(r);
  for(int i=0;i<8;i++) results[i]*=dens;
}

void GlauberDecayGridThick::intGlauberCosTheta(const double costheta, complex<double>* results, va_list ap){
  double r = va_arg(ap,double);
  int proton = va_arg(ap,int);
  double *pphiestimate = va_arg(ap,double*);  
  
  double sintheta = sqrt(1.-costheta*costheta);
  rombergerN(this,&GlauberDecayGridThick::intGlauberPhi,0.,2.*PI,8,results,PREC,3,5,pphiestimate,r,costheta,sintheta,proton);
  //cout << r << " " << acos(costheta)*RADTODEGR << " " << getFsiCorrelator()->getRindex() << endl;
  double src=getFsiCorrelator().getCorrGrid_interp(costheta,proton);
  results[1]*=src;
  results[3]*=src;
  results[5]*=src;
  results[7]*=src;
}
  
void GlauberDecayGridThick::intGlauberPhi(const double phi, complex<double>* results, va_list ap){
  
  double r = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  int proton = va_arg(ap,int);
  
  double cosphi, sinphi;
  sincos(phi,&sinphi,&cosphi);
  results[0]=results[2]=results[4]=results[6]=1.;
  
  for(size_t it=0;it<getParticles().size();it++){
    double zmom=getParticles()[it].calcZ(r,costheta,sintheta,cosphi,sinphi);
    bool check=(zmom>getParticles()[it].getHitz());
    if(getParticles()[it].getIncoming()) check=!check;
    if(check){
      complex<double> temp=getParticles()[it].getScatterfront(proton,getPnucleusthick())
			    *exp(-getParticles()[it].getBdist(r,costheta,sintheta,cosphi,sinphi,zmom)
			          /(2.*getParticles()[it].getBetasq(proton,getPnucleusthick())));
      results[0]*=1.-temp;
      results[2]*=1.-temp*getParticles()[it].getDecay_sigma(zmom,proton,getPnucleusthick());
      results[4]*=1.-temp*getParticles()[it].getCTsigma(zmom);
      results[6]*=1.-temp*getParticles()[it].getDecay_sigma(zmom,proton,getPnucleusthick())*getParticles()[it].getCTsigma(zmom);
    }
  }
  double gr=getFsiCorrelator().correlation(normr(r,costheta,sintheta,cosphi,sinphi,r_hit,costheta_hit,sintheta_hit,cosphi_hit,sinphi_hit));
  results[1]=results[0]*gr;
  results[3]=results[2]*gr;
  results[5]=results[4]*gr;
  results[7]=results[6]*gr;
}

void GlauberDecayGridThick::intGlauberRCT(const double r, complex<double> *results, va_list ap){
  int proton = va_arg(ap,int);
  double *pthetaestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  getFsiCorrelator().setRinterp(r);
  rombergerN(this,&GlauberDecayGridThick::intGlauberCosThetaCT,-1.,1.,4,results,PREC,3,8,pthetaestimate,r, proton,pphiestimate);
  double dens=getPnucleusthick()->getDensity(r,proton);
  for(int i=0;i<4;i++) results[i]*=dens;
}

void GlauberDecayGridThick::intGlauberCosThetaCT(const double costheta, complex<double>* results, va_list ap){
  double r = va_arg(ap,double);
  int proton = va_arg(ap,int);
  double *pphiestimate = va_arg(ap,double*);  
  
  double sintheta = sqrt(1.-costheta*costheta);
  rombergerN(this,&GlauberDecayGridThick::intGlauberPhiCT,0.,2.*PI,4,results,PREC,3,5,pphiestimate,r,costheta,sintheta,proton);
  double src=getFsiCorrelator().getCorrGrid_interp(costheta,proton);
  results[1]*=src;
  results[3]*=src;
}
  
void GlauberDecayGridThick::intGlauberPhiCT(const double phi, complex<double>* results, va_list ap){
  
  double r = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  int proton = va_arg(ap,int);
  
  double cosphi, sinphi;
  sincos(phi,&sinphi,&cosphi);
  results[0]=results[2]=1.;
  
  for(size_t it=0;it<getParticles().size();it++){
    double zmom=getParticles()[it].calcZ(r,costheta,sintheta,cosphi,sinphi);
    bool check=(zmom>getParticles()[it].getHitz());
    if(getParticles()[it].getIncoming()) check=!check;
    if(check){
      complex<double> temp=getParticles()[it].getScatterfront(proton,getPnucleusthick())
			    *exp(-getParticles()[it].getBdist(r,costheta,sintheta,cosphi,sinphi,zmom)
			          /(2.*getParticles()[it].getBetasq(proton,getPnucleusthick())));
      results[0]*=1.-temp*getParticles()[it].getCTsigma(zmom);
      results[2]*=1.-temp*getParticles()[it].getDecay_sigma(zmom,proton,getPnucleusthick())*getParticles()[it].getCTsigma(zmom);
    }
  }
  double gr=getFsiCorrelator().correlation(normr(r,costheta,sintheta,cosphi,sinphi,r_hit,costheta_hit,sintheta_hit,cosphi_hit,sinphi_hit));
  results[1]=results[0]*gr;
  results[3]=results[2]*gr;
  
}

	 

