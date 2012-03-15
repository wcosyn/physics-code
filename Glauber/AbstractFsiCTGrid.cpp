#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include "AbstractFsiCTGrid.hpp"
#include "constants.hpp"
#include <Utilfunctions.hpp>

AbstractFsiCTGrid::AbstractFsiCTGrid(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleus *pnucl, string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,homedir),filledctgrid(0){
}  
  
AbstractFsiCTGrid::~AbstractFsiCTGrid(){
  //cout << "Deleting FSI object" << endl;
  
}

//get corrfilename
const string AbstractFsiCTGrid::getFsi_Ct_Filename() const{
  return fsi_ct_filename;
}

//interpolation functions
complex<double> AbstractFsiCTGrid::getFsiCtGridFull_interpvec(const TVector3 &rvec){
  return getFsiCtGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiCTGrid::getFsiCtGridFull_interp3(const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiCtGridFull_interp();
}

complex<double> AbstractFsiCTGrid::getFsiCtGridFull_interp2(const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiCtGridFull_interp();
}


complex<double> AbstractFsiCTGrid::getFsiCtGridFull_interp1(const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiCtGridFull_interp();
}

//set filenames
void AbstractFsiCTGrid::setFilenames(string homedir){
  AbstractFsiGrid::setFilenames(homedir);
  //cout << "bla " << getFsi_Filename() << endl;
  fsi_ct_filename=getFsi_Filename();
  fsi_ct_filename.insert(fsi_ct_filename.size()-4,".CT");
  for(size_t i=0;i<getParticles().size();i++){
    fsi_ct_filename.insert(fsi_ct_filename.size()-4,("."+to_string(int(getParticles()[i].getHardScale()*10))));
  }
  //cout << fsi_ct_filename << endl;
}
  
//fills the grids, uses pure virtual functions!!!!
void AbstractFsiCTGrid::fillGrids(){
  AbstractFsiGrid::fillGrids();
  setFilenames(getDir());
  if(filledallgrid){
    filledctgrid=1;
    ofstream outfile2(fsi_ct_filename.c_str(),ios::out|ios::binary);
    if(outfile2.is_open()){
      //cout << "Writing out FSI+CT grid: " << fsi_ct_filename << endl;
      writeoutFsiCtGrid(outfile2);
      outfile2.close();
      return;
    }
    else{
      cerr << "could not open file for writing FSICT grid output: " << fsi_ct_filename << endl;
    }
  }
  
  else{
    ifstream infile2(fsi_ct_filename.c_str(),ios::in|ios::binary);
    //check if object has been created sometime earlier and read it in
    if(infile2.is_open()){
      cout << "Reading in FSI+CT grid from memory: " << fsi_ct_filename << endl;
      readinFsiCtGrid(infile2);
      filledctgrid=1;
      infile2.close();
    }
    else{
      cout << "Constructing FSI+CT grid" << endl;
      constructCtGrid();
      filledctgrid=1;
      ofstream outfile(fsi_ct_filename.c_str(),ios::out|ios::binary);
      if(outfile.is_open()){
	//cout << "Writing out FSI+CT grid: " << fsi_ct_filename << endl;
	writeoutFsiCtGrid(outfile);
	outfile.close();
	return;
      }
      else{
	cerr << "could not open file for writing corrgrid output: " << fsi_ct_filename << endl;
      }    
    }
  }
}
  