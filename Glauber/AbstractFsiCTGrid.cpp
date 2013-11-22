#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include "AbstractFsiCTGrid.hpp"
#include "constants.hpp"
#include <Utilfunctions.hpp>

AbstractFsiCTGrid::AbstractFsiCTGrid(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleus *pnucl, 
				     double prec, int integrator, string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir),filledctgrid(0){
  number_of_grids=2;
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
  setinterp(r,costheta,phi);
  return getFsiCtGridFull_interp();
}

complex<double> AbstractFsiCTGrid::getFsiCtGridFull_interp2(const double costheta, const double phi){
  setinterp(costheta,phi);
  return getFsiCtGridFull_interp();
}


complex<double> AbstractFsiCTGrid::getFsiCtGridFull_interp1(const double phi){
  setinterp(phi);
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
  AbstractFsiGrid::fillGrids(); // fill regular FSI (RMSGA) grid, this function is implemented in AbstractFsiGrid but it uses pure virtual functions
  setFilenames(getDir());
  if(filledallgrid){
    //ct grid was filled while constructing RMSGA grid
    filledctgrid=1;
    //write to file
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
  //grid is still empty
  else{
    ifstream infile2(fsi_ct_filename.c_str(),ios::in|ios::binary);
    //check if object has been created sometime earlier and read it in
    if(infile2.is_open()){
     // cout << "Reading in FSI+CT grid from memory: " << fsi_ct_filename << endl;
      readinFsiCtGrid(infile2);
      filledctgrid=filledallgrid=1;
      infile2.close();
    }
    //we have to calc the grid
    else{
      cout << "Constructing FSI+CT grid" << endl;
      constructCtGrid(); // ! pure virtual
      filledctgrid=filledallgrid=1;
      //write it out now too
      ofstream outfile(fsi_ct_filename.c_str(),ios::out|ios::binary);
      if(outfile.is_open()){
        //cout << "Writing out FSI+CT grid: " << fsi_ct_filename << endl;
        writeoutFsiCtGrid(outfile);
        outfile.close();
        return;
      } else {
        cerr << "could not open file for writing corrgrid output: " << fsi_ct_filename << endl;
      }
    }
  }
}

//updates the grids, uses pure virtual functions!!!!
void AbstractFsiCTGrid::updateGrids(){
  AbstractFsiGrid::updateGrids(); // implemented virtual function in AbstractFsiGrid but uses pure virtual functions.
                                    // the class that will eventually call this will have inherited the virtual function
                                    // updateGrids(), which uses the [former] pure virtual functions that are now implemented
                                    // within this class and hence will be called, as virtual functions are called using the 
                                    // object rather then handle type
  string old_fsi_ct_filename = fsi_ct_filename;
  setFilenames(getDir());
  //check is something is different so we needto update the grids
  if(old_fsi_ct_filename.compare(fsi_ct_filename)){
//     cout << "fsict grid not equal " << endl << fsi_ct_filename << endl << old_fsi_ct_filename << endl;
    ifstream infile2(fsi_ct_filename.c_str(),ios::in|ios::binary);
    //check if object has been created sometime earlier and read it in
    if(infile2.is_open()){
      //cout << "Reading in FSI+CT grid from memory: " << fsi_ct_filename << endl;
      readinFsiCtGrid(infile2);
      filledctgrid=filledallgrid=1;
      infile2.close();
    }
    else{
      cout << "Constructing FSI+CT grid" << endl;
      constructCtGrid();
      filledctgrid=filledallgrid=1;
      ofstream outfile(fsi_ct_filename.c_str(),ios::out|ios::binary);
      if(outfile.is_open()){
        //cout << "Writing out FSI+CT grid: " << fsi_ct_filename << endl;
        writeoutFsiCtGrid(outfile);
        outfile.close();
        return;
      } else {
        cerr << "could not open file for writing corrgrid output: " << fsi_ct_filename << endl;
      }
    }
  }
//   else cout << "fsict grid equal to the earlier one, doing nothing" << endl << fsi_ct_filename << endl << old_fsi_ct_filename << endl;
  
}
