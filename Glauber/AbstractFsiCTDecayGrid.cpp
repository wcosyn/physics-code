#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include "AbstractFsiCTDecayGrid.hpp"
#include "constants.hpp"
#include <Utilfunctions.hpp>

AbstractFsiCTDecayGrid::AbstractFsiCTDecayGrid(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleus *pnucl, 
					       double prec, int integrator, string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir),
AbstractFsiCTGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir){
  number_of_grids=4;
}  
  
AbstractFsiCTDecayGrid::~AbstractFsiCTDecayGrid(){
  //cout << "Deleting FSI object" << endl;
  
}

//interpolation functions
complex<double> AbstractFsiCTDecayGrid::getFsiDecayGridFull_interpvec(const TVector3 &rvec){
  return getFsiDecayGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiCTDecayGrid::getFsiDecayGridFull_interp3(const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiDecayGridFull_interp();
}

complex<double> AbstractFsiCTDecayGrid::getFsiDecayGridFull_interp2(const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiDecayGridFull_interp();
}


complex<double> AbstractFsiCTDecayGrid::getFsiDecayGridFull_interp1(const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiDecayGridFull_interp();
}

complex<double> AbstractFsiCTDecayGrid::getFsiCtDecayGridFull_interpvec(const TVector3 &rvec){
  return getFsiCtDecayGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiCTDecayGrid::getFsiCtDecayGridFull_interp3(const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiCtDecayGridFull_interp();
}

complex<double> AbstractFsiCTDecayGrid::getFsiCtDecayGridFull_interp2(const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiCtDecayGridFull_interp();
}


complex<double> AbstractFsiCTDecayGrid::getFsiCtDecayGridFull_interp1(const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiCtDecayGridFull_interp();
}

//set filenames
void AbstractFsiCTDecayGrid::setFilenames(string homedir){
  AbstractFsiCTGrid::setFilenames(homedir);
  fsi_ct_filename.insert(fsi_ct_filename.size()-4,".Decay");
  fsi_filename.insert(fsi_filename.size()-4,".Decay");
}
  
  