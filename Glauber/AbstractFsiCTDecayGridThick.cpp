#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include "AbstractFsiCTDecayGridThick.hpp"
#include <constants.hpp>
#include <Utilfunctions.hpp>

AbstractFsiCTDecayGridThick::AbstractFsiCTDecayGridThick(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleusThick *pnucl, 
					      string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,homedir),
AbstractFsiGridThick(r_grid,cth_grid,phi_grid,pnucl,homedir), 
AbstractFsiDecayGridThick(r_grid,cth_grid,phi_grid,pnucl,homedir), 
AbstractFsiCTGrid(r_grid,cth_grid,phi_grid,pnucl,homedir),
AbstractFsiCTDecayGrid(r_grid,cth_grid,phi_grid,pnucl,homedir){
  number_of_grids=8;
  
}  

AbstractFsiCTDecayGridThick::~AbstractFsiCTDecayGridThick(){
  
}

//interpolation functions
complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtGridFull_interpvec(const TVector3 &rvec){
  return getFsiSrcCtGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtGridFull_interp3(const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtGridFull_interp();
}

complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtGridFull_interp2(const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtGridFull_interp();
}


complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtGridFull_interp1(const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtGridFull_interp();
}


//interpolation functions
complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtDecayGridFull_interpvec(const TVector3 &rvec){
  return getFsiSrcCtDecayGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtDecayGridFull_interp3(const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtDecayGridFull_interp();
}

complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtDecayGridFull_interp2(const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtDecayGridFull_interp();
}


complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtDecayGridFull_interp1(const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtDecayGridFull_interp();
}


void AbstractFsiCTDecayGridThick::setFilenames(string homedir){
  AbstractFsiCTDecayGrid::setFilenames(homedir+"Thick.");
}
  