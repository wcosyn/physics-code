#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include "AbstractFsiCTDecayGridThick.hpp"
#include <constants.hpp>
#include <Utilfunctions.hpp>

AbstractFsiCTDecayGridThick::AbstractFsiCTDecayGridThick(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleusThick *pnucl, 
					      double prec, int integrator, string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir),
AbstractFsiGridThick(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir), 
AbstractFsiDecayGridThick(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir), 
AbstractFsiCTGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir),
AbstractFsiCTDecayGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir){
  number_of_grids=8;
  
}

AbstractFsiCTDecayGridThick::~AbstractFsiCTDecayGridThick(){
  
}

//interpolation functions
complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtGridFull_interpvec(const TVector3 &rvec){
  return getFsiSrcCtGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtGridFull_interp3(const double r, const double costheta, const double phi){
  setinterp(r,costheta,phi);
  return getFsiSrcCtGridFull_interp();
}

complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtGridFull_interp2(const double costheta, const double phi){
  setinterp(costheta,phi);
  return getFsiSrcCtGridFull_interp();
}


complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtGridFull_interp1(const double phi){
  setinterp(phi);
  return getFsiSrcCtGridFull_interp();
}

//interpolation functions
complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtDecayGridFull_interpvec(const TVector3 &rvec){
  return getFsiSrcCtDecayGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtDecayGridFull_interp3(const double r, const double costheta, const double phi){
  setinterp(r,costheta,phi);
  return getFsiSrcCtDecayGridFull_interp();
}

complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtDecayGridFull_interp2(const double costheta, const double phi){
  setinterp(costheta,phi);
  return getFsiSrcCtDecayGridFull_interp();
}

complex<double> AbstractFsiCTDecayGridThick::getFsiSrcCtDecayGridFull_interp1(const double phi){
  setinterp(phi);
  return getFsiSrcCtDecayGridFull_interp();
}

void AbstractFsiCTDecayGridThick::setFilenames(string homedir){
  AbstractFsiCTDecayGrid::setFilenames(homedir+"Thick.");
}
