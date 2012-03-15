#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include "AbstractFsiCTGridThick.hpp"
#include <constants.hpp>
#include <Utilfunctions.hpp>

AbstractFsiCTGridThick::AbstractFsiCTGridThick(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleusThick *pnucl, 
					      string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,homedir),
AbstractFsiGridThick(r_grid,cth_grid,phi_grid,pnucl,homedir), 
AbstractFsiCTGrid(r_grid,cth_grid,phi_grid,pnucl,homedir){
  
}  

AbstractFsiCTGridThick::~AbstractFsiCTGridThick(){
  
}

//interpolation functions
complex<double> AbstractFsiCTGridThick::getFsiSrcCtGridFull_interpvec(const TVector3 &rvec){
  return getFsiSrcCtGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiCTGridThick::getFsiSrcCtGridFull_interp3(const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtGridFull_interp();
}

complex<double> AbstractFsiCTGridThick::getFsiSrcCtGridFull_interp2(const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtGridFull_interp();
}


complex<double> AbstractFsiCTGridThick::getFsiSrcCtGridFull_interp1(const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtGridFull_interp();
}


void AbstractFsiCTGridThick::setFilenames(string homedir){
  AbstractFsiCTGrid::setFilenames(homedir+"Thick.");
}
  