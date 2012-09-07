#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include "AbstractFsiGridThick.hpp"
#include <constants.hpp>
#include <Utilfunctions.hpp>

AbstractFsiGridThick::AbstractFsiGridThick(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleusThick *pnucl, 
					   double prec, int integrator, string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator, homedir),
pnucleusthick(pnucl),
fsicorrelator(pnucl,90,90,homedir){
  number_of_grids=2;
  
  //fsicorrelator.printCorrGridAll();
}  

AbstractFsiGridThick::~AbstractFsiGridThick(){
  //cout << "Deleting FSI thickness object" << endl;
  
}

//get pnucleus
MeanFieldNucleusThick *AbstractFsiGridThick::getPnucleusthick() const{
  return pnucleusthick;
}

//get fsicorrelator
const FsiCorrelator& AbstractFsiGridThick::getFsiCorrelator() const{
  return fsicorrelator;
}

//get fsicorrelator
FsiCorrelator& AbstractFsiGridThick::getFsiCorrelator(){
  return fsicorrelator;
}
//interpolation functions
complex<double> AbstractFsiGridThick::getFsiSrcGridFull_interpvec(const TVector3 &rvec){
  return getFsiSrcGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiGridThick::getFsiSrcGridFull_interp3(const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcGridFull_interp();
}

complex<double> AbstractFsiGridThick::getFsiSrcGridFull_interp2(const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcGridFull_interp();
}

complex<double> AbstractFsiGridThick::getFsiSrcGridFull_interp1(const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcGridFull_interp();
}


void AbstractFsiGridThick::setFilenames(string homedir){
  AbstractFsiGrid::setFilenames(homedir+"Thick.");
}
  
  