#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include "AbstractFsiDecayGridThick.hpp"
#include <constants.hpp>
#include <Utilfunctions.hpp>

AbstractFsiDecayGridThick::AbstractFsiDecayGridThick(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleusThick *pnucl, 
					   double prec, int integrator, string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator, homedir),
AbstractFsiGridThick(r_grid,cth_grid,phi_grid,pnucl,prec,integrator, homedir)
{
  number_of_grids=4;
}

AbstractFsiDecayGridThick::~AbstractFsiDecayGridThick(){
  //cout << "Deleting FSI thickness object" << endl;
}

//interpolation functions
complex<double> AbstractFsiDecayGridThick::getFsiDecayGridFull_interpvec(const TVector3 &rvec){
  return getFsiDecayGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiDecayGridThick::getFsiDecayGridFull_interp3(const double r, const double costheta, const double phi){
  setinterp(r,costheta,phi);
  return getFsiDecayGridFull_interp();
}

complex<double> AbstractFsiDecayGridThick::getFsiDecayGridFull_interp2(const double costheta, const double phi){
  setinterp(costheta,phi);
  return getFsiDecayGridFull_interp();
}

complex<double> AbstractFsiDecayGridThick::getFsiDecayGridFull_interp1(const double phi){
  setinterp(phi);
  return getFsiDecayGridFull_interp();
}

//interpolation functions
complex<double> AbstractFsiDecayGridThick::getFsiSrcDecayGridFull_interpvec(const TVector3 &rvec){
  return getFsiSrcDecayGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiDecayGridThick::getFsiSrcDecayGridFull_interp3(const double r, const double costheta, const double phi){
  setinterp(r,costheta,phi);
  return getFsiSrcDecayGridFull_interp();
}

complex<double> AbstractFsiDecayGridThick::getFsiSrcDecayGridFull_interp2(const double costheta, const double phi){
  setinterp(costheta,phi);
  return getFsiSrcDecayGridFull_interp();
}

complex<double> AbstractFsiDecayGridThick::getFsiSrcDecayGridFull_interp1(const double phi){
  setinterp(phi);
  return getFsiSrcDecayGridFull_interp();
}

void AbstractFsiDecayGridThick::setFilenames(string homedir){
  AbstractFsiGridThick::setFilenames(homedir);
  fsi_filename.insert(fsi_filename.size()-4,".Decay");
}
