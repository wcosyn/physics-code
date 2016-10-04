#include <cstdlib>
#include <iostream>
#include <fstream>
#include <complex>
#include <cassert>

using namespace std;

#include "ROMEAGrid.hpp"
#include <Utilfunctions.hpp>
#include <AuxFunction.hpp>

//optical potentials from most recent Cooper, Hama parametrizations (Phys.Rev. C80 (2009) 034605)
extern"C"{
    void global_(int *ifit,double *tplab, double *aa, double r[1000], double vv[1000], double wv[1000], 
		 double vs[1000], double ws[1000], int *nstep);
}


//constructor, calls abstractfsigrid's constructor
ROMEAGrid::ROMEAGrid(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleus * pnucl,
                         double prec, int integrator, fit_type opticalfit, string dir)
  :AbstractFsiGrid(r_grid, cth_grid, phi_grid, pnucl, prec, integrator, dir),
  opticalsize(1000),opticalstep(pnucl->getRange()/1000.),fit(opticalfit),fsi_grid(NULL),treshold(NULL){

  //mem reservation for the grids, index [rgrid][thgrid][phgrid]
  //fsi grid
  if (fsi_grid == NULL) {
    //cout << "Mem reservation for fsi_grid" << endl;
    fsi_grid = new complex < double >**[getRgrid() + 1];
    for (int k = 0; k < (getRgrid() + 1); k++) {
      fsi_grid[k] = new complex < double >*[getCthgrid() + 1];
      for (int l = 0; l < (getCthgrid() + 1); l++) {
	fsi_grid[k][l] = new complex < double >[getPhigrid() + 1];
      }
    }
  }
  else {
    cerr << "fsi_grid should not be initialized more than once!!!!"<< endl;
    exit(1);
  }
}

//destructor
ROMEAGrid::~ROMEAGrid() {
  if (fsi_grid != NULL) {
    for (int k = 0; k < (getRgrid() + 1); k++) {
      for (int l = 0; l < (getCthgrid() + 1); l++) {
	delete[]fsi_grid[k][l];
      }
      delete[]fsi_grid[k];
    }
    delete[]fsi_grid;
  }
}

//interpolation of the grid after r,cth,phi have been set
complex < double >ROMEAGrid::getFsiGridFull_interp() {
  if (!filledgrid) {
    cerr << "You have to fill the grids first!" << endl;
    exit(1);
  }
  //sanity check to see that the right amount of knocked out nucleons end up in the final state!
  //if((getTotalProtonOut()==getProtonKnockout())&&(getTotalNeutronOut()==getNeutronKnockout())){
  complex < double >result = getInterp(fsi_grid);
  return result;
//   }
//   else{
//     cerr << "Mismatch between number of knocked out nucleons end number in final state!!!" << endl;
//     exit(1);
//   }
}

complex<double> ROMEAGrid::getFsiGridN_interp(const int grid){
   if (grid==0) return getFsiGridFull_interp();
  else {
    cerr << "gridindex out of range" << endl;
    assert(1==0);
  }
 
  
}
void ROMEAGrid::printFsi_grid() {
  cout << "Printing ROMEAarray for " << getFsi_Filename() << endl;
  printKnockout();
  for (int i = 0; i <= getRgrid(); i++) {
    double r = float (i) * getPnucleus()->getRange() / getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for (int j = 0; j <= getCthgrid(); j++) {
      double costheta = -2. * j / getCthgrid() + 1.;
      for (int k = 0; k <= getPhigrid(); k++) {
        double phi = (getAllinplane()? 1. : 2.) * PI * k / getPhigrid();
        if (std::isnan(phi))
          phi = 0.;
        setinterp(r,costheta,phi);
        complex < double >value = getFsiGridFull_interp();
        cout << r << " " << costheta << " " << phi * RADTODEGR << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }
}

void ROMEAGrid::print_grid(int gridindex){
  if (gridindex==0) printFsi_grid();
  else {
    cerr << "gridindex out of range" << endl;
    assert(1==0);
  }
}

//set filenames, don't forget to add ending "/" to dir if necessary!!!
void ROMEAGrid::setFilenames(string dir) {
  AbstractFsiGrid::setFilenames(dir + "ROMEA.");
  fsi_filename.insert(fsi_filename.size()-4,".fit"+to_string(fit));
}

//calc both fsi and fsi+ct grid
void ROMEAGrid::constructAllGrids() {
  //fill the opticalpotential array
  opticalpotential = new complex<double>*[getParticles().size()];
  double r[opticalsize]; double VV[opticalsize];double VS[opticalsize]; double WV[opticalsize]; double WS[opticalsize];
  for(int i=0;i<opticalsize;i++) r[i]=getPnucleus()->getRange()/opticalsize*i;
  for(size_t it=0;it<getParticles().size();it++){
    if(getParticles()[it].getParticletype()>1){
      cerr << "In ROMEA only ISI/FSI for nucleons!!!" <<endl;
      assert(1==0);
    }
    opticalpotential[it]=new complex<double>[opticalsize];
    double p=getParticles()[it].getP();
    double nucleons=getPnucleus()->getA()-getProtonKnockout()-getNeutronKnockout();
//     cout << nucleons << " " << fit << endl;
    global_(&fit,&p,&nucleons,r, VV, WV,VS,WS, &opticalsize);
    for(int i=0;i<opticalsize;i++) opticalpotential[it][i]=VS[i]+I_UNIT*WS[i]+getParticles()[it].getE()/getParticles()[it].getMass()*
      (VV[i]+I_UNIT*WV[i])+(-VV[i]*VV[i]+VS[i]*VS[i]+WV[i]*WV[i]-WS[i]*WS[i]+2.*I_UNIT*(VS[i]*WS[i]-VV[i]*WV[i]))
      /2./getParticles()[it].getMass();
  }
//    for(int i=0;i<opticalsize;i++) cout << r[i] << " " << opticalpotential[0][i] << " " << VV[i] << endl;
  //mem reservation for treshold array
  if (treshold == NULL) {
    //cout << "Mem reservation for treshold array" << endl;
    treshold = new int *[getCthgrid() + 1];
    for (int i = 0; i <= getCthgrid(); i++) {
      treshold[i] = new int[getPhigrid() + 1];
      for (int j = 0; j <= getPhigrid(); j++)
        treshold[i][j] = 0;
    }
  }
  else {
    cerr << "Treshold should only be initialized once!" << endl;
    exit(1);
  }
  
  //fill the arrays!
  for (int i = 0; i <= getRgrid(); i++) {
    r_hit = float (i) * getPnucleus()->getRange() / getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for (int j = 0; j <= getCthgrid(); j++) {
      costheta_hit = -2. * j / getCthgrid() + 1.;
      sintheta_hit = sqrt(1. - costheta_hit * costheta_hit);
      for (int k = 0; k <= getPhigrid(); k++) {
        phi_hit = (getAllinplane()? 1. : 2.) * PI * k / getPhigrid();
        if (std::isnan(phi_hit))
          phi_hit = 0.;
        sincos(phi_hit, &sinphi_hit, &cosphi_hit);
        if (treshold[j][k] == 0) {  //if treshold is true, we're close enough to one and we don't have to compute everything again and again
          for (size_t it = 0; it < getParticles().size(); it++) {
            getParticles()[it].setHitcoord(r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit);
          }
          calcROMEAphase(i, j, k);

          if (abs(fsi_grid[i][j][k]) > 0.999999) {
            treshold[j][k] = 1;
            if (j == 0 || j == getCthgrid()) {
              for (int ll = 1; ll <= getPhigrid(); ll++)
                treshold[j][ll] = 1;
            }
          }
          //r=0 symmetry shortcut
          if (i == 0) {
            for (j = 0; j <= getCthgrid(); j++) {
              for (k = 0; k <= getPhigrid(); k++) {
                fsi_grid[i][j][k] = fsi_grid[0][0][0];
                //cout << r_hit << " " << theta_hit*RADTODEGR << " " << phi_hit*RADTODEGR << " " << fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] << endl;
              }
            }
          }
          //theta=0 or Pi symmetry shortcut
          else if (j == 0 || j == getCthgrid()) {
            for (k = 1; k <= getPhigrid(); k++) {
              fsi_grid[i][j][k] = fsi_grid[i][j][0];
              //cout << r_hit << " " << theta_hit*RADTODEGR << " " << phi_hit*RADTODEGR << " " << fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] << endl;
            }
          }
        }
        //treshold reached, fill with 1s
        else {
	  fsi_grid[i][j][k] = 1.;
          //cout << r_hit << " " << theta_hit*RADTODEGR << " " << phi_hit*RADTODEGR << " " << fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] << endl;
        }
      }
    }
  }
  //mem cleanup
  for (int i = 0; i <= getCthgrid(); i++)
    delete[]treshold[i];
  delete[]treshold;
  treshold = NULL;
  for(size_t it=0;it<getParticles().size();it++) delete [] opticalpotential[it];
  delete[] opticalpotential;
}

//readin fsi grid
void ROMEAGrid::readinFsiGrid(ifstream & infile) {
  if (fsi_grid != NULL) {
    for (int k = 0; k < (getRgrid() + 1); k++) {
      for (int l = 0; l < (getCthgrid() + 1); l++) {
	for (int mm = 0; mm < (getPhigrid() + 1); mm++) {
	  infile.read(reinterpret_cast < char *>(&fsi_grid[k][l][mm]), sizeof(complex < double >));
	}
      }
    }
  }
  else {
    cerr << "fsi_grid not initialized!!!!" << endl;
    exit(1);
  }
}


//write fsi grid to file
void ROMEAGrid::writeoutFsiGrid(ofstream & outfile) {
  for (int k = 0; k < (getRgrid() + 1); k++) {
    for (int l = 0; l < (getCthgrid() + 1); l++) {
      for (int mm = 0; mm < (getPhigrid() + 1); mm++) {
	outfile.write(reinterpret_cast < char *>(&fsi_grid[k][l][mm]), sizeof(complex < double >));
      }
    }
  }
}




//calc glauberphases for one gridpoint
void ROMEAGrid::calcROMEAphase(const int i, const int j, const int k) {
  int res = 90;
  unsigned count = 0;
  double deeserror = 0.;
  if (integrator == 1 || integrator == 2) {
    fsi_grid[i][j][k] = 1.;
    for (size_t it = 0; it < getParticles().size(); it++) {	    
      numint::array < double, 1 > lower = { {getParticles()[it].getHitz()} };
      numint::array < double, 1 > upper = { {getPnucleus()->getRange()} };
      if(getParticles()[it].getIncoming()) {
	lower[0]=-getPnucleus()->getRange();
	upper[0]=getParticles()[it].getHitz();
      }
      ROMEAGrid::Ftor_one F;
      F.grid = this;
      F.it=it;
      numint::mdfunction < numint::vector_z, 1 > mdf;
      mdf.func = &Ftor_one::exec;
      mdf.param = &F;
      vector<complex<double> > ret(1, 0.);
      F.f = klaas_one_bound;
      if (integrator == 1)
	res = numint::cube_romb(mdf, lower, upper, 1.E-12, prec, ret, count, 0);
      else
	res = numint::cube_adaptive(mdf, lower, upper, 1.E-12, prec, 2E02,2E06, ret, count, 0);
      cout << i << " " << j << " " << k << " " << real(ret[0]) << " " << imag(ret[0]) << endl;
      fsi_grid[i][j][k] *= exp(-I_UNIT*getParticles()[it].getMass()/getParticles()[it].getP()*ret[0]/HBARC);
    }
    
  }
  else {cerr  << "integrator type not implemented" << endl; exit(1);}
//       cout << i << " " << j << " " << k << " " << level << " " << mm << fsi_grid[level][mm][i][j][k] << " " << fsi_ct_grid[level][mm][i][j][k] <<
//       " " << res << " " << count << endl;
}


void ROMEAGrid::klaas_one_bound(numint::vector_z & results, double z, ROMEAGrid & grid, size_t it) {
  results = numint::vector_z(1, 0.);
  double r=sqrt(z*z+grid.getParticles()[it].getHitbnorm()*grid.getParticles()[it].getHitbnorm());
  if (r > grid.getPnucleus()->getRange())
    return;
  results[0]= interpolate(grid.getOpticalpotential()[it],r,grid.getOpticalstep(),grid.getOpticalsize(),0);
}

