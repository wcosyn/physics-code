#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include "AbstractFsiGrid.hpp"
#include <constants.hpp>
#include <Utilfunctions.hpp>

AbstractFsiGrid::AbstractFsiGrid(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleus *pnucl, string homedir):
filledgrid(0),filledallgrid(0),dir(homedir+"/grids/"),rgrid(r_grid),cthgrid(cth_grid),phigrid(phi_grid),allinplane(0),totalprotonout(0), totalneutronout(0),
protonknockout(0), neutronknockout(0),pnucleus(pnucl),number_of_grids(1){
  
  invrstep=rgrid/pnucleus->getRange();
  invcthstep=cthgrid/2.;
  invphistep=0.5*phigrid/PI;
  
  
}  
  
AbstractFsiGrid::~AbstractFsiGrid(){
  //cout << "Deleting FSI object" << endl;
  
}

//get corrfilename
const string AbstractFsiGrid::getFsi_Filename() const{
  return fsi_filename;
}


//get dir
const string AbstractFsiGrid::getDir() const{
  return dir;
}

//get rgrid
int AbstractFsiGrid::getRgrid() const{
  return rgrid;
}

//get thgrid
int AbstractFsiGrid::getCthgrid() const{
  return cthgrid;
}

//get phigrid
int AbstractFsiGrid::getPhigrid() const{
  return phigrid;
}

//get invrstep
double AbstractFsiGrid::getInvRstep() const{
  return invrstep;
}

//get invthstep
double AbstractFsiGrid::getInvCthstep() const{
  return invcthstep;
}

//get invphistep
double AbstractFsiGrid::getInvPhistep() const{
  return invphistep;
}

//get s_interp
double AbstractFsiGrid::getS_interp() const{
  return s_interp;
}

//get t_interp
double AbstractFsiGrid::getT_interp() const{
  return t_interp;
}

//get u_interp
double AbstractFsiGrid::getU_interp() const{
  return u_interp;
}

//get comp_s_interp
double AbstractFsiGrid::getComp_s_interp() const{
  return comp_s_interp;
}

//get comp_t_interp
double AbstractFsiGrid::getComp_t_interp() const{
  return comp_t_interp;
}

//get comp_u_interp
double AbstractFsiGrid::getComp_u_interp() const{
  return comp_u_interp;
}

//get rindex
int AbstractFsiGrid::getRindex() const{
  return rindex;
}

//get thindex
int AbstractFsiGrid::getCthindex() const{
  return cthindex;
}

//get phiindex
int AbstractFsiGrid::getPhiindex() const{
  return phiindex;
}

//get pnucleus
MeanFieldNucleus *AbstractFsiGrid::getPnucleus() const{
  return pnucleus;
}

//made explicitly inline
//get particles vector
// vector<FastParticle*> AbstractFsiGrid::getParticles() const{
//  return particles; 
// }

//get knockoutlevels vector
vector<int> AbstractFsiGrid::getKnockoutLevels() const{
  return knockoutlevels;
}

//get knockoutm vector
vector<int> AbstractFsiGrid::getKnockoutM() const{
  return knockoutm;
}

//get allinpline
bool AbstractFsiGrid::getAllinplane() const{
  return allinplane;
}

  
//add particle to the particles vector
void AbstractFsiGrid::addParticle(FastParticle& newparticle){
  if(!filledgrid){
    if (particles.empty()){
      particles.push_back(newparticle);
      if(newparticle.getParticletype()==0){
	if(newparticle.getIncoming()) totalprotonout--;
	else totalprotonout++;
      }
      if(newparticle.getParticletype()==1){
	if(newparticle.getIncoming()) totalneutronout--;
	else totalneutronout++;
      }    
    }
    else{
      if(newparticle.getIncoming()){
	if(particles[0].getIncoming()){
	  cerr << "Not more than one incoming particle possible in the particle list!!!" << endl <<
	  "Replacing old incoming particle with this one" << endl;
	  if(particles[0].getParticletype()==0) totalprotonout++;
	  if(particles[0].getParticletype()==1) totalneutronout++;
	  particles[0]=newparticle;
	  if(particles[0].getParticletype()==0) totalprotonout--;
	  if(particles[0].getParticletype()==1) totalneutronout--;	
	}
	else particles.insert(particles.begin(),newparticle);
	if(particles[0].getParticletype()==0) totalprotonout--;
	if(particles[0].getParticletype()==1) totalneutronout--;	
      }
      else particles.push_back(newparticle);
      if(newparticle.getParticletype()==0) totalprotonout--;
      if(newparticle.getParticletype()==1) totalneutronout--;	    
    }
  }
  else{
    cerr << "You shouldn't add more particles after filling the fsi grids!!!" << endl;
    exit(1);
  }
}

//print the particles in the particles vector
void AbstractFsiGrid::printParticles() const{
  //cout << endl << endl << "Printing contents of particle vector in FSI grid" << endl;
  for(size_t i=0;i<particles.size();i++){
    cout <<"Particle " << i+1 << endl;
    particles[i].printParticle();
  }
}

//add a nucleon that gets knocked out
void AbstractFsiGrid::addKnockout(const int level, const int m){
  if(level<0||level>pnucleus->getTotalLevels()){
    cerr<<"invalid knockout level: " << level << endl;
    exit(1);
  }
  if(abs(m)>pnucleus->getJ_array()[level]){
    cerr<<"invalid m level: " <<pnucleus->getJ_array()[level] << " " << m << endl;
    exit(1);
  }
  for(size_t i=0;i<knockoutlevels.size();i++){
    if((level==knockoutlevels[i])&&(m==knockoutm[i])) {
      cerr<<"knockoutlevel is already present!!!" << endl << i+1 << " level: " << knockoutlevels[i] << ", m: " << knockoutm[i] << endl;
      exit(1);
    }
  }
  if(level<pnucleus->getPLevels()) protonknockout++;
  else neutronknockout++;
  knockoutlevels.push_back(level);
  knockoutm.push_back(m);
  return;
}

//empty knockout vectors
void AbstractFsiGrid::clearKnockout(){
 knockoutlevels.clear();
 knockoutm.clear();
 protonknockout=0;
 neutronknockout=0;
 return;
}

//print knockout vectors
void AbstractFsiGrid::printKnockout() const{
  cout << endl << endl << "Printing contents of knockout vector in FSI grid" << endl;
  for(size_t i=0;i<knockoutlevels.size();i++){
    cout <<"knockout " << i+1 << endl << "level: " << knockoutlevels[i] << ", m: " << knockoutm[i] << endl;
    
  }
}

int AbstractFsiGrid::getTotalProtonOut() const{
  return totalprotonout;
}

int AbstractFsiGrid::getTotalNeutronOut() const{
  return totalneutronout;
}

int AbstractFsiGrid::getProtonKnockout() const{
  return protonknockout;
}

int AbstractFsiGrid::getNeutronKnockout() const{
  return neutronknockout;
}



//set up inteprolation help variables for r
void AbstractFsiGrid::setRinterp(const double r){
  if(r>pnucleus->getRange()){
    cerr << "r out of range in : AbstractFsiGrid::setRinterp" << r << endl;
    exit(1);
  }
  rindex = int(floor(r*getInvRstep()));
  if(rindex==getRgrid()) rindex-=1;
  s_interp = (r*getInvRstep() - (rindex));
  comp_s_interp=1.-s_interp;
}

//set up interpolation help variables for theta
void AbstractFsiGrid::setCthinterp(const double costheta){
  cthindex = int(floor((-costheta+1.)*getInvCthstep()));
  if(cthindex==getCthgrid()) cthindex-=1;
  t_interp = ((-costheta+1.)*getInvCthstep() - (cthindex));  
  comp_t_interp=1.-t_interp;
}

//set up interpolation help variables for phi
void AbstractFsiGrid::setPhiinterp(const double phi){
  phiindex = int(floor(phi*getInvPhistep()));
  if(phiindex==getPhigrid()) phiindex-=1;
  u_interp = (phi*getInvPhistep() - (phiindex));  
  comp_u_interp=1.-u_interp;
}

//interpolation functions
complex<double> AbstractFsiGrid::getFsiGridFull_interpvec(const TVector3 &rvec){
  return getFsiGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiGrid::getFsiGridFull_interp3(const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiGridFull_interp();
}

complex<double> AbstractFsiGrid::getFsiGridFull_interp2(const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiGridFull_interp();
}

complex<double> AbstractFsiGrid::getFsiGridFull_interp1(const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiGridFull_interp();
}


//interpolation functions
complex<double> AbstractFsiGrid::getFsiGridN_interpvec(const int grid, const TVector3 &rvec){
  return getFsiGridN_interp3(grid, rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiGrid::getFsiGridN_interp3(const int grid, const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiGridN_interp(grid);
}

complex<double> AbstractFsiGrid::getFsiGridN_interp2(const int grid, const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiGridN_interp(grid);
}

complex<double> AbstractFsiGrid::getFsiGridN_interp1(const int grid, const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiGridN_interp(grid);
}


//3d interpolation of array
complex<double> AbstractFsiGrid::getInterp(complex<double> ***grid){
//   return getComp_u_interp()*(getComp_t_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()][getPhiindex()]
// 	+ getS_interp()*grid[getRindex()+1][getCthindex()][getPhiindex()]) 
// 	+ getT_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()+1][getPhiindex()] 
// 	+ getS_interp()*grid[getRindex()+1][getCthindex()+1][getPhiindex()]))
// 	+  getU_interp()*(getComp_t_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()][getPhiindex()+1]
// 	+ getS_interp()*grid[getRindex()+1][getCthindex()][getPhiindex()+1]) 
// 	+ getT_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()+1][getPhiindex()+1]
// 	+ getS_interp()*grid[getRindex()+1][getCthindex()+1][getPhiindex()+1]));
  return Interp3d(grid, getS_interp(), getT_interp(), getU_interp(), 
		  getComp_s_interp(), getComp_t_interp(), getComp_u_interp(), 
		  getRindex(), getCthindex(), getPhiindex());
}




//set filenames
void AbstractFsiGrid::setFilenames(string homedir){
  fsi_filename=homedir+pnucleus->getNucleusName();
  for(size_t i=0;i<particles.size();i++){
    if(particles[i].getIncoming()) fsi_filename+=".inc";
    fsi_filename+=".Ptype"+to_string(particles[i].getParticletype())+".p"+to_string(int(particles[i].getP()/10))
		  +".th"+to_string(int(round(particles[i].getCosTheta()*10)))+".ph"+to_string(int(round(particles[i].getPhi()*10)));
  }
  fsi_filename+=".r"+to_string(getRgrid())+".cth"+to_string(getCthgrid())+".phi"+to_string(getPhigrid())+".dat";
  //cout << fsi_filename << endl;
}
  
  
  
//fills the grids, uses pure virtual functions!!!!
void AbstractFsiGrid::fillGrids(){
  for(size_t i=0;i<getParticles().size();i++){
    allinplane=1;
    if((getParticles()[i].getPhi()!=0.)&&(getParticles()[i].getPhi()!=PI)){      
      allinplane=0;
      break;
    }
  }
  if(allinplane){
    //cout <<"All particles lie in the phi=0 or phi=Pi plane, special case activated" << endl;
    invphistep*=2.;
  }
  setFilenames(dir);
  ifstream infile(fsi_filename.c_str(),ios::in|ios::binary);
  //check if object has been created sometime earlier and read it in
  if(infile.is_open()){
    //cout << "Reading in FSI grid from memory: " << fsi_filename << endl;
    readinFsiGrid(infile);
    filledgrid=1;
    infile.close();
  }
  else{
    cout << "Constructing all grids" << endl;
    constructAllGrids();
    filledgrid=filledallgrid=1;
    ofstream outfile(fsi_filename.c_str(),ios::out|ios::binary);
    if(outfile.is_open()){
      //cout << "Writing out FSI grid: " << fsi_filename << endl;
      writeoutFsiGrid(outfile);
      outfile.close();
    }
    else{
      cerr << "could not open file for writing fsi grid output: " << fsi_filename << endl;
    }
  } 
}
  
