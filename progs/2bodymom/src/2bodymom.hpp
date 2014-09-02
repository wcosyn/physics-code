// the file that combines all the different thingies.
// for center of mass prob dens. for two nucleon knockout
// we have (1) RWPIA + GLAUBER for WIM's kinematics (symmetric kinematics)
// we have (2) RWPIA + GLAUBER for MAARTENs kinematics (note_kinematics) (missing momentum along z-axis)

// Some includes

#ifndef TWOBODYMOM_HPP
#define TWOBODYMOM_HPP

#include <TMFSpinor.hpp>
#include <TSpinor.h>
#include <FourVector.h>
#include <constants.hpp>
#include <MeanFieldNucleus.hpp>
#include <MeanFieldNucleusThick.hpp>
#include <GlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include "GlauberGridThick_SCX.hpp"
#include <FsiCorrelator.hpp>
#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ofstream;
#include <numint/numint.hpp>
#include <complex>
using std::complex;
#include <iomanip>
using std::setw;
#include <vector>
using std::vector;
using std::pair;
#include <fstream>
using std::ofstream;
#include <cmath>
#include <sstream>
using std::stringstream;
#include <cassert>
#include "event.hpp"
#include <cstdlib>

// some preprocessor directives

//#define SHAREDIR "/home/camille/Code/share" // for local
//#define SHAREDIR "/home/ccolle/Code/share" // for trillian
//#define SHAREDIR "/user/home/gent/vsc407/vsc40761/scratch/Code/share" // for gengar
#ifndef SHAREDIR
	#define SHAREDIR getShareDir()
#endif // SHAREDIR

inline std::string getShareDir(){
	std::stringstream ss;
	ss << getenv("HOME") << "/Code/share" ;
	return ss.str();
}

// structs
struct F;
struct F_Glauber;

// work function
void  dist2bodymom        (MeanFieldNucleusThick& nuc,std::vector<struct Event>& events,std::vector<double>& data,std::vector<double>& error,int fsi);
void  dist2bodymom_glauber(MeanFieldNucleusThick& nuc,struct Event& e,double& res, double& err);
void  dist2bodymom_rpwia  (MeanFieldNucleusThick& nuc,struct Event& e,double& res, double& err);
void  dist2bodymom_SCX    (MeanFieldNucleusThick& nuc,struct Event& e,double& res, double& err);

void integrandum_rwpia(complex<double>& res, const numint::array<double,3>& x, TVector3& P, TSpinor* u_pm, TSpinor* u_ps, MeanFieldNucleusThick* nuc,int shellindex1, int m1, int shellindex2,int m2);
void integrandum_glauber(complex<double>& res, const numint::array<double,3>& x, TVector3& P, TSpinor* u_pm, TSpinor* u_ps, MeanFieldNucleusThick* nuc, int shellindex1, int m1, int shellindex2,int m2, GlauberGridThick* grid);
void integrandum_SCX(complex<double>& res, const numint::array<double,3>& x, TVector3& P, TSpinor* u_pm, TSpinor* u_ps, MeanFieldNucleusThick* nuc, int shellindex1, int m1, int shellindex2,int m2, GlauberGridThick_SCX* grid);

// some util print functions
string vecString(TVector3& v);
void printArray(MeanFieldNucleusThick&, const int*);
void printNucleusInfo(MeanFieldNucleusThick&);


#endif // TWOBODYMOM_HPP
