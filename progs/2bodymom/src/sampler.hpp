#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include <vector>
#include "event.hpp"
#include "MeanFieldNucleusThick.hpp"
#include <cstdlib>
#include <cstdio>

/** This is a header containing declarations of different sampling routines.
 *  
 *  Only declare here the function that you want to expose to other routines,
 *  ideally this is only generateKinematics(...)
 *
 *  For a new sampling routine declare here a new generateKinematics(...) function
 *  and implement it in yourfilenameofchoice.cpp
 *
 *  preferrably keep everything in a seperate namespace as it is more
 *  clear to name the different sampling function that sample the same
 *  kinematic variables the same name...
 *  
 */


/** generate samples according to kinematics in arXiv:1401.6138 **/
namespace HALLA {
	/** possible values for KIN_SETTING are 500,625,750 **/
	const int KIN_SETTING = 625; // ~global variable, only visible in HALLA namespace though, so hopefully you forgive me
	void generateKinematics(std::vector<struct Event>& events, unsigned n, MeanFieldNucleusThick& nuc,int seed);
}

/** generate samples according to kinematics in PRL99, 072501 **/
namespace PRL99_072501 {
	/** possible values for KIN_SETTING are 350 450 550 **/
	const int KIN_SETTING = 550;
	void generateKinematics(std::vector<struct Event>& events, unsigned n, MeanFieldNucleusThick& nuc,int seed);
}

/** generate samples according to slightly modified version of PLB:722, 63. modified according to mails and stuff from CLAS collaboration. Don't know if there's a paper yet except for our own massDependence in the making describing these kinematics...**/
namespace HALLB {
	void generateKinematics(std::vector<struct Event>& events, unsigned n, MeanFieldNucleusThick& nuc,int type1,int type2, int seed);
}


#endif // SAMPLER_HPP
