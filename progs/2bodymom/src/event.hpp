#ifndef __EVENT_HPP__
#define __EVENT_HPP__

#include <TVector3.h>
#include <cstdio>
#include <string>
using std::string;
#include <sstream>
using std::stringstream;

/** A struct that will hold all kinematics
 *  for a calculation
 */
struct Event {
	
	enum STATUS { SURVIVED = 1, OUTSIDECUTS = 0, UNPHYSICAL = -1 };

	Event();
	double xB;
	double Q2;
	double omega;
	double mass1,mass2; // masses
	int type1,type2; // types of knockout particles (0 proton, 1 neutron), is in principle deducible from shellindex[1,2]
	int shellindex1,shellindex2;
	TVector3 k1,q,p1,p2;
	int status; /**< 1 means it survived the cuts, 0 means it is cut away by detector cuts, -1 means it is unphysical (energy conservation cannot be met in any way **/
	void setStatus( STATUS s) { status = s; }
};

string vecString(TVector3& v);

/** print out all the kinematics of the event
 *  note that xB, Q2 and omega are only calculated
 *  if calc_xB(nucleus,event e) has been called.
 *  This is because the shellindices and the excitation
 *  energies have to be known in order to calculate
 *  xB, Q2 and omega.
 *  @param e the struct to print the members of
 */
void print_event(struct Event& e);

#endif // __EVENT_HPP__
