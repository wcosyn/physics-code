#include "event.hpp"


Event::Event() : 
	xB(0.),
	Q2(0.),
	omega(0.),
	mass1(0.),
	mass2(0.),
	type1(-1),
	type2(-1),
	shellindex1(-1),
	shellindex2(-1),
	k1 (TVector3(0.,0.,0.)),
	q  (TVector3(0.,0.,0.)),
	p1 (TVector3(0.,0.,0.)),
	p2 (TVector3(0.,0.,0.)),
	status(666) // should never stay 666 //
{}

       

string vecString(TVector3& v){
	std::stringstream ss;
	ss << "[" << v.X() << ", " << v.Y() << ", " << v.Z() << "]";
	return ss.str();
}

void print_event(struct Event& e) {
	printf("k1          = %s \n",vecString(e.k1).c_str());
	printf("p2          = %s \n",vecString(e.p2).c_str());
	printf("p1          = %s \n",vecString(e.p1).c_str());
	printf("q           = %s \n",vecString(e.q).c_str());
	printf("xB          = %f \n",e.xB);
	printf("Q2          = %f \n",e.Q2);
	printf("omega       = %f \n",e.omega);
	printf("shellindex1 = %d \n",e.shellindex1);
	printf("shellindex2 = %d \n",e.shellindex2);
	printf("type1       = %d \n",e.type1);
	printf("type2       = %d \n",e.type2);
	printf("mass1       = %f \n",e.mass1);
	printf("mass2       = %f \n",e.mass2);
	printf("status      = %d \n",e.status);
}

