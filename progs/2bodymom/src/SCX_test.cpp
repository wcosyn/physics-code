#include "2bodymom.hpp"
#include "GlauberGridThick_SCX.hpp"
#include "MeanFieldNucleusThick.hpp"
#include "FastParticle.hpp"
#ifndef SHAREDIR
#define SHAREDIR "/home/ccolle/Code/share"
#endif
int main(){
	FastParticle fp(8,0,100,0.,0.,0.,0.,SHAREDIR); // type, beam part, momentum (MeV), ptheta,pphi,hard scale (CT),gamma (Decay width)
	MeanFieldNucleusThick nuc(MeanFieldNucleusThick::C,SHAREDIR);
	GlauberGridThick_SCX g(&nuc,fp,20,15);
	g.printGrid();
	//g.printDensityGrid();
	return 0;
}
