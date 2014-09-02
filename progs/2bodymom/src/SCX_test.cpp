#include "2bodymom.hpp"
#include "GlauberGridThick_SCX.hpp"
#include "MeanFieldNucleusThick.hpp"
#include "FastParticle.hpp"
#ifndef SHAREDIR
#define SHAREDIR "/home/ccolle/Code/share"
#endif
int main(){
	//FastParticle fp(8,0,100,0.,0.,0.,0.,SHAREDIR); // type, beam part, momentum (MeV), ptheta,pphi,hard scale (CT),gamma (Decay width)
	//FastParticle fp(8,0,100,M_PI/4.,0.,0.,0.,SHAREDIR);
	FastParticle fp(8,0,100,0.,0.,0.,0.,SHAREDIR);
	MeanFieldNucleusThick nuc(MeanFieldNucleusThick::C,SHAREDIR);
	GlauberGridThick_SCX g(&nuc,fp,20,15);
	g.addKnockoutParticle(nuc.getPLevels()+1);
	g.addKnockoutParticle(nuc.getPLevels()-1);
	g.constructGlauberGrid();
	//g.printGrid();
	for (double r = 0.; r<nuc.getRange(); r+=0.5) {
		for (double th = 0.0; th<=M_PI+0.05; th+=0.05) {
			double x[] = {r,cos(th),0.};
			complex<double> interp = g.getInterp(x);
			std::cout << x[0] << "\t" << x[1] << "\t" << interp.real() << "\t" << interp.imag() << std::endl;
		}
	}
	//g.printDensityGrid();
	return 0;
}
