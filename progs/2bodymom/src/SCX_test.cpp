#include "2bodymom.hpp"
#include "GlauberGridThick_SCX.hpp"
#include "MeanFieldNucleusThick.hpp"
#include "FastParticle.hpp"
#ifndef SHAREDIR
#define SHAREDIR "/home/ccolle/Code/share"
#endif

void testgrid();
void scatterFrontScaling();

int main(){
	//scatterFrontScaling();
	testgrid();
	return 0;
}
void testgrid(){
	//FastParticle fp(8,0,100,0.,0.,0.,0.,SHAREDIR); // type, beam part, momentum (MeV), ptheta,pphi,hard scale (CT),gamma (Decay width)
	//FastParticle fp(8,0,100,M_PI/4.,0.,0.,0.,SHAREDIR);
	FastParticle fp(8,0,1397.,M_PI/2.,3.*M_PI/4.,0.,0.,SHAREDIR);
	MeanFieldNucleusThick nuc(MeanFieldNucleusThick::C,SHAREDIR);
	GlauberGridThick_SCX g(&nuc,fp,15,10);
	g.addKnockoutParticle(nuc.getPLevels()+1);
	g.addKnockoutParticle(nuc.getPLevels()-1);
	//g.constructGlauberGrid();
	//g.printGrid();

	for (double x=-nuc.getRange(); x<=nuc.getRange(); x+=1){
	for (double y=-nuc.getRange(); y<=nuc.getRange(); y+=1){
	for (double z=-nuc.getRange(); z<=nuc.getRange(); z+=1){
		double r     = sqrt(x*x+y*y+z*z);
		if (r < nuc.getRange()){
			double costh = (r>0.)? z/r : 0.;
			double phi   = atan2(y,x);
			double pos[] = {r,costh,phi};
			complex<double> d = g.getInterp(pos);
			//std::cout << "r costheta phi " << r << ", " << costh << ", " << phi << std::endl;
			std::cout << x << "\t" << y << "\t" << z << "\t" << d.real() << "\t" << d.imag() << std::endl;
		}
	}
	}
	}

	/*for (double r = 0.; r<nuc.getRange(); r+=0.5) {
		for (double th = 0.0; th<=M_PI+0.05; th+=0.05) {
			double x[] = {r,cos(th),0.};
			complex<double> interp = g.getInterp(x);
			std::cout << x[0] << "\t" << x[1] << "\t" << interp.real() << "\t" << interp.imag() << std::endl;
		}
	}*/
	//g.printDensityGrid();
}

void scatterFrontScaling(){
	for (double p=0; p<1e5; p+=100){
		FastParticle fp(8,0,p,0.,0.,0.,0.,SHAREDIR);
		std::cout << p << "\t" <<fp.getScatterfront(0).real() << "\t" << fp.getScatterfront(0).imag() << std::endl;
	}
}
