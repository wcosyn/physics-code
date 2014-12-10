#include "2bodymom.hpp"
#include "GlauberGridThick_SCX.hpp"
#include "GlauberGridThick_SEL.hpp"
#include "MeanFieldNucleusThick.hpp"
#include "FastParticle.hpp"
#include "event.hpp"
#ifndef SHAREDIR
#define SHAREDIR "/home/ccolle/Code/share"
#endif

void testgrid();
void testgrid_el();
void scatterFrontScaling();
void createGrids();
void testIntegral(); // test integration routine thing for phase space unbiased stuff

/** a struct for single charge exchange nucleon nucleon rescattering **/
struct F_SCXm {
	static void exec(const numint::array<double,3> &x, void* param, complex<double> &ret){
		F_SCXm &p = * (F_SCXm *) param;
		p.f(ret,x,p.u_pm,p.u_ps,p.nuc,p.shellindex1,p.m1,p.shellindex2,p.m2,p.grid);
	}
	MeanFieldNucleusThick* nuc;
	TSpinor *u_pm,*u_ps;
	int shellindex1,m1,shellindex2,m2;
	GlauberGridThick_SCX* grid;

	void (*f)(complex<double>&, const numint::array<double,3>&,TSpinor*,TSpinor*,MeanFieldNucleusThick*,int,int,int,int,GlauberGridThick_SCX* grid);
};


void integrandum(complex<double>& res, const numint::array<double,3>& x, TSpinor* u_pm, TSpinor* u_ps, MeanFieldNucleusThick* nuc, int shellindex1, int m1, int shellindex2,int m2,GlauberGridThick_SCX* grid)
{
	res = std::norm(TSpinor::Bar(*u_pm)*TMFSpinor(*nuc,shellindex1,m1,x[0],x[1],x[2])*TSpinor::Bar(*u_ps)*TMFSpinor(*nuc,shellindex2,m2,x[0],x[1],x[2])*grid->getInterp(x)); // Wavefunctions already defined with factor r, so no r^2 that would appear from jacobian, no sin(theta) because we integrate over d(cos(theta))
	//TMFSpinor phi1(*nuc,shellindex1,m1,x[0],x[1],x[2]);
	//TMFSpinor phi2(*nuc,shellindex2,m2,x[0],x[1],x[2]);
	//res = phi1.H()*phi2;//*grid->getInterp(x)); // Wavefunctions already defined with factor r, so no r^2 that would appear from jacobian, no sin(theta) because we integrate over d(cos(theta))
}

void mdist2bodymom_SCX(MeanFieldNucleusThick& nuc,Event& e,double& res,double& error){
	double intgrl_res = 0.; // result of the integral we will return at the end
	double intgrl_err = 0.;
	// GLAUBER stuff -------------------------
	//FastParticle proton_pm(0,0,pm,0.,0.,SHAREDIR); // 1: 0=proton, 2: 0=not a beam particle, 3: pm= 3vector of particle, 4: hard scale, 5: decay width, 6 share dir
	
	FastParticle particle_p1(e.type1==0? 8 : 9,0,e.p1,0.,0.,SHAREDIR); // change elastic scattering nucleon (0,1) to SCX nucleon (8,9)
	FastParticle particle_p2(e.type2==0? 8 : 9,0,e.p2,0.,0.,SHAREDIR); // change elastic scattering nucleon (0,1) to SCX nucleon (8,9)
	cout << "# Type and mass of fast particle " << e.mass1 << " " << e.type1 << endl;
	cout << "# Type and mass of slow particle " << e.mass2 << " " << e.type2 << endl;
	bool SCX_LEADING = true;
	GlauberGridThick_SCX grid(&nuc,SCX_LEADING ? particle_p1 : particle_p2 ,25,20); // the particle you pass here is the one SCX FSIs are calculated upon
	grid.setArbitraryPhase(1.5*M_PI);
	grid.addKnockoutParticle(e.shellindex1); // add knockouts so density for SCX is adjusted 
	grid.addKnockoutParticle(e.shellindex2); // add knockouts so density for SCX is adjusted
	grid.constructGlauberGrid();             // call this after you set the correct knockout particles and the arbitrary phase! 
	// GLAUBER stuff END ---------------------------
	
	
	// init fourvectors to construct spinors
	// approx is neglecting the lower components
	// this is achieved by setting vector of fourvector to zero
	assert( e.mass1 > 0. && e.mass2 > 0.); // make sure masses are set correctly
	FourVector<double> f_k1 = FourVector<double>(sqrt(e.mass1*e.mass1 + e.k1.Mag2()),0.,0.,0.); // set lower component manually to zero
	FourVector<double> f_p2 = FourVector<double>(sqrt(e.mass2*e.mass2 + e.p2.Mag2()),0.,0.,0.); // set lower component manually to zero
	
	// prepare integration
	numint::array<double,3> lowerb = {{0.,-1.,0.}}; // r, cos(theta), phi
	numint::array<double,3> upperb = {{nuc.getRange(),1,2.*PI}}; // r, cos(theta), phi
	numint::mdfunction<complex<double>,3> mdf;

	// make struct and set struct parameters
	F_SCXm f;
	f.nuc = &nuc;
	f.f = integrandum;
	mdf.func = &F_SCXm::exec;
	mdf.param = &f;
	f.shellindex1 = e.shellindex1;
	f.shellindex2 = e.shellindex2;
	f.grid = &grid;
	TSpinor::Polarization::State pol[2] = {TSpinor::Polarization::kDown,TSpinor::Polarization::kUp}; // make a list to use my beloved basic forloops
	for (int s1=0; s1<2; s1++) // sum m_s1
	{
		TSpinor u_k1 = TSpinor(f_k1,e.mass1,TSpinor::Polarization(0.,0.,pol[s1]),TSpinor::kUnity);
		f.u_pm = &u_k1;
		for (int s2=0; s2<2; s2++) // sum m_s2
		{
			TSpinor u_p2 = TSpinor(f_p2,e.mass2,TSpinor::Polarization(0.,0.,pol[s2]),TSpinor::kUnity);
			f.u_ps = &u_p2;
			for (int m1 = -nuc.getJ_array()[e.shellindex1]; m1 <= nuc.getJ_array()[e.shellindex1]; m1+=2) // sum m1
			{
				// if shellindices are equal, m2 > m1 else we will have double counting!!!
				for (int m2 = (e.shellindex1 == e.shellindex2)? m1+2 : -nuc.getJ_array()[e.shellindex2] ; m2 <= nuc.getJ_array()[e.shellindex2]; m2+=2) // Pauli principle -- never forget --
				{
					if ((e.shellindex1==e.shellindex2) && (m1==m2)) m2+=2; // pauli!!!
					if (m2 > nuc.getJ_array()[e.shellindex2]) break; // we tried to change m2 because pauli but m2 is now too large
					//cout << "shellind 1 : " << setw(2) << e.shellindex1 << " shellind2 " << setw(2) <<e.shellindex2 <<" m1 : " << setw(2) << m1 << " , m2 : " << setw(2) << m2 << " pol1 : " << setw(2) << pol[s1] << " pol2: " << setw(2) << pol[s2] << endl;
					f.m1 = m1;
					f.m2 = m2;
					int minEval =   10000;
					int maxEval =  500000;
					/*
					unsigned count;
					double prec=1e-06;
					complex<double> ret;
					int succes = numint::cube_adaptive(mdf,lowerb,upperb,1e-06,prec,minEval,maxEval,ret,count,0);
					*/
					//complex<double> err(0.,0.);
					
					//--------- CUHRE---------------//
					/*	
					int nregions,neval,fail;
					complex<double> ret,err,prob;
					double epsabs = 1e-8;
					double epsrel = 1e-4;
					cuhre(mdf,lowerb,upperb,1,epsrel,epsabs,0x00,minEval,maxEval,11,NULL,nregions,neval,fail,ret,err,prob);
					*/
					
					//-------- Suave ----------------//
					int nregions,neval,fail;
					complex<double> ret,err,prob;
					double epsabs = 1e-8;
					double epsrel = 1e-4;
					int seed = 113337;
					int nnew = 5000;
					double flatness = 0.2;
					suave(mdf,lowerb,upperb,1,epsrel,epsabs,0x00,seed,minEval,maxEval,nnew,flatness,NULL,nregions,neval,fail,ret,err,prob);
					
					std::cout << "shell1 m1 = " << e.shellindex1 << " " << m1 << " shell2 m2 = " << e.shellindex2 << " " << m2 << std::endl;	
					std::cout << " ret is " << ret << "  real part is " << real(ret) << "  complex part is " << imag(ret) << std::endl;
					intgrl_res += abs(ret); // norm is abs val squared
					intgrl_err += 2.*fabs(real( ret*conj(err) )); // basic error propagation ignoring term of err^2
					//cout << succes << "  --  " << count << " -- " << ret << "    " << intgrl_res << endl;
				} // m2 loop
			} // m1 loop
		} // s2 loop
	} // s1 loop
	//cout << "res " << intgrl_res << "  err  " << intgrl_err << endl;
	res   = intgrl_res;
	error = intgrl_err;
}

void testIntegrand(){
	Event e;
	e.type1 = 0;
	e.type2 = 0;
	e.mass1 = 938.5;
	e.mass2 = 938.5;
	e.p1    = TVector3(0.,0.,1000.);
	e.p2    = TVector3(0.,0.,1000.);
	
	MeanFieldNucleusThick nuc(MeanFieldNucleus::Fe,SHAREDIR);
	double result = 0.;
	double error  = 0.;
	
	for (int shellindex1 = 0; shellindex1 < nuc.getPLevels(); shellindex1++){
		for (int shellindex2 = shellindex1; shellindex2 < nuc.getPLevels(); shellindex2++){
			double res,err;
			e.shellindex1 = shellindex1;
			e.shellindex2 = shellindex2;
			mdist2bodymom_SCX(nuc,e,res,err);
			result += res;
			error  += err;
		}
	}
	std::cout << result << "    " << error << std::endl;
}

int main(){
	//scatterFrontScaling();
	testgrid_el();
	//testgrid();
	//createGrids();
	//testIntegrand();
	return 0;
}

void createGrids(){
	FastParticle fpp(8,0,1000.,0.,0.,0.,0.,SHAREDIR);
	FastParticle fpn(9,0,1000.,0.,0.,0.,0.,SHAREDIR);
	MeanFieldNucleusThick C(MeanFieldNucleusThick::C,SHAREDIR);
	MeanFieldNucleusThick Al(MeanFieldNucleusThick::Al,SHAREDIR);
	MeanFieldNucleusThick Fe(MeanFieldNucleusThick::Fe,SHAREDIR);
	MeanFieldNucleusThick Pb(MeanFieldNucleusThick::Pb,SHAREDIR);
	GlauberGridThick_SCX gCp(&C,fpp,25,20);
	GlauberGridThick_SCX gCn(&C,fpn,25,20);
	GlauberGridThick_SCX gAlp(&Al,fpp,25,20);
	GlauberGridThick_SCX gAln(&Al,fpn,25,20);
	GlauberGridThick_SCX gFep(&Fe,fpp,25,20);
	GlauberGridThick_SCX gFen(&Fe,fpn,25,20);
	GlauberGridThick_SCX gPbp(&Pb,fpp,25,20);
	GlauberGridThick_SCX gPbn(&Pb,fpn,25,20);
}

void testgrid_el(){
	FastParticle fp(1,0,250.,0.,0.,0.,0.,SHAREDIR); // type, beam part, momentum (MeV), ptheta,pphi,hard scale (CT),gamma (Decay width)
	//FastParticle fp(8,0,100,M_PI/4.,0.,0.,0.,SHAREDIR);
	//FastParticle fp(8,0,1397.,M_PI/2.,3.*M_PI/4.,0.,0.,SHAREDIR);
	MeanFieldNucleusThick nuc(MeanFieldNucleusThick::Fe,SHAREDIR);
	GlauberGridThick_SEL g(&nuc,fp,25,20);
	g.addKnockoutParticle(nuc.getPLevels()+1);
	g.addKnockoutParticle(nuc.getPLevels()-1);
	g.constructGlauberGrid();
	g.printGrid();
	/*
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
	}*/

	/*for (double r = 0.; r<nuc.getRange(); r+=0.5) {
		for (double th = 0.0; th<=M_PI+0.05; th+=0.05) {
			double x[] = {r,cos(th),0.};
			complex<double> interp = g.getInterp(x);
			std::cout << x[0] << "\t" << x[1] << "\t" << interp.real() << "\t" << interp.imag() << std::endl;
		}
	}*/
	//g.printDensityGrid();
}
void testgrid(){
	FastParticle fp(8,0,250.,0.,0.,0.,0.,SHAREDIR); // type, beam part, momentum (MeV), ptheta,pphi,hard scale (CT),gamma (Decay width)
	//FastParticle fp(8,0,100,M_PI/4.,0.,0.,0.,SHAREDIR);
	//FastParticle fp(8,0,1397.,M_PI/2.,3.*M_PI/4.,0.,0.,SHAREDIR);
	MeanFieldNucleusThick nuc(MeanFieldNucleusThick::Fe,SHAREDIR);
	GlauberGridThick_SCX g(&nuc,fp,25,20);
	g.addKnockoutParticle(nuc.getPLevels()+1);
	g.addKnockoutParticle(nuc.getPLevels()-1);
	g.setArbitraryPhase(0.5*M_PI);
	g.constructGlauberGrid();
	g.printGrid();
	/*
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
	}*/

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
	for (double p=350; p<2000; p+=10){ // p is in MeV
		FastParticle fp(8,0,p,0.,0.,0.,0.,SHAREDIR);
		FastParticle fpel(0,0,p,0.,0.,0.,0.,SHAREDIR);
		std::cout << p << "\t" <<fp.getScatterfront(0).real() << "\t" << fp.getScatterfront(0).imag() << "\t";
		std::cout << fpel.getScatterfront(0).real() << "\t" << fpel.getScatterfront(0).imag() << "\t";
		std::cout << fpel.getScatterfront(1).real() << "\t" << fpel.getScatterfront(1).imag() << std::endl;
	}
	std::ofstream beta("beta_dep.dat");
	for (double p=350; p<2000; p+=10){
		FastParticle fp(8,0,p,0.,0.,0.,0.,SHAREDIR);
		FastParticle fpel(0,0,p,0.,0.,0.,0.,SHAREDIR);
		beta << p << "\t" << fp.getBetasq(0) << "\t" << fp.getBetasq(1) << "\t";
		beta << fpel.getBetasq(0) << "\t" << fpel.getBetasq(1) << "\t" << std::endl;
	}
	beta.close();
}
