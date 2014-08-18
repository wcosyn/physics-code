#include "2bodymom.hpp"
#include "event.hpp"
#include "parser.hpp"
#include "sampler.hpp"
#include <cassert>

int main(int argc, char** argv){
	std::vector<struct Event> events;
	MeanFieldNucleusThick nuc(MeanFieldNucleus::Al,SHAREDIR);
	printf("\n***********************************************************\n");
	printf("Assuming nucleus to be %s check this please!\n",nuc.getNucleusName().c_str());
	printf("************************************************************\n");
	/** input is of the format
	 *  [./exec] [input file] [output file] [fsi]
	 *  or
	 *  [./exec] [input file] [output file] [fsi] [start] [stop]
	 */
	if (argc == 4){ // second argument is file containing events
		EventParser::read_events(argv[1],events,nuc);
	} else if (argc == 6) { // read in only a part of the kinematics file. specify start and stop line. Also id for unique output names./exec [kinfile] [outfile] [fsi] [start] [stop]	
		EventParser::read_events(argv[1],events,atoi(argv[4]),atoi(argv[5]),nuc);
	} else {
		fprintf(stderr,"Unsupported number of cmd line arguments.\n");
		fprintf(stderr,"Expected input of the form: \n");
		fprintf(stderr,"  $> [./exec] [input file] [output file] [fsi] \n");
		fprintf(stderr,"or\n");
		fprintf(stderr,"  $> [./exec] [input file] [output file] [fsi] [start] [stop] \n\n");
		exit(-1);
	}
	int fsi = atoi(argv[3]);
	std::vector<double> data;
	std::vector<double> err;
	// calculation
	dist2bodymom(nuc,events,data,err,fsi);
	//----------------------------
	//
	//output
	EventParser::write_data(argv[2],data,err);
	//---------------------------
}

// integration structs
struct F {
	static void exec(const numint::array<double,3> &x, void *param, complex<double> &ret){
		F &p = * (F *) param; // cast void pointer to an struct F pointer
		p.f(ret,x,p.P,p.u_pm,p.u_ps,p.nuc,p.shellindex1,p.m1,p.shellindex2,p.m2);
	}
	MeanFieldNucleusThick* nuc;
	TSpinor *u_pm,*u_ps;
	int shellindex1,m1,shellindex2,m2;
	TVector3 P;

	void (*f)(complex<double>& , const numint::array<double,3>&, TVector3& P,TSpinor*, TSpinor*, MeanFieldNucleusThick*, int, int,int,int);
};

struct F_Glauber : F { 
	static void exec(const numint::array<double,3> &x, void *param, complex<double> &ret){
		F_Glauber &p = * (F_Glauber *) param; // cast void pointer to an struct F pointer
		p.f(ret,x,p.P,p.u_pm,p.u_ps,p.nuc,p.shellindex1,p.m1,p.shellindex2,p.m2,p.grid);
	}
	GlauberGridThick* grid;
	
	void (*f)(complex<double>& , const numint::array<double,3>&, TVector3& P,TSpinor*, TSpinor*, MeanFieldNucleusThick*,int, int,int,int,GlauberGridThick*);
};
/** expose dist2bodymom function so that i can use it in python or whatever **/
/*extern "C" double dist2bodymom(void* nuc,int stat,
					 double xB,
					 double Q2,
					 double omega,
					 double mass1,
					 double mass2,
					 int type1,
					 int type2,
					 int shellindex1,
					 int shellindex2,
					 double k1_x,double k1_y,double k1_z,
					 double q_x ,double q_y ,double q_z,
					 double p1_x,double p1_y,double p1_z,
					 double p2_x,double p2_y,double p2_z) {
	MeanFieldNucleusThick
*/

// MAIN FUNCTION YOU SHOULD ONLY CALL THIS ONE
void dist2bodymom(MeanFieldNucleusThick& nuc,std::vector<struct Event>& events, std::vector<double>& data,std::vector<double>& error,bool fsi){
	data.resize(events.size());
	error.resize(events.size());
	for (unsigned i=0; i<events.size();i++){
		assert( nuc.getTotalLevels() > events[i].shellindex1 ); // hopefully these will catch things when you are using the wrong nucleus
		assert( nuc.getTotalLevels() > events[i].shellindex2 ); // idem, only works if you selected the nucleus too small
		if (events[i].status == Event::SURVIVED) { // only calculate survived events
			(fsi) ? dist2bodymom_glauber(nuc,events[i],data[i],error[i]) : dist2bodymom_rpwia(nuc,events[i],data[i],error[i]);
		} else {
			data[i]  = 0.;
			error[i] = 0.;
		}
	}
}

void dist2bodymom_glauber(MeanFieldNucleusThick& nuc,Event& e,double& res,double& error){
	double intgrl_res = 0.; // result of the integral we will return at the end
	double intgrl_err = 0.;
	// GLAUBER stuff -------------------------
	//FastParticle proton_pm(0,0,pm,0.,0.,SHAREDIR); // 1: 0=proton, 2: 0=not a beam particle, 3: pm= 3vector of particle, 4: hard scale, 5: decay width, 6 share dir
	FastParticle particle_p1(e.type1,0,e.p1,0.,0.,SHAREDIR);
	FastParticle particle_p2(e.type2,0,e.p2,0.,0.,SHAREDIR);
	cout << "Sigmap of fast particle    " << particle_p1.getSigmap() << endl;
	cout << "Sigman of fast particle    " << particle_p1.getSigman() << endl;
	cout << "Sigmap of recoil particle "  << particle_p2.getSigmap() << endl;
	cout << "Sigman of recoil particle "  << particle_p2.getSigman() << endl;
	GlauberGridThick grid(60,20,5,&nuc,1e-03,2,SHAREDIR); // 1: 60=rgrid, 2: 20=thetagrid, 3: 10=phigrid, 4: nucleus, 5: 1e-6= prec, 5: [int]=integrator, 6=share dir
	grid.clearParticles();
	grid.addParticle(particle_p1);
	grid.addParticle(particle_p2);
	cout << "Calculating grids... "; cout.flush();
	grid.updateGrids();
	cout << "  [DONE] " << endl; cout.flush();
	grid.clearKnockout();
	// GLAUBER stuff END ---------------------------
	
	
	// init fourvectors to construct spinors
	// approx is neglecting the lower components
	// this is achieved by setting vector of fourvector to zero
	assert( e.mass1 > 0. && e.mass2 > 0.); // make sure masses are set correctly
	FourVector<double> f_k1 = FourVector<double>(sqrt(e.mass1*e.mass1 + e.k1.Mag2()),0.,0.,0.); // set lower component manually to zero
	FourVector<double> f_p2 = FourVector<double>(sqrt(e.mass2*e.mass2 + e.p2.Mag2()),0.,0.,0.); // set lower component manually to zero
	
	// prepare integration
	numint::array<double,3> lowerb = {{0.,1.,0.}}; // r, cos(theta), phi
	numint::array<double,3> upperb = {{nuc.getRange(),-1,2.*PI}}; // r, cos(theta), phi
	numint::mdfunction<complex<double>,3> mdf;

	// make struct and set struct parameters
	F_Glauber f;
	f.nuc = &nuc;
	f.f = integrandum_glauber;
	mdf.func = &F_Glauber::exec;
	mdf.param = &f;
	f.shellindex1 = e.shellindex1;
	f.shellindex2 = e.shellindex2;
	f.P = e.k1 + e.p2; // c.m. momentum
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
					cout << "shellind 1 : " << setw(2) << e.shellindex1 << " shellind2 " << setw(2) <<e.shellindex2 <<" m1 : " << setw(2) << m1 << " , m2 : " << setw(2) << m2 << " pol1 : " << setw(2) << pol[s1] << " pol2: " << setw(2) << pol[s2] << endl;
					grid.clearKnockout();
					grid.addKnockout(e.shellindex1,m1);
					grid.addKnockout(e.shellindex2,m2);
					f.m1 = m1;
					f.m2 = m2;
					int minEval =  20000;
					int maxEval = 100000;
					//unsigned count;
					//int succes = numint::cube_adaptive(mdf,lowerb,upperb,1e-06,prec,minEval,maxEval,ret,count,0);
					//complex<double> err(0.,0.);
					
					//--------- CUHRE---------------//
					int nregions,neval,fail;
					complex<double> ret,err,prob;
					double epsabs = 1e-8;
					double epsrel = 1e-4;
					cuhre(mdf,lowerb,upperb,1,epsrel,epsabs,0x00,minEval,maxEval,11,NULL,nregions,neval,fail,ret,err,prob);
					intgrl_res += norm(ret); // norm is abs val squared
					intgrl_err += 2.*fabs(real( ret*conj(err) )); // basic error propagation ignoring term of err^2
					//cout << succes << "  --  " << count << " -- " << ret << "    " << intgrl_res << endl;
				} // m2 loop
			} // m1 loop
		} // s2 loop
	} // s1 loop
	cout << "res " << intgrl_res << "  err  " << intgrl_err << endl;
	res   = (1./pow(2.*PI,3.))*intgrl_res;
	error = (1./pow(2.*PI,3.))*intgrl_err;
}

void dist2bodymom_rpwia(MeanFieldNucleusThick& nuc,Event& e,double& res,double& error){
	double intgrl_res = 0.; // result of the integral we will return at the end
	double intgrl_err = 0.;
	cout << "running rpwia calculations for event: " << endl;
	print_event(e);
	
	// init fourvectors to construct spinors
	assert( e.mass1 > 0. && e.mass2 > 0.);
	FourVector<double> f_k1 = FourVector<double>(sqrt(e.mass1*e.mass1+e.k1*e.k1), 0, 0, 0);
	FourVector<double> f_p2 = FourVector<double>(sqrt(e.mass2*e.mass2+e.p2*e.p2), 0, 0, 0);

	// prepare integration
	numint::array<double,3> lowerb = {{0.,1.,0.}}; // r, cos(theta), phi
	numint::array<double,3> upperb = {{nuc.getRange(),-1,2.*PI}}; // r, cos(theta), phi
	numint::mdfunction<complex<double>,3> mdf;

	// make struct and set struct parameters
	F f;
	f.nuc = &nuc;
	f.f = integrandum_rwpia;
	mdf.func = &F::exec;
	mdf.param = &f;
	f.shellindex1 = e.shellindex1;
	f.shellindex2 = e.shellindex2;
	f.P = e.k1 + e.p2 ; // c.m. momentum P = k1 + p2
	//cout << vecString(f.P) << "\t Mag: " << f.P.Mag() << "\t";
	TSpinor::Polarization::State pol[2] = {TSpinor::Polarization::kDown,TSpinor::Polarization::kUp}; // make a list to use my beloved basic forloops
	for (int s1=0; s1<2; s1++) // sum s_1
	{
		TSpinor u_k1 = TSpinor(f_k1,e.mass1,TSpinor::Polarization(0.,0.,pol[s1]),TSpinor::kUnity);
		f.u_pm = &u_k1;
		for (int s2=0; s2<2; s2++) // sum s_2
		{
			TSpinor u_p2 = TSpinor(f_p2,e.mass2,TSpinor::Polarization(0.,0.,pol[s2]),TSpinor::kUnity);
			f.u_ps = &u_p2;
			for (int m1 = -nuc.getJ_array()[e.shellindex1]; m1 <= nuc.getJ_array()[e.shellindex1]; m1+=2) // sum m1
			{
				// if shellindices are equal, m2 > m1 else we will have double counting!!!
				for (int m2 = (e.shellindex1 == e.shellindex2)? m1+2 : -nuc.getJ_array()[e.shellindex2] ; m2 <= nuc.getJ_array()[e.shellindex2]; m2+=2) // Pauli principle -- never forget --
				{
					if ((e.shellindex1==e.shellindex2) && (m1==m2)) m2+=2; // pauli!!!
					if (m2 > nuc.getJ_array()[e.shellindex2]) break; // we tried to change m2 because pauli but m2 is now 2 large
					cout << "shellind 1 : " << setw(2) << e.shellindex1 << " shellind2 " << setw(2) << e.shellindex2 <<" m1 : " << setw(2) << m1 << " , m2 : " << setw(2) << m2 << " pol1 : " << setw(2) << pol[s1] << " pol2: " << setw(2) << pol[s2] << endl;
					//cout << "nuc.getJarray()[" <<  shellindex1 << "] = " << nuc.getJ_array()[shellindex1] << endl;
					//cout << "nuc.getJarray()[" <<  shellindex2 << "] = " << nuc.getJ_array()[shellindex2] << endl;
					//cout << vecString(f.P) << endl;
					//Matrix<4,1> pmmat = u_pm;
					//cout << "Spinor pm: [" << pmmat(0,0) << ", " << pmmat(1,0) << ", " << pmmat(2,0) << ", " << pmmat(3,0) << "]" << endl;
					f.m1 = m1;
					f.m2 = m2;
					//double prec=1e-08;
					//complex<double> ret;
					int minEval =  20000;
					int maxEval = 100000;
					//unsigned count;
					//int succes = numint::cube_adaptive(mdf,lowerb,upperb,1e-06,prec,minEval,maxEval,ret,count,0);
					//complex<double> err(0.,0.);
					
					//--------- CUHRE---------------//
					int nregions,neval,fail;
					complex<double> ret,err,prob;
					double epsabs = 1e-8;
					double epsrel = 1e-4;
					cuhre(mdf,lowerb,upperb,1,epsrel,epsabs,0x00,minEval,maxEval,11,NULL,nregions,neval,fail,ret,err,prob);
					
					/*	
					int nnew = 10000;
					double flatness = 20.;
					complex<double> ret, err, prob;
					double epsabs = 1e-8;
					double epsrel = 1e-3;
					int nregions,neval,fail;
					suave(mdf,lowerb,upperb,1,epsrel,epsabs,0x00,1234,minEval,maxEval,nnew,flatness,NULL,nregions,neval,fail,ret,err,prob);
					*/
					//cout << " s1,s2,m1,m2 " << s1 << s2 << m1 << m2 << " result " << ret << " error " << err << endl;	
					//cout << " error on norm " << 2.*fabs(real( ret*conj(err) )) << endl;
					intgrl_res += norm(ret); // norm is abs val squared
					intgrl_err += 2.*fabs(real( ret*conj(err) )); // basic error propagation ignoring term of err^2
					//cout << "neval " << neval << endl;
					//assert(succes==0);
				//	cout << succes << "  --  " << count << " -- " << ret << "    " << intgrl_res << endl;
				} // m2 loop
			} // m1 loop
		} // s2 loop
	} // s1 loop
	cout << "res " << intgrl_res << "  err  " << intgrl_err << endl;
	res   = (1./pow(2.*PI,3.))*intgrl_res;
	error = (1./pow(2.*PI,3.))*intgrl_err;
	//return intgrl_res;
}

void integrandum_rwpia(complex<double>& res, const numint::array<double,3>& x, TVector3& P, TSpinor* u_pm, TSpinor* u_ps, MeanFieldNucleusThick* nuc, int shellindex1, int m1, int shellindex2,int m2)
{
	double sinth = sqrt(1.-x[1]*x[1]);
	double PR = P*TVector3(x[0]*sinth*cos(x[2]),x[0]*sinth*sin(x[2]),x[0]*x[1]); // dot product of vector P with vector R
	res = exp(-I_UNIT*INVHBARC*PR)*TSpinor::Bar(*u_pm)*TMFSpinor(*nuc,shellindex1,m1,x[0],x[1],x[2])*TSpinor::Bar(*u_ps)*TMFSpinor(*nuc,shellindex2,m2,x[0],x[1],x[2]); // Wavefunctions already defined with factor r, so no r^2 that would appear from jacobian, no sin(theta) because we integrate over d(cos(theta))
}

void integrandum_glauber(complex<double>& res, const numint::array<double,3>& x, TVector3& P, TSpinor* u_pm, TSpinor* u_ps, MeanFieldNucleusThick* nuc, int shellindex1, int m1, int shellindex2,int m2, GlauberGridThick* grid)
{
	complex<double> rwpia;
	integrandum_rwpia(rwpia,x,P,u_pm,u_ps,nuc,shellindex1,m1,shellindex2,m2);
	res = rwpia*grid->getFsiGridFull_interp3(x[0],x[1],x[2]);
}
