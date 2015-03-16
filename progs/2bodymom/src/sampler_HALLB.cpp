#include "sampler.hpp"
#include <fstream>
#include <string>
#include <cstdlib>
#include <cassert>
#include <random>
#include <iostream>

/** axis convention z // k1
 *  q in x-z plane
 */
/** note that as in sampler_HALLA you request n shellIndependent events. This means that you will receive at LEAST n events and at MOST n*shellcombos events! **/

/** forward declarations here instead of header sampler.hpp to prevent pollution of sampler.hpp **/
namespace HALLB {
	bool cut(Event& e); // check if a event survives our cuts
	bool get_q(Event& e, const MeanFieldNucleus& nuc, const TVector3& Pcm,const double qmag); // get angle of q
	//bool get_k1(Event& e,const MeanFieldNucleus& nuc, const TVector3& Pcm); // get k1, returns false if no physical solution is possible, else true
	void get_Pcm(TVector3& Pcm); // get a possible c.m. momentum vector

	// nested class
	// a xB-Q2 sampler to get xB-Q2 from exp distribution
	class xBQ2_sampler {
		public:
			xBQ2_sampler(const char* inputfile,unsigned,unsigned,double,double,double,double);
			~xBQ2_sampler();
			void sample(double& xB,double &Q2);
			unsigned tries() { return _tries;}
			unsigned accept(){ return _accept;}
		private:
			void init2d(); // init 2d array
			void del2d (); // delete 2d array
			double max2d (); // find maximum of 2d array
			void norm  (); // normalize 2d array

			double** _H;
			double* _nval,* _mval;
			double _nstep,_mstep;
			unsigned _n,_m;

			unsigned _tries,_accept;
	};
}

/** n is the number of events per shell. So the actual number of events that will be generated will be 
 *  number of shell combo's * n
 */
void HALLB::generateKinematics(std::vector<struct Event>& events, unsigned n, MeanFieldNucleusThick& nuc, int type1, int type2, int seed) {
	srand(seed);
	std::default_random_engine generator;
	generator.seed(seed);
	std::uniform_int_distribution<int> shell1_distribution;
	std::uniform_int_distribution<int> shell2_distribution;
	std::cout << "[Info] Particle types are " << type1 << " and " << type2 << std::endl;
	if (type1==0) // proton knock out
		shell1_distribution = std::uniform_int_distribution<int>(0,nuc.getPLevels()-1); // shell index generator, bounds inclusive
	else  // neutron knock out
		shell1_distribution = std::uniform_int_distribution<int>(nuc.getPLevels(),nuc.getTotalLevels()-1);
	if (type2==0) // proton knock out
		shell2_distribution = std::uniform_int_distribution<int>(0,nuc.getPLevels()-1); // shell index generator, bounds inclusive
	else  // neutron knock out
		shell2_distribution = std::uniform_int_distribution<int>(nuc.getPLevels(),nuc.getTotalLevels()-1);
	
		
	events.clear();
	std::string path = "../HALLB/data/";
	std::string file;
	printf("Nucleus name is %s\n",nuc.getNucleusName().c_str());
	if (nuc.getNucleusName().compare("He")==0)
		file = "C_xB_Q2.txt"; // use carbon as we do not have He xB Q2 distribution
	else if (nuc.getNucleusName().compare("C")==0)
		file = "C_xB_Q2.txt";
	else if (nuc.getNucleusName().compare("Al")==0)
		file = "Al_xB_Q2.txt";
	else if (nuc.getNucleusName().compare("Fe")==0)
		file = "Fe_xB_Q2.txt";
	else if(nuc.getNucleusName().compare("Pb")==0)
		file = "Pb_xB_Q2.txt";
	else {
		fprintf(stderr,"Warning: I could not check the validity of the xB Q2 distribution input for this nucleus\n");
		exit(-1);
	}
	std::string inputfile = path+file;
	
	xBQ2_sampler s(inputfile.c_str(),100,400,1.,2.,1.e6,6.e6); // 100 rows, 400 columns, xB rows in range 1 to 2, Q2 columns in range 1 to 6 [GeV]^2 of 1e6 to 6e6 [MeV]^2


	struct sampler_logger {
		sampler_logger() : tries(0),cuts(0),invalid(0),accept(0) {}
		unsigned tries;
		unsigned cuts;
		unsigned invalid;	
		unsigned accept;
	};
	struct sampler_logger logger;

	while (logger.accept < n){
		/** shell independent quantities **/
		/** xB,Q2,omega **/
		double xB,Q2;
		s.sample(xB,Q2);
		double omega = Q2/(2.*MASSP*xB);
		/** |q| **/
		double qmag = sqrt( omega*omega + Q2 );
		/** Pcm **/
		TVector3 Pcm = TVector3();	
		get_Pcm(Pcm);
		TVector3 k1  = TVector3(0.,0.,300. + 300.*rand()/RAND_MAX );
		/** shellindices **/
		int shellindex1 = shell1_distribution(generator); // get a random shellindex1
		int shellindex2 = shell2_distribution(generator); // get a random shellindex2
		if (type1==type2){
			while (shellindex1 > shellindex2){ // if you didn't do this the weights of different shellcombos would be twice as large as same shellcombos because 1,2 == 2,1 and hence twice prob to prob equiv pair if different shells
				shellindex1 = shell1_distribution(generator);
				shellindex2 = shell2_distribution(generator);
			}
		}	
		logger.tries++;
		Event e;
		e.xB          = xB;
		e.Q2          = Q2;
		e.omega       = omega;
		
		e.k1          = k1;
		e.p2          = Pcm-k1;
		
		e.shellindex1 = shellindex1;
		e.shellindex2 = shellindex2;
		e.mass1       = (type1==0)? MASSP : MASSN;
		e.mass2       = (type2==0)? MASSP : MASSN;
		e.type1       = type1; // 0 is for proton
		e.type2       = type2; // 0 is for proton
			
		/** q **/
		bool valid = get_q(e,nuc,Pcm,qmag);
			
		if (valid) { // valid, q has physical value
			e.p1 = e.k1 + e.q; // set the final piece of the kinematics
			/** check energy conservation just to be sure **/
			double M_Amin2 = nuc.getMassA_min_pp()+nuc.getExcitation()[shellindex1]+nuc.getExcitation()[shellindex2];
			double initE   = omega + nuc.getMassA();
			double finalE  = sqrt( M_Amin2*M_Amin2 + Pcm.Mag2()) + sqrt( e.mass1*e.mass1 + e.p1.Mag2()) + sqrt( e.mass2*e.mass2 + e.p2.Mag2());
			assert(fabs(initE - finalE) < 1e-9); // allow diff of 1 meV as units are MeV
			if (cut(e)) { // valid && survived cuts
				logger.accept++;
				e.status = Event::SURVIVED; // 1 means survival, 0 means cut away
			} else { // valid but did not survive cuts
				e.status = Event::OUTSIDECUTS;
				logger.cuts++;
			}
		} else { // invalid
			e.status = Event::UNPHYSICAL;
			logger.invalid++;
		}
		events.push_back(e); // save the event, even though it is most probably not going to make it trough our cuts
		if (events.size() % 1000 == 0)
			printf("Sampler has sampled %8d number of events of which %8d are valid and survived cuts.\n",events.size(),logger.accept);
	} 

	char fname[128];
	sprintf(fname,"HALLB_%s_N%d.sampler",nuc.getNucleusName().c_str(),n);
	FILE *fp = fopen(fname,"w");
	fprintf(fp,"============================================\n");
	fprintf(fp,"               SAMPLER REPORT               \n");
	fprintf(fp,"============================================\n\n");
	fprintf(fp,"** xB-Q2 rejection sampling ** \n");
	fprintf(fp,"------------------------------ \n");
	fprintf(fp," Total number of tries        : %d\n",s.tries());
	fprintf(fp," Total number of accepts      : %d\n",s.accept());
	fprintf(fp," Total xB-Q2 sampler rejected : %d\n",s.tries()-s.accept());
	fprintf(fp," xB-Q2 sampling efficiency    : %.2f %%\n", 100.*s.accept()/s.tries());
	fprintf(fp,"\n\n");
	fprintf(fp,"** q - Pcm sampling **\n");
	fprintf(fp,"-----------------------\n");
	fprintf(fp,"Total number of tries                       : % d\n",logger.tries);
	fprintf(fp,"Total number of invalids (|cos(theta)| > 1.): % d\n",logger.invalid);
	fprintf(fp,"Total number of physical valid events       : % d\n",logger.accept+logger.cuts);
	fprintf(fp,"Total number of cuts                        : % d\n",logger.cuts);
	fprintf(fp,"Total number accepted                       : % d\n",logger.accept);
	fprintf(fp,"\n\n");
	fclose(fp);
}
			
				
bool HALLB::cut(Event& e){
	return     ( e.xB > 1.2)
		&& ( e.k1.Mag() > 300.)
		&& ( e.p2.Mag() > 350.)
		&& ( e.p1.Angle(e.q) < 25*DEGRTORAD)
		&& ( e.p1.Mag() / e.q.Mag() >= 0.60);
		//&& ( e.p1.Mag() / e.q.Mag() < 0.96);
}

/** Warning it is criticial that the following fields of Event are set!
 * Event::shellindex1
 * Event::shellindex2
 * Event::mass1
 * Event::mass2
 * Event::k1
 * Event::p2
 *
 * it is assumed that particle "1" is the leading proton
 * it is assumed that particle "2" is the recoiling proton
 */
bool HALLB::get_q(Event& e, const MeanFieldNucleus& nuc, const TVector3& Pcm, const double qmag){
	double M_Amin2 = nuc.getMassA_min_pp() + nuc.getExcitation()[e.shellindex1] + nuc.getExcitation()[e.shellindex2];
	double E_Amin2 = sqrt( M_Amin2*M_Amin2 + Pcm.Mag2() );
	double E_2     = sqrt( e.mass2*e.mass2 + e.p2.Mag2());
	double A       = nuc.getMassA() + e.omega - E_Amin2 - E_2;

	double cosThetaq = ( A*A - e.mass1*e.mass1 - e.k1.Mag2() - qmag*qmag )/(2.*e.k1.Mag()*qmag);
	
	if ( fabs(cosThetaq) <= 1.0 ){ // valid!
		e.q = qmag*TVector3(sqrt(1.-cosThetaq*cosThetaq),0.,cosThetaq);
		return true;
	} else { // invalid
		e.q = TVector3(0.,0.,0.);
		return false;
	}
}
/** Warning it is critical that the following fields of Event are set!
 *  Event::shellindex1
 *  Event::shellindex2
 *  Event::mass1
 *  Event::mass2
 *  Event::omega
 *  Event::q
 *
bool HALLB::get_k1(Event& e,const MeanFieldNucleus& nuc, const TVector3& Pcm){
	double M_Amin2 = nuc.getMassA_min_pp() + nuc.getExcitation()[e.shellindex1] + nuc.getExcitation()[e.shellindex2];
	double E_Amin2 = sqrt( M_Amin2*M_Amin2 + Pcm.Mag2() );
	double deltaM2 = e.mass1*e.mass1 - e.mass2*e.mass2;
	
	double A = nuc.getMassA() + e.omega - E_Amin2;
	double B = e.q.Z()+Pcm.Z();
	double C = A*A + deltaM2 + e.q.Mag2() - Pcm.Mag2();
	
	double a = 4.*(B*B-A*A);
	double b = 4.*(B*C-2.*A*A*e.q.Z());
	double c = C*C - 4.*A*A*(e.mass1*e.mass1 + e.q.Mag2());

	double D = b*b - 4.*a*c;
	if ( D < 0. ) { // no physical solutions
		return false;
	} else {
		double k1_zm = (-b - sqrt(D))/2./a;
		double k1_zp = (-b + sqrt(D))/2./a;
		double k1_mag = 0.;
		if (k1_zm > 0.){ // if k1_zm is pos so will k1_zp be, accept lowest value though as mom distribution falls of exponentially
			k1_mag = k1_zm;
		} else if (k1_zp > 0.){
			k1_mag = k1_zp;
		} else if ( k1_zm < 0 && k1_zp < 0){
			printf("Both roots negative!!! %f and %f\n",k1_zm,k1_zp);
			return false;
		} else {
			printf("ERROR YOU SHOULD NEVER END UP HERE\n");
			exit(-1);
		}
		e.k1 = TVector3(0.,0.,k1_mag);
	}
	return true;
}
*/
void HALLB::get_Pcm(TVector3& Pcm){
	Pcm.SetX(-600.0 + 1200.*(double) rand()/RAND_MAX);
	Pcm.SetY(-600.0 + 1200.*(double) rand()/RAND_MAX);
	Pcm.SetZ(-600.0 + 1200.*(double) rand()/RAND_MAX);
}

/***********************************************
 *                                             *
 *  NESTED CLASS FOR XB Q2 REJECTION SAMPLING  *
 *                                             *
 ***********************************************/

HALLB::xBQ2_sampler::xBQ2_sampler(const char* inputfile,unsigned n, unsigned m,double nstart,double nstop, double mstart, double mstop) : _n(n), _m(m), _tries(0), _accept(0) {
       init2d();
	_nval = new double[n]; // contains LEFT bin edges!
 	_mval = new double[m]; // contains LEFT bin edges!
	_nstep = (nstop-nstart)/n;
	_mstep = (mstop-mstart)/m;
	
	for (unsigned i=0; i<n;i++){
		_nval[i] = nstart + (double)i*_nstep;
	}
	for (unsigned i=0; i<m;i++){
		_mval[i] = mstart + (double)i*_mstep;
	}
	std::ifstream data(inputfile);
	if (data.is_open()){
		for (unsigned x=0; x<n; x++){
			for (unsigned y=0; y<m; y++) {
				data >> _H[x][y];
			}
		}
		norm(); // normalise the freshly read in _H[][]
	} else {
		printf("Error could not open file %s\n",inputfile);
	}
}

void HALLB::xBQ2_sampler::init2d(){
	_H = new double*[_n];
	for (unsigned i=0; i<_n;i++)
		_H[i] = new double[_m];
}

void HALLB::xBQ2_sampler::del2d(){
	for (unsigned i=0; i<_n; i++)
		delete[] _H[i];
	delete[] _H;
}

HALLB::xBQ2_sampler::~xBQ2_sampler(){
	del2d();
	delete[] _nval;
	delete[] _mval;
}

double HALLB::xBQ2_sampler::max2d(){
	double max = -1e99;
	for (unsigned i=0;i<_n;i++){
		for (unsigned j=0; j<_m; j++) {
			if (_H[i][j] > max)
				max = _H[i][j];
		}
	}
	return max;
}

void HALLB::xBQ2_sampler::norm(){
	double max = max2d();
	for (unsigned i=0;i<_n;i++){
		for (unsigned j=0; j<_m; j++) {
			_H[i][j]/=max;
		}
	}
}

// request a xB-Q2 sample
void HALLB::xBQ2_sampler::sample(double& xB,double &Q2) {
	
	int i = rand() % _n;
	int j = rand() % _m;
	double xx = (double) rand()/RAND_MAX;

	while ( _H[i][j] < xx ){ // reject sample
		xx = (double) rand()/RAND_MAX;
		i  = rand() % _n;
		j  = rand() % _m;
		_tries++;
	}
	_accept++;	
	xB = _nval[i] + _nstep*rand()/RAND_MAX;
	Q2 = _mval[j] + _mstep*rand()/RAND_MAX;
}

