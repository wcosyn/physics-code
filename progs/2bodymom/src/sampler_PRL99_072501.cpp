#include "2bodymom.hpp"
#include "event.hpp"
#include <cstdlib>
#include <vector>
#include <cstdio>
#include <cassert>
#include "sampler.hpp"

/** pass a vector of events. Each shell is calculated for given shell (excitation energy)
 *  independent quantities. Only the shells that yield physical possible solutions are kept.
 *  Note that n is the number of events that will be calculated for each possible shell combinations
 *  Hence the total length of events will generally exceed the n you supplied. The maximum length of the
 *  vector is obviously n*shellcombos
 *
 *  axis convention is z-axis along incoming electron beam!
 *  
 *        x  ^
 *           |
 *           |     e_i
 *           |  --->
 *           |
 *    -------o----------> z
 *           |y
 *           |
 *           
 *  
 *  **/
namespace PRL99_072501 {
	// some forward declarations
	TVector3 get_e_p2(); // get random direction for p2
	bool validate_p2(double); // determine if |p2| is within detector acceptance
	TVector3 get_ef(); // get a random vector for e'
	TVector3 get_p1();  // get a random vector for p1
	
	/** n is number of events without loop over shellcombos, total
	 *  number of events generated will be at most n*shellcombos
	 */
	void generateKinematics(std::vector<struct Event>& events, unsigned n,MeanFieldNucleusThick& nuc,int seed) {
		
		srand(seed);

		events.clear();
		double Eei  = 4627.; // MeV
		TVector3 ei = TVector3(0.,0.,Eei); // @ ultra relativistic E = |p|c
		printf("Seed is : %d\n",seed);

		struct sampler_logger{ // a logger to check sampling efficiency
			sampler_logger() : tries(0), p2_choices(0), p2_invalid(0), p2_nosol(0), seed(0) {}
			unsigned tries; // total number of events generated
			unsigned p2_choices; // total number of times i had two valid solutions for p2_mag
			unsigned p2_invalid; // total number of times p2 was invalidated by detector cuts
			unsigned p2_nosol;   // total number of times D<0 and hence no physical solution for p2
			int seed;            // seed of random generator
		};
		struct sampler_logger logger;
		logger.seed = seed; // set the seed for logger

		unsigned int shellIndep_counter = 0; // don't count events from same shells
		while (shellIndep_counter < n) {
			//   ------------------------   //
			// shell independent quantities //
			TVector3 ef     = get_ef();
			double Eef      = ef.Mag(); // @ ultra relativistic E = |p|c
			double omega    = Eei - Eef;
			TVector3 q      = ei - ef;
			double Q2       = q.Mag2() - omega*omega; // Q^2 = q^2 - w^2
			double xB       = Q2/(2.*MASSP*omega);	
			TVector3 e_p2   = get_e_p2(); // random direction for e_p2
			TVector3 p1     = get_p1();
			double E1       = sqrt( MASSP*MASSP + p1.Mag2() );
			TVector3 k1     = p1-q;
			double ek1_ep2  = k1.Unit().Dot(e_p2); // (e_{k1} \cdot e_{p2})
			//   -----------------------   //
			//
			unsigned oldNumberOfEvents = events.size(); // copy number of events before loop starts
			for (int shellindex1 = 0; shellindex1 < nuc.getPLevels(); shellindex1++){
				for (int shellindex2 = shellindex1; shellindex2 < nuc.getPLevels(); shellindex2++) {
					logger.tries++;
					struct Event ev;
					ev.xB          = xB;
					ev.Q2          = Q2;
					ev.omega       = omega;
					ev.mass1       = MASSP;
					ev.mass2       = MASSP;
					ev.k1          = k1;
					ev.p1          = p1;
					ev.q           = q;
					ev.shellindex1 = shellindex1;
					ev.shellindex2 = shellindex2;

					double Eexc    = nuc.getExcitation()[shellindex1] + nuc.getExcitation()[shellindex2];
					double M_Amin2 = nuc.getMassA_min_pp() + Eexc;
						
					double B = (omega + nuc.getMassA() - E1);
					double A = 0.5*(M_Amin2*M_Amin2 + ev.mass2*ev.mass2 + k1.Mag2() - B*B);

					double a = (ev.mass2*ev.mass2 + M_Amin2*M_Amin2 + k1.Mag2()*( 1. - ek1_ep2*ek1_ep2) - 2.*A);
					double b = 2.*k1.Mag()*ek1_ep2*( ev.mass2*ev.mass2 - A);
					double c = M_Amin2*M_Amin2*ev.mass2*ev.mass2 + k1.Mag2()*ev.mass2*ev.mass2 - A*A;
					
					double D = b*b - 4.*a*c;
					
					if ( D >= 0){ // found two solutions
						double p2_1 = ( -b - sqrt(D) )/(2.*a);
						double p2_2 = ( -b + sqrt(D) )/(2.*a);
						double p2_mag = 0.;
						bool accept = true;	
						if ( validate_p2(p2_1) && !validate_p2(p2_2)) { // only first solution is within acceptance
							p2_mag = p2_1;
						} else if (!validate_p2(p2_1) && validate_p2(p2_2) ) { // only second solution is within acceptance
							p2_mag = p2_2;
						} else if (validate_p2(p2_1) && validate_p2(p2_2) ) { // both solutions are valid, make choice
							p2_mag = std::min(p2_1,p2_2); // choose the minimum one, ideally from exp decay distr... 
							logger.p2_choices++;
						} else { // no solution within detector acceptance
							accept = false;
							logger.p2_invalid++;
						}
						if (accept){
							ev.p2  = p2_mag*e_p2;
							TVector3 Pcm = ev.k1 + ev.p2;
							double initE  = ev.omega + nuc.getMassA();
							double finalE = sqrt( M_Amin2*M_Amin2 + Pcm.Mag2()) + sqrt(ev.mass1*ev.mass1 + ev.p1.Mag2()) + sqrt(ev.mass2*ev.mass2 + ev.p2.Mag2());
							assert(fabs(initE-finalE) < 1e-9); // allow energy discrepancy of 1 meV as units are MeV
							events.push_back(ev);   // store the event
							//printf("Difference in init and final energy: %e \n",fabs(initE-finalE));
						}
					} else { logger.p2_nosol++;}
				} // shellindex2
			} // shellindex1
			if (events.size() > oldNumberOfEvents){ // at least one shellcombo made it trough and was appended to events.
				shellIndep_counter++;
			}
		} // while ( shellIndep_counter < n)

		// write out logger info
		char fname[128];
		sprintf(fname,"PRL99_072501_%s_%d_N%d.sampler",nuc.getNucleusName().c_str(),KIN_SETTING,n);
		FILE* fp = fopen(fname,"w");
		fprintf(fp,"\n===========================================================\n");
		fprintf(fp,"                 PRL99, 072501 SAMPLER REPORT                \n");
		fprintf(fp,"===========================================================\n\n");
		fprintf(fp,"Seed : %d\n",logger.seed);
		fprintf(fp,"Kinematic setting : %d\n",KIN_SETTING);
		fprintf(fp,"Total number of shell combinations                    : %d \n",nuc.getPLevels()*(nuc.getPLevels()+1)/2);
		fprintf(fp,"Total number of events requested/accepted             : %d / %d \n",n,events.size());
		fprintf(fp,"Total number of events accepted/generated             : %d / %d \n",events.size(),logger.tries);
		fprintf(fp,"Total number of events trashed                        : %d \n",logger.tries-events.size());
		fprintf(fp,"Total number of events that had 2 valid sols for |p2| : %d \n",logger.p2_choices);
		fprintf(fp,"Total number of events invalidated by |p2| cut        : %d \n",logger.p2_invalid);
		fprintf(fp,"Total number of events w.o. phys. solution (D<0)      : %d \n",logger.p2_nosol);
		fprintf(fp,"\n\n=========================================================\n\n");
		fclose(fp);
	}

	/** is magnitude of recoiling proton within detector acceptance? **/
	bool validate_p2(double p2){
		return ( p2 > 250. && p2 < 900.);
	}

	TVector3 get_e_p2(){
		double Omega = 0.096; // in sr
		double costhMin = 1. - Omega/(2.*PI); // maximum deviation angle of top cone of ps, minimum cosine of that angle
		
		double u = (double) rand()/RAND_MAX; // u in range [0;1]
		//double v = 0.5*(1.+cosThMax) + (double) rand()/RAND_MAX* 0.5*(1.-cosThMax); // v in range [0.5*(1+cos(thmax));1]
		double costh = costhMin + (double) rand()/RAND_MAX * (1. - costhMin);
		
		double th    = acos(costh); // acos returns always in [0,pi], this is what we want
		double phi   = 2.*PI*u; 
		
		TVector3 p2   = TVector3();
		p2.SetMagThetaPhi(1.,th,phi); // unit vector, mag = 1, angles
		if (KIN_SETTING == 350 || KIN_SETTING == 450 || KIN_SETTING == 550 )
			p2.RotateY(-99.*DEGRTORAD); //
		else
			fprintf(stderr,"KIN SETTING %d NOT SUPPORTED!\n",KIN_SETTING);
		assert( std::fabs(1.-p2.Mag()) < 1e-10); // check if e_p2 is indeed unit vector, can be removed if found to never fail
		return p2;
	}

	TVector3 get_ef(){
		double hacc = 0.03; // in rad
		double vacc = 0.06; // in rad
		TVector3 ef = TVector3(0.,0.,1.); // create // at z-axis with unit magnitude
		
		double u = (double) rand()/RAND_MAX;
		double v = (double) rand()/RAND_MAX;
		
		ef.RotateY((1.-2.*u)*hacc); // (1-2*u) is random number between -1 and 1, so rotation is between -h and h rad
		ef.RotateX((1.-2.*v)*vacc);
		//ef *= 3762. + (2.*w - 1)*0.04*3762. ; // ef varied by +- 4%
		ef *= 3762.; // fixed ef
		ef.RotateY(DEGRTORAD*19.5); // rotate to 20.3 degrees, now pointing in right mean direction
		return ef;
	}

	TVector3 get_p1(){
		double hacc = 0.03; // in rad
		double vacc = 0.06; // in rad
		TVector3 p = TVector3(0.,0.,1.); // create at z-axis with unit magnitude
		
		double u = (double) rand()/RAND_MAX;
		double v = (double) rand()/RAND_MAX;
		double w = (double) rand()/RAND_MAX;
		
		// this is not exactly uniform but for small accuracy angles good approx
		p.RotateY((1.-2.*u)*hacc); // (1-2*u) is random number between -1 and 1, so rotation is between -h and h rad
		p.RotateX((1.-2.*v)*vacc);
		
		if (KIN_SETTING==350) {     
			p *= 1450 + (2.*w-1.)*0.04*1450; // 4 % spread on p
			p.RotateY(DEGRTORAD*(-40.1)); //  now pointing in right mean direction
		} else if (KIN_SETTING==450) {  
			p *= 1420 + (2.*w-1.)*0.04*1420; // 4 % spread on p
			p.RotateY(DEGRTORAD*(-35.8)); //  now pointing in right mean direction
		} else if (KIN_SETTING==550) {  // 750 MEV SETTINGS
			p*= 1360  + (2.*w-1.)*0.04*1360; // 4 % spread on p
			p.RotateY(DEGRTORAD*(-32.)); //  now pointing in the right mean direction
		} else {
			fprintf(stderr,"KIN SETTING %d NOT SUPPORTED\n",KIN_SETTING);
		}
		return p;
	}
} // NAMESPACE HALLA

