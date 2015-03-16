#include "event.hpp"
#include "parser.hpp"
#include "2bodymom.hpp"
#include <numint/numint.hpp>
#include "ClassGridThick_SCX.hpp"
#include <FastParticle.hpp>
#include <vector>
#include <complex>

void calculate_CX_probabilities(std::vector<struct Event>& events,MeanFieldNucleusThick& nuc,std::vector<double*>&,std::vector<double*>&);
void calculate_CX_probabilities(struct Event&,MeanFieldNucleusThick&,double*,double*);

struct F_SCX{
	static int exec( const int* ndim, const double xu[], const int* ncomp, double f[],void *param){
		struct F_SCX *p = (struct F_SCX* ) param;
		/** first of all map unit cube coordinates to r,cos theta,phi **/
		double x[3];
		x[0] = xu[0]*p->nuc->getRange(); // map [0,1] -> [0,R]
		x[1] = xu[1]*2.-1; // map[0,1] -> [-1,1]
		x[2] = xu[2]*2.*M_PI; // map[0,1] -> [0,2PI]
		/** construct single particle the wave functions ! these include a factor r! dim is [fm^{-1/2}] instead of [fm^{-3/2]! **/
		std::complex<double> wvf_1[4];
		std::complex<double> wvf_2[4];
		p->nuc->getWaveFunction(wvf_1,p->event->shellindex1,p->m1,x[0],x[1],x[2]);
		p->nuc->getWaveFunction(wvf_2,p->event->shellindex2,p->m2,x[0],x[1],x[2]);
		/** construct single particle densities **/
		double rho_1=0.;
		double rho_2=0.;
		for (unsigned i=0;i<2;i++){ // should run to 4, but we are neglecting lower components!
			rho_1 += std::norm(wvf_1[i]);
			rho_2 += std::norm(wvf_2[i]);
		}
		/** construct FSI factor **/
		double fsi = std::norm(p->fsi_grid->getFsiGridFull_interp3(x[0],x[1],x[2]) );
		/** construct different probability factors **/
		double p1 = p->cx_grid_1->getInterp(x);
		double p2 = p->cx_grid_2->getInterp(x);
		/** fill different components of f, with some cheap tricks to minimize operations **/
		/** [0] = no cx scattering, 
		 *  [1] is only first particle scattering, 
		 *  [2] is only second particle scattering, 
		 *  [3] is both scattering **/
		f[0] = rho_1*rho_2*fsi/(x[0]*x[0]);// divide by r^2. Each density has already factor r^2 included. "\int dr r^2 rho r^2 rho" is one r^2 too many
		f[1] = f[0]*p1;
		f[2] = f[0]*p2;
		f[3] = f[1]*p2;
		return 0; // always "succes"
	}
	Event *event;
	ClassGridThick_SCX *cx_grid_1,*cx_grid_2;
	GlauberGridThick *fsi_grid;
	MeanFieldNucleusThick* nuc;
	int m1,m2;
};

int main(int argc, char** argv){
	std::vector<struct Event> events;
	/** input is of the format
	 *  [./exec] [nucleus] [input file] [output file] [fsi]
	 *  or
	 *  [./exec] [nucleus] [input file] [output file] [fsi] [start] [stop]
	 */
	if (!(argc==4 || argc==6)) {
		fprintf(stderr,"Unsupported number of cmd line arguments.\n");
		fprintf(stderr,"Expected input of the form: \n");
		fprintf(stderr,"  $> [./exec] [nucleus] [input file] [output file] \n");
		fprintf(stderr,"or\n");
		fprintf(stderr,"  $> [./exec] [nucleus] [input file] [output file] [start] [stop] \n\n");
		exit(-1);
	}
	MeanFieldNucleusThick nuc(MeanFieldNucleus::TypeNames.at(std::string(argv[1])),SHAREDIR); // if this throws out of range you probably supplied incorrect nucleus name!
	if (argc == 4) // second argument is file containing events
		EventParser::read_events(argv[2],events,nuc);
	else // read in only a part of the kinematics file. specify start and stop line. Also id for unique output names./exec [kinfile] [outfile] [fsi] [start] [stop]	
		EventParser::read_events(argv[2],events,atoi(argv[4]),atoi(argv[5]),nuc);
	
	/** make vectors to store data and errors **/
	std::vector<double*> res,err;
	res.resize(events.size());
	err.resize(events.size());
	for (unsigned i=0;i<events.size();i++){
		res[i] = new double[3]; // why 3? (one particle scatters, other particle scatters, both scatter, makes 3
		err[i] = new double[3]; // idem
	}
	/** DO ALL THE WORK! **/
	calculate_CX_probabilities(events,nuc,res,err);
	
	/** the results now stored in res and err are not yet in the correct format (see notes) **/
	FILE* fout;
	fout = fopen(argv[3],"w");
	if (fout != NULL){ // file was succesfully opened
		fprintf(fout,"#[0]: no scattering [1]: leading cx [2]: recoil cx [3]: both cx");
		fprintf(fout,"[4]: no scattering error [5]: leading cx error [6]: recoil cx error [7]: both cx error\n");
		for (unsigned i=0;i<events.size();i++){
			double no_scx   = 1. - res[i][0] - res[i][1] + res[i][2];
			double scx_p1   = res[i][0] - res[i][2]; // particle 1 scatters, 2 does not
			double scx_p2   = res[i][1] - res[i][2]; // particle 2 scatters, 1 does not
			double scx_both = res[i][2]; // both scatter
			double no_scx_err   = err[i][0] + err[i][1] + err[i][2];
			double scx_p1_err   = err[i][0] + err[i][2];
			double scx_p2_err   = err[i][1] + err[i][2];
			double scx_both_err = err[i][2];
			fprintf(stdout,"%f   %f   %f   %f\n",no_scx,scx_p1,scx_p2,scx_both);
			fprintf(fout,"%f   %f   %f   %f      ",no_scx,scx_p1,scx_p2,scx_both);
			fprintf(fout,"%f   %f   %f   %f\n",no_scx_err,scx_p1_err,scx_p2_err,scx_both_err);
		}
	} else {
		std::cout << "[ERROR] Could not open file: " << argv[3] << std::endl;
		exit(-1);
	}
	fclose(fout);
	
	
	/** free the data and error vectors **/
	for (unsigned i=0;i<events.size();i++){
		delete[] res[i];
		delete[] err[i];
	}
	return 0;
}

void calculate_CX_probabilities(std::vector<struct Event>& events,MeanFieldNucleusThick& nuc,std::vector<double*>& res, std::vector<double*>& err){
	for (unsigned i=0;i<events.size();i++){
		std::cout << "[Info] processing event nr " << i+1 << " of " << events.size() << "  :: shellindices are " <<  events[i].shellindex1 << ", " << events[i].shellindex2 << std::endl;
		calculate_CX_probabilities(events[i],nuc,res[i],err[i]);
	}
}

void calculate_CX_probabilities(struct Event& event,MeanFieldNucleusThick& nuc,double* res, double* err){
	// make sure res and err are set to zero //
	res[0]=0.; res[1]=0.; res[2]=0.;
	err[0]=0.; err[1]=0.; err[2]=0.;
	/** make the fast particles required for the elastic FSI grid **/
	FastParticle fp1_fsi(event.type1,0,event.p1,0.,0.,SHAREDIR);
	FastParticle fp2_fsi(event.type2,0,event.p2,0.,0.,SHAREDIR);
	/** set up the Glauber FSI grid **/
	
	GlauberGridThick grid(60,20,5,&nuc,1e-03,2,SHAREDIR); // 1: 60=rgrid, 2: 20=thetagrid, 3: 10=phigrid, 4: nucleus, 5: 1e-6= prec, 5: [int]=integrator, 6=share dir
	grid.clearParticles();
	grid.addParticle(fp1_fsi);
	grid.addParticle(fp2_fsi);
	cout << "[GlauberGridThick::info] Calculating grids... "; cout.flush();
	grid.updateGrids();
	cout << "  [DONE] " << endl; cout.flush();
	grid.clearKnockout();
	
	/** make the fast particles required for the classical CX grids **/
	FastParticle fp1_cx((event.type1==0)? FastParticle::P_CLASS_SCX : FastParticle::N_CLASS_SCX,0,event.p1,0.,0.,SHAREDIR); // change elastic scattering nucleon 0 is proton 1 is neutron
	FastParticle fp2_cx((event.type2==0)? FastParticle::P_CLASS_SCX : FastParticle::N_CLASS_SCX,0,event.p2,0.,0.,SHAREDIR); // change elastic scattering nucleon 0 is proton 1 is neutron
	/** now make the classical CX grids, one for each particle **/
	ClassGridThick_SCX cx_grid_1(&nuc,&fp1_cx,25,25);
	ClassGridThick_SCX cx_grid_2(&nuc,&fp2_cx,25,25);
	cx_grid_1.addKnockoutParticle(event.shellindex1);
	cx_grid_1.addKnockoutParticle(event.shellindex2);
	cx_grid_1.constructGrid();
	cx_grid_2.addKnockoutParticle(event.shellindex1);
	cx_grid_2.addKnockoutParticle(event.shellindex2);
	cx_grid_2.constructGrid();
	/** set up the integration struct parameters **/
	struct F_SCX F;
	F.cx_grid_1 = &cx_grid_1;
	F.cx_grid_2 = &cx_grid_2;
	F.fsi_grid  = &grid;
	F.nuc       = &nuc;
	F.event     = &event;
	
	unsigned m_combos=0;
	for (int m1 = -nuc.getJ_array()[event.shellindex1]; m1 <= nuc.getJ_array()[event.shellindex1]; m1+=2){ // sum m1
		for (int m2 = (event.shellindex1 == event.shellindex2)? m1+2 : -nuc.getJ_array()[event.shellindex2] ; m2 <= nuc.getJ_array()[event.shellindex2]; m2+=2) { // Pauli principle -- never forget --
			if ((event.shellindex1==event.shellindex2) && (m1==m2)) m2+=2; // pauli!!!
			if (m2 > nuc.getJ_array()[event.shellindex2]) break; // we tried to change m2 because pauli but m2 is now too large
				/** set parameters to new m values **/
				F.m1=m1;
				F.m2=m2;
				/** add correct m knockouts to the FSI grid, not needed for classical grids **/
				grid.clearKnockout();
				grid.addKnockout(event.shellindex1,m1);
				grid.addKnockout(event.shellindex2,m2);
				
				/** general integration parameters **/
				int ndim = 3; // r,costh,phi
				int ncomp = 4; // see notes, 1, P_1, P_2 and P_1*P_2
				integrand_t integr = F_SCX::exec;
				void* userdata = (void*) &F; // cast the integrator struct to void pointer
				int nvec = 1; // unless you now what SIMD is leave this to 1
				double epsrel = 1e-6;
				double epsabs = 1e-6;
				int flags = 0x00;
				int mineval = 10000;
				int maxeval = 10000000;
				char* statefile = NULL;
				int neval,fail;
				double integral[ncomp];
				double error[ncomp];
				double prob[ncomp];
				double J = 4.*M_PI*nuc.getRange();
				/** VEGAS specific arguments **/
				int seed=1234;
				int nstart=10000;
				int nincrease=20000;
				int nbatch=10000;
				int gridno=0;
				
				Vegas(ndim,ncomp,integr,userdata,nvec,
					epsrel,epsabs,flags,seed,
					mineval,maxeval,nstart,nincrease,
					nbatch,gridno,statefile,&neval,
					&fail,integral,error,prob);
				/** vegas done still have to take care of jacobian!**/
				for (int i=0;i<ncomp;i++){
					integral[i]*=J;
					error[i]*=J;
				}
				//std::cout << "# integration result [ " << integral[0] << ", " << integral[1] << ", " << integral[2] << ", " << integral[3] << " ]" << std::endl;
				//std::cout << "# error       result [ " << error[0] << ", " << error[1] << ", " << error[2] << ", " << error[3] << " ]" << std::endl;
				/** now normalise the probabilities so we have proper probs
				 *  for the two shells, shellindex1,m1 and shellindex2, m2
				 *  equiv to n_1, \kappa_1, m_1 and n_2, \kappa_2, m_2
				 * */
				for (int i=1;i<4;i++){ // yes start @ 1, normalisation N is in integral[0]
					double rel_err = (error[0]/integral[0] + error[i]/integral[i]);
					integral[i] /= integral[0]; // normalise
					error[i] = integral[i]*rel_err; // propagate errors
				}
				m_combos++; // keep track of total m1,m2 combinations
				for (int i=0;i<3;i++){ // copy the 3 probs, P_alpha, P_beta, P_alpha*P_beta, see notes
					res[i] += integral[i+1];
					err[i] += error[i+1];
				}
		} // m2 sum
	} // m1 sum
	//std::cout << "# sum        result [ " << res[0] << ", " << res[1] << ", " << res[2] << " ]" << std::endl;
	//std::cout << "# total m_combos is " << m_combos << " dividing result by it " << std::endl;
	// average over different m combos
	for (int i=0;i<3;i++){
		res[i] /= m_combos;
		err[i] /= m_combos;
	}
	//std::cout << "# sum        result [ " << res[0] << ", " << res[1] << ", " << res[2] << " ]" << std::endl;
}
