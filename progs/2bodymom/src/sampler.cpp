#include "sampler.hpp"
#include "event.hpp"
#include <vector>
#include "MeanFieldNucleusThick.hpp"
#include "parser.hpp"
#include "2bodymom.hpp"

#ifndef SHAREDIR
#define SHAREDIR "/home/camille/Code/share"
#endif


void halla_kinematics(int, char**);
void hallb_kinematics(int, char**);
void prl99_072501_kinematics(int, char**);

int main(int argc, char* argv[]){
	//halla_kinematics(argc,argv);
	hallb_kinematics(argc,argv);
	//prl99_072501_kinematics(argc,argv);
	return 0;
}


void halla_kinematics(int argc, char* argv[]){
	/** cmd line argument parsing **/
	if (argc != 3){
		fprintf(stderr,"expected ./executable [seed] [number of events]\n");
		exit(-1);
	}
	int seed    = atoi(argv[1]); // avoid special case "1" which resets the generator, causing seed 0 and 1 to generate the same set of random numbers
	unsigned nevents = (unsigned) atoi(argv[2]);
	
	std::vector<struct Event> events;
	MeanFieldNucleusThick nuc(MeanFieldNucleus::He,SHAREDIR);
	HALLA::generateKinematics(events,nevents,nuc,seed);

	char fname[128];
	printf("HALL kin setting is %d\n",HALLA::KIN_SETTING);
	sprintf(fname,"HALLA_%s_%d_N%d.kin",nuc.getNucleusName().c_str(),HALLA::KIN_SETTING,nevents);
	EventParser::write_events(fname,events);
}

void prl99_072501_kinematics(int argc, char* argv[]){
	if (argc != 3){
		fprintf(stderr,"expeced [./exec] [seed] [number of events]\n");
		exit(-1);
	}
	int seed = atoi(argv[1]);
	unsigned nevents = (unsigned) atoi(argv[2]);

	std::vector<struct Event> events;
	MeanFieldNucleusThick nuc(MeanFieldNucleus::C,SHAREDIR);
	PRL99_072501::generateKinematics(events,nevents,nuc,seed);

	char fname[128];
	printf("PRL 072501 KIN setting is %d\n",PRL99_072501::KIN_SETTING);
	sprintf(fname,"PRL99_072501_%s_%d_N%d.kin",nuc.getNucleusName().c_str(),PRL99_072501::KIN_SETTING,nevents);
	EventParser::write_events(fname,events);
}

void hallb_kinematics(int argc, char* argv[]){
	
	if (argc != 3){
		fprintf(stderr,"expected ./executable [seed] [number of events]\n");
		exit(-1);
	}	
	int seed=atoi(argv[1]);
	unsigned nevents = (unsigned) atoi(argv[2]);
	
	std::vector<struct Event> events;
	MeanFieldNucleusThick nuc(MeanFieldNucleus::C,SHAREDIR);
	HALLB::generateKinematics(events,nevents,nuc,0,1,seed);

	char fname[128];
	sprintf(fname,"HALLB_%s_N%d.kin",nuc.getNucleusName().c_str(),nevents);
	EventParser::write_events(fname,events);
}
