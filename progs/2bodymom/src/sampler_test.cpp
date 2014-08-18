#include "sampler.hpp"
#include "event.hpp"
#include <vector>
#include "MeanFieldNucleusThick.hpp"
#include "parser.hpp"
#include "2bodymom.hpp"

#ifndef SHAREDIR
#define SHAREDIR "/home/camille/Code/share"
#endif

int main(int argc, char* argv[]){
	/** cmd line argument parsing **/
	if (argc != 3){
		printf("expected ./executable [id (used for seed)] [number of events]\n");
		exit(-1);
	}
	int seed    = atoi(argv[1])+17; // avoid special case "1" which resets the generator, causing seed 0 and 1 to generate the same set of random numbers
	unsigned nevents = (unsigned) atoi(argv[2]);
	

	std::vector<struct Event> events;
	MeanFieldNucleusThick nuc(MeanFieldNucleus::He,SHAREDIR);
	HALLA::generateKinematics(events,nevents,nuc,seed);
	EventParser::write_events("parse_writetest.txt",events);

	std::vector<struct Event> readEvents;
	EventParser::read_events("parse_writetest.txt",readEvents,10,15);
	EventParser::write_events("parse_readwritetest.txt",readEvents);

	return 0;
}
