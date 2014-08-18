#include "parser.hpp"
#include <constants.hpp>

void EventParser::write_events(char* fname, std::vector< struct Event>& events) {
	FILE* fp;
	fp = fopen(fname,"w");
	
	// indicate which columns contain what and how many columns they span, so that my python script will now which columns are what
	fprintf(fp,"#:: stat[1] cut[1] xB[1] Q2[1] omega[1] q[3] k1[3] p1[3] p2[3] shellindex1[1] shellindex2[1] ::#\n"); 	
	// print column labels
	fprintf(fp,"#%4s    ","stat");
	fprintf(fp,"%8s    %8s    %8s     ","xB","Q2","omega");
	fprintf(fp,"%8s    %8s    %8s    ","q_x","q_y","q_z");
	fprintf(fp,"%8s    %8s    %8s    ","k1_x","k1_y","k1_z");
	fprintf(fp,"%8s    %8s    %8s    ","p1_x","p1_y","p1_z");
	fprintf(fp,"%8s    %8s    %8s    ","p2_x","p2_y","p2_z");
	fprintf(fp,"%8s    %8s","shell1","shell2");
	fprintf(fp,"\n");

	for(unsigned i=0; i<events.size(); i++){ 
		fprintf(fp,"%4d    ",events[i].status);
		fprintf(fp,"%8.3f    %8.1f    %8.2f    ",events[i].xB,events[i].Q2,events[i].omega);
		fprintf(fp,"%8.2f    %8.2f    %8.2f    ",events[i].q.X() ,events[i].q.Y() ,events[i].q.Z());
		fprintf(fp,"%8.2f    %8.2f    %8.2f    ",events[i].k1.X(),events[i].k1.Y(),events[i].k1.Z());
		fprintf(fp,"%8.2f    %8.2f    %8.2f    ",events[i].p1.X(),events[i].p1.Y(),events[i].p1.Z());
		fprintf(fp,"%8.2f    %8.2f    %8.2f    ",events[i].p2.X(),events[i].p2.Y(),events[i].p2.Z());
		fprintf(fp,"%8d    %8d    ",events[i].shellindex1,events[i].shellindex2);
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void EventParser::write_data(char* fname, std::vector<double>& data) {
	FILE* fp;
	fp = fopen(fname,"w");
	if (fp){
		for (unsigned i=0; i<data.size(); i++)
			fprintf(fp,"%8.4e\n",data[i]);
		fclose(fp);
	} else {
		fprintf(stderr,"Error: could not open file %s. \n",fname);
	}
}

void EventParser::write_data(char* fname, std::vector<double>& data, std::vector<double>& err) {
	FILE* fp;
	fp = fopen(fname,"w");
	if (fp){
		for (unsigned i=0; i<data.size(); i++)
			fprintf(fp,"%8.4e    %8.4e\n",data[i],err[i]);
		fclose(fp);
	} else {
		fprintf(stderr,"Error: could not open file %s. \n",fname);
	}
}

/** no cast checks are done at all, you better pass the right variable type or terrible things will happen
 *  unused for now but in the future if fields change place or new fields are introduced this might become
 *  handy to keep compatibility across versions **/
void EventParser::setattr(Event& e,char* member,void* value){
	if ( strcmp(member,"xB") == 0){
		e.xB    = *(double*) value;
	} else if (strcmp(member,"Q2"   ) == 0){
		e.Q2    = *(double*)   value;
	} else if (strcmp(member,"omega") == 0){
		e.omega = *(double*)   value;
	} else if (strcmp(member,"mass1") == 0){
		e.mass1 = *(double*)    value;
	} else if (strcmp(member,"mass2") == 0){
		e.mass2 = *(double*)    value;
	} else if (strcmp(member,"q"    ) == 0){
		e.q     = *(TVector3*) value;
	} else if (strcmp(member,"k1"   ) == 0){
		e.k1    = *(TVector3*) value;
	} else if (strcmp(member,"p1"   ) == 0){
		e.p1    = *(TVector3*) value;
	} else if (strcmp(member,"p2"   ) == 0){
		e.p2    = *(TVector3*) value;
	} else if (strcmp(member,"shellindex1") == 0) {
		e.shellindex1 = *(int*) value;
	} else if (strcmp(member,"shellindex2") == 0) {
		e.shellindex2 = *(int*) value;
	} else if (strcmp(member,"status") == 0) {
		e.status = *(int*) value;
	} else {
		fprintf(stderr,"Error: EventParser::setattr called on non-supported member field %s of Event\n",member);
		exit(-1);
	}
}

void EventParser::read_events(char* fname, std::vector< struct Event>& events,MeanFieldNucleus& nuc) {
	events.clear();

	std::ifstream file(fname);
	std::string line;
	if (file.is_open()){
		while(getline(file,line)){
			std::vector<double> tokens = splitToDoubles(line,' ');
			if (tokens.size() > 0) { // this is not a comment line
				Event e;
				e.status = (int) tokens[0];
				e.xB = tokens[1];
				e.Q2 = tokens[2];
				e.omega = tokens[3];
				e.q = TVector3(tokens[4],tokens[5],tokens[6]);
				e.k1 = TVector3(tokens[7],tokens[8],tokens[9]);
				e.p1 = TVector3(tokens[10],tokens[11],tokens[12]);
				e.p2 = TVector3(tokens[13],tokens[14],tokens[15]);
				e.shellindex1 = (int) tokens[16];
				e.shellindex2 = (int) tokens[17];

				// assume we are always dealing with protons for now...
				e.type1 = (e.shellindex1 < nuc.getPLevels())? 0 : 1; // 0 is proton,
				e.type2 = (e.shellindex2 < nuc.getPLevels())? 0 : 1; // 0 is proton
				e.mass1 = (e.type1==0)? MASSP : MASSN;
				e.mass2 = (e.type2==0)? MASSP : MASSN;
				events.push_back(e);
			}
		}
		printf("Parser read in %d events.\n",events.size());
	} else {
		fprintf(stderr,"Error: EventParser::read_events: could not open file %s\n",fname);
	}	
}

void EventParser::read_events(char* fname, std::vector< struct Event>& trimmed_events, unsigned int start,unsigned int stop,MeanFieldNucleus& nuc){
	std::vector<struct Event> all_events;
	read_events(fname,all_events,nuc);
	printf("read in all events with size %d \n",all_events.size());
	assert(stop>start && start < all_events.size());// && stop <= all_events.size());
	if (stop >= all_events.size())
		stop = all_events.size();
	trimmed_events.resize(stop-start);

	for (unsigned i=0; i<stop-start;i++){
		printf("setting element %d of trimmed events with %d of all_events.\n",i,start+i);
		trimmed_events[i] = all_events[start+i];
	}
}

std::vector<double> EventParser::splitToDoubles(const string &s, char delim){
	std::vector<double> elems;
	std::vector<string> svec = split(s,delim);
	for (unsigned i=0; i<svec.size(); i++) {
		elems.push_back(atof(svec[i].c_str()));
	}
	return elems;
}

std::vector<string> EventParser::split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	std::string item;
	// use stdlib to tokenize the string
	std::stringstream ss(s);
	while (getline(ss, item, delim))
		if(!item.empty())
			elems.push_back(item);
		
	if (elems[0].c_str()[0] == '#'){ // this is a comment line, clear this shizzle
		elems.clear();
	}
	return elems;
}
