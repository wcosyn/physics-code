#ifndef PARSER_HPP
#define PARSER_HPP

#include "event.hpp"
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cassert>
#include "MeanFieldNucleus.hpp"

/** write and read in two nucleon knockout events **/
namespace EventParser {
		void write_events(char*,std::vector< struct Event >& );
		void read_events (char*,std::vector< struct Event >&,MeanFieldNucleus& );
		void read_events (char*,std::vector< struct Event >&,unsigned int,unsigned int,MeanFieldNucleus& );
		void write_data  (char*,std::vector< double >&);
		void write_data  (char*,std::vector< double >&, std::vector< double >&);
		void setattr(Event&,char*,void*); // pass void pointer to cast to correct value

		std::vector<double> splitToDoubles(const std::string&, char delim);
		std::vector<string> split(const std::string&, char delim);
};


#endif // parser.hpp
