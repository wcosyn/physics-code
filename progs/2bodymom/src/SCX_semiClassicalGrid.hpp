
#ifndef SCX_SEMICLASSICALGRID_HPP
#define SCX_SEMICLASSICALGRID_HPP

#include "MeanFieldNucleus.hpp"
#include "event.hpp"
#include <numint/numint.hpp>

class SCX_semiClassicalGrid{
	public:
		SCX_semiClassicalGrid(MeanFieldNucleusThick* nuc,Event& e,int bpoints, int zpoints);
		SCX_semiClassicalGrid();
		
	private:
		MeanFieldNucleusThick* _nuc;
		Event _event;

};


#endif // SCX_SEMICLASSICALGRID_HPP
