/** A simple C-interface for MeanFieldNucleus for use with for 
 *  example ctypes-Python.
 *  It is not the cleanest interface as a lot of void* pointer
 *  passing is done but it results in relatively short and easy
 *  to read "glue"/export headers.
 *
 *  @author Camille Colle
 */

#include "MeanFieldNucleus.hpp"
#define SHAREDIR "/home/camille/Code/share"

extern "C" {
	void* get_MeanFieldNucleus(int nuctype) {
		return new MeanFieldNucleus(nuctype,SHAREDIR);
	}
	
	int getZ(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getZ();
	}

	int getA(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getA();
	}

	int getPLevels(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getPLevels();
	}
	
	int getNLevels(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getNLevels();
	}
	
	int getTotalLevels(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getTotalLevels();
	}
	
	int getFinalMProton(void* nuc){
		return ((MeanFieldNucleus*) nuc)->getFinalMProton();
	}
	
	int getFinalMNeutron(void* nuc){
		return ((MeanFieldNucleus*) nuc)->getFinalMNeutron();
	}
	
	const double* getExcitation(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getExcitation();
	}
	
	double getMassA(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getMassA();
	}
	
	double getMassA_min_pp(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getMassA_min_pp();
	}
	
	double getMassA_min_pn(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getMassA_min_pn();
	}
	
	double getMassA_min_nn(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getMassA_min_nn();
	}
	
	const int* getN_array(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getN_array();
	}
	
	const int* getJ_array(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getJ_array();
	}
 
	const int* getL_array(void* nuc) {
		return ((MeanFieldNucleus*) nuc)->getL_array();
	}
	void del_MeanFieldNucleus(void* nuc) {
		delete (MeanFieldNucleus*) nuc;
	}

}
