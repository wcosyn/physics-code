/*
 * TMPI.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: February,25 2009
 *
 */

#ifdef STRANGEDEUTERON_MPI
#include <mpi.h>
#endif
#include "TMPI.h"
#include <cstdlib>
#include <cstdio>
#include <iostream>

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TMPI                                                                  //
//                                                                       //
// Global object coordinating multi-process calculations                 //
//
// TMPI makes sure the strangecalc-deuteron library works in singleton
// and multi-process mode.
//
// The user is advised to create a single instance of TMPI at the start
// of main(). When this object goes out-of-scope, the MPI environment
// will be finalized correctly.
//
// Attempts to create a TMPI object when strangecalc-deuteron wasn't
// compiled with mpicc will result in a runtime error.
//
// TMPI provides some functionality that may be useful in a multi-process
// situation:
//  - Cout()
//  - Cerr()
//  - IsMaster()
//  - IsSlave()
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TMPI)

//_____________________________________________________________________
std::ofstream TMPI::fgDevNull("/dev/null");

//_____________________________________________________________________
MPI_Op TMPI::MPI_COMBINE_RESULTS;

//_____________________________________________________________________
bool TMPI::fgMpiInitialized = false;

//_____________________________________________________________________
int TMPI::fgRank = 0;

//_____________________________________________________________________
int TMPI::fgNproc = 1;
MPI_Datatype TMPI::datatype = 0;


//_____________________________________________________________________
TMPI::TMPI(int* argc, char*** argv)
  : TObject()
{
  // Constructor
  //
  // Creating a TMPI object will initialize the MPI environment. 
  // Do not call MPI_init after calling this constructor.
  //
  // For obvious reasons, users are not allowed to create more than
  // one instance of TMPI.
#ifndef STRANGEDEUTERON_MPI
  Warning("TMPI(int*,char***)","strangecalc-deuteron wasn't compiled for MPI.");

#else
  // Make sure this constructor is only called once per process
  if( fgMpiInitialized )
    Cerr() << "TMPI::TMPI(int*,char***): Only one instance of TMPI allowed."
	   << std::endl;
  
  MPI_Init(argc,argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &fgRank);
  MPI_Comm_size(MPI_COMM_WORLD, &fgNproc);

  MPI_Op_create(&CombineResults, true, &MPI_COMBINE_RESULTS); 

  fgMpiInitialized = true;
  int count =2;
  int lengths[2] = {NROFRES,1};
  MPI_Aint offsets[2] = {0,NROFRES*sizeof(double)};
  MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
  MPI_Type_struct(count, lengths, offsets, types, &datatype);
  MPI_Type_commit(&datatype);
  //std::cout << datatype << std::endl;
#endif
}

//_____________________________________________________________________
TMPI::~TMPI()
{
  // Destructor
  //
  // A call to MPI_Finalize will be made. The user doesn't need to worry about this.
#ifdef STRANGEDEUTERON_MPI
cout << TMPI::Rank() << " finalize" << endl; 
  MPI_Finalize();
#endif
}

//_____________________________________________________________________
TMPI::TMPI(const TMPI&)
{
  // Copy constructor is meaningless
  std::exit(1);
}

//_____________________________________________________________________
TMPI&
TMPI::operator=(const TMPI&)
{
  // Assignment operator is meaningless
  std::exit(1);
}

//_____________________________________________________________________
std::ostream&
TMPI::Cout()
{
  // Redirect standard output to /dev/null in case of slave processes.
  if( IsSlave() ) return fgDevNull;
  return std::cout;
}

//_____________________________________________________________________
std::ostream&
TMPI::Cerr()
{
  // Redirect standard error output to /dev/null in case of slave processes.
  if( IsSlave() ) return fgDevNull;
  return std::cerr;
}

//_____________________________________________________________________
void TMPI::SilenceSlaves(bool off)
{
  // Turns off the output of all slaves.
  // This is done by redirecting cout and cerr to /dev/null.
  //
  // This function can be used to suppress multiple output from 
  // functions from different libraries. Don't forget to turn 
  // the output back on at the end.
  if( IsMaster() ) return;
  
  // C++ Stream buffers
  static std::streambuf *devNullBuffer = fgDevNull.rdbuf();
  static std::streambuf *coutBuffer = std::cout.rdbuf();
  static std::streambuf *cerrBuffer = std::cerr.rdbuf();

  // C file descriptor and position
  static int fd;
  static fpos_t pos;
  
  std::cout.flush();
  fflush(stdout);

  if( off ) {
    std::cout.rdbuf( devNullBuffer );
    std::cerr.rdbuf( devNullBuffer );
    fgetpos(stdout, &pos);
    fd = dup(fileno(stdout));
    freopen("/dev/null", "w", stdout);
  }
  else {
    std::cout.rdbuf( coutBuffer );
    std::cerr.rdbuf( cerrBuffer );
    dup2(fd, fileno(stdout));
    close(fd);
    clearerr(stdout);
    fsetpos(stdout, &pos);
  }
}

//_____________________________________________________________________
void
TMPI::GatherResults(int length, double **valueArray)
{
  // If a for loop is parallelized as follows:
  // for(int i=TMPI::Rank(); i<max; i += TMPI::NumberOfProcesses() ) { x[i]=...}
  // Than a call to TMPI::GatherResults(max,x), will make sure
  // that the array x of the master process holds all the results
  // from the slaves.
  //
  // Obviously this function doesn't do anything if strangecalc-deuteron 
  // isn't compiled with MPI or when the app is run with a single CPU.
#ifdef STRANGEDEUTERON_MPI
  if( NumberOfProcesses() > 1 ) {
    // We store the results of each process in a struct of type
    // ValueSave_t. The 'save' member specifies whether the result
    // should be witheld in the final result.
    ValueSave_t *valueSaveArray_in = new ValueSave_t[length];
    for(int i=0; i<length; ++i) {
      if( TMPI::Rank()==(i%NumberOfProcesses()) ) {
	for(int j=0;j<NROFRES;j++){
	  valueSaveArray_in[i].value[j] = valueArray[i][j];
	}
	valueSaveArray_in[i].save = 1;
      } else {
	for(int j=0;j<NROFRES;j++) valueSaveArray_in[i].value[j] = 0.;
	valueSaveArray_in[i].save = 0;
      }
    }
    ValueSave_t *valueSaveArray_out = new ValueSave_t[length];

    MPI_Reduce( valueSaveArray_in, valueSaveArray_out, length, 
		datatype, MPI_COMBINE_RESULTS, 
		kMasterRank, MPI_COMM_WORLD);
  
    if( IsMaster() )
      for(int i=0; i<length; ++i)
	for(int j=0; j<NROFRES; j++) {
	  valueArray[i][j] = valueSaveArray_out[i].value[j];
	}

    delete[] valueSaveArray_in;
    delete[] valueSaveArray_out;
  }
#endif
}

//_____________________________________________________________________
void
TMPI::CombineResults(void* invec, void* inoutvec,
		     int* len, MPI_Datatype* datatype)
{
    
  // This is our custom MPI operator that is used in the MPI_Reduce
  // call in TMPI::GatherResults(int,double*). Look at the documentation
  // of that function to understand what happens here.
  //
  // The MPI operator is initialized in TMPI's constructor.
  for(int i=0; i<(*len); ++i) {
    if( ((ValueSave_t*)invec)[i].save ) {
      ((ValueSave_t*)inoutvec)[i].save = ((ValueSave_t*)invec)[i].save;
      for(int j=0;j<NROFRES;j++) ((ValueSave_t*)inoutvec)[i].value[j] = ((ValueSave_t*)invec)[i].value[j];
    }
  }
}
