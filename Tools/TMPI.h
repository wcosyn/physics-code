/*
 * TMPI.h
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: February,25 2009
 *
 */

#ifndef TMPI_H
#define TMPI_H
#define STRANGEDEUTERON_MPI

#define NROFRES 9

#if defined(STRANGEDEUTERON_MPI) && !defined(__MAKECINT__)
#include <mpi.h>
#else
typedef int MPI_Op;
typedef int MPI_Datatype;
#endif
#include <fstream>
#include <TObject.h>

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TMPI                                                                  //
//                                                                       //
// Global object coordinating multi-process calculations                 //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

class TMPI : public TObject /* FINAL */
{
 public:
  // Store results to be gathered in TMPI::GatherResults.
  typedef struct { 
    double value[NROFRES];
    int    save; 
  } ValueSave_t;

 public:
  TMPI(int*,char***); // can only be called once.
  ~TMPI();
  
  static bool IsMaster() { return fgRank==kMasterRank; }
  static bool IsSlave() { return fgRank!=kMasterRank; }
  static int Rank() { return fgRank; }
  static int NumberOfProcesses() { return fgNproc; }
  static int RankOfMaster() { return kMasterRank; }
 
  static std::ostream& Cout();
  static std::ostream& Cerr();
  static void SilenceSlaves(bool=true);

  static void GatherResults(int,double**);
      
 private:
  TMPI(const TMPI&);
  TMPI& operator=(const TMPI&);
  static void CombineResults(void*,void*,int*,MPI_Datatype*);

 private:
  static const int kMasterRank=0; // Rank of the master process
  static int fgRank;  // The rank of the current process
  static int fgNproc; // The total number of processes
  static bool fgMpiInitialized; // Has a TMPI object been created?
  static std::ofstream fgDevNull; // stream to /dev/null
  static MPI_Op MPI_COMBINE_RESULTS; // Our custom MPI operator
  static MPI_Datatype datatype;
  ClassDef(TMPI,0); // Global object coordinating multi-process calculations

}; // class TMPI

#endif // TMPI_H
