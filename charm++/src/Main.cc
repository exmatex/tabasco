#include "TabaSCo.decl.h"
#include "Main.h"
#include "CoarseScaleModel.h"
#include "FineScaleModel.h"
#include "NearestNeighborSearch.h"
#include "Interpolate.h"
#include "DBInterface.h"
#include "input.h"

#include <cstring>

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_CoarseScaleModel coarseScaleArray;
/* readonly */ CProxy_FineScaleModel fineScaleArray;
/* readonly */ CProxy_NearestNeighborSearch nnsArray;
/* readonly */ CProxy_Interpolate interpolateArray;
/* readonly */ CProxy_DBInterface DBArray;

/*readonly*/ int coarseType;
/*readonly*/ int coarseCount;
/*readonly*/ int NBR_LIMIT;
/*readonly*/ bool useAdaptiveSampling;
/*readonly*/ Real_t stopTime;

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg)
 {
   doneCount = 0;

   // Display info about this run
   
  CkPrintf("**************************************************\n");
  CkPrintf("**************************************************\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**  Running \"Charm++ 3D TabaSCo %d processors   **\n", CkNumPes());
  CkPrintf("**                                              **\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**************************************************\n");
  CkPrintf("**************************************************\n");

  // Set the mainProxy readonly to point to a
  // proxy for the Main chare object (this
  // chare object).
  mainProxy = thisProxy;

  // Get input values from json file
  Input in;
  char input_file[1024];
  strcpy(input_file, msg->argv[1]);
  parse_input((string)input_file, &in);

  // Get coarse model type
  coarseType = in.coarseType;

  // Get number of coarse model chares
  coarseCount = in.coarseCount;

  // Check for adaptive sampling flag
  useAdaptiveSampling = (in.useAdaptiveSampling == 1) ? true : false; 

  // Get simulation stop time
  stopTime = in.stopTime;

  // Neighbor limit
  NBR_LIMIT = 10;


  // Set other parameters

  // Display of some basic information
  // Elements breakdown
  // Chares breakdown [number of chares per processor]
  CkPrintf("Lulesh (Charm++)\n"
           "  Coarse type: %d\n"
           "  Coarse count: %d\n"
           "  UseAdaptiveSampling: %d\n"
           "  Simulation stop time: %e\n", 
          coarseType, coarseCount, 
          ((useAdaptiveSampling == true) ? 1 : 0), stopTime );

  // Setup chare array size
  CkArrayOptions opts(coarseCount);
  //Setup round robin map
  CProxy_RRMap rrMap = CProxy_RRMap::ckNew();
  opts.setMap(rrMap);
  
  // Create coarse scale model, Lulesh
  coarseScaleArray = CProxy_CoarseScaleModel::ckNew(opts);

/*
  // Create DB interfaces
  DBArray = CProxy_DBInterface::ckNew(opts);

  // Create Nearest Neighbor Searches
  nnsArray = CProxy_NearestNeighborSearch::ckNew(opts);

  // Create interpolates
  interpolateArray = CProxy_Interpolate::ckNew(opts);
*/

  // Create fine scale models
  fineScaleArray = CProxy_FineScaleModel::ckNew();  
  fineScaleArray.doneInserting();

  // Start simulation
  run(coarseCount, useAdaptiveSampling, stopTime);

}
  
// Constructor needed for chare object migration (ignore for now)
// // NOTE: This constructor does not need to appear in the ".ci" file
Main::Main(CkMigrateMessage* msg)
{
}

#include "TabaSCo.def.h"
