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
/*readonly*/ int fineType;
/*readonly*/ int nnsType;
/*readonly*/ int pointDim;
/*readonly*/ int numTrees;

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

  // Get coarse model parameters - type, count, adptive sampling flag
  coarseType = in.coarseType;
  coarseCount = in.coarseCount;
  useAdaptiveSampling = (in.useAdaptiveSampling == 1) ? true : false; 

  // Get fine model parameters - type
  fineType = in.fineType;

  // Get Nearest Neighbor Search parameters
  nnsType = in.nnsType;
  pointDim = in.pointDim;
  numTrees = in.numTrees;

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
           "  Simulation stop time: %e\n"
           "  Fine type: %d\n"
           "  NNS type: %d\n"
           "  NNS point dimension: %d\n"
           "  NNS number of trees: %d\n", 
          coarseType, coarseCount, 
          ((useAdaptiveSampling == true) ? 1 : 0), stopTime,
          fineType, nnsType, pointDim, numTrees);

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
*/

  // Create Nearest Neighbor Searches
  nnsArray = CProxy_NearestNeighborSearch::ckNew(opts);
  nnsArray.initialize(nnsType, pointDim, numTrees);

/*
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
