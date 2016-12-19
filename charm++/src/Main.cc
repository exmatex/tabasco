#include "TabaSCo.decl.h"
#include "Main.h"
#include "CoarseScaleModel.h"
#include "FineScaleModel.h"
#include "NearestNeighborSearch.h"
#include "Interpolate.h"
#include "Evaluate.h"
#include "DBInterface.h"
#include "DBMap.h"
#include "DBVecMessage.h"
#include "input.h"
#include "ModelDB_Enums.h"
#include "NodeMaps.h"

#include <cstring>
#include <vector>

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_CoarseScaleModel coarseScaleArray;
/* readonly */ CProxy_FineScaleModel fineScaleArray;
/* readonly */ CProxy_NearestNeighborSearch nnsArray;
/* readonly */ CProxy_Interpolate interpolateArray;
/* readonly */ CProxy_Evaluate evaluateArray;
/* readonly */ CProxy_DBInterface DBArray;
/* readonly */ CProxy_KaswellMap KaswellMapProxy;
/* readonly */ CProxy_HaswellMap HaswellMapProxy;

/*readonly*/ int coarseType;
/*readonly*/ int coarseCount;
/*readonly*/ int NBR_LIMIT;
/*readonly*/ bool useAdaptiveSampling;
/*readonly*/ Real_t stopTime;
/*readonly*/ int maxSteps;
/*readonly*/ int fineType;
/*readonly*/ int nnsType;
/*readonly*/ int nnsCount;
/*readonly*/ int pointDim;
/*readonly*/ int numTrees;
/*readonly*/ int interpType;
/*readonly*/ int interpCount;
/*readonly*/ int evalType;
/*readonly*/ int evalCount;
/*readonly*/ int dbType;
/*readonly*/ int dbCount;
/*readonly*/ bool dbRemote; 

/*readonly*/ int file_parts; 
/*readonly*/ int visit_data_interval; 
/*readonly*/ int edgeElems; 
/*readonly*/ int heightElems;
/*readonly*/ int timerRate;

///TODO: Make these less static
const int NUM_Kaswell_RANKS = 8;


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

  if(msg->argc < 2){
    std::cerr << "Missing argument (json file)" << std::endl;
    CkExit();
  }
  // Get input values from json file
  Input in;
  char input_file[1024];
  strcpy(input_file, msg->argv[1]);
  parse_input((string)input_file, &in);

  // Get coarse model parameters - type, count, adptive sampling flag
  coarseType = in.coarseType;
  coarseCount = in.coarseCount;
  useAdaptiveSampling = (in.useAdaptiveSampling == 1) ? true : false;
  stopTime = in.stopTime;
  maxSteps = in.maxSteps;
  edgeElems = in.edgeElems;
  heightElems = in.heightElems;
  timerRate = in.timerRate;

  // Get fine model parameters - type
  fineType = in.fineType;

  // Get Nearest Neighbor Search parameters
  nnsType = in.nnsType;
  nnsCount= in.nnsCount;
  pointDim = in.pointDim;
  numTrees = in.numTrees;

  // Get Interpolate parameters
  interpType = in.interpType;
  interpCount = in.interpCount;

  // Get Evaluate parameters
  evalType = in.evalType;
  evalCount = in.evalCount;

  // Get DBInterface parameters
  int tempDBCount;
  tempDBCount = in.dbCount;
  dbRemote = (in.dbRemote != 0) ? true : false;
  if(dbRemote == false and ( in.dbType == 0 or in.dbType == 4))
  {
	//Specifically ensure we can't have local redis
	CkError("WARNING: Non-Remote Redis is Not Supported\n");
	dbType = 1;
  }
  else
  {
	  dbType = in.dbType;
  }

  // Neighbor limit
  NBR_LIMIT = 10;

  // Get SILO parameters
  file_parts = (in.file_parts != 0) ? in.file_parts : 1;
  visit_data_interval = in.visit_data_interval;

  //Set DBCount before moving forward for output purposes
  // Create DB interfaces with desired backend, if we need them
  ///TODO: Probably should merge the db and kaswell map because we may not want redis on a kaswell...
  if(dbRemote == true)
  {
    int numDBs = tempDBCount;
    if(CmiCpuTopologyEnabled())
    {
        numDBs = std::min(CmiNumPhysicalNodes(), tempDBCount);
    }
    //CkPrintf("numDBs = %d\n", numDBs);
    CProxy_DBMap DBMapProxy = CProxy_DBMap::ckNew(numDBs);
    //Create array options for DB Interface
    CkArrayOptions dbMapOptions(numDBs);
    dbMapOptions.setMap(DBMapProxy);

    DBArray = CProxy_DBInterface::ckNew(dbType, dbMapOptions);
    dbCount = numDBs;
  }
  else
  {
    dbCount = 0;
  }

  // Set other parameters

  // Build Kaswell map
	std::vector<int> kaswellVec;
	std::vector<int> haswellVec;
  for(int node = 0; node < CmiNumPhysicalNodes(); node++)
  {
    if(this->isNodeKaswell(node) == true)
    {
      //Push to kaswell list
	  int * peList;
	  int nPEs;
	  CmiGetPesOnPhysicalNode(node, &peList, &nPEs);
	  for(int i = 0; i < nPEs; i++)
	  {
        kaswellVec.push_back(peList[i]);
	  }
    }
	else
	{
	  int * peList;
	  int nPEs;
	  CmiGetPesOnPhysicalNode(node, &peList, &nPEs);
	  for(int i = 0; i < nPEs; i++)
	  {
        haswellVec.push_back(peList[i]);
	  }
	}
  }
  //Add safety so that this can work even if not on hybrid platform
  //If no kaswells, everything is a kaswell
  if(kaswellVec.size() == 0)
  {
	for(int node = 0; node < CmiNumPhysicalNodes(); node++)
	{
	  int * peList;
	  int nPEs;
	  CmiGetPesOnPhysicalNode(node, &peList, &nPEs);
	  for(int i = 0; i < nPEs; i++)
	  {
        kaswellVec.push_back(peList[i]);
	  }
	}
  }
  //And if no haswells, everything is a haswell
  if(haswellVec.size() == 0)
  {
	for(int node = 0; node < CmiNumPhysicalNodes(); node++)
	{
	  int * peList;
	  int nPEs;
	  CmiGetPesOnPhysicalNode(node, &peList, &nPEs);
	  for(int i = 0; i < nPEs; i++)
	  {
        haswellVec.push_back(peList[i]);
	  }
	}
  }
  //Build maps
	KaswellMapProxy = CProxy_KaswellMap::ckNew(kaswellVec);
	HaswellMapProxy = CProxy_HaswellMap::ckNew(haswellVec);

  // Display of some basic information
  // Elements breakdown
  // Chares breakdown [number of chares per processor]
  CkPrintf("Lulesh (Charm++)\n"
           "  Coarse type: %d\n"
           "  Coarse count: %d\n"
           "  UseAdaptiveSampling: %d\n"
           "  Simulation stop time: %e\n"
           "  Simulations max steps: %d\n"
           "  Simulation Height Elems: %d\n"
           "  Simulation Edge Elems: %d\n"
           "  Fine type: %d\n"
           "  NNS type: %d\n"
           "  NNS count: %d\n"
           "  NNS point dimension: %d\n"
           "  NNS number of trees: %d\n"
           "  Interpolate type: %d\n"
           "  Interpolate count: %d\n"
           "  Evaluate type: %d\n"
           "  Evaluate count: %d\n"
           "  DBInterface type: %s\n"
           "  DBInterface count: %d\n"
           "  Use Remote DB %d\n"
           "  Number of SILO Files for Single Domain: %d\n"
           "  Visit Data Interval: %d\n"
           "  Adaptive Timer Sampling Rate: %d\n"
		   "  Number of Kaswell PEs: %d\n"
		   "  Number of Haswell PEs: %d\n",
          coarseType, coarseCount, 
          ((useAdaptiveSampling == true) ? 1 : 0), stopTime, maxSteps, heightElems, 
          edgeElems, fineType, nnsType, nnsCount, pointDim, numTrees,
          interpType, interpCount, evalType, evalCount, 
          SingletonDBBackendStrings[dbType], dbCount, dbRemote, file_parts, 
          visit_data_interval, timerRate, kaswellVec.size(), haswellVec.size());

  //Setup round robin map
  CProxy_RRMap rrMap = CProxy_RRMap::ckNew();
  
  // Setup chare array size
/*
  CkArrayOptions coarseOpts(coarseCount);
  coarseOpts.setMap(rrMap);
  coarseScaleArray = CProxy_CoarseScaleModel::ckNew(coarseOpts);
*/

  // Evenly distribute CoarseScaleModels across PEs
  CkArrayOptions coarseScaleOptions(coarseCount);
  coarseScaleOptions.setMap(HaswellMapProxy);
  coarseScaleArray = CProxy_CoarseScaleModel::ckNew(haswellVec.data(), haswellVec.size(), coarseCount);
  /*
  int inc = CkNumPes() / coarseCount;
  int whichPE = 0;
  for (int i = 0; i < coarseCount; i++)
  {
    coarseScaleArray(i).insert(whichPE);
    whichPE += inc;
  }
  */
  /*
  for(int i = 0; i < coarseCount; i++)
  {
	int whichPE = haswellVec[i % haswellVec.size()];
	coarseScaleArray(i).insert(haswellVec.data(), haswellVec.size(), whichPE);
  }
  */
  coarseScaleArray.doneInserting();
  
#ifdef SILO
  coarseScaleArray.setSiloParams(file_parts, visit_data_interval);
#endif
  coarseScaleArray.setRemoteDB(dbRemote);

  // Create Nearest Neighbor Searches
#ifdef NNS_AS_CHARE
  CkArrayOptions nnsOpts(nnsCount);
  nnsOpts.setMap(HaswellMapProxy);
  nnsArray = CProxy_NearestNeighborSearch::ckNew(nnsOpts);
  nnsArray.initialize(nnsType, pointDim, numTrees);
#endif

/*
  // Create interpolates
  interpolateArray = CProxy_Interpolate::ckNew(opts);
*/

  // Create Evaluates
#ifdef EVAL_AS_CHARE
  CkArrayOptions evalOpts(evalCount);
  evalOpts.setMap(KaswellMapProxy);
  evaluateArray = CProxy_Evaluate::ckNew(evalOpts);
  evaluateArray.initialize(evalType);
#endif

  // Create fine scale models
  fineScaleArray = CProxy_FineScaleModel::ckNew();  
  fineScaleArray.doneInserting();

  // Start simulation
  run(coarseCount, nnsCount, interpCount, evalCount, dbCount, useAdaptiveSampling, stopTime, maxSteps, timerRate);

}
  
// Constructor needed for chare object migration (ignore for now)
// // NOTE: This constructor does not need to appear in the ".ci" file
Main::Main(CkMigrateMessage* msg)
{
}

bool Main::isNodeKaswell(int physNodeID)
{
	///TODO: Consider a more robust check to determine if this is a ``kaswell'' node
	int nRanks = CmiNumPesOnPhysicalNode(physNodeID);
	if(nRanks == NUM_Kaswell_RANKS)
	{
		return true;
	}
	else
	{
		return false;
	}
}

#include "TabaSCo.def.h"
