#include "CoM4.decl.h"
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

/*readonly*/ int elemDimX;
/*readonly*/ int elemDimY;
/*readonly*/ int elemDimZ;
/*readonly*/ int blockDimX;
/*readonly*/ int blockDimY;
/*readonly*/ int blockDimZ;
/*readonly*/ int ghostDimX;
/*readonly*/ int ghostDimY;
/*readonly*/ int ghostDimZ;
/*readonly*/ int chareDimX;
/*readonly*/ int chareDimY;
/*readonly*/ int chareDimZ;
/*readonly*/ double charesPerPE;
/*readonly*/ int numElems;
/*readonly*/ int numNodes;
/*readonly*/ int ghostElems;
/*readonly*/ int NBR_LIMIT;

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg)
 {
   doneCount = 0;

   // Display info about this run
   
  CkPrintf("**************************************************\n");
  CkPrintf("**************************************************\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**   Running \"Charm++ 3D CoM4 %d processors     **\n", CkNumPes());
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

  // Get element and block dimension sizes
  elemDimX = in.elemDimX;
  elemDimY = in.elemDimY;
  elemDimZ = in.elemDimZ;

  blockDimX = in.blockDimX;
  blockDimY = in.blockDimY;
  blockDimZ = in.blockDimZ;

  NBR_LIMIT = 10;

  if ((elemDimX < 1) || (elemDimY < 1) || (elemDimZ < 1) ||
      (blockDimZ < 1) || (blockDimY < 1) || (blockDimZ < 1)) {
    CkPrintf("Usage: elemDim and blockDim must take positive values\n");
    CkExit();
  }
  if ((elemDimX % blockDimX) || (elemDimY % blockDimY) || (elemDimZ % blockDimZ)) {
    CkPrintf("Usage: elemDim must be divided evenly by blockDim\n");
    CkExit();
  }

  // Set other parameters
  ghostDimX = blockDimX+1;
  ghostDimY = blockDimY+1;
  ghostDimZ = blockDimZ+1;
  chareDimX = elemDimX/blockDimX;
  chareDimY = elemDimY/blockDimY;
  chareDimZ = elemDimZ/blockDimZ;
  charesPerPE = (double)(chareDimX*chareDimY*chareDimZ)/CkNumPes();
  if(charesPerPE < 1) { charesPerPE = 1; }
  numElems = blockDimX*blockDimY*blockDimZ;
  numNodes = ghostDimX*ghostDimY*ghostDimZ;
  ghostElems = numElems + 2*blockDimY*blockDimZ
               + 2*blockDimX*blockDimZ + 2*blockDimX*blockDimY;

  // Display of some basic information
  // Elements breakdown
  // Chares breakdown [number of chares per processor]
  CkPrintf("Lulesh (Charm++)\n"
           "  Elements: %d (%d x %d x %d)\n"
           "  Chares: %d [%.2g] (%d x %d x %d)\n"
           "  Num Elements = %d\n",
           elemDimX*elemDimY*elemDimZ,
           elemDimX, elemDimY, elemDimZ,
           chareDimX*chareDimY*chareDimZ, charesPerPE,
           chareDimX, chareDimY, chareDimZ,
           numElems);

  totalChares = chareDimX*chareDimY*chareDimZ;
  
  // Setup chare array size
  CkArrayOptions opts(chareDimX, chareDimY, chareDimZ);
  //Setup round robin map
  CProxy_RRMap rrMap = CProxy_RRMap::ckNew();
  opts.setMap(rrMap);
  
  // Create coarse scale model
  coarseScaleArray = CProxy_CoarseScaleModel::ckNew(opts);

  // Create DB interfaces
  DBArray = CProxy_DBInterface::ckNew(opts);

  // Create Nearest Neighbor Searches
  nnsArray = CProxy_NearestNeighborSearch::ckNew(opts);

  // Create interpolates
  interpolateArray = CProxy_Interpolate::ckNew(opts);

  // Create fine scale models
  fineScaleArray = CProxy_FineScaleModel::ckNew();  
  fineScaleArray.doneInserting();

  // Start coarse scale models
  coarseScaleArray.run(in.maxTimesteps, numElems);

}
  
// Constructor needed for chare object migration (ignore for now)
// // NOTE: This constructor does not need to appear in the ".ci" file
Main::Main(CkMigrateMessage* msg)
{
}

// End processing
void Main::done()
{
  doneCount++;
  if (doneCount >= totalChares)
    CkExit();
}

#include "CoM4.def.h"
