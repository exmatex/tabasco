#include "CoPPM.decl.h"
#include "Main.hpp"
#include "Domain.hpp"
#include "CoarseScaleModel.hpp"
#include "FineScaleModel.hpp"
#include "DBInterface.hpp"
#include "CoPPM.hpp"
#include "input.hpp"

#include <cstring>

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Domain domainArray;
/* readonly */ CProxy_CoarseScaleModel coarseScaleArray;
/* readonly */ CProxy_FineScaleModel fineScaleArray;
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

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg)
 {
   doneCount = 0;

   // Display info about this run
   
  CkPrintf("**************************************************\n");
  CkPrintf("**************************************************\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**   Running \"Charm++ 3D CoPPM %d processors     **\n", CkNumPes());
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
           "  Chares: %d [%.2g] (%d x %d x %d)\n",
           elemDimX*elemDimY*elemDimZ,
           elemDimX, elemDimY, elemDimZ,
           chareDimX*chareDimY*chareDimZ, charesPerPE,
           chareDimX, chareDimY, chareDimZ);

  totalChares = chareDimX*chareDimY*chareDimZ;
  
  // Setup chare array size
  CkArrayOptions opts(chareDimX, chareDimY, chareDimZ);
  //Setup round robin map
  CProxy_RRMap rrMap = CProxy_RRMap::ckNew();
  opts.setMap(rrMap);
  
  // Create domain chare array
  domainArray = CProxy_Domain::ckNew(opts);
  
  // Create coarse scale models, 1 per domain
  coarseScaleArray = CProxy_CoarseScaleModel::ckNew(opts);

  // Create DB interfaces, 1 per domain
  DBArray = CProxy_DBInterface::ckNew(opts);

  // Create fine scale models
  fineScaleArray = CProxy_FineScaleModel::ckNew();  

  // For Round Robin insertion
  int numPes = CkNumPes();
  int currPE = -1;

  for (int i = 0; i < chareDimX; i++)
  {
    for (int j = 0; j < chareDimY; j++)
    {
      for (int k = 0; k < chareDimZ; k++)
      {
        for (int l = 0; l < numElems; l++)
        {
          fineScaleArray(i, j, k, l).insert( (++currPE) % numPes);
        }
      }
    }
  }

  fineScaleArray.doneInserting();

  // Start coarse scale models
  coarseScaleArray.run(in.maxTimesteps, numElems);

  // Start main program
  //mainProxy.go(in);
  // Print exit msg after sucessful run
  //CkPrintf("Exiting charm++\n");
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
// Executes main of CoPPM
// and exits charm after successful run
void Main::go(Input in)
{
  main_CoPPM(in, domainArray);

  CkExit();
}

#include "CoPPM.def.h"
