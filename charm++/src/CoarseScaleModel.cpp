#include "CoM4.decl.h"
#include "CoarseScaleModel.hpp"
#include "DBInterface.hpp"
#include "FineScaleModel.hpp"

extern CProxy_Main mainProxy;
extern CProxy_FineScaleModel fineScaleArray;

CoarseScaleModel::CoarseScaleModel()
{
  
  printf("CoarseScaleModel created on PE %d Index %d %d %d\n", 
      CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.z);
}

CoarseScaleModel::CoarseScaleModel(CkMigrateMessage *msg)
{

}

CoarseScaleModel::~CoarseScaleModel()
{

}

void CoarseScaleModel::pup(PUP::er &p)
{

}

void CoarseScaleModel::run(int ntimesteps, int nelems)
{
  maxTimesteps = ntimesteps;
  numElems = nelems;
  //printf("CoarseScaleModel running\n");

  for (int i = 0; i < maxTimesteps;i++)
  {
    for (int j = 0; j < numElems; j++)
    {
      fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, j).run(i);
    }
  }

  // This chare is done
  //printf("CoarseScaleModel is done\n");
  mainProxy.done();
}
