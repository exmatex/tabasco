#include "CoM4.decl.h"
#include "FineScaleModel.hpp"
#include "Interpolate.hpp"
#include "DBInterface.hpp"

extern CProxy_Main mainProxy;
extern CProxy_Interpolate interpolateArray;
extern CProxy_DBInterface DBArray;

FineScaleModel::FineScaleModel()
{
  // Ordering for a 4D array chare is w, x, y, z 
  printf("FineScaleModel created on PE %d Index %d %d %d %d\n", 
      CkMyPe(), thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z);
}

FineScaleModel::FineScaleModel(CkMigrateMessage *msg)
{

}

FineScaleModel::~FineScaleModel()
{

}

void FineScaleModel::pup(PUP::er &p)
{

}

void FineScaleModel::run(int iter)
{
  printf("FineScaleModel %d %d %d %d running iter %d\n",
    thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, iter);

  interpolateArray(thisIndex.w, thisIndex.x, thisIndex.y).run(iter);

  DBArray(thisIndex.w, thisIndex.x, thisIndex.y).run();
}
