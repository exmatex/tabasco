#include "CoM4.decl.h"
#include "FineScaleModel.hpp"

extern CProxy_Main mainProxy;

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
}
