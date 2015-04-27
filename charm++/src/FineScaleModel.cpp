#include "CoM4.decl.h"
#include "FineScaleModel.hpp"
#include "NearestNeighborSearch.hpp"
#include "Interpolate.hpp"
#include "DBInterface.hpp"

extern CProxy_Main mainProxy;
extern CProxy_CoarseScaleModel coarseScaleArray;
extern CProxy_NearestNeighborSearch nnsArray;
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
  CBase_FineScaleModel::pup(p);
  p|newPt;
  p|currentIter;
}

void FineScaleModel::query2(int iter)
{
  currentIter = iter;

  printf("FineScaleModel %d %d %d %d query iter %d\n",
    thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, iter);

  // Check for neighbors
  nnsArray(thisIndex.w, thisIndex.x, thisIndex.y).getNeighbors();

  // There are enough neighbors, do interpolation
  interpolateArray(thisIndex.w, thisIndex.x, thisIndex.y).run();

  // There are not enough neighbors, evaluate
  // Or interpolation did not converge
  evaluate();

  newPt = thisIndex.z;

  // Save new point to DB
  DBArray.put();

  // Send back new point to Coarse Scale model
  coarseScaleArray(thisIndex.w, thisIndex.x, thisIndex.y).receiveNewPoint(thisIndex.z, currentIter, newPt);
}

void FineScaleModel::evaluate()
{
  printf("FineScaleModel evaluate\n");
}
