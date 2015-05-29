#include "TabaSCo.decl.h"
#include "FineScaleModel.h"
#include "NearestNeighborSearch.h"
#include "Interpolate.h"
#include "DBInterface.h"

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
  p|nbrCount;
  p|nbrData;
}

void FineScaleModel::query2(int iter, int qPt)
{
  currentIter = iter;

  printf("FineScaleModel %d %d %d %d query iter %d\n",
    thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, iter);

  // Check for neighbors
  nnsArray(thisIndex.w, thisIndex.x, thisIndex.y).getNeighbors(thisIndex.z, qPt);

  // There are enough neighbors, do interpolation
  interpolateArray(thisIndex.w, thisIndex.x, thisIndex.y).run(thisIndex.z, nbrCount, nbrData, qPt);

  // There are not enough neighbors, evaluate
  // Or interpolation did not converge
  evaluate(qPt);

  newPt = thisIndex.z;

  // Save new point to DB
  DBArray(thisIndex.w, thisIndex.x, thisIndex.y).put(1, newPt);

  // Send back new point to Coarse Scale model
  //coarseScaleArray(thisIndex.w, thisIndex.x, thisIndex.y).receiveNewPoint(thisIndex.z, currentIter, newPt);
}

void FineScaleModel::evaluate(int qPt)
{
  printf("FineScaleModel evaluate\n");
 
  newPt = thisIndex.z;

  receivePoint(newPt);
}

void FineScaleModel::requestNeighbors(int qPt)
{
  printf("FineScaleModel requestNeighbors\n");
 
  nnsArray(thisIndex.w, thisIndex.x, thisIndex.y).getNeighbors(thisIndex.z, qPt);
}

void FineScaleModel::requestInterpolation(int nbrCount, int nbrData, int qPt)
{
  printf("FineScaleModel requestInterpolation\n");

    interpolateArray(thisIndex.w, thisIndex.x, thisIndex.y).run(thisIndex.z, nbrCount, nbrData, qPt);
}

void FineScaleModel::requestDBStore(int cPt)
{
  printf("FineScaleModel requestDBStore\n");
 
  newPt = cPt;
  DBArray(thisIndex.w, thisIndex.x, thisIndex.y).put(1, newPt);
}

void FineScaleModel::sendNewPoint2Coarse(int elnum, int iter, int cPt)
{
  printf("FineScaleModel sendNewPoint\n");

  newPt = cPt;
  //coarseScaleArray(thisIndex.w, thisIndex.x, thisIndex.y).receiveNewPoint(elnum, iter, newPt);
}
