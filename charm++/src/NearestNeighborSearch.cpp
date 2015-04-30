#include "CoM4.decl.h"
#include "NearestNeighborSearch.hpp"
#include "DBInterface.hpp"

extern CProxy_Main mainProxy;
extern CProxy_DBInterface DBArray;
extern CProxy_FineScaleModel fineScaleArray;
extern int NBR_LIMIT;

NearestNeighborSearch::NearestNeighborSearch()
{
  
  printf("NearestNeighborSearch created on PE %d Index %d %d %d\n", 
      CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.z);
}

NearestNeighborSearch::NearestNeighborSearch(CkMigrateMessage *msg)
{

}

NearestNeighborSearch::~NearestNeighborSearch()
{

}

void NearestNeighborSearch::pup(PUP::er &p)
{
  CBase_NearestNeighborSearch::pup(p);
  p|nbrCount;
  p|nbrIndex;
  p|nbrData;
}

void NearestNeighborSearch::getIndex(int qPt, int nbrCount, int nbrIndex)
{
  printf("NearestNeighborSearch getIndex\n");
}

void NearestNeighborSearch::putIndex(int qPt, int nbrIndex)
{
  printf("NearestNeighborSearch putIndex\n");
}

void NearestNeighborSearch::requestDBGet(int nbrCount, int nbrIndex, int nbrData)
{
  printf("NearestNeighborSearch requestDBGet\n");

  if (nbrCount == NBR_LIMIT)
    DBArray(thisIndex.x, thisIndex.y, thisIndex.z).get(nbrCount, nbrIndex, nbrData);
}

void NearestNeighborSearch::sendNeighbors(int elnum, int nbrCount, int nbrData)
{
  printf("NearestNeighborSearch sendNeighbors\n");

  fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, elnum).receiveNeighbors(nbrCount, nbrData);
}
