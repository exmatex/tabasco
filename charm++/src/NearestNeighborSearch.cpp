#include "CoM4.decl.h"
#include "NearestNeighborSearch.hpp"
#include "DBInterface.hpp"

extern CProxy_Main mainProxy;
extern CProxy_DBInterface DBArray;

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

}

void NearestNeighborSearch::getIndex()
{
  printf("NearestNeighborSearch getIndex\n");
}

void NearestNeighborSearch::putIndex()
{
  printf("NearestNeighborSearch putIndex\n");
}

void NearestNeighborSearch::getNeighbors()
{
  printf("NearestNeighborSearch getNeighbors\n");

  // Check for neighbor indices
  getIndex();
  
  // Get data from DB for each neighbor
  DBArray.get();
}

void NearestNeighborSearch::get()
{
  printf("NearestNeighborSearch get\n");
}

