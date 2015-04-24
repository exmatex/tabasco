#include "CoM4.decl.h"
#include "NearestNeighborSearch.hpp"

extern CProxy_Main mainProxy;

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

void NearestNeighborSearch::run(int iter)
{
  printf("NearestNeighborSearch running iter %d\n", iter);
}
