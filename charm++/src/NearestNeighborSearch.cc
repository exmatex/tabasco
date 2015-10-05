#include "TabaSCo.decl.h"
#include "NearestNeighborSearch.h"
#include "DBInterface.h"

#include "adaptive_sampling/interpolation_database/kriging_database/ApproxNearestNeighborsFLANN.h"

#include <vector>

extern CProxy_Main mainProxy;
extern CProxy_DBInterface DBArray;
extern CProxy_FineScaleModel fineScaleArray;
extern int NBR_LIMIT;

NearestNeighborSearch::NearestNeighborSearch()
{
  printf("NearestNeighborSearch created on PE %d Index %d\n",
      CkMyPe(), thisIndex);

}

void NearestNeighborSearch::initialize(int ntype, int dim, int ntrees)
{
  int nchecks = 20;

  if (ntype == 0) 
  {
    //ann = (ApproxNearestNeighbors*)(new ApproxNearestNeighborsFLANN(dim, ntrees, nchecks));
    printf("FLANN NNS not added yet\n");
  }
  else
  {
    printf("MTree NNS not added yet\n");
  }

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

int NearestNeighborSearch::insert(std::vector<double>& point, uint128_t& key)
{
  ann->insert(point, key);
}

void NearestNeighborSearch::insert(std::vector<uint128_t> const& keys)
{
  ann->insert(keys);
}

void NearestNeighborSearch::remove(int id)
{
  ann->remove(id);
}

void NearestNeighborSearch::knn(std::vector<double> const& x,
                     int k,
                     std::vector<int> &ids,
                     std::vector<uint128_t> &keys,
                     std::vector<double> &dists)
{
  ann->knn(x, k, ids, keys, dists);
}

uint128_t NearestNeighborSearch::getKey(int id)
{
  uint128_t rkey = ann->getKey(id);
  return rkey; 
}

/*
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

  //if (nbrCount == NBR_LIMIT)
    DBArray(thisIndex.x, thisIndex.y, thisIndex.z).get(nbrCount, nbrIndex, nbrData);
}

void NearestNeighborSearch::sendNeighbors(int elnum, int nbrCount, int nbrData)
{
  printf("NearestNeighborSearch sendNeighbors\n");

  nbrCount = NBR_LIMIT;
  //fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, elnum).receiveNeighbors(nbrCount, nbrData);
}
*/
