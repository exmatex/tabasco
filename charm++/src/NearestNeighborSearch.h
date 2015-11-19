#ifndef _NEARESTNEIGHBORSEARCH_H_
#define _NEARESTNEIGHBORSEARCH_H_

#include "TabaSCo.decl.h"

#ifdef FLANN
#include "ApproxNearestNeighborsFLANN.h"
#else
#include "ApproxNearestNeighborsMTree.h"
#endif

class NearestNeighborSearch : public CBase_NearestNeighborSearch {
  private:
    NearestNeighborSearch_SDAG_CODE;
    int nbrCount;
    int nbrIndex;   
    int nbrData;   

    ApproxNearestNeighbors* ann;

  public:
  
  NearestNeighborSearch();
  NearestNeighborSearch(CkMigrateMessage *msg);
  ~NearestNeighborSearch();
  void pup(PUP::er &p);

  // Entry methods
  void initialize(int ntype, int dim, int ntrees);
  int insert(std::vector<double>& point, uint128_t& key);
  void insert(std::vector<uint128_t> const& keys);
  void remove(int id);
  void knn(std::vector<double> const& x,
                     int k,
                     std::vector<int> &ids,
                     std::vector<uint128_t> &keys,
                     std::vector<double> &dists);
  uint128_t getKey(int id);

};

#endif
