#ifndef _NEARESTNEIGHBORSEARCH_H_
#define _NEARESTNEIGHBORSEARCH_H_

#include "TabaSCo.decl.h"

#ifdef FLANN
#include "ApproxNearestNeighborsFLANN.h"
#else
#include "ApproxNearestNeighborsMTree.h"
#endif

#include "types.h"

#include <vector>

class int_message : public CMessage_int_message {
  public:
    int value;
    int_message() {}
    int_message(int val) : value(val) {}
};

class uint128_message : public CMessage_uint128_message {
  public:
    uint128_t value;
    uint128_message() {}
    uint128_message(uint128_t val) : value(val) {}
};

class knn_message : public CMessage_knn_message {
  public:
    int k;
    int size;
    int *ids;
    uint128_t *keys;
    double *dists;

    knn_message() {}
    knn_message(int k, std::vector<int>& ids, std::vector<uint128_t>& keys, std::vector<double>& dists);
    knn_message(char* buf, int k);
    ~knn_message();
    static void *pack(knn_message *);
    static knn_message *unpack(void *);
    void copyTo(std::vector<int>& ids, std::vector<uint128_t>& keys, std::vector<double>& dists);
};

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
  void insertAsync(std::vector<double> point, uint128_t key);
  void insert(std::vector<uint128_t> keys);
  void remove(int id);

  void insert(std::vector<double> point, uint128_t key, const CkCallback &cb);
  void knn(std::vector<double> x, int k, const CkCallback &cb);
  void getKey(int id, const CkCallback &cb);

};

#endif
