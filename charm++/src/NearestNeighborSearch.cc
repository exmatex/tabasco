#include "TabaSCo.decl.h"
#include "NearestNeighborSearch.h"
#include "DBInterface.h"

// Approximate nearest neighbor search options
#include "ApproxNearestNeighborsFLANN.h"
#include "ApproxNearestNeighborsMTree.h"

extern CProxy_Main mainProxy;
extern CProxy_DBInterface DBArray;
extern CProxy_FineScaleModel fineScaleArray;
extern CProxy_NearestNeighborSearch nnsArray;
extern int NBR_LIMIT;

knn_message::knn_message(int k, std::vector<int>& vids, std::vector<uint128_t>& vkeys, std::vector<double>& vdists) : k(k), size(vids.size())
{
  ids = new int[size];
  keys = new uint128_t[size];
  dists = new double[size];

  for (int i = 0; i < size; i++)
  {
    ids[i] = vids[i];
  }
  for (int i = 0; i < size; i++)
  {
    keys[i] = vkeys[i];
  }
  for (int i  = 0; i < size; i++)
  {
    dists[i] = vdists[i];
  }
}

knn_message::knn_message(char* buf, int k) : k(k)
{
  char* p = buf;

  memcpy(&k, p, sizeof(int));
  p += sizeof(int);
  memcpy(&size, p, sizeof(int));
  p += sizeof(int);

  ids = new int[size];  
  keys = new uint128_t[size];
  dists = new double[size];

  memcpy(ids, p, size*sizeof(int));
  p += size*sizeof(int);
  memcpy(keys, p, size*sizeof(uint128_t));
  p += size*sizeof(uint128_t);
  memcpy(dists, p, size*sizeof(double));
}

void *
knn_message::pack(knn_message* m)
{
  int msize = 2 * sizeof(int)+m->size*sizeof(int)+m->size*sizeof(double)+m->size*sizeof(uint128_t);
  char *p = (char*)CkAllocBuffer(m, msize);
  char *buf = p;

  memcpy(buf, &m->k, sizeof(int));
  buf += sizeof(int);
  memcpy(buf, &m->size, sizeof(int));
  buf += sizeof(int);
  memcpy(buf, m->ids, m->size * sizeof(int));
  buf += m->size * sizeof(int);
  memcpy(buf, m->keys, m->size * sizeof(uint128_t));
  buf += m->size * sizeof(uint128_t);
  memcpy(buf, m->dists, m->size * sizeof(double));

  CkFreeMsg(m);
  return (void*) p;
}

knn_message *
knn_message::unpack(void *inbuf)
{
   char *buf = (char *) inbuf;
   int k;
   memcpy(&k, buf, sizeof(int));

   knn_message* pmsg = (knn_message*)CkAllocBuffer(inbuf, sizeof(knn_message));
   pmsg = new ((void*)pmsg) knn_message(buf, k);

   CkFreeMsg(buf);
   return pmsg;
}

void
knn_message::copyTo(std::vector<int>& vids, std::vector<uint128_t>& vkeys, std::vector<double>& vdists)
{
  for (int i = 0; i < size; i++)
  {
    vids.push_back(ids[i]);
  }

  for (int i = 0; i < size; i++)
  {
    vkeys.push_back(keys[i]);
  }

  for (int i = 0; i < size; i++)
  {
    vdists.push_back(dists[i]);
  }
}

knn_message::~knn_message()
{
  delete [] ids;
  delete [] keys;
  delete [] dists;
}

NearestNeighborSearch::NearestNeighborSearch()
{
/*
  printf("NearestNeighborSearch created on PE %d Index %d\n",
      CkMyPe(), thisIndex);
*/

}

void NearestNeighborSearch::initialize(int ntype, int dim, int ntrees)
{
  if (ntype == 0) 
  {
#ifdef FLANN
    int nchecks = 20;
    ann = (ApproxNearestNeighbors*)(new ApproxNearestNeighborsFLANN(dim, ntrees, nchecks));
    //printf("FLANN NNS creatd.\n");
#endif
  }
  else
  {
    std::string mtreeDirectoryName = ".";
    ann = (ApproxNearestNeighbors*)(new ApproxNearestNeighborsMTree(dim,
                                                                    "kriging_model_database",
                                                                    mtreeDirectoryName,
                                                                    &(std::cout),
                                                                    false));
    //printf("MTree NNS created.\n");
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

void NearestNeighborSearch::insert(std::vector<double> point, uint128_t key, const CkCallback &cb)
{
  int id = ann->insert(point, key);

  int_message* msg = new int_message(id);

  cb.send(msg);
}

void NearestNeighborSearch::insertAsync(std::vector<double> point, uint128_t key)
{
  ann->insert(point, key);
}

void NearestNeighborSearch::insert(std::vector<uint128_t> keys)
{
  ann->insert(keys);
}

void NearestNeighborSearch::remove(int id)
{
  ann->remove(id);
}

void NearestNeighborSearch::knn(std::vector<double> x, int k, const CkCallback &cb)
{

  std::vector<int> ids;
  std::vector<uint128_t> keys;
  std::vector<double> dists;

  ann->knn(x, k, ids, keys, dists);

  knn_message* msg = new knn_message(k, ids, keys, dists);

  cb.send(msg);
}

void NearestNeighborSearch::getKey(int id, const CkCallback &cb)
{
  uint128_t rkey = ann->getKey(id);

  uint128_message *msg = new uint128_message(rkey);
  
  cb.send(msg);
}
