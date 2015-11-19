#include "TabaSCo.decl.h"
#include "NearestNeighborSearch.h"
#include "DBInterface.h"

extern CProxy_Main mainProxy;
extern CProxy_DBInterface DBArray;
extern CProxy_FineScaleModel fineScaleArray;
extern CProxy_NearestNeighborSearch nnsArray;
extern int NBR_LIMIT;

knn_message::knn_message(int k, std::vector<int>& vids, std::vector<uint128_t>& vkeys, std::vector<double>& vdists) : k(k)
{
  ids = new int[k];
  keys = new uint128_t[k];
  dists = new double[k];

  for (int i =0; i < k; i++)
  {
    ids[i] =  vids[i];
    keys[i] = vkeys[i];
    dists[i] = vdists[i];
  }
}

knn_message::knn_message(char* buf, int k) : k(k)
{
  ids = new int[k];
  keys = new uint128_t[k];
  dists = new double[k];

  char* p = buf;

  memcpy(ids, p, k*sizeof(int));
  p += k*sizeof(int);
  memcpy(keys, p, k*sizeof(uint128_t));
  p += k*sizeof(uint128_t);
  memcpy(dists, p, k*sizeof(double));
}

void *
knn_message::pack(knn_message* m)
{
  int msize = sizeof(int)+m->k*sizeof(int)+m->k*sizeof(double)+m->k*sizeof(uint128_t);
  char *p = (char*)CkAllocBuffer(m, msize);
  char *buf = p;

  memcpy(buf, &m->k, sizeof(int));
  buf += sizeof(int);
  memcpy(buf, m->ids, m->k * sizeof(int));
  buf += m->k * sizeof(int);
  memcpy(buf, m->keys, m->k * sizeof(uint128_t));
  buf += m->k * sizeof(uint128_t);
  memcpy(buf, m->dists, m->k * sizeof(double));

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

NearestNeighborSearch::NearestNeighborSearch()
{
  printf("NearestNeighborSearch created on PE %d Index %d\n",
      CkMyPe(), thisIndex);

}

void NearestNeighborSearch::initialize(int ntype, int dim, int ntrees)
{
  if (ntype == 0) 
  {
#ifdef FLANN
    int nchecks = 20;
    ann = (ApproxNearestNeighbors*)(new ApproxNearestNeighborsFLANN(dim, ntrees, nchecks));
    printf("FLANN NNS creatd.\n");
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
    printf("MTree NNS created.\n");
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

int_message* NearestNeighborSearch::insert(std::vector<double> point, uint128_t key)
{
  int id = ann->insert(point, key);

  int_message* msg = new int_message(id);
  return msg;
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

knn_message* NearestNeighborSearch::knn(std::vector<double> x, int k)
//NearestNeighborSearch::knn(std::vector<double> x, int k, const CkCallback &cb)
{
  std::vector<int> ids(k);
  std::vector<uint128_t> keys(k);
  std::vector<double> dists(k);

  ann->knn(x, k, ids, keys, dists);

  knn_message* msg = new knn_message(k, ids, keys, dists);
  return msg;

  //cb.send((void*&)msg);
}

uint128_message* NearestNeighborSearch::getKey(int id)
{
  uint128_t rkey = ann->getKey(id);

  uint128_message *msg = new uint128_message(rkey);
  return msg;
}

