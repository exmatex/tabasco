#include "CoM4.decl.h"
#include "CoarseScaleModel.h"
#include "DBInterface.h"
#include "FineScaleModel.h"

extern CProxy_Main mainProxy;
extern CProxy_FineScaleModel fineScaleArray;

void *
Msg::pack(Msg* m)
{
  int *p = (int *) CkAllocBuffer(m, 3*sizeof(int));
  int *t = p;

  *t = m->whichEl; t++;
  *t = m->whichIter; t++;
  *t = m->newPt; t++;
  CkFreeMsg(m);
  return(p);
}

Msg *
Msg::unpack(void *buf)
{
   int *in = (int *) buf;
   Msg *t = new (CkAllocBuffer(in, sizeof(Msg))) Msg;
   t->whichEl = in[0];
   t->whichIter = in[1];
   t->newPt = in[2];
   CkFreeMsg(buf);
   return t;
}

CoarseScaleModel::CoarseScaleModel()
{
  
  printf("CoarseScaleModel created on PE %d Index %d %d %d\n", 
      CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.z);
}

CoarseScaleModel::CoarseScaleModel(CkMigrateMessage *msg)
{

}

CoarseScaleModel::~CoarseScaleModel()
{

}

void CoarseScaleModel::pup(PUP::er &p)
{
  CBase_CoarseScaleModel::pup(p);
  p|numElemGhosts;
  p|numNodeGhosts;
  p|ghostNodeCount;
  p|ghostElemCount;
  p|maxTimesteps;
  p|numElems;
  p|nstep;
  p|tstep;
  p|e;
  p|currentPt;
  p|count;
}

void CoarseScaleModel::startElementFineScaleQuery(int step, int nelems)
{
  nstep = step;

  //CkCallback cb(CkIndex_CoarseScaleModel::receiveNewPoint((Msg*)NULL), thisProxy);
  for (int j = 0; j < nelems; j++)
  {
    //fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, j).query(j, nstep, currentPt, cb);
    fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, j).query(j, nstep, currentPt);
  }

}

void CoarseScaleModel::updateElement(int whichEl, int whichIter, int newPt)
{
  int currentPt = newPt;
  printf("Iter %d Coarse %d %d %d Element %d update newPt %d\n", whichIter, thisIndex.x, thisIndex.y, thisIndex.z, whichEl, newPt);
}

void CoarseScaleModel::haloExchange()
{
  printf("halo exchange\n");
}
