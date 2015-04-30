#include "CoM4.decl.h"
#include "CoarseScaleModel.hpp"
#include "DBInterface.hpp"
#include "FineScaleModel.hpp"

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
  p|maxTimesteps;
  p|numElems;
  p|nstep;
  p|tstep;
  p|e;
  p|currentPt;
}

void CoarseScaleModel::startElementFineScaleQuery(int step, int nelems)
{
  nstep = step;
  numElems = nelems;
  printf("CoarseScaleModel startElementFineScaleQuery\n");

  CkCallback cb(CkIndex_CoarseScaleModel::receiveNewPoint((Msg*)NULL), thisProxy);
  for (int j = 0; j < numElems; j++)
  {
    fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, j).query(j, nstep, currentPt, cb);
  }

}

void CoarseScaleModel::updateElement(int whichEl, int whichIter, int newPt)
{
  int currentPt = newPt;
  printf("Iter %d Element %d %d %d %d update newPt %d\n", whichIter, thisIndex.x, thisIndex.y, thisIndex.z, whichEl, newPt);
}
