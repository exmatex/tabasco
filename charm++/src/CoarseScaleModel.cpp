#include "CoM4.decl.h"
#include "CoarseScaleModel.hpp"
#include "DBInterface.hpp"
#include "FineScaleModel.hpp"

extern CProxy_Main mainProxy;
extern CProxy_FineScaleModel fineScaleArray;

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

  for (int j = 0; j < numElems; j++)
  {
    fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, j).query(nstep);
  }

  // This chare is done
  //printf("CoarseScaleModel is done\n");
  //mainProxy.done();
}

void CoarseScaleModel::updateElement(int whichEl, int whichIter, int newPt)
{
  int currentPt = newPt;
  printf("Iter %d Element %d %d %d %d update newPt %d\n", whichIter, thisIndex.x, thisIndex.y, thisIndex.z, whichEl, newPt);
}
