#include "CoM4.decl.h"
#include "Interpolate.hpp"

extern CProxy_Main mainProxy;
extern CProxy_FineScaleModel fineScaleArray;

Interpolate::Interpolate()
{
  
  printf("Interpolate created on PE %d Index %d %d %d\n", 
      CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.z);
}

Interpolate::Interpolate(CkMigrateMessage *msg)
{

}

Interpolate::~Interpolate()
{

}

void Interpolate::pup(PUP::er &p)
{
  CBase_Interpolate::pup(p);
  p|newPt;
}

void Interpolate::run(int elnum, int nbrCount, int nbrData, int qPt)
{
  printf("Interpolate run\n");

  newPt = thisIndex.z;

  int converged = thisIndex.x;

  // Send back new point, if converged
  // else call evaluate
  if (converged == 1)
  {
    newPt = thisIndex.z;
    fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, elnum).receivePoint(newPt);
  }
  else
    fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, elnum).evaluate(qPt);
}
