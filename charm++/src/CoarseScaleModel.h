#ifndef _COARSE_SCALE_MODEL_H_
#define _COARSE_SCALE_MODEL_H_

#include "TabaSCo.decl.h"

#include "LULESH/lulesh.h"

class Msg : public CMessage_Msg
{
  public:

    int whichEl;
    int whichIter;
    int newPt;

    static void *pack(Msg *);
    static Msg *unpack(void *);
};

class CoarseScaleModel : public CBase_CoarseScaleModel {
  private:
    CoarseScaleModel_SDAG_CODE

    Lulesh *lulesh;

    int maxTimesteps;
    int numElems;
    int nstep;
    int tstep;
    int e;
    int currentPt;
    int count;
     
  public:
  
  CoarseScaleModel();
  CoarseScaleModel(CkMigrateMessage *msg);
  ~CoarseScaleModel();
  void pup(PUP::er &p);

  // Entry methods
  void startElementFineScaleQuery(int step, int nelems);
  void updateElement(int whichEl, int whichIter, int newPt);
  void go(int myRank, int numRanks);
};

#endif
