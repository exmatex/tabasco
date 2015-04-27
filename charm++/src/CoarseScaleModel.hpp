#ifndef _COARSE_SCALE_MODEL_H_
#define _COARSE_SCALE_MODEL_H_

#include "CoM4.decl.h"

class CoarseScaleModel : public CBase_CoarseScaleModel {
  private:
    CoarseScaleModel_SDAG_CODE
    int maxTimesteps;
    int numElems;
    int nstep;
    int tstep;
    int e;
    int currentPt;
     
  public:
  
  CoarseScaleModel();
  CoarseScaleModel(CkMigrateMessage *msg);
  ~CoarseScaleModel();
  void pup(PUP::er &p);

  // Entry methods
  void startElementFineScaleQuery(int step, int nelems);
  void updateElement(int whichEl, int whichIter, int newPt);
};

#endif
