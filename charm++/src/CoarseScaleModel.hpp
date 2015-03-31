#ifndef _COARSE_SCALE_MODEL_H_
#define _COARSE_SCALE_MODEL_H_

#include "CoPPM.decl.h"

class CoarseScaleModel : public CBase_CoarseScaleModel {
  private:
    int maxTimesteps;
    int numElems;
     
  public:
  
  CoarseScaleModel();
  CoarseScaleModel(CkMigrateMessage *msg);
  ~CoarseScaleModel();
  void pup(PUP::er &p);

  // Entry methods
  void run(int ntimesteps, int nelems);
};

#endif
