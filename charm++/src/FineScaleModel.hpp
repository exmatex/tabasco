#ifndef _FINE_SCALE_MODEL_H_
#define _FINE_SCALE_MODEL_H_

#include "CoM4.decl.h"

class FineScaleModel : public CBase_FineScaleModel {
  private:
    FineScaleModel_SDAG_CODE
    int newPt;
    int currentIter;
       
  public:
  
  FineScaleModel();
  FineScaleModel(CkMigrateMessage *msg);
  ~FineScaleModel();
  void pup(PUP::er &p);

  // Entry methods
  void evaluate();
  void query2(int iter);
};

#endif
