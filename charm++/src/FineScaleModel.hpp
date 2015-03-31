#ifndef _FINE_SCALE_MODEL_H_
#define _FINE_SCALE_MODEL_H_

#include "CoPPM.decl.h"

class FineScaleModel : public CBase_FineScaleModel {
  private:
       
  public:
  
  FineScaleModel();
  FineScaleModel(CkMigrateMessage *msg);
  ~FineScaleModel();
  void pup(PUP::er &p);

  // Entry methods
  void run(int iter);
};

#endif
