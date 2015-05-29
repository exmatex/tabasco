#ifndef _FINE_SCALE_MODEL_H_
#define _FINE_SCALE_MODEL_H_

#include "TabaSCo.decl.h"

class FineScaleModel : public CBase_FineScaleModel {
  private:
    FineScaleModel_SDAG_CODE
    int newPt;
    int currentIter;
    int nbrCount;
    int nbrData;
       
  public:
  
  FineScaleModel();
  FineScaleModel(CkMigrateMessage *msg);
  ~FineScaleModel();
  void pup(PUP::er &p);

  // Entry methods
  void evaluate(int qPt);
  void query2(int iter, int qPt);
  void requestNeighbors(int qPt);
  void requestInterpolation(int nbrCount, int nbrData, int qPt);
  void requestDBStore(int cPt);
  void sendNewPoint2Coarse(int elnum, int iter, int cPt);

};

#endif
