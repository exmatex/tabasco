#ifndef _FINE_SCALE_MODEL_H_
#define _FINE_SCALE_MODEL_H_

#include "TabaSCo.decl.h"

#include "ElastoViscoPlasticity.h"

class FineScaleModel : public CBase_FineScaleModel {
  private:
    FineScaleModel_SDAG_CODE
    int newPt;
    int currentIter;
    int nbrCount;
    int nbrData;
       
  public:
  
    Constitutive* cm; 

  FineScaleModel();
  FineScaleModel(bool use_adaptive_sampling);
  FineScaleModel(CkMigrateMessage *msg);
  ~FineScaleModel();
  void pup(PUP::er &p);

/*
  void evaluate(int qPt);
  void query2(int iter, int qPt);
  void requestNeighbors(int qPt);
  void requestInterpolation(int nbrCount, int nbrData, int qPt);
  void requestDBStore(int cPt);
  void sendNewPoint2Coarse(int elnum, int iter, int cPt);
*/

  // Entry methods
  void advance();

};

#endif
