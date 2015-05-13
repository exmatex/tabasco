#ifndef _NEARESTNEIGHBORSEARCH_H_
#define _NEARESTNEIGHBORSEARCH_H_

#include "CoM4.decl.h"

class NearestNeighborSearch : public CBase_NearestNeighborSearch {
  private:
    NearestNeighborSearch_SDAG_CODE;
    int nbrCount;
    int nbrIndex;   
    int nbrData;   

  public:
  
  NearestNeighborSearch();
  NearestNeighborSearch(CkMigrateMessage *msg);
  ~NearestNeighborSearch();
  void pup(PUP::er &p);

  // Entry methods
  void getIndex(int qPt, int nbrCount, int nbrIndex);
  void putIndex(int qPt, int nbrIndex);
  void requestDBGet(int nbrCount, int nbrIndex, int nbrData);
  void sendNeighbors(int elnum, int nbrCount, int nbrData);
};

#endif
