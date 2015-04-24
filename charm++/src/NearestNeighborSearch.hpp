#ifndef _NEARESTNEIGHBORSEARCH_H_
#define _NEARESTNEIGHBORSEARCH_H_

#include "CoM4.decl.h"

class NearestNeighborSearch : public CBase_NearestNeighborSearch {
  private:
       
  public:
  
  NearestNeighborSearch();
  NearestNeighborSearch(CkMigrateMessage *msg);
  ~NearestNeighborSearch();
  void pup(PUP::er &p);

  // Entry methods
  void run(int iter);
};

#endif
