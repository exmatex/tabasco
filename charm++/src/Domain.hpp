#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "CoPPM.decl.h"

class Domain : public CBase_Domain {
  private:
       
  // Charm communication variables
  int numElemGhosts;
  int numNodeGhosts;
  int ghostNodeCount;
  int ghostElemCount;

  public:
  
  Domain();
  Domain(CkMigrateMessage *msg);
  ~Domain();
  void pup(PUP::er &p);

  // Entry methods
  void run();
};

#endif
