#ifndef _INTERPOLATE_H_
#define _INTERPOLATE_H_

#include "CoM4.decl.h"

class Interpolate : public CBase_Interpolate {
  private:
       
  public:
  
  Interpolate();
  Interpolate(CkMigrateMessage *msg);
  ~Interpolate();
  void pup(PUP::er &p);

  // Entry methods
  void run();
};

#endif
