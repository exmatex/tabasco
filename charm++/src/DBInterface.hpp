#ifndef _DB_INTERFACE_H_
#define _DB_INTERFACE_H_

#include "CoM4.decl.h"

class DBInterface : public CBase_DBInterface {
  private:
       
  public:
  
  DBInterface();
  DBInterface(CkMigrateMessage *msg);
  ~DBInterface();
  void pup(PUP::er &p);

  // Entry methods
  void run();
};

#endif
