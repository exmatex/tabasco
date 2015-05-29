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
  void get(int pCount, int pIndex, int pData);
  void put(int pCount, int pData);
};

#endif
