#ifndef _DB_INTERFACE_H_
#define _DB_INTERFACE_H_

#include "TabaSCo.decl.h"

#include "DBVecMessage.h"
#include "SingletonDB.h"


class DBInterface : public CBase_DBInterface {
  private:
      SingletonDB & dbRef;
  public:

  DBInterface(int backType=0);
  DBInterface(CkMigrateMessage *msg);
  ~DBInterface();
  void pup(PUP::er &p);

  // Entry methods
  ///TODO: is vector safe to pass? Charm++ can auto-pup it, so probably?
  void push(uint128_t key, std::vector<double> buf, unsigned long key_length);
  DBVecMessage * pull(uint128_t key);
  DBVecMessage * pull_key(uint128_t key);
  void erase(uint128_t key);

};

#endif
