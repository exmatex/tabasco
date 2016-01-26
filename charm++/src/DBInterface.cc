#include "TabaSCo.decl.h"
#include "DBInterface.h"

extern CProxy_Main mainProxy;
extern CProxy_NearestNeighborSearch nnsArray;

DBInterface::DBInterface(int backType)
: dbRef(SingletonDB::getInstance(static_cast<SingletonDBBackendEnum>(backType))) ///TODO: Find a MUCH better way to do this
{
  CkPrintf("DBInterface created on PE %d Index %d\n",
      CkMyPe(), thisIndex);
  //Tell charm++ to not migrate this chare
  setMigratable(false);
}

DBInterface::DBInterface(CkMigrateMessage *msg)
: dbRef(SingletonDB::getInstance())
{

}

DBInterface::~DBInterface()
{

}

void DBInterface::pup(PUP::er &p)
{
    //This should never be called... Probably
    ///TODO: Verify this is only called if migrating, not context switching, and then uncomment error
    //CkPrintf("ERROR: DBInterface %d:%d Migrating\n", CkMyPe(), thisIndex.x);
}

void DBInterface::push(uint128_t key, std::vector<double> buf, unsigned long key_length)
{
    //Simple, just directly pass the arguments
    this->dbRef.push(key, buf, key_length);
}

DBVecMessage * DBInterface::pull(uint128_t key)
{
    std::vector<double> retVec = this->dbRef.pull(key);
    DBVecMessage * msg = new DBVecMessage(retVec);
    return msg;
}

DBVecMessage * DBInterface::pull_key(uint128_t key)
{
    std::vector<double> retVec = this->dbRef.pull_key(key);
    DBVecMessage * msg = new DBVecMessage(retVec);
    return msg;
}

void DBInterface::erase(uint128_t key)
{
    //Simple, just directly pass the arguments
    this->dbRef.erase(key);
}
