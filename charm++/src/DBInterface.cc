#include "TabaSCo.decl.h"
#include "DBInterface.h"
#include "SingletonDB.h"

#include "SingletonDB_HashMap.h"

///TODO: Refactor the dummy backend stuff
#ifdef REDIS
#include "SingletonDB_Redis.h"
#else
#include "SingletonDB_Dummy.h"
typedef SingletonDB_Dummy SingletonDB_Redis;
#endif // REDIS



extern CProxy_Main mainProxy;
extern CProxy_NearestNeighborSearch nnsArray;

DBInterface::DBInterface(int numDBNodes, int backType)
{
  //Get the backType as an enum for easy use
  SingletonDBBackendEnum eBackType = static_cast<SingletonDBBackendEnum>(backType);
  bool isSpawner = thisIndex < numDBNodes;

  //Iterate over backend types
	if(eBackType == REDIS_DB)
	{
		this->dbRef = new SingletonDB_Redis(2, false, !isSpawner);
	}
	else if(backType == DIST_REDIS_DB)
	{
		this->dbRef = new SingletonDB_Redis(2, true, !isSpawner);
	}
	else if(backType == HASHMAP_DB)
	{
		///TODO: Consider a warning if isSpawner is false?
		this->dbRef = new SingletonDB_HashMap();
	}
	else
	{
		CkError("WARNING: Invalid DB Backend\n");
		CkExit();
	}
  CkPrintf("DBInterface created on PE %d Index %d\n",
      CkMyPe(), thisIndex);
  //Tell charm++ to not migrate this chare
  setMigratable(false);
}

DBInterface::DBInterface(CkMigrateMessage *msg)
{

}

DBInterface::~DBInterface()
{

}

void DBInterface::pup(PUP::er &p)
{
    //This should never be called... Probably
    ///TODO: Verify this is only called if migrating, not context switching, and then uncomment error
    //CkError("ERROR: DBInterface %d:%d Migrating\n", CkMyPe(), thisIndex.x);
}

void DBInterface::push(uint128_t key, std::vector<double> buf, unsigned long key_length)
{
    //Simple, just directly pass the arguments
    this->dbRef->push(key, buf, key_length);
}

DBVecMessage * DBInterface::pull(uint128_t key)
{
    std::vector<double> retVec = this->dbRef->pull(key);
    DBVecMessage * msg = new DBVecMessage(retVec);
    return msg;
}

DBVecMessage * DBInterface::pull_key(uint128_t key)
{
    std::vector<double> retVec = this->dbRef->pull_key(key);
    DBVecMessage * msg = new DBVecMessage(retVec);
    return msg;
}

void DBInterface::erase(uint128_t key)
{
    //Simple, just directly pass the arguments
    this->dbRef->erase(key);
}
