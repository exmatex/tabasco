#include "TabaSCo.decl.h"
#include "DBInterface.h"

extern CProxy_Main mainProxy;
extern CProxy_NearestNeighborSearch nnsArray;

DBInterface::DBInterface()
{
  
  printf("DBInterface created on PE %d Index %d %d %d\n", 
      CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.z);
}

DBInterface::DBInterface(CkMigrateMessage *msg)
{

}

DBInterface::~DBInterface()
{

}

void DBInterface::pup(PUP::er &p)
{

}

void DBInterface::get(int pCount, int pIndex, int pData)
{
  printf("DBInterface get\n");

  // Call to DB

  nnsArray(thisIndex.x, thisIndex.y, thisIndex.z).receiveData(pCount, pData);
}

void DBInterface::put(int pCount, int pData)
{
  printf("DBInterface put\n");
}
