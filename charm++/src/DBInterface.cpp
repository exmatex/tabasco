#include "CoPPM.decl.h"
#include "DBInterface.hpp"

extern CProxy_Main mainProxy;

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

void DBInterface::run()
{
  printf("DBInterface running\n");
}
