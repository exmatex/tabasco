#include "CoM4.decl.h"
#include "Domain.hpp"

extern CProxy_Main mainProxy;

Domain::Domain()
{
  
  printf("domain created on PE %d Index %d %d %d\n", CkMyPe(), thisIndex.x, thisIndex.
         y, thisIndex.z);
}

Domain::Domain(CkMigrateMessage *msg)
{

}

Domain::~Domain()
{

}

void Domain::pup(PUP::er &p)
{

}

void Domain::run()
{

}
