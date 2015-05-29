#include "TabaSCo.decl.h"
#include "Domain.h"

extern CProxy_Main mainProxy;

Domain::Domain()
{
  
  printf("Domain created on PE %d Index %d %d %d\n", 
         CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.z);
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

void Domain::haloExchange()
{

}

