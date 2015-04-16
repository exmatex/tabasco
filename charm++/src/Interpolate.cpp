#include "CoM4.decl.h"
#include "Interpolate.hpp"

extern CProxy_Main mainProxy;

Interpolate::Interpolate()
{
  
  printf("Interpolate created on PE %d Index %d %d %d\n", 
      CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.z);
}

Interpolate::Interpolate(CkMigrateMessage *msg)
{

}

Interpolate::~Interpolate()
{

}

void Interpolate::pup(PUP::er &p)
{

}

void Interpolate::run(int iter)
{
  printf("Interpolate running iter %d\n", iter);
}
