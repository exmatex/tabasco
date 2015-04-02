/****************C-STANDARDS**************/
#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <cfloat>
#include <math.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <ctime>
/****************CPP-STANDARDS************/
#include <map>
#include <vector>
#include <list>
/****************CHARM++******************/
#include "Main.hpp"
#include "CoM4.decl.h"
#define printf CkPrintf
/****************C-STUFF******************/
#include "types.hpp"
/****************CPP-STUFF***************/
#include "CoM4.hpp"
#include "Domain.hpp"

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Domain domainArray;
extern /* readonly */ CProxy_CoarseScaleModel coarseScaleArray;
extern /* readonly */ CProxy_FineScaleModel fineScaleArray;
extern /* readonly */ CProxy_DBInterface DBArray;

void main_CoM4(Input in, CProxy_Domain domainArray)
{
  printf("In main_CoM4\n");

}
