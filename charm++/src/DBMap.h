#ifndef _DB_MAP_H_
#define _DB_MAP_H_

#include "TabaSCo.decl.h"

#include <vector>

class DBMap : public CkArrayMap {
  private:

  int dbNodeCount;
  bool useNodeList;
  std::vector<int> nodeList;

  public:

  DBMap(int dbNodeCount);
  DBMap(std::vector<int> nodeList);

  inline int procNum(int, const CkArrayIndex &iIndex);
};

#endif
