#ifndef _MAIN_H_
#define _MAIN_H_

#include "TabaSCo.decl.h"
#include "types.h"

class Main : public CBase_Main {

  private:
    Main_SDAG_CODE
    int doneCount;
    int totalChares;

    bool isNodeKaswell(int physNodeID);
  public:

  /// Constructors ///
  Main(CkArgMsg* msg);
  Main(CkMigrateMessage* msg);

};
#endif //__MAIN_H__
