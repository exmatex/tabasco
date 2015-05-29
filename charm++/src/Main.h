#ifndef _MAIN_H_
#define _MAIN_H_

#include "TabaSCo.decl.h"
#include "types.h"

class Main : public CBase_Main {

  private:
    int doneCount;
    int totalChares;

  public:

  /// Constructors ///
  Main(CkArgMsg* msg);
  Main(CkMigrateMessage* msg);

  /// Entry Methods ///
  void done();
};
#endif //__MAIN_H__
