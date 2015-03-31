#ifndef _MAIN_H_
#define _MAIN_H_

#include "CoPPM.decl.h"
#include "types.hpp"

class Main : public CBase_Main {

  private:
    int doneCount;
    int totalChares;

  public:

  /// Constructors ///
  Main(CkArgMsg* msg);
  Main(CkMigrateMessage* msg);

  /// Entry Methods ///
  void go(Input in);
  void done();
};
#endif //__MAIN_H__
