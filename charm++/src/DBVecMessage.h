#ifndef _DB_VEC_MESSAGE_H_
#define _DB_VEC_MESSAGE_H_

#include "TabaSCo.decl.h"

#include <vector>

class DBVecMessage : public CMessage_DBVecMessage {
  private:
  public:
	int length;
	double * data;
	DBVecMessage(std::vector<double> vec);
	DBVecMessage(double * vec, int dim);
	static void *pack(DBVecMessage *);
	static DBVecMessage *unpack(void *);
};

#endif
