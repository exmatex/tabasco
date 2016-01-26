#ifndef _EVALUATE_H_
#define _EVALUATE_H_

#include "TabaSCo.decl.h"
#include "Taylor.h"
#include "tensor.h"

#include <vector>

class eval_message : public CMessage_eval_message {
  public:
    int size;
    double* value;

    eval_message() {}
    eval_message(std::vector<double>& val);
    eval_message(char* buf);
    ~eval_message();
    static void *pack(eval_message *);
    static eval_message *unpack(void *);
    void copyTo(std::vector<double>& val);
};

class evaln_message : public CMessage_evaln_message {
  public:
    int vsize;
    int dsize;
    double* value;
    double* deriv;

    evaln_message() {}
    evaln_message(int vsize, double* out_val, int dsize, double*  out_deriv);
    evaln_message(Tensor2Sym& out_val, Tensor4LSym& out_deriv);
    evaln_message(char* buf);
    ~evaln_message();
    static void *pack(evaln_message *);
    static evaln_message *unpack(void *);
    void copyTo(double* out_val, double* out_deriv);
    void copyTo(Tensor2Sym& out_val, Tensor4LSym& out_deriv);
};


class Evaluate : public CBase_Evaluate {
  private:
    Evaluate_SDAG_CODE;
    int evalType;     

  public:
    Plasticity *pm;

  Evaluate();
  Evaluate(CkMigrateMessage *msg);
  ~Evaluate();
  void pup(PUP::er &p);

  // Entry methods
  void initialize(int etype);
  void eval(std::vector<double> point, const CkCallback &cb);
  void evalNative(Tensor2Sym in, const CkCallback &cb);
};

#endif
