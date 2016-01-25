#include "TabaSCo.decl.h"
#include "FineScaleModel.h"
#include "Evaluate.h"
#include "Taylor.h"

extern CProxy_Main mainProxy;
extern CProxy_FineScaleModel fineScaleArray;

eval_message::eval_message(std::vector<double>& val) : size(val.size())
{
  value = new double[size];

  for (int i = 0; i < size; i++) 
  {
    value[i] = val[i];
  }
}

eval_message::eval_message(char* buf)
{
  char* p = buf;

  memcpy(&size, p, sizeof(int));
  p += sizeof(int);

  value = new double[size];

  memcpy(value, p, size*sizeof(double));
}

void *
eval_message::pack(eval_message* m) 
{
  int msize = sizeof(int) + m->size*sizeof(double);
  char *p = (char*)CkAllocBuffer(m, msize);
  char *buf = p;

  memcpy(buf, &m->size, sizeof(int));
  buf += sizeof(int);
  memcpy(buf, m->value, m->size * sizeof(double));

  CkFreeMsg(m);
  return (void*) p;
}

eval_message *
eval_message::unpack(void *inbuf)
{
   char *buf = (char *) inbuf;

   eval_message* emsg = (eval_message*)CkAllocBuffer(inbuf, sizeof(eval_message));
   emsg = new ((void*)emsg) eval_message(buf);

   CkFreeMsg(buf);
   return emsg;
}

void
eval_message::copyTo(std::vector<double>& val)
{
  for (int i = 0; i < size; i++)
  {
    val.push_back(value[i]);
  }
}

eval_message::~eval_message()
{
  delete [] value;
}

evaln_message::evaln_message(int vsize, double* out_val, int dsize, double* out_deriv) : vsize(vsize), dsize(dsize)
{
  value = new double[vsize];
  deriv = new double[dsize];

  memcpy(value, out_val, vsize*sizeof(double));
  memcpy(deriv, out_deriv, dsize*sizeof(double));
}

evaln_message::evaln_message(Tensor2Sym& out_val, Tensor4LSym& out_deriv) : vsize(6), dsize(36)
{
  value = new double[vsize];
  deriv = new double[dsize];

  memcpy(value, out_val.a, vsize*sizeof(double));
  memcpy(deriv, out_deriv.a, dsize*sizeof(double));
}

evaln_message::evaln_message(char* buf)
{
  char* p = buf;

  memcpy(&vsize, p, sizeof(int));
  p += sizeof(int);
  memcpy(&dsize, p, sizeof(int));
  p += sizeof(int);

  value = new double[vsize];
  deriv = new double[dsize];

  memcpy(value, p, vsize*sizeof(double));
  p += vsize*sizeof(double);
  memcpy(deriv, p, dsize*sizeof(double));
}

void *
evaln_message::pack(evaln_message* m)
{
  int msize = 2 * sizeof(int) + m->vsize*sizeof(double) + m->dsize*sizeof(double);
  char *p = (char*)CkAllocBuffer(m, msize);
  char *buf = p;

  memcpy(buf, &m->vsize, sizeof(int));
  buf += sizeof(int);
  memcpy(buf, &m->dsize, sizeof(int));
  buf += sizeof(int);
  memcpy(buf, m->value, m->vsize * sizeof(double));
  buf += m->vsize * sizeof(double);
  memcpy(buf, m->deriv, m->dsize * sizeof(double));

  CkFreeMsg(m);
  return (void*) p;
}

evaln_message *
evaln_message::unpack(void *inbuf)
{
   char *buf = (char *) inbuf;

   evaln_message* emsg = (evaln_message*)CkAllocBuffer(inbuf, sizeof(evaln_message));
   emsg = new ((void*)emsg) evaln_message(buf);

   CkFreeMsg(buf);
   return emsg;
}

void
evaln_message::copyTo(double* out_val, double* out_deriv)
{
  memcpy(out_val, value, vsize*sizeof(double));
  memcpy(out_deriv, deriv, dsize*sizeof(double));
}

void
evaln_message::copyTo(Tensor2Sym& out_val, Tensor4LSym& out_deriv)
{
  memcpy(out_val.a, value, vsize*sizeof(double));
  memcpy(out_deriv.a, deriv, dsize*sizeof(double));
}

evaln_message::~evaln_message()
{
  delete [] value;
  delete [] deriv;
}


Evaluate::Evaluate()
{
  pm = NULL;

}

Evaluate::Evaluate(CkMigrateMessage *msg)
{

}

Evaluate::~Evaluate()
{

}

void Evaluate::pup(PUP::er &p)
{
  CBase_Evaluate::pup(p);
}

void Evaluate::initialize(int etype)
{
  evalType = etype;

  // Taylor plasticity model
  if (evalType == 0) {
      double D_0 = 1.e-2;
      double m = 1./20.;
      double g = 2.e-3; // (Mbar)

      if (pm == NULL) {
        pm = (Plasticity*)(new Taylor(D_0, m, g));
      }
  }

  // VPSC plasticity model
  else {

  }

}

void Evaluate::eval(std::vector<double> point, const CkCallback &cb)
{
  // Taylor plasticity model
  if (evalType == 0) {
    std::vector<double> value(pm->valueAndDerivativeDimension());

    pm->evaluate(point, value);

    eval_message* msg = new eval_message(value);

    cb.send(msg);
  }
  // VPSC plasticity model
  else {

  }
}

void Evaluate::evalNative(Tensor2Sym in, const CkCallback &cb)
{
  // Taylor plasticity model
  if (evalType == 0) {
    Tensor2Sym out_value;
    Tensor4LSym out_derivative;

    pm->evaluateNative(in, out_value, out_derivative);
    evaln_message* msg = new evaln_message(out_value, out_derivative);

    cb.send(msg);
  }
  else {

  }
}
