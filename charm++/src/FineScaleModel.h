#ifndef _FINE_SCALE_MODEL_H_
#define _FINE_SCALE_MODEL_H_

#include "TabaSCo.decl.h"

#include "ElastoViscoPlasticity.h"
#include "tensor.h"

// PUP operator for Tensor2Sym
inline void operator|(PUP::er &p, Tensor2Sym &tensor)
{
  PUParray(p, tensor.a, 6);
}

// PUP operator for Tensor2Gen
inline void operator|(PUP::er &p, Tensor2Gen &tensor)
{
  PUParray(p, tensor.a, 9);
}

// PUP operator for ConstitutiveData
inline void operator|(PUP::er &p, ConstitutiveData &cdata)
{
  p|cdata.sigma_prime;
  p|cdata.num_models;
  p|cdata.num_point_value_pairs;
  p|cdata.num_Newton_iters;

}

class FineScaleModel : public CBase_FineScaleModel {
  private:
    FineScaleModel_SDAG_CODE
    int newPt;
    int currentIter;
    int nbrCount;
    int nbrData;

    size_t stateSize;
       
  public:
  
  Constitutive* cm; 

  FineScaleModel();
  FineScaleModel(int state_size, bool use_adaptive_sampling);
  //void initialize(int state_size, bool use_adaptive_sampling);
  FineScaleModel(CkMigrateMessage *msg);
  ~FineScaleModel();
  void pup(PUP::er &p);

/*
  void evaluate(int qPt);
  void query2(int iter, int qPt);
  void requestNeighbors(int qPt);
  void requestInterpolation(int nbrCount, int nbrData, int qPt);
  void requestDBStore(int cPt);
  void sendNewPoint2Coarse(int elnum, int iter, int cPt);
*/

  // Entry methods
  void advance(const double delta_t, const Tensor2Gen& L_new, const double volume_change, int ssize, char* state);

};

#endif
