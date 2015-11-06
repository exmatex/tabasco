#include "TabaSCo.decl.h"
#include "FineScaleModel.h"
#include "NearestNeighborSearch.h"
#include "Interpolate.h"
#include "DBInterface.h"

#include "Constitutive.h"
#include "Taylor.h"
#include "MieGruneisen.h"

extern CProxy_Main mainProxy;
extern CProxy_CoarseScaleModel coarseScaleArray;
extern CProxy_NearestNeighborSearch nnsArray;
extern CProxy_Interpolate interpolateArray;
extern CProxy_DBInterface DBArray;

FineScaleModel::FineScaleModel()
{}

FineScaleModel::FineScaleModel(int state_size, bool use_adaptive_sampling)
//void FineScaleModel::initialize(int state_size, bool use_adaptive_sampling)
{
  // Ordering for a 2D array chare is x, y 
/*
  printf("FineScaleModel created on PE %d Index %d %d state size = %d\n", 
      CkMyPe(), thisIndex.x, thisIndex.y, state_size);
*/

  stateSize = state_size;

  ConstitutiveGlobal cm_global;

  // Construct the fine-scale plasticity model
  double D_0 = 1.e-2;
  double m = 1./20.;
  double g = 2.e-3; // (Mbar)
  Plasticity* plasticity_model = (Plasticity*)(new Taylor(D_0, m, g));

  // Construct the equation of state
  EOS* eos_model;
  /* 
   From Table 1 (converted from GPa to Mbar) in P. J. Maudlin et al.,
   "On the modeling of the Taylor cylinder impact test for orthotropic
   textured materials: experiments and simulations", Inter. J.
   Plasticity 15 (1999), pp. 139-166.
  */
  double k1 = 1.968;  // Mbar
  double k2 = 2.598;  // Mbar
  double k3 = 2.566;  // Mbar
  double Gamma = 1.60;  // dimensionless
  eos_model = (EOS*)(new MieGruneisen(k1, k2, k3, Gamma));

  // Construct the constitutive model
  double bulk_modulus = 1.94; // Tantallum (Mbar)
  double shear_modulus = 6.9e-1; // Tantallum (Mbar)
 
  // Make empty Tensor
  Tensor2Gen L;

  cm = (Constitutive*)(new ElastoViscoPlasticity(cm_global, L, bulk_modulus, shear_modulus, eos_model,
                                                 plasticity_model, use_adaptive_sampling, stateSize));
}

FineScaleModel::FineScaleModel(CkMigrateMessage *msg)
{
  printf("Migrate Fine Scale chare %d %d\n", thisIndex.x, thisIndex.y);
}

FineScaleModel::~FineScaleModel()
{

}

void FineScaleModel::pup(PUP::er &p)
{
  CBase_FineScaleModel::pup(p);
  p|newPt;
  p|currentIter;
  p|nbrCount;
  p|nbrData;
  p|stateSize;
}

/*
void FineScaleModel::query2(int iter, int qPt)
{
  currentIter = iter;

  printf("FineScaleModel %d %d %d %d query iter %d\n",
    thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, iter);

  // Check for neighbors
  nnsArray(thisIndex.w, thisIndex.x, thisIndex.y).getNeighbors(thisIndex.z, qPt);

  // There are enough neighbors, do interpolation
  interpolateArray(thisIndex.w, thisIndex.x, thisIndex.y).run(thisIndex.z, nbrCount, nbrData, qPt);

  // There are not enough neighbors, evaluate
  // Or interpolation did not converge
  evaluate(qPt);

  newPt = thisIndex.z;

  // Save new point to DB
  DBArray(thisIndex.w, thisIndex.x, thisIndex.y).put(1, newPt);

  // Send back new point to Coarse Scale model
  //coarseScaleArray(thisIndex.w, thisIndex.x, thisIndex.y).receiveNewPoint(thisIndex.z, currentIter, newPt);
}

void FineScaleModel::evaluate(int qPt)
{
  printf("FineScaleModel evaluate\n");
 
  newPt = thisIndex.z;

  receivePoint(newPt);
}

void FineScaleModel::requestNeighbors(int qPt)
{
  printf("FineScaleModel requestNeighbors\n");
 
  nnsArray(thisIndex.w, thisIndex.x, thisIndex.y).getNeighbors(thisIndex.z, qPt);
}

void FineScaleModel::requestInterpolation(int nbrCount, int nbrData, int qPt)
{
  printf("FineScaleModel requestInterpolation\n");

    interpolateArray(thisIndex.w, thisIndex.x, thisIndex.y).run(thisIndex.z, nbrCount, nbrData, qPt);
}

void FineScaleModel::requestDBStore(int cPt)
{
  printf("FineScaleModel requestDBStore\n");
 
  newPt = cPt;
  DBArray(thisIndex.w, thisIndex.x, thisIndex.y).put(1, newPt);
}

void FineScaleModel::sendNewPoint2Coarse(int elnum, int iter, int cPt)
{
  printf("FineScaleModel sendNewPoint\n");

  newPt = cPt;
  //coarseScaleArray(thisIndex.w, thisIndex.x, thisIndex.y).receiveNewPoint(elnum, iter, newPt);
}
*/

// Execute advance for element and return results to calling coarse model
void FineScaleModel::advance(const double delta_t, const Tensor2Gen& L_new, const double volume_change, int ssize, char* state)
{
  ConstitutiveData cm_data = cm->advance(delta_t, L_new, volume_change, state);
  int num_samples = 0;
  int num_interpolations = 0;
  if (cm->adaptiveSamplingEnabled())
  {
    num_samples = cm->getNumSamples();
    num_interpolations = cm->getNumSuccessfulInterpolations();
  }
  coarseScaleArray(thisIndex.x).receiveAdvanceResults(thisIndex.y, cm_data, ssize, state, num_samples, num_interpolations);
}
