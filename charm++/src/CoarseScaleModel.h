#ifndef _COARSE_SCALE_MODEL_H_
#define _COARSE_SCALE_MODEL_H_

#include "TabaSCo.decl.h"

#include "LULESH/lulesh.h"

class Msg : public CMessage_Msg
{
  public:

    int whichEl;
    int whichIter;
    int newPt;

    static void *pack(Msg *);
    static Msg *unpack(void *);
};

// Reduction Operations
void registerMinReal(void);
CkReductionMsg *minReal(int nMsg, CkReductionMsg **msgs);

class CoarseScaleModel : public CBase_CoarseScaleModel {
  private:
    CoarseScaleModel_SDAG_CODE

    Lulesh *lulesh;

    CkCallback *cbTime;
    CkCallback *cbIters;
    CkCallback *cbTotal;
    CkCallback *cbFdone;

    int maxTimesteps;
    int numElems;
    int nstep;
    int tstep;
    int e;
    int currentPt;
    int count;

    int total_samples;
    int total_interpolations;

    int max_nonlinear_iters;
    int max_local_newton_iters;

    int use_adaptive_sampling;

    size_t* state_size;

    int file_parts;
    int visit_data_interval;
     
  public:
  
  CoarseScaleModel();
  CoarseScaleModel(CkMigrateMessage *msg);
  ~CoarseScaleModel();
  void pup(PUP::er &p);
  int* calcRange(int chareCount);

  // Entry methods
  void initialize(int numRanks, bool useAdaptiveSampling, Real_t stopTime);
  void ConstructFineScaleModel(bool useAdaptiveSampling);
  void LagrangeNodal1();
  void LagrangeNodal2();
  void LagrangeElements();
  void LagrangeElements2();
  void CalcTimeConstraintsForElems();
  void TimeIncrement();
  void updateTimeIncrement(Real_t reducedt);
  void TimeIncrement2();
  void UpdateStressForElems();
  void updateAdvanceResults(int elemNum, ConstitutiveData cm_data, int ssize, char state[], int num_samples, int num_interpolations);
  void afterAdvance();
  void UpdateStressForElems2(int reducedIters);

  void sendDataNodes(int xferFields, Real_t **fieldData);
  void sendDataElems(int xferFields, Real_t **fieldData);

  void receiveDataNodes(int msgType, int size, int xferFields, Real_t **fieldData, Real_t rdata[]);
  void receiveDataElems(int msgType, int size, int xferFields, Real_t rdata[]);

  void sendNodalMass();
  void updateNodalMass(int msgType, int rsize, Real_t rdata[]);

  void sendForce();
  void updateForce(int msgType, int rsize, Real_t rdata[]);
  
  void sendVelocityGrad();
  void updateVelocityGrad(int msgType, int rsize, Real_t rdata[]);

  void sendPositionVelocity();
  void updatePositionVelocity(int msgType, int rsize, Real_t rdata[]);

  void makeADump(bool forceDump);
  void setSiloParams(int numParts, int dataInterval);
};

#endif
