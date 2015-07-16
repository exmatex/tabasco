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

    int maxTimesteps;
    int numElems;
    int nstep;
    int tstep;
    int e;
    int currentPt;
    int count;
     
  public:
  
  CoarseScaleModel();
  CoarseScaleModel(CkMigrateMessage *msg);
  ~CoarseScaleModel();
  void pup(PUP::er &p);

  // Entry methods
  void initialize(int numRanks, bool useAdaptiveSampling);
  void LagrangeNodal1();
  void LagrangeNodal2();
  void LagrangeElements();
  void CalcTimeConstraintsForElems();
  void TimeIncrement();
  void updateTimeIncrement(Real_t reducedt);
  void TimeIncrement2();
  void UpdateStressForElems();

  void sendData(int xferFields, Real_t **fieldData);
  void receiveData(int msgType, int size, int xferFields, Real_t **fieldData, Real_t rdata[]);

  void sendNodalMass();
  void updateNodalMass(int msgType, int rsize, Real_t rdata[]);
  //void receiveNodalMass(int msgType, int size, Real_t rdata[]);

  void sendForce();
  void updateForce(int msgType, int rsize, Real_t rdata[]);
  //void receiveForce(int msgType, int size, Real_t rdata[]);
  
  void sendPositionVelocity();
  void updatePositionVelocity(int msgType, int rsize, Real_t rdata[]);
  //void receivePositionVelocity(int msgType, int size, Real_t rdata[]);

  void startElementFineScaleQuery(int step, int nelems);
  void updateElement(int whichEl, int whichIter, int newPt);
  void go(int numRanks, bool useAdaptiveSampling);
};

#endif
