#include "TabaSCo.decl.h"
#include "CoarseScaleModel.h"
#include "DBInterface.h"
#include "FineScaleModel.h"

#include "Constitutive.h"
#include "tensor.h"

#ifdef SILO
#include "LULESH/siloDump.h"
#endif

extern CProxy_Main mainProxy;
extern CProxy_CoarseScaleModel coarseScaleArray;
extern CProxy_FineScaleModel fineScaleArray;

extern int coarseCount;
extern int nnsType;
extern int nnsCount;
extern int interpType;
extern int interpCount;
extern int evalType;
extern int evalCount;
extern int dbType;
extern int dbCount;

void *
Msg::pack(Msg* m)
{
  int *p = (int *) CkAllocBuffer(m, 3*sizeof(int));
  int *t = p;

  *t = m->whichEl; t++;
  *t = m->whichIter; t++;
  *t = m->newPt; t++;
  CkFreeMsg(m);
  return(p);
}

Msg *
Msg::unpack(void *buf)
{
   int *in = (int *) buf;
   Msg *t = new (CkAllocBuffer(in, sizeof(Msg))) Msg;
   t->whichEl = in[0];
   t->whichIter = in[1];
   t->newPt = in[2];
   CkFreeMsg(buf);
   return t;
}

void CoarseScaleModel::updateTimeIncrement(Real_t reducedt)
{

  //printf("%d In timeReduce %e\n", thisIndex, reducedt);

  Real_t newdt = reducedt;

  Real_t olddt = lulesh->domain.deltatime() ;

  Real_t ratio = newdt / olddt ;
  if (ratio >= Real_t(1.0)) {
    if (ratio < lulesh->domain.deltatimemultlb()) {
      newdt = olddt ;
    }
    else if (ratio > lulesh->domain.deltatimemultub()) {
      newdt = olddt*lulesh->domain.deltatimemultub() ;
    }
  }

  if (newdt > lulesh->domain.dtmax()) {
    newdt = lulesh->domain.dtmax() ;
  }
  lulesh->domain.deltatime() = newdt ;

}

CoarseScaleModel::CoarseScaleModel()
{
  
  printf("CoarseScaleModel created on PE %d Index %d\n", 
      CkMyPe(), thisIndex);

  lulesh = new Lulesh();

  // Set time reduction and iters reduction
  cbTime = new CkCallback(CkReductionTarget(CoarseScaleModel, reduceTimeIncrement), thisProxy);
  cbIters = new CkCallback(CkReductionTarget(CoarseScaleModel, reduceIters), thisProxy);
  cbTotal = new CkCallback(CkReductionTarget(CoarseScaleModel, reduceTotal), thisProxy);

  total_samples = 0;
  total_interpolations = 0;

}

CoarseScaleModel::CoarseScaleModel(CkMigrateMessage *msg)
{

}

CoarseScaleModel::~CoarseScaleModel()
{
 free(state_size);
 free(lulesh);
}

void CoarseScaleModel::pup(PUP::er &p)
{
  CBase_CoarseScaleModel::pup(p);
  p|maxTimesteps;
  p|numElems;
  p|nstep;
  p|tstep;
  p|e;
  p|currentPt;
  p|total_samples;
  p|total_interpolations;
  p|count;
  p|max_nonlinear_iters;
  p|max_local_newton_iters;
  p|use_adaptive_sampling;
}

void CoarseScaleModel::initialize(int numRanks, int nnsCount, int interpCount, int evalCount, int dbCount, bool useAdaptiveSampling, int edgeElems, int heightElems, Real_t stopTime)
{
  lulesh->Initialize(thisIndex, numRanks, edgeElems, heightElems, stopTime);

  lulesh->domain.stoptime() = stopTime;

  numElems = lulesh->domain.numElem();
  //printf("numElems = %d\n", numElems);

  state_size = new size_t[numElems];
  use_adaptive_sampling = (useAdaptiveSampling) ? 1 : 0;
  ConstructFineScaleModel(numRanks, nnsCount, interpCount, evalCount, dbCount, useAdaptiveSampling);

}

void CoarseScaleModel::ConstructFineScaleModel(int numRanks, int nnsCount, int interpCount, int evalCount, int dbCount, bool useAdaptiveSampling)
{
  ModelDatabase * global_modelDB = nullptr;
  ApproxNearestNeighbors* global_ann = nullptr;
  int flanning = 0;
  int flann_n_trees=1;
  int flann_n_checks=20; 
  int global_ns=0;
  // Create Lulesh-local parts of fine scale models
  lulesh->ConstructFineScaleModel(useAdaptiveSampling,global_modelDB,global_ann,flanning,flann_n_trees,flann_n_checks,global_ns);

  // Now create 2-D fine scale model chares
  numElems = lulesh->domain.numElem();

  int* nnsRange = calcRange(nnsCount, numRanks);
  int* interpRange = calcRange(interpCount, numRanks);
  int* dbRange = calcRange(dbCount, numRanks);
  int* evalRange = calcRange(evalCount, numRanks);

  printf("CoarseScaleModel %d: nns %d %d\n", thisIndex, nnsRange[0], nnsRange[1]);
  printf("CoarseScaleModel %d: interpolate %d %d\n", thisIndex, interpRange[0], interpRange[1]);
  printf("CoarseScaleModel %d: dbinterface %d %d\n", thisIndex, dbRange[0], dbRange[1]);
  printf("CoarseScaleModel %d: evaluate %d %d\n", thisIndex, evalRange[0], evalRange[1]);

  int nnsIndex = nnsRange[0];
  int interpIndex = interpRange[0];
  int dbIndex = dbRange[0];
  int evalIndex = evalRange[0];

  for (Index_t i = 0; i < numElems; ++i) {
    state_size[i] = lulesh->domain.cm(i)->getStateSize();

    fineScaleArray(thisIndex, i).insert(state_size[i], useAdaptiveSampling, nnsIndex, interpIndex, dbIndex, remoteDB, evalIndex);
    
    nnsIndex++;
    if (nnsIndex >= nnsRange[1]) nnsIndex = nnsRange[0];
    interpIndex++;
    if (interpIndex >= interpRange[1]) interpIndex = interpRange[0];
    dbIndex++;
    if (dbIndex >= dbRange[1]) dbIndex = dbRange[0];
    evalIndex++;
    if (evalIndex >= evalRange[1]) evalIndex = evalRange[0];
  }
  printf("%d FineScaleModels created count = %d\n", thisIndex, numElems);
}

int* CoarseScaleModel::calcRange(int chareCount, int numRanks)
{
  int* crange = new int[2];

  int perCoarseChare = chareCount / numRanks;
  crange[0] = perCoarseChare * thisIndex;
  crange[1] = crange[0] + perCoarseChare;
  if (crange[1] > chareCount)
    crange[1] = chareCount; 
  if (thisIndex == numRanks-1 && crange[1] < chareCount)
    crange[1] = chareCount;
  return crange;
}

void CoarseScaleModel::LagrangeNodal1()
{
  lulesh->LagrangeNodal1();
}

void CoarseScaleModel::LagrangeNodal2()
{
  lulesh->LagrangeNodal2();
}

void CoarseScaleModel::LagrangeElements()
{
  lulesh->LagrangeElements();
}

void CoarseScaleModel::LagrangeElements2()
{
  lulesh->LagrangeElements2();
}

void CoarseScaleModel::CalcTimeConstraintsForElems()
{
  lulesh->CalcTimeConstraintsForElems();
}

void CoarseScaleModel::TimeIncrement()
{

  if ((lulesh->domain.dtfixed() <= Real_t(0.0)) && (lulesh->domain.cycle() != Int_t(0))) {

    Real_t gnewdt = Real_t(1.0e+20) ;
    if (lulesh->domain.dtcourant() < gnewdt) {
       gnewdt = lulesh->domain.dtcourant() / Real_t(2.0) ;
    }
    if (lulesh->domain.dthydro() < gnewdt) {
       gnewdt = lulesh->domain.dthydro() * Real_t(2.0) / Real_t(3.0) ;
    }

    //printf("%d fs_dt_mod = %e\n", thisIndex, lulesh->finescale_dt_modifier);
    gnewdt *= lulesh->finescale_dt_modifier;

    // Reduction
    if (lulesh->domain.numSlices() == 1)
      updateTimeIncrement(gnewdt);
    else
      contribute(sizeof(double), (void *)&gnewdt, CkReduction::min_double, *cbTime);

  }

}

void CoarseScaleModel::TimeIncrement2()
{
  Real_t targetdt = lulesh->domain.stoptime() - lulesh->domain.time() ;

 // * TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE *
  if ((targetdt > lulesh->domain.deltatime()) &&
      (targetdt < (Real_t(4.0) * lulesh->domain.deltatime() / Real_t(3.0))) ) {
    targetdt = Real_t(2.0) * lulesh->domain.deltatime() / Real_t(3.0) ;
  }

  if (targetdt < lulesh->domain.deltatime()) {
    lulesh->domain.deltatime() = targetdt ;
  }

  lulesh->domain.time() += lulesh->domain.deltatime() ;

  ++lulesh->domain.cycle() ;

}


void CoarseScaleModel::UpdateStressForElems()
{
  // For this Coarse model's elements
  total_samples = 0;
  total_interpolations = 0;
  max_nonlinear_iters = 0;
  max_local_newton_iters = 0;
  numElems = lulesh->domain.numElem();

#ifdef _OPENMP
#pragma omp for
#endif
  for (Index_t k=0; k < numElems; ++k) {

    // advance constitutive model
    fineScaleArray(thisIndex, k).advance(lulesh->domain.deltatime(),
                                        lulesh->domain.cm_vel_grad(k),
                                        lulesh->domain.cm_vol_chng(k),
                                        state_size[k],
                                        (char*)lulesh->domain.cm_state(k));
  }
}


void CoarseScaleModel::updateAdvanceResults(int elemNum, ConstitutiveData cm_data, int ssize, char* state, int num_samples, int num_interpolations)
{
  // set new state
  lulesh->domain.cm(elemNum)->setState((void *)state);

  int num_iters = cm_data.num_Newton_iters;
  if (num_iters > max_local_newton_iters) max_local_newton_iters = num_iters;

  const Tensor2Sym& sigma_prime = cm_data.sigma_prime;

  Real_t sx  = lulesh->domain.sx(elemNum) = sigma_prime(1,1);
  Real_t sy  = lulesh->domain.sy(elemNum) = sigma_prime(2,2);
  Real_t sz  = - sx - sy;
  Real_t txy = lulesh->domain.txy(elemNum) = sigma_prime(2,1);
  Real_t txz = lulesh->domain.txz(elemNum) = sigma_prime(3,1);
  Real_t tyz = lulesh->domain.tyz(elemNum) = sigma_prime(3,2);
  
  lulesh->domain.mises(elemNum) = SQRT( Real_t(0.5) * ( (sy - sz)*(sy - sz) + (sz - sx)*(sz - sx) + (sx - sy)*(sx - sy) )
                    + Real_t(3.0) * ( txy*txy + txz*txz + tyz*tyz) );

  if (max_local_newton_iters > max_nonlinear_iters)
  {
    max_nonlinear_iters = max_local_newton_iters;
  }

  lulesh->domain.cm(elemNum)->getState(lulesh->domain.cm_state(elemNum));

  total_samples += num_samples;
  total_interpolations += num_interpolations;
}

void CoarseScaleModel::afterAdvance()
{
  // Reduction
  if (lulesh->domain.numSlices() == 1)
    UpdateStressForElems2(max_nonlinear_iters);
  else
  {
    contribute(sizeof(int), (void *)&max_nonlinear_iters, CkReduction::max_int, *cbIters);
    int total[2];
    total[0] = total_samples;
    total[1] = total_interpolations;
    contribute(2*sizeof(int), (void *)total, CkReduction::sum_int, *cbTotal);
  }

}

/* previous version
void CoarseScaleModel::UpdateStressForElems()
{
  
  int new_iters = lulesh->UpdateStressForElems();

  if (lulesh->domain.numSlices() == 1) {
    lulesh->UpdateStressForElems2(new_iters);
    return;
  }

  // Reduction on iters
  //printf("%d local iters = %d\n", thisIndex, new_iters);
  contribute(sizeof(int), (void *)&new_iters, CkReduction::max_int, *cbIters);

}
*/

void CoarseScaleModel::UpdateStressForElems2(int reducedIters)
{

  //printf("%d In UpdateStressForElems2 reducedIters = %d\n", thisIndex, reducedIters);

  lulesh->UpdateStressForElems2(reducedIters);

}

void CoarseScaleModel::sendNodalMass()
{
  int xferFields = 1;

  Real_t *fd = &lulesh->domain.nodalMass(0);
  Real_t **fieldData = &fd;

  // Send nodal mass to neighbors
  sendDataNodes(xferFields, fieldData);
}

#define NBR_M1 11
#define NBR_P1 22
void CoarseScaleModel::sendDataNodes(int xferFields, Real_t **fieldData)
{
  if (lulesh->domain.numSlices() == 1) return ;

  int size = lulesh->domain.commNodes();
  int offset = lulesh->domain.sliceHeight();
  int *iset = lulesh->domain.planeNodeIds;

  Real_t *sdataM1 = new Real_t[xferFields*size];
  Real_t *sdataP1 = new Real_t[xferFields*size];

  int myRank = lulesh->domain.sliceLoc() ;

  Real_t *destAddr ;
  bool planeMin, planeMax ;

  /* assume communication to 2 neighbors by default */
  planeMin = planeMax = true ;
  if (myRank == 0) {
    planeMin = false ;
  }
  if (myRank == (lulesh->domain.numSlices()-1)) {
    planeMax = false ;
  }

  /* send to myRank-1 and myRank+1 neighbors */

  /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
  if (planeMin) {
    destAddr = sdataM1;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *srcAddr = fieldData[fi] ;
      for (Index_t ii=0; ii<size; ++ii) {
        destAddr[ii] = srcAddr[iset[ii]] ;
        //printf("%d to M1 %d %d = %e\n", lulesh->domain.sliceLoc(), ii, iset[ii], srcAddr[iset[ii]]);
      }
      destAddr += size ;
    }
    destAddr -= xferFields*size ;

    // Send to myRank-1 neighbor
    if (xferFields == 1)
      thisProxy(thisIndex-1).receiveNodalMass(NBR_M1, xferFields*size, sdataM1);
    else if (xferFields == 3)
      thisProxy(thisIndex-1).receiveForce(NBR_M1, xferFields*size, sdataM1);
    else if (xferFields == 6)
      thisProxy(thisIndex-1).receivePositionVelocity(NBR_M1, xferFields*size, sdataM1);
    else
      printf("No where to send!\n");

  }

  if (planeMax && xferFields < 6) {
    destAddr = sdataP1;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *srcAddr = &fieldData[fi][offset] ;
      for (Index_t ii=0; ii<size; ++ii) {
        destAddr[ii] = srcAddr[iset[ii]] ;
        //printf("%d to P1 %d %d = %e\n", lulesh->domain.sliceLoc(), ii, iset[ii], srcAddr[iset[ii]]);
      }
      destAddr += size ;
    }
    destAddr -= xferFields*size ;

    // Send to myRank+1 neighbor
    if (xferFields == 1)
      thisProxy(thisIndex+1).receiveNodalMass(NBR_P1, xferFields*size, sdataP1);
    else if (xferFields == 3)
      thisProxy(thisIndex+1).receiveForce(NBR_P1, xferFields*size, sdataP1);
/*
    else if (xferFields == 6)
      thisProxy(thisIndex+1).receivePositionVelocity(NBR_P1, xferFields*size, sdataP1);
*/
    else
      printf("No send P1!\n");
  }

  delete [] sdataM1;
  delete [] sdataP1;
}

void CoarseScaleModel::sendDataElems(int xferFields, Real_t **fieldData)
{
  if (lulesh->domain.numSlices() == 1) return ;

  int size = lulesh->domain.commElems();
  int offset = lulesh->domain.sliceHeight() - 1;
  int *iset = lulesh->domain.planeElemIds;

  Real_t *sdataM1 = new Real_t[xferFields*size];
  Real_t *sdataP1 = new Real_t[xferFields*size];

  int myRank = lulesh->domain.sliceLoc() ;

  Real_t *destAddr ;
  bool planeMin, planeMax ;

  /* assume communication to 2 neighbors by default */
  planeMin = planeMax = true ;
  if (myRank == 0) {
    planeMin = false ;
  }
  if (myRank == (lulesh->domain.numSlices()-1)) {
    planeMax = false ;
  }

  /* send to myRank-1 and myRank+1 neighbors */

  /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
  if (planeMin) {
    destAddr = sdataM1;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *srcAddr = fieldData[fi] ;
      for (Index_t ii=0; ii<size; ++ii) {
        destAddr[ii] = srcAddr[iset[ii]] ;
        //printf("%d to M1 %d %d = %e\n", lulesh->domain.sliceLoc(), ii, iset[ii], srcAddr[iset[ii]]);
      }
      destAddr += size ;
    }
    destAddr -= xferFields*size ;

    // Send to myRank+1 neighbor
    if (xferFields == 1)
      thisProxy(thisIndex-1).receiveVelocityGrad(NBR_M1, xferFields*size, sdataM1);
    else
      printf("No where to send!\n");

  }

  if (planeMax) {
    destAddr = sdataP1;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *srcAddr = &fieldData[fi][offset] ;
      for (Index_t ii=0; ii<size; ++ii) {
        destAddr[ii] = srcAddr[iset[ii]] ;
        //printf("%d to P1 %d %d = %e\n", lulesh->domain.sliceLoc(), ii, iset[ii], srcAddr[iset[ii]]);
      }
      destAddr += size ;
    }
    destAddr -= xferFields*size ;

    // Send to myRank+1 neighbor
    if (xferFields == 1)
      thisProxy(thisIndex+1).receiveVelocityGrad(NBR_P1, xferFields*size, sdataP1);
    else
      printf("No send P1!\n");
  }

  delete [] sdataM1;
  delete [] sdataP1;
}

void CoarseScaleModel::updateNodalMass(int msgType, int rsize, Real_t rdata[])
{
  //printf("%d In updateNodalMass\n", thisIndex);

  int xferFields = 1;
  int size = lulesh->domain.commNodes();

  Real_t *fd = &(lulesh->domain.nodalMass(0)) ;
  Real_t **fieldData = &fd;

  // Receive nodal mass data from L/R neighbor
  receiveDataNodes(msgType, size, xferFields, fieldData, rdata);
}

void CoarseScaleModel::receiveDataNodes(int msgType, int size, int xferFields,
   Real_t **fieldData, Real_t rdata[])
{
  if (lulesh->domain.numSlices() == 1) return ;

  int offset = lulesh->domain.sliceHeight();
  int *iset = lulesh->domain.planeNodeIds;
 
  /* summation order should be from smallest value to largest */
  /* or we could try out kahan summation! */

  /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
  Real_t *srcAddr;
  if (msgType == NBR_P1) {
    //printf("%d received from rank-1 size = %d xferfields = %d\n", lulesh->domain.sliceLoc(), size, xferFields);
    srcAddr = rdata;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *destAddr = fieldData[fi] ;
      for (Index_t i=0; i<size; ++i) {
        if (xferFields < 6)
        destAddr[iset[i]] += srcAddr[i] ;
        else
        destAddr[iset[i]] = srcAddr[i] ;
        //printf("%d from M1 %d %d = %e\n", lulesh->domain.sliceLoc(), i, iset[i], srcAddr[i]);
      }
      srcAddr += size ;
    }
  }
  else if (msgType == NBR_M1) {
    //printf("%d received from rank+1  size = %d xferfields = %d\n", lulesh->domain.sliceLoc(), size, xferFields);
    srcAddr = rdata;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *destAddr = &(fieldData[fi][offset]) ;
      for (Index_t i=0; i<size; ++i) {
        if (xferFields < 6)
        destAddr[iset[i]] += srcAddr[i] ;
        else
        destAddr[iset[i]] = srcAddr[i] ; 
        //printf("%d from P1 %d %d = %e\n", lulesh->domain.sliceLoc(), i, iset[i], srcAddr[i]);
      }
      srcAddr += size ;
    }
  }

  //delete [] rdata;
}

void CoarseScaleModel::receiveDataElems(int msgType, int size, int xferFields,
   Real_t rdata[])
{
  if (lulesh->domain.numSlices() == 1) return ;

  int offset = lulesh->domain.sliceHeight() - 1;
  int *iset = lulesh->domain.planeElemIds;

  Real_t *fieldData[1] ;
  /* point into ghost data area */
   fieldData[0] = &(lulesh->domain.delv_xi(lulesh->domain.numElem())) ;

  /* summation order should be from smallest value to largest */
  /* or we could try out kahan summation! */

  /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
  Real_t *srcAddr;
  if (msgType == NBR_P1) {
    //printf("%d received from rank-1 size = %d xferfields = %d\n", lulesh->domain.sliceLoc(), size, xferFields);
    srcAddr = rdata;
    //fieldData[0] += size;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *destAddr = fieldData[fi] ;
      for (Index_t i=0; i<size; ++i) {
        destAddr[i] = srcAddr[i] ;
        //destAddr[iset[i]] = srcAddr[i] ;
        //printf("%d from P1 %d %d = %e\n", lulesh->domain.sliceLoc(), i, iset[i], srcAddr[i]);
      }
      srcAddr += size ;
      fieldData[fi] += size;
    }
  }
  else if (msgType == NBR_M1) {
    //printf("%d received from rank+1  size = %d xferfields = %d\n", lulesh->domain.sliceLoc(), size, xferFields);
    srcAddr = rdata;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *destAddr = fieldData[fi] ;
      //Real_t *destAddr = &(fieldData[fi][offset]) ;
      for (Index_t i=0; i<size; ++i) {
        destAddr[i] = srcAddr[i] ;
        //destAddr[iset[i]] = srcAddr[i] ;
        //printf("%d from P1 %d %d = %e\n", lulesh->domain.sliceLoc(), i, iset[i], srcAddr[i]);
      }
      srcAddr += size ;
      fieldData[fi] += size ;
    }
  }

}

void CoarseScaleModel::sendForce()
{
  int xferFields = 3;

  Real_t *fieldData[3] ;

  fieldData[0] = &(lulesh->domain.fx(0)) ;
  fieldData[1] = &(lulesh->domain.fy(0)) ;
  fieldData[2] = &(lulesh->domain.fz(0)) ;

  sendDataNodes(xferFields, fieldData);
}

void CoarseScaleModel::updateForce(int msgType, int rsize, Real_t rdata[])
{
  //printf("%d In updateForce\n", thisIndex);

  int xferFields = 3;
  int size = lulesh->domain.commNodes();

  Real_t *fieldData[3] ;

  fieldData[0] = &(lulesh->domain.fx(0)) ;
  fieldData[1] = &(lulesh->domain.fy(0)) ;
  fieldData[2] = &(lulesh->domain.fz(0)) ;

  // Receive force data from L/R neighbor
  receiveDataNodes(msgType, size, xferFields, fieldData, rdata);

}

void CoarseScaleModel::sendVelocityGrad()
{
  int xferFields = 1;

  Real_t *fieldData[1] ;

  fieldData[0] = &(lulesh->domain.delv_xi(0)) ;

  sendDataElems(xferFields, fieldData);
}

void CoarseScaleModel::updateVelocityGrad(int msgType, int rsize, Real_t rdata[])
{
  //printf("%d In updateVelocityGrad\n", thisIndex);

  int xferFields = 1;
  int size = lulesh->domain.commElems();

  Real_t *fieldData[1] ;

  fieldData[0] = &(lulesh->domain.delv_xi(0)) ;

  // Receive force data from L/R neighbor
  receiveDataElems(msgType, size, xferFields, rdata);
}
  
void CoarseScaleModel::sendPositionVelocity()
{
  int xferFields = 6;

  Real_t *fieldData[6] ;

  fieldData[0] = &(lulesh->domain.x(0)) ;
  fieldData[1] = &(lulesh->domain.y(0)) ;
  fieldData[2] = &(lulesh->domain.z(0)) ;
  fieldData[3] = &(lulesh->domain.xd(0)) ;
  fieldData[4] = &(lulesh->domain.yd(0)) ;
  fieldData[5] = &(lulesh->domain.zd(0)) ;

  sendDataNodes(xferFields, fieldData);
}

void CoarseScaleModel::updatePositionVelocity(int msgType, int rsize, Real_t rdata[])
{ 
  //printf("%d In updatPositionVelocity\n", thisIndex);

  int xferFields = 6;
  int size = lulesh->domain.commNodes();

  Real_t *fieldData[6] ;

  fieldData[0] = &(lulesh->domain.x(0)) ;
  fieldData[1] = &(lulesh->domain.y(0)) ;
  fieldData[2] = &(lulesh->domain.z(0)) ;
  fieldData[3] = &(lulesh->domain.xd(0)) ;
  fieldData[4] = &(lulesh->domain.yd(0)) ;
  fieldData[5] = &(lulesh->domain.zd(0)) ;

  // Receive force data from L/R neighbor
  receiveDataNodes(msgType, size, xferFields, fieldData, rdata);    
}

void CoarseScaleModel::makeADump(bool forceDump)
{
#ifdef SILO
   int debug_topology = 0; //Eventually pull this out
   int sampling = this->use_adaptive_sampling; //Verify these are the same
   //int sampling = 1;
   if ((this->visit_data_interval != 0) and ( (lulesh->domain.cycle() % this->visit_data_interval == 0) or forceDump)) {
      DumpDomain(&(lulesh->domain), lulesh->domain.sliceLoc(), lulesh->domain.numSlices(), 
                 ((lulesh->domain.numSlices() == 1) ? this->file_parts : 0), sampling, debug_topology ) ;
   }
#endif
}


void CoarseScaleModel::setSiloParams(int numParts, int dataInterval)
{
  this->file_parts = numParts;
  this->visit_data_interval = dataInterval;
}

void CoarseScaleModel::setRemoteDB(bool remoteDB)
{
  this->remoteDB = remoteDB;
}
