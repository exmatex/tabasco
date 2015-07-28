#include "TabaSCo.decl.h"
#include "CoarseScaleModel.h"
#include "DBInterface.h"
#include "FineScaleModel.h"

extern CProxy_Main mainProxy;
extern CProxy_CoarseScaleModel coarseScaleArray;
extern CProxy_FineScaleModel fineScaleArray;

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
  cbTime = new CkCallback(CkReductionTarget(CoarseScaleModel, reduceTimeIncrement), coarseScaleArray);
  cbIters = new CkCallback(CkReductionTarget(CoarseScaleModel, reduceIters), coarseScaleArray);

}

CoarseScaleModel::CoarseScaleModel(CkMigrateMessage *msg)
{

}

CoarseScaleModel::~CoarseScaleModel()
{
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
  p|count;
}

void CoarseScaleModel::startElementFineScaleQuery(int step, int nelems)
{
  nstep = step;

  //CkCallback cb(CkIndex_CoarseScaleModel::receiveNewPoint((Msg*)NULL), thisProxy);
/*
  for (int j = 0; j < nelems; j++)
  {
    //fineScaleArray(thisIndex.x, thisIndex.y, thisIndex.z, j).query(j, nstep, currentPt, cb);
    fineScaleArray(thisIndex.x, j).query(j, nstep, currentPt);
  }
*/
}

void CoarseScaleModel::updateElement(int whichEl, int whichIter, int newPt)
{
  int currentPt = newPt;
	  printf("Iter %d Coarse %d Element %d update newPt %d\n", whichIter, thisIndex, whichEl, newPt);
}

void CoarseScaleModel::initialize(int numRanks, bool useAdaptiveSampling)
{
  lulesh->Initialize(thisIndex, numRanks);
  lulesh->ConstructFineScaleModel(useAdaptiveSampling);
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
  
  int new_iters = lulesh->UpdateStressForElems();

  if (lulesh->domain.numSlices() == 1) {
    lulesh->UpdateStressForElems2(new_iters);
    return;
  }

  // Reduction on iters
  //printf("%d local iters = %d\n", thisIndex, new_iters);
  contribute(sizeof(int), (void *)&new_iters, CkReduction::max_int, *cbIters);

}

void CoarseScaleModel::UpdateStressForElems2(int reducedIters)
{

  //printf("%d iters = %d\n", thisIndex, reducedIters);

  lulesh->UpdateStressForElems2(reducedIters);

}

void CoarseScaleModel::go(int numRanks, bool useAdaptiveSampling)
{
  lulesh->go(thisIndex, numRanks, useAdaptiveSampling);
}

void CoarseScaleModel::sendNodalMass()
{
  int xferFields = 1;

  Real_t *fd = &lulesh->domain.nodalMass(0);
  Real_t **fieldData = &fd;

  // Send nodal mass to neighbors
  sendDataNodes(xferFields, fieldData);
}

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
