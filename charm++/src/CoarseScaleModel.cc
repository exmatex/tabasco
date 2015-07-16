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
  //printf("%d In timeReduce\n", thisIndex);

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

void CoarseScaleModel::CalcTimeConstraintsForElems()
{
  lulesh->CalcTimeConstraintsForElems();
}

void CoarseScaleModel::TimeIncrement()
{
  if (lulesh->domain.numSlices() == 1) {
    lulesh->TimeIncrement();
    --lulesh->domain.cycle();
    return;
  }

  if ((lulesh->domain.dtfixed() <= Real_t(0.0)) && (lulesh->domain.cycle() != Int_t(0))) {

    //Real_t olddt = lulesh->domain.deltatime() ;

    Real_t gnewdt = Real_t(1.0e+20) ;
    //Real_t newdt;
    if (lulesh->domain.dtcourant() < gnewdt) {
       gnewdt = lulesh->domain.dtcourant() / Real_t(2.0) ;
    }
    if (lulesh->domain.dthydro() < gnewdt) {
       gnewdt = lulesh->domain.dthydro() * Real_t(2.0) / Real_t(3.0) ;
    }

    gnewdt *= lulesh->finescale_dt_modifier;

    // Reduction
    contribute(sizeof(double), (void *)&gnewdt, CkReduction::min_double);
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
  lulesh->UpdateStressForElems();
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
  sendData(xferFields, fieldData);
}

void CoarseScaleModel::sendData(int xferFields, Real_t **fieldData)
{
  int size = lulesh->domain.commNodes();
  int offset = lulesh->domain.sliceHeight();
  int *iset = lulesh->domain.planeNodeIds;

  Real_t *sdataM1 = new Real_t[xferFields*size];
  Real_t *sdataP1 = new Real_t[xferFields*size];

  if (lulesh->domain.numSlices() == 1) return ;

  int myRank = lulesh->domain.sliceLoc() ;

  Real_t *destAddr ;
  bool planeMin, planeMax ;

  /* assume communication to 2 neighbors by default */
  planeMin = planeMax = true ;
  if (lulesh->domain.sliceLoc() == 0) {
    planeMin = false ;
  }
  if (lulesh->domain.sliceLoc() == (lulesh->domain.numSlices()-1)) {
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

void CoarseScaleModel::updateNodalMass(int msgType, int rsize, Real_t rdata[])
{
  int xferFields = 1;
  int size = lulesh->domain.commNodes();

  Real_t *fd = &(lulesh->domain.nodalMass(0)) ;
  Real_t **fieldData = &fd;

  // Receive nodal mass data from L/R neighbor
  receiveData(msgType, size, xferFields, fieldData, rdata);
}

void CoarseScaleModel::receiveData(int msgType, int size, int xferFields,
   Real_t **fieldData, Real_t rdata[])
{
  int offset = lulesh->domain.sliceHeight();
  int *iset = lulesh->domain.planeNodeIds;
 
  if (lulesh->domain.numSlices() == 1) return ;
  /* summation order should be from smallest value to largest */
  /* or we could try out kahan summation! */

  /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
  Real_t *srcAddr;
  if (msgType == NBR_P1) {
    srcAddr = rdata;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *destAddr = fieldData[fi] ;
      for (Index_t i=0; i<size; ++i) {
        destAddr[iset[i]] += srcAddr[i] ;
        //printf("d %d %d = %e\n", i, iset[i], destAddr[iset[i]]);
      }
      srcAddr += size ;
    }
  }
  else if (msgType == NBR_M1) {
    srcAddr = rdata;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
      Real_t *destAddr = &(fieldData[fi][offset]) ;
      for (Index_t i=0; i<size; ++i) {
        destAddr[iset[i]] += srcAddr[i] ;
        //printf("d %d %d = %e\n", i, iset[i], destAddr[iset[i]]);
      }
      srcAddr += size ;
    }
  }

}

/*
void CoarseScaleModel::receiveNodalMass(int msgType, int size, Real_t rdata[])
{

}
*/

void CoarseScaleModel::sendForce()
{
  int xferFields = 3;

  Real_t *fieldData[3] ;

  fieldData[0] = &(lulesh->domain.fx(0)) ;
  fieldData[1] = &(lulesh->domain.fy(0)) ;
  fieldData[2] = &(lulesh->domain.fz(0)) ;

  sendData(xferFields, fieldData);
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
  receiveData(msgType, size, xferFields, fieldData, rdata);
}

/*
void CoarseScaleModel::receiveForce(int msgType, int size, Real_t rdata[])
{

}
*/

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

  sendData(xferFields, fieldData);
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
  receiveData(msgType, size, xferFields, fieldData, rdata);    
} 

/*
void CoarseScaleModel::receivePositionVelocity(int msgType, int size, Real_t rdata[])
{

}
*/
