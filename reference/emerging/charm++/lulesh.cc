/**
 * Original port of LULESH by Felix Wang (6/2012)
 * 
 * Ghost Data Notation/Orientation:
 *     Z (up)
 *     ^
 *     |  Standard 3D Cartesian Coordinate System
 *     |
 *     X - - > Y (right)
 *  (front)
 * Chare indexing is: [z][y][x]
 *   given (x,y,z) we have idx = z*dimY*dimZ + y*dimX + x
 */

#include <vector>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pup_stl.h"
#include "lulesh.h"

//
// Charm++ Bookkeeping Definitions
//
#define DIR_SELF    0
#define DIR_FRONT   1
#define DIR_BACK    2
#define DIR_RIGHT   3
#define DIR_LEFT    4
#define DIR_UP      5
#define DIR_DOWN    6
#define CORNER_FRU  7
#define CORNER_FRD  8
#define CORNER_FLU  9
#define CORNER_FLD  10
#define CORNER_BRU  11
#define CORNER_BRD  12
#define CORNER_BLU  13
#define CORNER_BLD  14
#define EDGE_PXPY   15
#define EDGE_PXMY   16
#define EDGE_MXPY   17
#define EDGE_MXMY   18
#define EDGE_PXPZ   19
#define EDGE_PXMZ   20
#define EDGE_MXPZ   21
#define EDGE_MXMZ   22
#define EDGE_PYPZ   23
#define EDGE_PYMZ   24
#define EDGE_MYPZ   25
#define EDGE_MYMZ   26
#define FRONT_OFFSET(x,y,z) 0
#define BACK_OFFSET(x,y,z)  y*z
#define RIGHT_OFFSET(x,y,z) 2*y*z
#define LEFT_OFFSET(x,y,z)  2*y*z + x*z
#define UP_OFFSET(x,y,z)    2*y*z + 2*x*z
#define DOWN_OFFSET(x,y,z)  2*y*z + 2*x*z + x*y

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int elemDimX;
/*readonly*/ int elemDimY;
/*readonly*/ int elemDimZ;
/*readonly*/ int blockDimX;
/*readonly*/ int blockDimY;
/*readonly*/ int blockDimZ;
/*readonly*/ int ghostDimX;
/*readonly*/ int ghostDimY;
/*readonly*/ int ghostDimZ;
/*readonly*/ int chareDimX;
/*readonly*/ int chareDimY;
/*readonly*/ int chareDimZ;
/*readonly*/ double charesPerPE;
/*readonly*/ int numElems;
/*readonly*/ int numNodes;
/*readonly*/ int ghostElems;
/*readonly*/ int frontOffset;
/*readonly*/ int backOffset;
/*readonly*/ int rightOffset;
/*readonly*/ int leftOffset;
/*readonly*/ int upOffset;
/*readonly*/ int downOffset;
/*readonly*/ Real_t stopTime;
/*readonly*/ int lbFrequency;
/*readonly*/ int cpFrequency;


/*******************************************
 *
 * Reduction
 *
 *******************************************/
CkReduction::reducerType min_real;
/*initnode*/
void registerMinReal(void) {
  min_real = CkReduction::addReducer(minReal);
}
CkReductionMsg *minReal(int nMsg, CkReductionMsg **msgs) {
  Real_t ret = DT_MAX;
  //CkPrintf("number of messages: %d\n", nMsg);
  for (int i = 0; i < nMsg; i++) {
    // Sanity check:
    CkAssert(msgs[i]->getSize() == sizeof(Real_t));
    // Extract data and reduce
    Real_t m = *(Real_t *)msgs[i]->getData();
    ret = std::min(ret,m);
  }
  return CkReductionMsg::buildNew(sizeof(Real_t),&ret);
}

CkReduction::reducerType ave_time;
/*initnode*/
void registerAveTime(void) {
  ave_time = CkReduction::addReducer(aveTime);
}
CkReductionMsg *aveTime(int nMsg, CkReductionMsg **msgs) {
  // There are as many return values as there are timers
  // in the code.
  double ret[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (int i = 0; i < nMsg; i++) {
    // Sanity check
    CkAssert(msgs[i]->getSize() == 11*sizeof(double));
    // Extract data and reduce
    double *m = (double *)msgs[i]->getData();
    ret[0] = m[0];
    ret[1] += m[1]/nMsg;
    ret[2] += m[2]/nMsg;
    ret[3] += m[3]/nMsg;
    ret[4] += m[4]/nMsg;
    ret[5] += m[5]/nMsg;
    ret[6] += m[6]/nMsg;
    ret[7] += m[7]/nMsg;
    ret[8] += m[8]/nMsg;
    ret[9] += m[9]/nMsg;
    ret[10] += m[10]/nMsg;
  }
  return CkReductionMsg::buildNew(11*sizeof(double),ret);
}


/*******************************************
 *
 * Main
 *
 *******************************************/
Main::Main(CkArgMsg *msg) {
  // Process command line arguments
#if WRITE_TO_FILE
  if ((msg->argc != 4) && (msg->argc != 8)) {
    CkPrintf("Usage: Takes exactly 3 or 7 command line arguments\n"
             "  %s [elemDim] [blockDim] [output file]\n"
             "  %s [elemDimX] [elemDimY] [elemDimZ] [blockDimX] [blockDimY] [blockDimZ] [output file]\n",
             msg->argv[0], msg->argv[0]);
#else
  if ((msg->argc != 3) && (msg->argc != 7)) {
    CkPrintf("Usage: Takes exactly 2 or 6 command line arguments\n"
             "  %s [elemDim] [blockDim]\n"
             "  %s [elemDimX] [elemDimY] [elemDimZ] [blockDimX] [blockDimY] [blockDimZ]\n",
             msg->argv[0], msg->argv[0]);
#endif
    CkExit();
  }
#if WRITE_TO_FILE
  else if (msg->argc == 4) {
#else
  else if (msg->argc == 3) {
#endif
    elemDimX = elemDimY = elemDimZ = atoi(msg->argv[1]);
    blockDimX = blockDimY = blockDimZ = atoi(msg->argv[2]);
#if WRITE_TO_FILE
    fileOut = msg->argv[3];
#endif
  }
  else {
    elemDimX = atoi(msg->argv[1]);
    elemDimY = atoi(msg->argv[2]);
    elemDimZ = atoi(msg->argv[3]);
    blockDimX = atoi(msg->argv[4]);
    blockDimY = atoi(msg->argv[5]);
    blockDimZ = atoi(msg->argv[6]);
#if WRITE_TO_FILE
    fileOut = msg->argv[7];
#endif
  }
  if ((elemDimX < 1) || (elemDimY < 1) || (elemDimZ < 1) ||
      (blockDimZ < 1) || (blockDimY < 1) || (blockDimZ < 1)) {
    CkPrintf("Usage: elemDim and blockDim must take positive values\n");
    CkExit();
  }
  if ((elemDimX % blockDimX) || (elemDimY % blockDimY) || (elemDimZ % blockDimZ)) {
    CkPrintf("Usage: elemDim must be divided evenly by blockDim\n");
    CkExit();
  }
  delete msg;
  // Set additional readonly parameters
  ghostDimX = blockDimX+1;
  ghostDimY = blockDimY+1;
  ghostDimZ = blockDimZ+1;
  chareDimX = elemDimX/blockDimX;
  chareDimY = elemDimY/blockDimY;
  chareDimZ = elemDimZ/blockDimZ;
  charesPerPE = (double)(chareDimX*chareDimY*chareDimZ)/CkNumPes();
  if(charesPerPE < 1) { charesPerPE = 1; }
  numElems = blockDimX*blockDimY*blockDimZ;
  numNodes = ghostDimX*ghostDimY*ghostDimZ;
  ghostElems = numElems + 2*blockDimY*blockDimZ
               + 2*blockDimX*blockDimZ + 2*blockDimX*blockDimY;
  // Offsets and other bookkeeping
  frontOffset = numElems + FRONT_OFFSET(blockDimX,blockDimY,blockDimZ);
  backOffset  = numElems + BACK_OFFSET(blockDimX,blockDimY,blockDimZ);
  rightOffset = numElems + RIGHT_OFFSET(blockDimX,blockDimY,blockDimZ);
  leftOffset  = numElems + LEFT_OFFSET(blockDimX,blockDimY,blockDimZ);
  upOffset    = numElems + UP_OFFSET(blockDimX,blockDimY,blockDimZ);
  downOffset  = numElems + DOWN_OFFSET(blockDimX,blockDimY,blockDimZ);
  stopTime = STOP_TIME;
  lbFrequency = LB_FREQUENCY;
  cpFrequency = CP_FREQUENCY;
  mainProxy = thisProxy;

  // Display of some basic information
  // Elements breakdown
  // Chares breakdown [number of chares per processor]
  CkPrintf("Lulesh (Charm++)\n"
           "  Elements: %d (%d x %d x %d)\n"
           "  Chares: %d [%.2g] (%d x %d x %d)\n",
           elemDimX*elemDimY*elemDimZ,
           elemDimX, elemDimY, elemDimZ,
           chareDimX*chareDimY*chareDimZ, charesPerPE,
           chareDimX, chareDimY, chareDimZ);

  // Timing Code
  startSimTime = CmiWallTimer();
  // Initialize the domains and set reduction client
  domains = CProxy_Domain::ckNew(chareDimX,chareDimY,chareDimZ);
  CkCallback *cb = new CkCallback(CkReductionTarget(Domain, updateTimeIncrement), domains);
  domains.ckSetReductionClient(cb);
}

Main::Main(CkMigrateMessage *msg): CBase_Main(msg) {
  delete msg;
}

void
Main::pup(PUP::er &p) {
  CBase_Main::pup(p);
  p|domains;
  p|fileOut;
  // Timing Variables
  p|startSimTime;
}

void
Main::averageTimeCheckin(CkReductionMsg *msg) {
  double *timers = (double *)(msg->getData());
  delete msg;
  // Compensate for approximately how many
  // Chares there are per PE
  iterTimers[0] = timers[0];                // Iterations
  iterTimers[1] = timers[1];                // Itertime
#if ITER_TIMING
  iterTimers[2] = timers[2]*charesPerPE;    // Iterticks
  iterTimers[3] = timers[3]*charesPerPE;    // Stress
  iterTimers[4] = timers[4]*charesPerPE;    // Hourglass
  iterTimers[5] = timers[5]*charesPerPE;    // Positions
  iterTimers[6] = timers[6]*charesPerPE;    // Kinematics
  iterTimers[7] = timers[7]*charesPerPE;    // Viscocity
  iterTimers[8] = timers[8]*charesPerPE;    // Materials
  iterTimers[9] = timers[9]*charesPerPE;    // Volume
  iterTimers[10] = timers[10]*charesPerPE;  // Time Constraints
#endif
  // Timing Code
  endSimTime = CmiWallTimer();
  simTime = endSimTime - startSimTime;
  // Print summary
  totalIterations = (int)iterTimers[0];
  CkPrintf("Lulesh (Charm++) Completed\n"
           "  #Cores: %d #Iterations: %d Total Time: %.4e (s)\n",
           CkNumPes(), totalIterations, iterTimers[1]);
  // Finish by writing to file
#if WRITE_TO_FILE
  writeToFile();
#endif
  CkExit();
}

void
Main::initCheckin(CkReductionMsg *msg) {
  delete msg;
  domains.run();
}

void
Main::writeToFile() {
  // Open file
  FILE *pData;
	pData = fopen(fileOut.c_str(),"a");
  if (pData != NULL) {
#if ITER_TIMING
    // Print more information for when we gather more information
    double physicsTime = iterTimers[3] + iterTimers[4] + iterTimers[5]
                       + iterTimers[6] + iterTimers[7] + iterTimers[8]
                       + iterTimers[9] + iterTimers[10];
    // blockDim, version, PEs, iterations, simTime, iterTime, time/iter
    fprintf(pData, "%d Charm %d %d %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
        blockDimX, CkNumPes(), totalIterations,
        simTime, iterTimers[1], iterTimers[1]/totalIterations,
        iterTimers[2]/totalIterations, physicsTime/totalIterations,
        iterTimers[3]/totalIterations, iterTimers[4]/totalIterations,
        iterTimers[5]/totalIterations, iterTimers[6]/totalIterations,
        iterTimers[7]/totalIterations, iterTimers[8]/totalIterations,
        iterTimers[9]/totalIterations, iterTimers[10]/totalIterations);
#else
    // elemDim, version, blockDim, PEs, iterations, simTime, time/iter
    fprintf(pData, "%d Charm %d %d %.12e %.12e\n",
        blockDimX, CkNumPes(), totalIterations,
        iterTimers[1], iterTimers[1]/totalIterations);
#endif
	  fclose(pData);
  }
  else {
    fprintf(stderr, "Could not open file for writing\n");
  }
  // Print correctness
  domains(0,0,0).printEnergy();
}

/*******************************************
 *
 * Domain
 *
 *******************************************/
Domain::Domain() {
  __sdag_init();
  usesAtSync = CmiTrue;
#if ITER_TIMING
  // Initialize Timers
  stressTime = 0;
  hgTime = 0;
  posTime = 0;
  kinTime = 0;
  viscTime = 0;
  matTime = 0;
  volTime = 0;
  consTime = 0;
#endif 
  // Initialize Lulesh Data
  InitializeLuleshData();

  // Send Nodal Mass across domains
  sendNodalMass();
}

Domain::Domain(CkMigrateMessage *msg): CBase_Domain(msg) {
  __sdag_init();
  delete msg;
}

Domain::~Domain() {
  // Delete various items that may be allocated
  // Ghosts
  domain->ng_front.clear();
  domain->ng_back.clear();
  domain->ng_right.clear();
  domain->ng_left.clear();
  domain->ng_up.clear();
  domain->ng_down.clear();
  // Node Persist
  domain->m_x.clear();
  domain->m_y.clear();
  domain->m_z.clear();
  domain->m_xd.clear();
  domain->m_yd.clear();
  domain->m_zd.clear();
  domain->m_xdd.clear();
  domain->m_ydd.clear();
  domain->m_zdd.clear();
  domain->m_fx.clear();
  domain->m_fy.clear();
  domain->m_fz.clear();
  domain->m_nodalMass.clear();
  // Nodesets
  domain->m_symmX.clear();
  domain->m_symmY.clear();
  domain->m_symmZ.clear();
  // Elem Persistent
  domain->m_matElemlist.clear();
  domain->m_nodelist.clear();
  domain->m_lxim.clear();
  domain->m_lxip.clear();
  domain->m_letam.clear();
  domain->m_letap.clear();
  domain->m_lzetam.clear();
  domain->m_lzetap.clear();
  domain->m_elemBC.clear();
  domain->m_e.clear();
  domain->m_p.clear();
  domain->m_q.clear();
  domain->m_ql.clear();
  domain->m_qq.clear();
  domain->m_v.clear();
  domain->m_volo.clear();
  domain->m_delv.clear();
  domain->m_vdov.clear();
  domain->m_arealg.clear();
  domain->m_ss.clear();
  domain->m_elemMass.clear();
  // Elem Temporary
  domain->m_dxx.clear();
  domain->m_dyy.clear();
  domain->m_dzz.clear();
  domain->m_delv_xi.clear();
  domain->m_delv_eta.clear();
  domain->m_delv_zeta.clear();
  domain->m_delx_xi.clear();
  domain->m_delx_eta.clear();
  domain->m_delx_zeta.clear();
  domain->m_vnew.clear();
  domain->m_determ.clear();
  // Everything else
  delete domain;
}

// Pack and Unpack routine for migrating Chares
//
void
Domain::pup(PUP::er &p) {
  CBase_Domain::pup(p);
  __sdag_pup(p);
  // Clocking Variables
  p|iterations;
  p|startIterTime;
  p|endIterTime;
  p|iterTime;
#if ITER_TIMING
  // Timing Variables
  p|startIterTicks;
  p|iterTicks;
  p|tempTime;
  p|stressTime;
  p|hgTime;
  p|posTime;
  p|kinTime;
  p|viscTime;
  p|matTime;
  p|volTime;
  p|consTime;
#endif
  // Charm Communication
  p|numElemGhosts;
  p|numNodeGhosts;
  p|ghostNodeCount;
  p|ghostElemCount;
  // Domain
  //
  if(p.isUnpacking()) {
    domain = new Dom;
    domain->ng_front.resize(3*ghostDimY*ghostDimZ,0.0);
    domain->ng_back.resize(3*ghostDimY*ghostDimZ,0.0);
    domain->ng_right.resize(3*ghostDimX*ghostDimZ,0.0);
    domain->ng_left.resize(3*ghostDimX*ghostDimZ,0.0);
    domain->ng_up.resize(3*ghostDimX*ghostDimY,0.0);
    domain->ng_down.resize(3*ghostDimX*ghostDimY,0.0);
    domain->m_determ.resize(numElems, 0.0);
  }
  // Simulation Related 
  p|domain->deltaTime;
  p|domain->totalTime;
  // TODO: Figure out which of these don't need
  //       to be transferred and simply need to
  //       be reallocated on the other side
  // Node centered persistent
  p|domain->m_x;
  p|domain->m_y;
  p|domain->m_z;
  p|domain->m_xd;
  p|domain->m_yd;
  p|domain->m_zd;
  p|domain->m_xdd;
  p|domain->m_ydd;
  p|domain->m_zdd;
  p|domain->m_fx;
  p|domain->m_fy;
  p|domain->m_fz;
  p|domain->m_nodalMass;
  // Node centered nodesets
  p|domain->m_symmX;
  p|domain->m_symmY;
  p|domain->m_symmZ;
  // Elem centered persistent
  p|domain->m_matElemlist;
  p|domain->m_nodelist;
  p|domain->m_lxim;
  p|domain->m_lxip;
  p|domain->m_letam;
  p|domain->m_letap;
  p|domain->m_lzetam;
  p|domain->m_lzetap;
  p|domain->m_elemBC;
  p|domain->m_e;
  p|domain->m_p;
  p|domain->m_q;
  p|domain->m_ql;
  p|domain->m_qq;
  p|domain->m_v;
  p|domain->m_volo;
  p|domain->m_delv;
  p|domain->m_vdov;
  p|domain->m_arealg;
  p|domain->m_ss;
  p|domain->m_elemMass;
  // Elem centered temporary
  p|domain->m_dxx;
  p|domain->m_dyy;
  p|domain->m_dzz;
  p|domain->m_delv_xi;
  p|domain->m_delv_eta;
  p|domain->m_delv_zeta;
  p|domain->m_delx_xi;
  p|domain->m_delx_eta;
  p|domain->m_delx_zeta;
  p|domain->m_vnew;
  // Parameters
  p|domain->u_cut;
  p|domain->hgcoef;
  p|domain->qstop;
  p|domain->monoq_max_slope;
  p|domain->monoq_limiter_mult;
  p|domain->e_cut;
  p|domain->p_cut;
  p|domain->ss4o3;
  p|domain->q_cut;
  p|domain->v_cut;
  p|domain->qlc_monoq;
  p|domain->qqc_monoq;
  p|domain->qqc;
  p|domain->eosvmax;
  p|domain->eosvmin;
  p|domain->pmin;
  p|domain->emin;
  p|domain->dvovmax;
  p|domain->refdens;
  p|domain->m_dtcourant;
  p|domain->m_dthydro;
}

// Load Balancing Routines
//
void
Domain::startLB() {
  // Save variables here so that we don't
  // need to check back in with the mainchare
  AtSync();
}

void
Domain::ResumeFromSync() {
  beginIteration();
}

// -----------------------------------------
//
// Charm++ Program Flow

void
Domain::beginIteration() {
#if ITER_TIMING
  tempTime = getticks();
#endif
  // Computation
  // Calculate Forces of Lagrange Leap Frog
  LagrangeNodal_Forces();
#if ITER_TIMING
  hgTime += elapsed(getticks(),tempTime);
#endif

#if CHECKPOINT_DEBUG
  if (domain->totalTime >= 4.0e-5 && domain->totalTime <= 4.2e-5 && thisIndex.x == 1 && thisIndex.y == 0 &&  thisIndex.z == 1) {
    if(CkHasCheckpoints()) {
      CkDieNow();
    }
  }
#endif
}

void
Domain::resumeNodeIteration() {
  // Computation
  // Compute Node positions and Element Kinematics
#if ITER_TIMING
  tempTime = getticks();
#endif
  LagrangeNodal_Positions();
#if ITER_TIMING
  posTime += elapsed(getticks(),tempTime);
  tempTime = getticks();
#endif
  // Compute the kinematics and viscocity gradients
  LagrangeElements_KinVisc();
#if ITER_TIMING
  kinTime += elapsed(getticks(),tempTime);
#endif
}
void
Domain::resumeElemIteration() {
#if ITER_TIMING
  tempTime = getticks();
#endif
  // Computation
  // Compute Q for Elems, apply material
  //   properties, update volumes
  LagrangeElements_QEOS();
#if ITER_TIMING
  volTime += elapsed(getticks(),tempTime);
#endif

#if APPLY_TIME_CONSTRAINTS
#if ITER_TIMING
  tempTime = getticks();
#endif
  // Timestep Updates
  CalcTimeConstraintsForElems();
#if ITER_TIMING
  consTime += elapsed(getticks(),tempTime);
#endif
  // Perform Reduction to main chare
  contribute(sizeof(Real_t), &domain->deltaTime, min_real);
#else
  contribute(sizeof(Real_t), NULL, CkReduction::nop);
#endif
}

/*
void
Domain::updateTimeIncrement(Real_t dt) {
#if USE_FIXED_DT
  // Test for fixed dt
  domain->deltaTime = DT_FIXED;
#elif APPLY_TIME_CONSTRAINTS
  domain->deltaTime = dt;
#else
  domain->deltaTime = DT_FIXED;
#endif

#if ITER_TIMING
  if (iterations == 0) {
    startIterTime = CmiWallTimer();
    startIterTicks = getticks();
  }
#endif
  
  domain->totalTime += domain->deltaTime; 
  // Iteration controls
  if (domain->deltaTime == Real_t(0.0) || iterations >= MAX_ITERATIONS) {
    domain->totalTime -= domain->deltaTime;
    // Stop Condition
    //
    // Timing Code
    endIterTime = CmiWallTimer();
#if ITER_TIMING
    iterTicks = elapsed(getticks(),startIterTicks);
    iterTime = endIterTime - startIterTime;
#endif
#if WRITE_TO_FILE
    // Average the timer values across all domains
    averageTime();
#else
    CkExit();
#endif
  }
  else {
    // Continue Condition
    //
    iterations++;
#if LULESH_SHOW_PROGRESS
    // Display Timing
    CkPrintf("%d: time = %e, dt = %e\n",iterations, domain->totalTime, domain->deltaTime);
#endif
#if PERFORM_CHECKPOINT
    if (iterations % cpFrequency == 0) {
      CkPrintf("Triggering Checkpointing [%d] to:\n", iterations);
      CkCallback cpCB(CkIndex_Domain::beginIteration(),thisProxy);
      CkStartMemCheckpoint(cpCB);
    }
#if PERFORM_LOAD_BALANCE
    else if (iterations % lbFrequency == 0) {
      CkPrintf("Triggering Load Balancing [%d]...\n", iterations);
      startLB();
    }
#endif
    else {
      beginIteration();
    }
#else
#if PERFORM_LOAD_BALANCE
    if (iterations % lbFrequency == 0) {
      CkPrintf("Triggering Load Balancing [%d]...\n", iterations);
      startLB();
    }
    else {
      beginIteration();
    }
#else
    beginIteration();
#endif
#endif
  }
}
*/

void
Domain::averageTime() {
  double timers[11];
  timers[0] = (double)iterations;
  timers[1] = iterTime;
#if ITER_TIMING
  timers[2] = iterTicks;
  timers[3] = stressTime;
  timers[4] = hgTime;
  timers[5] = posTime;
  timers[6] = kinTime;
  timers[7] = viscTime;
  timers[8] = matTime;
  timers[9] = volTime;
  timers[10] = consTime;
#endif
  CkCallback *cb = new CkCallback(CkIndex_Main::averageTimeCheckin(NULL),mainProxy);
  contribute(11*sizeof(double), timers, ave_time, *cb);
}

// -----------------------------------------
//
// Charm++ Send and Receive Ghost Data

void
Domain::sendNodeGhosts() {
  Real_t *frontGhost = new Real_t[ghostDimY*ghostDimZ*3];
  Real_t *backGhost = new Real_t[ghostDimY*ghostDimZ*3];
  Real_t *rightGhost = new Real_t[ghostDimX*ghostDimZ*3];
  Real_t *leftGhost = new Real_t[ghostDimX*ghostDimZ*3];
  Real_t *upGhost = new Real_t[ghostDimX*ghostDimY*3];
  Real_t *downGhost = new Real_t[ghostDimX*ghostDimY*3];
  Real_t *xEdge = new Real_t[ghostDimX*3];
  Real_t *yEdge = new Real_t[ghostDimY*3];
  Real_t *zEdge = new Real_t[ghostDimZ*3];
  Real_t *cornerGhost = new Real_t[3];
  // Create Ghost data
  // yz plane
  for (Index_t k = 0; k < ghostDimZ; k++) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      frontGhost[(k*ghostDimY+j)*3]   = domain->m_fx[k*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)];
      frontGhost[(k*ghostDimY+j)*3+1] = domain->m_fy[k*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)];
      frontGhost[(k*ghostDimY+j)*3+2] = domain->m_fz[k*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)];
      backGhost[(k*ghostDimY+j)*3]   = domain->m_fx[k*ghostDimY*ghostDimX+j*ghostDimX];
      backGhost[(k*ghostDimY+j)*3+1] = domain->m_fy[k*ghostDimY*ghostDimX+j*ghostDimX];
      backGhost[(k*ghostDimY+j)*3+2] = domain->m_fz[k*ghostDimY*ghostDimX+j*ghostDimX];
    }
  }
  // xz plane
  for (Index_t k = 0; k < ghostDimZ; k++) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      rightGhost[(k*ghostDimX+i)*3]   = domain->m_fx[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i];
      rightGhost[(k*ghostDimX+i)*3+1] = domain->m_fy[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i];
      rightGhost[(k*ghostDimX+i)*3+2] = domain->m_fz[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i];
      leftGhost[(k*ghostDimX+i)*3]   = domain->m_fx[k*ghostDimY*ghostDimX+i];
      leftGhost[(k*ghostDimX+i)*3+1] = domain->m_fy[k*ghostDimY*ghostDimX+i];
      leftGhost[(k*ghostDimX+i)*3+2] = domain->m_fz[k*ghostDimY*ghostDimX+i];
    }
  }
  // xy plane
  for (Index_t j = 0; j < ghostDimY; j++) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      upGhost[(j*ghostDimX+i)*3] = domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+i];
      upGhost[(j*ghostDimX+i)*3+1] = domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+i];
      upGhost[(j*ghostDimX+i)*3+2] = domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+i];
      downGhost[(j*ghostDimX+i)*3] = domain->m_fx[j*ghostDimX+i];
      downGhost[(j*ghostDimX+i)*3+1] = domain->m_fy[j*ghostDimX+i];
      downGhost[(j*ghostDimX+i)*3+2] = domain->m_fz[j*ghostDimX+i];
    }
  }
  // corners
  if (thisIndex.x != chareDimX-1 && thisIndex.y != chareDimY-1 && thisIndex.z != chareDimZ-1) {
    cornerGhost[0] = domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimY-1)*ghostDimX +(ghostDimX-1)];
    cornerGhost[1] = domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimY-1)*ghostDimX +(ghostDimX-1)];
    cornerGhost[2] = domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimY-1)*ghostDimX +(ghostDimX-1)];
    thisProxy(thisIndex.x+1,thisIndex.y+1,thisIndex.z+1)
      .receiveNodeGhosts(CORNER_FRU,1,1,cornerGhost);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.y != chareDimY-1 && thisIndex.z != 0) {
    cornerGhost[0] = domain->m_fx[(ghostDimY-1)*ghostDimX +(ghostDimX-1)];
    cornerGhost[1] = domain->m_fy[(ghostDimY-1)*ghostDimX +(ghostDimX-1)];
    cornerGhost[2] = domain->m_fz[(ghostDimY-1)*ghostDimX +(ghostDimX-1)];
    thisProxy(thisIndex.x+1,thisIndex.y+1,thisIndex.z-1)
      .receiveNodeGhosts(CORNER_FRD,1,1,cornerGhost);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.y != 0 && thisIndex.z != chareDimZ-1) {
    cornerGhost[0] = domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimX-1)];
    cornerGhost[1] = domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimX-1)];
    cornerGhost[2] = domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimX-1)];
    thisProxy(thisIndex.x+1,thisIndex.y-1,thisIndex.z+1)
      .receiveNodeGhosts(CORNER_FLU,1,1,cornerGhost);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.y != 0 && thisIndex.z != 0) {
    cornerGhost[0] = domain->m_fx[(ghostDimX-1)];
    cornerGhost[1] = domain->m_fy[(ghostDimX-1)];
    cornerGhost[2] = domain->m_fz[(ghostDimX-1)];
    thisProxy(thisIndex.x+1,thisIndex.y-1,thisIndex.z-1)
      .receiveNodeGhosts(CORNER_FLD,1,1,cornerGhost);
  }
  if (thisIndex.x != 0 && thisIndex.y != chareDimY-1 && thisIndex.z != chareDimZ-1) {
    cornerGhost[0] = domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimY-1)*ghostDimX];
    cornerGhost[1] = domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimY-1)*ghostDimX];
    cornerGhost[2] = domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimY-1)*ghostDimX];
    thisProxy(thisIndex.x-1,thisIndex.y+1,thisIndex.z+1)
      .receiveNodeGhosts(CORNER_BRU,1,1,cornerGhost);
  }
  if (thisIndex.x != 0 && thisIndex.y != chareDimY-1 && thisIndex.z != 0) {
    cornerGhost[0] = domain->m_fx[(ghostDimY-1)*ghostDimX];
    cornerGhost[1] = domain->m_fy[(ghostDimY-1)*ghostDimX];
    cornerGhost[2] = domain->m_fz[(ghostDimY-1)*ghostDimX];
    thisProxy(thisIndex.x-1,thisIndex.y+1,thisIndex.z-1)
      .receiveNodeGhosts(CORNER_BRD,1,1,cornerGhost);
  }
  if (thisIndex.x != 0 && thisIndex.y != 0 && thisIndex.z != chareDimZ-1) {
    cornerGhost[0] = domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX];
    cornerGhost[1] = domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX];
    cornerGhost[2] = domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX];
    thisProxy(thisIndex.x-1,thisIndex.y-1,thisIndex.z+1)
      .receiveNodeGhosts(CORNER_BLU,1,1,cornerGhost);
  }
  if (thisIndex.x != 0 && thisIndex.y != 0 && thisIndex.z != 0) {
    cornerGhost[0] = domain->m_fx[0];
    cornerGhost[1] = domain->m_fy[0];
    cornerGhost[2] = domain->m_fz[0];
    thisProxy(thisIndex.x-1,thisIndex.y-1,thisIndex.z-1)
      .receiveNodeGhosts(CORNER_BLD,1,1,cornerGhost);
  }
  // z edges
  if (thisIndex.x != chareDimX-1 && thisIndex.y != chareDimY-1) {
    for (Index_t k = 0; k < ghostDimZ; k++) {
      zEdge[k*3]   = domain->m_fx[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+(ghostDimX-1)];
      zEdge[k*3+1] = domain->m_fy[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+(ghostDimX-1)];
      zEdge[k*3+2] = domain->m_fz[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+(ghostDimX-1)];
    }
    thisProxy(thisIndex.x+1, thisIndex.y+1, thisIndex.z)
      .receiveNodeGhosts(EDGE_PXPY,1,ghostDimZ,zEdge);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.y != 0) {
    for (Index_t k = 0; k < ghostDimZ; k++) {
      zEdge[k*3]   = domain->m_fx[k*ghostDimY*ghostDimX+(ghostDimX-1)];
      zEdge[k*3+1] = domain->m_fy[k*ghostDimY*ghostDimX+(ghostDimX-1)];
      zEdge[k*3+2] = domain->m_fz[k*ghostDimY*ghostDimX+(ghostDimX-1)];
    }
    thisProxy(thisIndex.x+1, thisIndex.y-1, thisIndex.z)
      .receiveNodeGhosts(EDGE_PXMY,1,ghostDimZ,zEdge);
  }
  if (thisIndex.x != 0 && thisIndex.y != chareDimY-1) {
    for (Index_t k = 0; k < ghostDimZ; k++) {
      zEdge[k*3]   = domain->m_fx[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX];
      zEdge[k*3+1] = domain->m_fy[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX];
      zEdge[k*3+2] = domain->m_fz[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX];
    }
    thisProxy(thisIndex.x-1, thisIndex.y+1, thisIndex.z)
      .receiveNodeGhosts(EDGE_MXPY,1,ghostDimZ,zEdge);
  }
  if (thisIndex.x != 0 && thisIndex.y != 0) {
    for (Index_t k = 0; k < ghostDimZ; k++) {
      zEdge[k*3]   = domain->m_fx[k*ghostDimY*ghostDimX];
      zEdge[k*3+1] = domain->m_fy[k*ghostDimY*ghostDimX];
      zEdge[k*3+2] = domain->m_fz[k*ghostDimY*ghostDimX];
    }
    thisProxy(thisIndex.x-1, thisIndex.y-1, thisIndex.z)
      .receiveNodeGhosts(EDGE_MXMY,1,ghostDimZ,zEdge);
  }
  // y edges
  if (thisIndex.x != chareDimX-1 && thisIndex.z != chareDimZ-1) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      yEdge[j*3]   = domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)];
      yEdge[j*3+1] = domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)];
      yEdge[j*3+2] = domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)];
    }
    thisProxy(thisIndex.x+1, thisIndex.y, thisIndex.z+1)
      .receiveNodeGhosts(EDGE_PXPZ,1,ghostDimY,yEdge);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.z != 0) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      yEdge[j*3]   = domain->m_fx[j*ghostDimX+(ghostDimX-1)];
      yEdge[j*3+1] = domain->m_fy[j*ghostDimX+(ghostDimX-1)];
      yEdge[j*3+2] = domain->m_fz[j*ghostDimX+(ghostDimX-1)];
    }
    thisProxy(thisIndex.x+1, thisIndex.y, thisIndex.z-1)
      .receiveNodeGhosts(EDGE_PXMZ,1,ghostDimY,yEdge);
  }
  if (thisIndex.x != 0 && thisIndex.z != chareDimZ-1) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      yEdge[j*3]   = domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX];
      yEdge[j*3+1] = domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX];
      yEdge[j*3+2] = domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX];
    }
    thisProxy(thisIndex.x-1, thisIndex.y, thisIndex.z+1)
      .receiveNodeGhosts(EDGE_MXPZ,1,ghostDimY,yEdge);
  }
  if (thisIndex.x != 0 && thisIndex.z != 0) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      yEdge[j*3]   = domain->m_fx[j*ghostDimX];
      yEdge[j*3+1] = domain->m_fy[j*ghostDimX];
      yEdge[j*3+2] = domain->m_fz[j*ghostDimX];
    }
    thisProxy(thisIndex.x-1, thisIndex.y, thisIndex.z-1)
      .receiveNodeGhosts(EDGE_MXMZ,1,ghostDimY,yEdge);
  }
  // x edges
  if (thisIndex.y != chareDimY-1 && thisIndex.z != chareDimZ-1) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      xEdge[i*3]   = domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i];
      xEdge[i*3+1] = domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i];
      xEdge[i*3+2] = domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i];
    }
    thisProxy(thisIndex.x, thisIndex.y+1, thisIndex.z+1)
      .receiveNodeGhosts(EDGE_PYPZ,1,ghostDimX,xEdge);
  }
  if (thisIndex.y != chareDimY-1 && thisIndex.z != 0) { 
    for (Index_t i = 0; i < ghostDimX; i++) {
      xEdge[i*3]   = domain->m_fx[(ghostDimY-1)*ghostDimX+i];
      xEdge[i*3+1] = domain->m_fy[(ghostDimY-1)*ghostDimX+i];
      xEdge[i*3+2] = domain->m_fz[(ghostDimY-1)*ghostDimX+i];
    }
    thisProxy(thisIndex.x, thisIndex.y+1, thisIndex.z-1)
      .receiveNodeGhosts(EDGE_PYMZ,1,ghostDimX,xEdge);
  }
  if (thisIndex.y != 0 && thisIndex.z != chareDimZ-1) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      xEdge[i*3]   = domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX+i];
      xEdge[i*3+1] = domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX+i];
      xEdge[i*3+2] = domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX+i];
    }
    thisProxy(thisIndex.x, thisIndex.y-1, thisIndex.z+1)
      .receiveNodeGhosts(EDGE_MYPZ,1,ghostDimX,xEdge);
  }
  if (thisIndex.y != 0 && thisIndex.z != 0) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      xEdge[i*3]   = domain->m_fx[i];
      xEdge[i*3+1] = domain->m_fy[i];
      xEdge[i*3+2] = domain->m_fz[i];
    }
    thisProxy(thisIndex.x, thisIndex.y-1, thisIndex.z-1)
      .receiveNodeGhosts(EDGE_MYMZ,1,ghostDimX,xEdge);
  }

  // Send Ghost Data faces (direction sent)
  // front (+x)
  if (thisIndex.x != chareDimX-1) {
    thisProxy(thisIndex.x+1,thisIndex.y,thisIndex.z)
      .receiveNodeGhosts(DIR_FRONT,ghostDimY,ghostDimZ,frontGhost);
  }
  // back  (-x)
  if (thisIndex.x != 0) {
    thisProxy(thisIndex.x-1,thisIndex.y,thisIndex.z)
      .receiveNodeGhosts(DIR_BACK,ghostDimY,ghostDimZ,backGhost);
  }
  // right (+y)
  if (thisIndex.y != chareDimY-1) {
    thisProxy(thisIndex.x,thisIndex.y+1,thisIndex.z)
      .receiveNodeGhosts(DIR_RIGHT,ghostDimX,ghostDimZ,rightGhost);
  }
  // left  (-y)
  if (thisIndex.y != 0) {
    thisProxy(thisIndex.x,thisIndex.y-1,thisIndex.z)
      .receiveNodeGhosts(DIR_LEFT,ghostDimX,ghostDimZ,leftGhost);
  }
  // up    (+z)
  if (thisIndex.z != chareDimZ-1) {
    thisProxy(thisIndex.x,thisIndex.y,thisIndex.z+1)
      .receiveNodeGhosts(DIR_UP,ghostDimX,ghostDimY,upGhost);
  }
  // down  (-z)
  if (thisIndex.z != 0) {
    thisProxy(thisIndex.x,thisIndex.y,thisIndex.z-1)
      .receiveNodeGhosts(DIR_DOWN,ghostDimX,ghostDimY,downGhost);
  }

#if CHARM_DEBUG
  CkPrintf("(%d,%d,%d) Node Ghosts Sent %d/%d\n",
           thisIndex.x, thisIndex.y, thisIndex.z, ghostNodeCount, numNodeGhosts);
#endif
  // Delete Ghost Data
  delete[] cornerGhost;
  delete[] zEdge;
  delete[] yEdge;
  delete[] xEdge;
  delete[] downGhost;
  delete[] upGhost;
  delete[] leftGhost;
  delete[] rightGhost;
  delete[] backGhost;
  delete[] frontGhost;
}

void
Domain::processNodeGhosts(int dir, int width, int height,
                      Real_t ghostData[]) {
  // Process Ghost Data (direction received)
  switch (dir) {
  // front (-x)
    case DIR_FRONT:
      for (Index_t j = 0; j < 3*height*width; j++) {
        domain->ng_front[j] += ghostData[j];
      }
      break;
  // back  (+x)
    case DIR_BACK:
      for (Index_t j = 0; j < 3*height*width; j++) {
        domain->ng_back[j] += ghostData[j];
      }
      break;
  // right (-y)
    case DIR_RIGHT:
      for (Index_t j = 0; j < 3*height*width; j++) {
        domain->ng_right[j] += ghostData[j];
      }
      break;
  // left  (+y)
    case DIR_LEFT:
      for (Index_t j = 0; j < 3*height*width; j++) {
        domain->ng_left[j] += ghostData[j];
      }
      break;
  // up    (-z)
    case DIR_UP:
      for (Index_t j = 0; j < 3*height*width; j++) {
        domain->ng_up[j] += ghostData[j];
      }
      break;
  // down  (+z)
    case DIR_DOWN:
      for (Index_t j = 0; j < 3*height*width; j++) {
        domain->ng_down[j] += ghostData[j];
      }
      break;
  // Corners
    case CORNER_FRU:
      domain->ng_front[0] += ghostData[0];
      domain->ng_front[1] += ghostData[1];
      domain->ng_front[2] += ghostData[2];
      break;
    case CORNER_FRD:
      domain->ng_front[((ghostDimZ-1)*ghostDimY)*3] += ghostData[0];
      domain->ng_front[((ghostDimZ-1)*ghostDimY)*3+1] += ghostData[1];
      domain->ng_front[((ghostDimZ-1)*ghostDimY)*3+2] += ghostData[2];
      break;
    case CORNER_FLU:
      domain->ng_front[((ghostDimY-1))*3] += ghostData[0];
      domain->ng_front[((ghostDimY-1))*3+1] += ghostData[1];
      domain->ng_front[((ghostDimY-1))*3+2] += ghostData[2];
      break;
    case CORNER_FLD:
      domain->ng_front[((ghostDimZ-1)*ghostDimY+(ghostDimY-1))*3] += ghostData[0];
      domain->ng_front[((ghostDimZ-1)*ghostDimY+(ghostDimY-1))*3+1] += ghostData[1];
      domain->ng_front[((ghostDimZ-1)*ghostDimY+(ghostDimY-1))*3+2] += ghostData[2];
      break;
    case CORNER_BRU:
      domain->ng_back[0] += ghostData[0];
      domain->ng_back[1] += ghostData[1];
      domain->ng_back[2] += ghostData[2];
      break;
    case CORNER_BRD:
      domain->ng_back[((ghostDimZ-1)*ghostDimY)*3] += ghostData[0];
      domain->ng_back[((ghostDimZ-1)*ghostDimY)*3+1] += ghostData[1];
      domain->ng_back[((ghostDimZ-1)*ghostDimY)*3+2] += ghostData[2];
      break;
    case CORNER_BLU:
      domain->ng_back[((ghostDimY-1))*3] += ghostData[0];
      domain->ng_back[((ghostDimY-1))*3+1] += ghostData[1];
      domain->ng_back[((ghostDimY-1))*3+2] += ghostData[2];
      break;
    case CORNER_BLD:
      domain->ng_back[((ghostDimZ-1)*ghostDimY+(ghostDimY-1))*3] += ghostData[0];
      domain->ng_back[((ghostDimZ-1)*ghostDimY+(ghostDimY-1))*3+1] += ghostData[1];
      domain->ng_back[((ghostDimZ-1)*ghostDimY+(ghostDimY-1))*3+2] += ghostData[2];
      break;
  // z edges
    case EDGE_PXPY:
      for (Index_t k = 0; k < ghostDimZ; k++) {
        domain->ng_front[(k*ghostDimY)*3] += ghostData[k*3];
        domain->ng_front[(k*ghostDimY)*3+1] += ghostData[k*3+1];
        domain->ng_front[(k*ghostDimY)*3+2] += ghostData[k*3+2];
      }
      break;
    case EDGE_PXMY:
      for (Index_t k = 0; k < ghostDimZ; k++) {
        domain->ng_front[(k*ghostDimY+(ghostDimY-1))*3] += ghostData[k*3];
        domain->ng_front[(k*ghostDimY+(ghostDimY-1))*3+1] += ghostData[k*3+1];
        domain->ng_front[(k*ghostDimY+(ghostDimY-1))*3+2] += ghostData[k*3+2];
      }
      break;
    case EDGE_MXPY:
      for (Index_t k = 0; k < ghostDimZ; k++) {
        domain->ng_back[(k*ghostDimY)*3] += ghostData[k*3];
        domain->ng_back[(k*ghostDimY)*3+1] += ghostData[k*3+1];
        domain->ng_back[(k*ghostDimY)*3+2] += ghostData[k*3+2];
      }
      break;
    case EDGE_MXMY:
      for (Index_t k = 0; k < ghostDimZ; k++) {
        domain->ng_back[(k*ghostDimY+(ghostDimY-1))*3] += ghostData[k*3];
        domain->ng_back[(k*ghostDimY+(ghostDimY-1))*3+1] += ghostData[k*3+1];
        domain->ng_back[(k*ghostDimY+(ghostDimY-1))*3+2] += ghostData[k*3+2];
      }
      break;
  // y edges
    case EDGE_PXPZ:
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->ng_front[j*3] += ghostData[j*3];
        domain->ng_front[j*3+1] += ghostData[j*3+1];
        domain->ng_front[j*3+2] += ghostData[j*3+2];
      }
      break;
    case EDGE_PXMZ:
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->ng_front[((ghostDimZ-1)*ghostDimY+j)*3] += ghostData[j*3];
        domain->ng_front[((ghostDimZ-1)*ghostDimY+j)*3+1] += ghostData[j*3+1];
        domain->ng_front[((ghostDimZ-1)*ghostDimY+j)*3+2] += ghostData[j*3+2];
      }
      break;
    case EDGE_MXPZ:
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->ng_back[j*3] += ghostData[j*3];
        domain->ng_back[j*3+1] += ghostData[j*3+1];
        domain->ng_back[j*3+2] += ghostData[j*3+2];
      }
      break;
    case EDGE_MXMZ:
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->ng_back[((ghostDimZ-1)*ghostDimY+j)*3] += ghostData[j*3];
        domain->ng_back[((ghostDimZ-1)*ghostDimY+j)*3+1] += ghostData[j*3+1];
        domain->ng_back[((ghostDimZ-1)*ghostDimY+j)*3+2] += ghostData[j*3+2];
      }
      break;
  // x edges
    case EDGE_PYPZ:
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->ng_right[i*3] += ghostData[i*3];
        domain->ng_right[i*3+1] += ghostData[i*3+1];
        domain->ng_right[i*3+2] += ghostData[i*3+2];
      }
      break;
    case EDGE_PYMZ:
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->ng_right[((ghostDimZ-1)*ghostDimX+i)*3] += ghostData[i*3];
        domain->ng_right[((ghostDimZ-1)*ghostDimX+i)*3+1] += ghostData[i*3+1];
        domain->ng_right[((ghostDimZ-1)*ghostDimX+i)*3+2] += ghostData[i*3+2];
      }
      break;
    case EDGE_MYPZ:
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->ng_left[i*3] += ghostData[i*3];
        domain->ng_left[i*3+1] += ghostData[i*3+1];
        domain->ng_left[i*3+2] += ghostData[i*3+2];
      }
      break;
    case EDGE_MYMZ:
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->ng_left[((ghostDimZ-1)*ghostDimX+i)*3] += ghostData[i*3];
        domain->ng_left[((ghostDimZ-1)*ghostDimX+i)*3+1] += ghostData[i*3+1];
        domain->ng_left[((ghostDimZ-1)*ghostDimX+i)*3+2] += ghostData[i*3+2];
      }
      break;
    default:
      CkAbort("Error: Node Ghost Direction\n");
  } 
#if CHARM_DEBUG
  CkPrintf("(%d,%d,%d): [%d] Nodes[%d] %d/%d\n",
           thisIndex.x, thisIndex.y, thisIndex.z, 
            thisIndex.z*chareDimY*chareDimX + thisIndex.y*chareDimX + thisIndex.x, dir,
           ghostNodeCount+1, numNodeGhosts);
#endif
}

void
Domain::updateForceGhosts() {
    // Add in temporary force values
    for (Index_t k = 0; k < ghostDimZ; k++) {
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->m_fx[k*ghostDimY*ghostDimX+j*ghostDimX] += domain->ng_front[(k*ghostDimY+j)*3];
        domain->m_fy[k*ghostDimY*ghostDimX+j*ghostDimX] += domain->ng_front[(k*ghostDimY+j)*3+1];
        domain->m_fz[k*ghostDimY*ghostDimX+j*ghostDimX] += domain->ng_front[(k*ghostDimY+j)*3+2];
        domain->m_fx[k*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)] += domain->ng_back[(k*ghostDimY+j)*3];
        domain->m_fy[k*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)] += domain->ng_back[(k*ghostDimY+j)*3+1];
        domain->m_fz[k*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)] += domain->ng_back[(k*ghostDimY+j)*3+2];
      }
    }
    for (Index_t k = 0; k < ghostDimZ; k++) {
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->m_fx[k*ghostDimY*ghostDimX+i] += domain->ng_right[(k*ghostDimX+i)*3];
        domain->m_fy[k*ghostDimY*ghostDimX+i] += domain->ng_right[(k*ghostDimX+i)*3+1];
        domain->m_fz[k*ghostDimY*ghostDimX+i] += domain->ng_right[(k*ghostDimX+i)*3+2];
        domain->m_fx[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i] += domain->ng_left[(k*ghostDimX+i)*3];
        domain->m_fy[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i] += domain->ng_left[(k*ghostDimX+i)*3+1];
        domain->m_fz[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i] += domain->ng_left[(k*ghostDimX+i)*3+2];
      }
    }
    for (Index_t j = 0; j < ghostDimY; j++) {
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->m_fx[j*ghostDimX+i] += domain->ng_up[(j*ghostDimX+i)*3];
        domain->m_fy[j*ghostDimX+i] += domain->ng_up[(j*ghostDimX+i)*3+1];
        domain->m_fz[j*ghostDimX+i] += domain->ng_up[(j*ghostDimX+i)*3+2];
        domain->m_fx[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+i] += domain->ng_down[(j*ghostDimX+i)*3];
        domain->m_fy[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+i] += domain->ng_down[(j*ghostDimX+i)*3+1];
        domain->m_fz[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+i] += domain->ng_down[(j*ghostDimX+i)*3+2];
      }
    }
    // Reset ghost count and node ghosts
    // for next iteration
    for (Index_t i = 0; i < 3*ghostDimZ*ghostDimY; i++) {
      domain->ng_front[i] = Real_t(0.0);
      domain->ng_back[i] = Real_t(0.0);
    }
    for (Index_t i = 0; i < 3*ghostDimZ*ghostDimX; i++) {
      domain->ng_right[i] = Real_t(0.0);
      domain->ng_left[i] = Real_t(0.0);
    }
    for (Index_t i = 0; i < 3*ghostDimY*ghostDimX; i++) {
      domain->ng_up[i] = Real_t(0.0);
      domain->ng_down[i] = Real_t(0.0);
    }
}

void
Domain::sendElemGhosts() {
  Real_t *frontGhost = new Real_t[blockDimY*blockDimZ*3];
  Real_t *backGhost = new Real_t[blockDimY*blockDimZ*3];
  Real_t *rightGhost = new Real_t[blockDimX*blockDimZ*3];
  Real_t *leftGhost = new Real_t[blockDimX*blockDimZ*3];
  Real_t *upGhost = new Real_t[blockDimX*blockDimY*3];
  Real_t *downGhost = new Real_t[blockDimX*blockDimY*3];
  // Create Elem Ghost data
  // yz plane
  for (Index_t k = 0; k < blockDimZ; k++) {
    for (Index_t j = 0; j < blockDimY; j++) {
      frontGhost[                        k*blockDimY + j] =
        domain->m_delv_xi[k*blockDimY*blockDimX + j*blockDimX + (blockDimX-1)];
      frontGhost[  blockDimZ*blockDimY + k*blockDimY + j] =
        domain->m_delv_eta[k*blockDimY*blockDimX + j*blockDimX + (blockDimX-1)];
      frontGhost[2*blockDimZ*blockDimY + k*blockDimY + j] =
        domain->m_delv_zeta[k*blockDimY*blockDimX + j*blockDimX + (blockDimX-1)];
      backGhost[                        k*blockDimY + j] =
        domain->m_delv_xi[k*blockDimY*blockDimX + j*blockDimX];
      backGhost[  blockDimZ*blockDimY + k*blockDimY + j] =
        domain->m_delv_eta[k*blockDimY*blockDimX + j*blockDimX];
      backGhost[2*blockDimZ*blockDimY + k*blockDimY + j] =
        domain->m_delv_zeta[k*blockDimY*blockDimX + j*blockDimX];
    }
  }
  // xz plane
  for (Index_t k = 0; k < blockDimZ; k++) {
    for (Index_t i = 0; i < blockDimX; i++) {
      rightGhost[                        k*blockDimX + i] =
        domain->m_delv_xi[k*blockDimY*blockDimX + (blockDimY-1)*blockDimX + i];
      rightGhost[  blockDimZ*blockDimX + k*blockDimX + i] =
        domain->m_delv_eta[k*blockDimY*blockDimX + (blockDimY-1)*blockDimX + i];
      rightGhost[2*blockDimZ*blockDimX + k*blockDimX + i] =
        domain->m_delv_zeta[k*blockDimY*blockDimX + (blockDimY-1)*blockDimX + i];
      leftGhost[                        k*blockDimX + i] =
        domain->m_delv_xi[k*blockDimY*blockDimX + i];
      leftGhost[  blockDimZ*blockDimX + k*blockDimX + i] =
        domain->m_delv_eta[k*blockDimY*blockDimX + i];
      leftGhost[2*blockDimZ*blockDimX + k*blockDimX + i] =
        domain->m_delv_zeta[k*blockDimY*blockDimX + i];
    }
  }
  // xy plane
  for (Index_t j = 0; j < blockDimY; j++) {
    for (Index_t i = 0; i < blockDimX; i++) {
      upGhost[                        j*blockDimX + i] =
        domain->m_delv_xi[(blockDimZ-1)*blockDimY*blockDimX + j*blockDimX + i];
      upGhost[  blockDimY*blockDimX + j*blockDimX + i] =
        domain->m_delv_eta[(blockDimZ-1)*blockDimY*blockDimX + j*blockDimX + i];
      upGhost[2*blockDimY*blockDimX + j*blockDimX + i] =
        domain->m_delv_zeta[(blockDimZ-1)*blockDimY*blockDimX + j*blockDimX + i];
      downGhost[                        j*blockDimX + i] =
        domain->m_delv_xi[j*blockDimX + i];
      downGhost[  blockDimY*blockDimX + j*blockDimX + i] =
        domain->m_delv_eta[j*blockDimX + i];
      downGhost[2*blockDimY*blockDimX + j*blockDimX + i] =
        domain->m_delv_zeta[j*blockDimX + i];
    }
  }
  // Send Ghost Data (direction sent)
  // front (+x)
  if (thisIndex.x != chareDimX-1) {
    thisProxy(thisIndex.x+1,thisIndex.y,thisIndex.z)
      .receiveElemGhosts(DIR_FRONT,blockDimY,blockDimZ,frontGhost);
  }
  // back  (-x)
  if (thisIndex.x != 0) {
    thisProxy(thisIndex.x-1,thisIndex.y,thisIndex.z)
      .receiveElemGhosts(DIR_BACK,blockDimY,blockDimZ,backGhost);
  }
  // right (+y)
  if (thisIndex.y != chareDimY-1) {
    thisProxy(thisIndex.x,thisIndex.y+1,thisIndex.z)
      .receiveElemGhosts(DIR_RIGHT,blockDimX,blockDimZ,rightGhost);
  }
  // left  (-y)
  if (thisIndex.y != 0) {
    thisProxy(thisIndex.x,thisIndex.y-1,thisIndex.z)
      .receiveElemGhosts(DIR_LEFT,blockDimX,blockDimZ,leftGhost);
  }
  // up    (+z)
  if (thisIndex.z != chareDimZ-1) {
    thisProxy(thisIndex.x,thisIndex.y,thisIndex.z+1)
      .receiveElemGhosts(DIR_UP,blockDimX,blockDimY,upGhost);
  }
  // down  (-z)
  if (thisIndex.z != 0) {
    thisProxy(thisIndex.x,thisIndex.y,thisIndex.z-1)
      .receiveElemGhosts(DIR_DOWN,blockDimX,blockDimY,downGhost);
  }

#if CHARM_DEBUG
  CkPrintf("(%d,%d,%d) Elem Ghosts Sent %d/%d\n",
           thisIndex.x, thisIndex.y, thisIndex.z, ghostElemCount, numElemGhosts);
#endif
  // Delete Ghost Data
  delete[] downGhost;
  delete[] upGhost;
  delete[] leftGhost;
  delete[] rightGhost;
  delete[] backGhost;
  delete[] frontGhost;
}

void
Domain::processElemGhosts(int dir, int width, int height,
    Real_t ghostData[]) {
  // Process Ghost Data (direction received)
  switch (dir) {
    // front (-x)
    case DIR_FRONT:
      for (Index_t k = 0; k < blockDimZ; k++) {
        for (Index_t j = 0; j < blockDimY; j++) {
          domain->m_delv_xi[backOffset + k*blockDimY+j] = ghostData[k*blockDimY+j];
          domain->m_delv_eta[backOffset + k*blockDimY+j] = ghostData[blockDimZ*blockDimY+k*blockDimY+j];
          domain->m_delv_zeta[backOffset + k*blockDimY+j] = ghostData[2*blockDimZ*blockDimY+k*blockDimY+j];
        }
      }
      break;
      // back  (+x)
    case DIR_BACK:
      for (Index_t k = 0; k < blockDimZ; k++) {
        for (Index_t j = 0; j < blockDimY; j++) {
          domain->m_delv_xi[frontOffset + k*blockDimY+j] = ghostData[k*blockDimY+j];
          domain->m_delv_eta[frontOffset + k*blockDimY+j] = ghostData[blockDimZ*blockDimY+k*blockDimY+j];
          domain->m_delv_zeta[frontOffset + k*blockDimY+j] = ghostData[2*blockDimZ*blockDimY+k*blockDimY+j];
        }
      }
      break;
      // right (-y)
    case DIR_RIGHT:
      for (Index_t k = 0; k < blockDimZ; k++) {
        for (Index_t i = 0; i < blockDimX; i++) {
          domain->m_delv_xi[leftOffset + k*blockDimX+i] = ghostData[k*blockDimX+i];
          domain->m_delv_eta[leftOffset + k*blockDimX+i] = ghostData[blockDimZ*blockDimX+k*blockDimX+i];
          domain->m_delv_zeta[leftOffset + k*blockDimX+i] = ghostData[2*blockDimZ*blockDimX+k*blockDimX+i];
        }
      }
      break;
      // left  (+y)
    case DIR_LEFT:
      for (Index_t k = 0; k < blockDimZ; k++) {
        for (Index_t i = 0; i < blockDimX; i++) {
          domain->m_delv_xi[rightOffset + k*blockDimX+i] = ghostData[k*blockDimX+i];
          domain->m_delv_eta[rightOffset + k*blockDimX+i] = ghostData[blockDimZ*blockDimX+k*blockDimX+i];
          domain->m_delv_zeta[rightOffset + k*blockDimX+i] = ghostData[2*blockDimZ*blockDimX+k*blockDimX+i];
        }
      }
      break;
      // up    (-z)
    case DIR_UP:
      for (Index_t j = 0; j < blockDimY; j++) {
        for (Index_t i = 0; i < blockDimX; i++) {
          domain->m_delv_xi[downOffset + j*blockDimX+i] = ghostData[j*blockDimX+i];
          domain->m_delv_eta[downOffset + j*blockDimX+i] = ghostData[blockDimY*blockDimX+j*blockDimX+i];
          domain->m_delv_zeta[downOffset + j*blockDimX+i] = ghostData[2*blockDimY*blockDimX+j*blockDimX+i];
        }
      }
      break;
      // down  (+z)
    case DIR_DOWN:
      for (Index_t j = 0; j < blockDimY; j++) {
        for (Index_t i = 0; i < blockDimX; i++) {
          domain->m_delv_xi[upOffset + j*blockDimX+i] = ghostData[j*blockDimX+i];
          domain->m_delv_eta[upOffset + j*blockDimX+i] = ghostData[blockDimY*blockDimX+j*blockDimX+i];
          domain->m_delv_zeta[upOffset + j*blockDimX+i] = ghostData[2*blockDimY*blockDimX+j*blockDimX+i];
        }
      }
      break;
    default:
      CkAbort("Error: Elem Ghost Direction\n");
  } 
#if CHARM_DEBUG
  CkPrintf("(%d,%d,%d): Elems[%d] %d/%d\n",
           thisIndex.x, thisIndex.y, thisIndex.z, dir,
           ghostElemCount+1, numElemGhosts);
#endif
}

void
Domain::sendNodalMass() {
  Real_t *frontGhost = new Real_t[ghostDimY*ghostDimZ];
  Real_t *backGhost = new Real_t[ghostDimY*ghostDimZ];
  Real_t *rightGhost = new Real_t[ghostDimX*ghostDimZ];
  Real_t *leftGhost = new Real_t[ghostDimX*ghostDimZ];
  Real_t *upGhost = new Real_t[ghostDimX*ghostDimY];
  Real_t *downGhost = new Real_t[ghostDimX*ghostDimY];
  Real_t *xEdge = new Real_t[ghostDimX];
  Real_t *yEdge = new Real_t[ghostDimY];
  Real_t *zEdge = new Real_t[ghostDimZ];
  Real_t *cornerGhost = new Real_t[1];
  // Create Ghost data
  // yz plane
  for (Index_t k = 0; k < ghostDimZ; k++) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      frontGhost[k*ghostDimY+j] = domain->m_nodalMass[k*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)];
      backGhost[k*ghostDimY+j] = domain->m_nodalMass[k*ghostDimY*ghostDimX+j*ghostDimX];
    }
  }
  // xz plane
  for (Index_t k = 0; k < ghostDimZ; k++) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      rightGhost[k*ghostDimX+i] = domain->m_nodalMass[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i];
      leftGhost[k*ghostDimX+i] = domain->m_nodalMass[k*ghostDimY*ghostDimX+i];
    }
  }
  // xy plane
  for (Index_t j = 0; j < ghostDimY; j++) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      upGhost[j*ghostDimX+i] = domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+i];
      downGhost[j*ghostDimX+i] = domain->m_nodalMass[j*ghostDimX+i];
    }
  }
  // corners
  if (thisIndex.x != chareDimX-1 && thisIndex.y != chareDimY-1 && thisIndex.z != chareDimZ-1) {
    cornerGhost[0] = domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimY-1)*ghostDimX +(ghostDimX-1)];
    thisProxy(thisIndex.x+1,thisIndex.y+1,thisIndex.z+1)
      .receiveNodalMass(CORNER_FRU,1,1,cornerGhost);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.y != chareDimY-1 && thisIndex.z != 0) {
    cornerGhost[0] = domain->m_nodalMass[(ghostDimY-1)*ghostDimX +(ghostDimX-1)];
    thisProxy(thisIndex.x+1,thisIndex.y+1,thisIndex.z-1)
      .receiveNodalMass(CORNER_FRD,1,1,cornerGhost);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.y != 0 && thisIndex.z != chareDimZ-1) {
    cornerGhost[0] = domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimX-1)];
    thisProxy(thisIndex.x+1,thisIndex.y-1,thisIndex.z+1)
      .receiveNodalMass(CORNER_FLU,1,1,cornerGhost);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.y != 0 && thisIndex.z != 0) {
    cornerGhost[0] = domain->m_nodalMass[(ghostDimX-1)];
    thisProxy(thisIndex.x+1,thisIndex.y-1,thisIndex.z-1)
      .receiveNodalMass(CORNER_FLD,1,1,cornerGhost);
  }
  if (thisIndex.x != 0 && thisIndex.y != chareDimY-1 && thisIndex.z != chareDimZ-1) {
    cornerGhost[0] = domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX +(ghostDimY-1)*ghostDimX];
    thisProxy(thisIndex.x-1,thisIndex.y+1,thisIndex.z+1)
      .receiveNodalMass(CORNER_BRU,1,1,cornerGhost);
  }
  if (thisIndex.x != 0 && thisIndex.y != chareDimY-1 && thisIndex.z != 0) {
    cornerGhost[0] = domain->m_nodalMass[(ghostDimY-1)*ghostDimX];
    thisProxy(thisIndex.x-1,thisIndex.y+1,thisIndex.z-1)
      .receiveNodalMass(CORNER_BRD,1,1,cornerGhost);
  }
  if (thisIndex.x != 0 && thisIndex.y != 0 && thisIndex.z != chareDimZ-1) {
    cornerGhost[0] = domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX];
    thisProxy(thisIndex.x-1,thisIndex.y-1,thisIndex.z+1)
      .receiveNodalMass(CORNER_BLU,1,1,cornerGhost);
  }
  if (thisIndex.x != 0 && thisIndex.y != 0 && thisIndex.z != 0) {
    cornerGhost[0] = domain->m_nodalMass[0];
    thisProxy(thisIndex.x-1,thisIndex.y-1,thisIndex.z-1)
      .receiveNodalMass(CORNER_BLD,1,1,cornerGhost);
  }
  // z edges
  if (thisIndex.x != chareDimX-1 && thisIndex.y != chareDimY-1) {
    for (Index_t k = 0; k < ghostDimZ; k++) {
      zEdge[k]   = domain->m_nodalMass[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+(ghostDimX-1)];
    }
    thisProxy(thisIndex.x+1, thisIndex.y+1, thisIndex.z)
      .receiveNodalMass(EDGE_PXPY,1,ghostDimZ,zEdge);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.y != 0) {
    for (Index_t k = 0; k < ghostDimZ; k++) {
      zEdge[k]   = domain->m_nodalMass[k*ghostDimY*ghostDimX+(ghostDimX-1)];
    }
    thisProxy(thisIndex.x+1, thisIndex.y-1, thisIndex.z)
      .receiveNodalMass(EDGE_PXMY,1,ghostDimZ,zEdge);
  }
  if (thisIndex.x != 0 && thisIndex.y != chareDimY-1) {
    for (Index_t k = 0; k < ghostDimZ; k++) {
      zEdge[k]   = domain->m_nodalMass[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX];
    }
    thisProxy(thisIndex.x-1, thisIndex.y+1, thisIndex.z)
      .receiveNodalMass(EDGE_MXPY,1,ghostDimZ,zEdge);
  }
  if (thisIndex.x != 0 && thisIndex.y != 0) {
    for (Index_t k = 0; k < ghostDimZ; k++) {
      zEdge[k]   = domain->m_nodalMass[k*ghostDimY*ghostDimX];
    }
    thisProxy(thisIndex.x-1, thisIndex.y-1, thisIndex.z)
      .receiveNodalMass(EDGE_MXMY,1,ghostDimZ,zEdge);
  }
  // y edges
  if (thisIndex.x != chareDimX-1 && thisIndex.z != chareDimZ-1) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      yEdge[j]   = domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)];
    }
    thisProxy(thisIndex.x+1, thisIndex.y, thisIndex.z+1)
      .receiveNodalMass(EDGE_PXPZ,1,ghostDimY,yEdge);
  }
  if (thisIndex.x != chareDimX-1 && thisIndex.z != 0) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      yEdge[j]   = domain->m_nodalMass[j*ghostDimX+(ghostDimX-1)];
    }
    thisProxy(thisIndex.x+1, thisIndex.y, thisIndex.z-1)
      .receiveNodalMass(EDGE_PXMZ,1,ghostDimY,yEdge);
  }
  if (thisIndex.x != 0 && thisIndex.z != chareDimZ-1) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      yEdge[j]   = domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX];
    }
    thisProxy(thisIndex.x-1, thisIndex.y, thisIndex.z+1)
      .receiveNodalMass(EDGE_MXPZ,1,ghostDimY,yEdge);
  }
  if (thisIndex.x != 0 && thisIndex.z != 0) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      yEdge[j]   = domain->m_nodalMass[j*ghostDimX];
    }
    thisProxy(thisIndex.x-1, thisIndex.y, thisIndex.z-1)
      .receiveNodalMass(EDGE_MXMZ,1,ghostDimY,yEdge);
  }
  // x edges
  if (thisIndex.y != chareDimY-1 && thisIndex.z != chareDimZ-1) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      xEdge[i]   = domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i];
    }
    thisProxy(thisIndex.x, thisIndex.y+1, thisIndex.z+1)
      .receiveNodalMass(EDGE_PYPZ,1,ghostDimX,xEdge);
  }
  if (thisIndex.y != chareDimY-1 && thisIndex.z != 0) { 
    for (Index_t i = 0; i < ghostDimX; i++) {
      xEdge[i]   = domain->m_nodalMass[(ghostDimY-1)*ghostDimX+i];
    }
    thisProxy(thisIndex.x, thisIndex.y+1, thisIndex.z-1)
      .receiveNodalMass(EDGE_PYMZ,1,ghostDimX,xEdge);
  }
  if (thisIndex.y != 0 && thisIndex.z != chareDimZ-1) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      xEdge[i]   = domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX+i];
    }
    thisProxy(thisIndex.x, thisIndex.y-1, thisIndex.z+1)
      .receiveNodalMass(EDGE_MYPZ,1,ghostDimX,xEdge);
  }
  if (thisIndex.y != 0 && thisIndex.z != 0) {
    for (Index_t i = 0; i < ghostDimX; i++) {
      xEdge[i]   = domain->m_nodalMass[i];
    }
    thisProxy(thisIndex.x, thisIndex.y-1, thisIndex.z-1)
      .receiveNodalMass(EDGE_MYMZ,1,ghostDimX,xEdge);
  }

  // Send Ghost Data faces (direction sent)
  // front (+x)
  if (thisIndex.x != chareDimX-1) {
    thisProxy(thisIndex.x+1,thisIndex.y,thisIndex.z)
      .receiveNodalMass(DIR_FRONT,ghostDimY,ghostDimZ,frontGhost);
  }
  // back  (-x)
  if (thisIndex.x != 0) {
    thisProxy(thisIndex.x-1,thisIndex.y,thisIndex.z)
      .receiveNodalMass(DIR_BACK,ghostDimY,ghostDimZ,backGhost);
  }
  // right (+y)
  if (thisIndex.y != chareDimY-1) {
    thisProxy(thisIndex.x,thisIndex.y+1,thisIndex.z)
      .receiveNodalMass(DIR_RIGHT,ghostDimX,ghostDimZ,rightGhost);
  }
  // left  (-y)
  if (thisIndex.y != 0) {
    thisProxy(thisIndex.x,thisIndex.y-1,thisIndex.z)
      .receiveNodalMass(DIR_LEFT,ghostDimX,ghostDimZ,leftGhost);
  }
  // up    (+z)
  if (thisIndex.z != chareDimZ-1) {
    thisProxy(thisIndex.x,thisIndex.y,thisIndex.z+1)
      .receiveNodalMass(DIR_UP,ghostDimX,ghostDimY,upGhost);
  }
  // down  (-z)
  if (thisIndex.z != 0) {
    thisProxy(thisIndex.x,thisIndex.y,thisIndex.z-1)
      .receiveNodalMass(DIR_DOWN,ghostDimX,ghostDimY,downGhost);
  }

  // Update itself
  thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)
    .receiveNodalMass(DIR_SELF,0,0,NULL);

#if CHARM_DEBUG
  CkPrintf("(%d,%d,%d) Node Ghosts Sent %d/%d\n",
           thisIndex.x, thisIndex.y, thisIndex.z, ghostNodeCount, numNodeGhosts);
#endif
  // Delete Ghost Data
  delete[] cornerGhost;
  delete[] zEdge;
  delete[] yEdge;
  delete[] xEdge;
  delete[] downGhost;
  delete[] upGhost;
  delete[] leftGhost;
  delete[] rightGhost;
  delete[] backGhost;
  delete[] frontGhost;
}

void
Domain::receiveNodalMass(int dir, int width, int height,
                      Real_t ghostData[]) {
  // Process Ghost Data (direction received)
  switch (dir) {
    case DIR_SELF:
      break;
  // front (-x)
    case DIR_FRONT:
      for (Index_t j = 0; j < height*width; j++) {
        domain->ng_front[j] += ghostData[j];
      }
      break;
  // back  (+x)
    case DIR_BACK:
      for (Index_t j = 0; j < height*width; j++) {
        domain->ng_back[j] += ghostData[j];
      }
      break;
  // right (-y)
    case DIR_RIGHT:
      for (Index_t j = 0; j < height*width; j++) {
        domain->ng_right[j] += ghostData[j];
      }
      break;
  // left  (+y)
    case DIR_LEFT:
      for (Index_t j = 0; j < height*width; j++) {
        domain->ng_left[j] += ghostData[j];
      }
      break;
  // up    (-z)
    case DIR_UP:
      for (Index_t j = 0; j < height*width; j++) {
        domain->ng_up[j] += ghostData[j];
      }
      break;
  // down  (+z)
    case DIR_DOWN:
      for (Index_t j = 0; j < height*width; j++) {
        domain->ng_down[j] += ghostData[j];
      }
      break;
  // Corners
    case CORNER_FRU:
      domain->ng_front[0] += ghostData[0];
      break;
    case CORNER_FRD:
      domain->ng_front[(ghostDimZ-1)*ghostDimY] += ghostData[0];
      break;
    case CORNER_FLU:
      domain->ng_front[(ghostDimY-1)] += ghostData[0];
      break;
    case CORNER_FLD:
      domain->ng_front[(ghostDimZ-1)*ghostDimY+(ghostDimY-1)] += ghostData[0];
      break;
    case CORNER_BRU:
      domain->ng_back[0] += ghostData[0];
      break;
    case CORNER_BRD:
      domain->ng_back[(ghostDimZ-1)*ghostDimY] += ghostData[0];
      break;
    case CORNER_BLU:
      domain->ng_back[(ghostDimY-1)] += ghostData[0];
      break;
    case CORNER_BLD:
      domain->ng_back[(ghostDimZ-1)*ghostDimY+(ghostDimY-1)] += ghostData[0];
      break;
  // z edges
    case EDGE_PXPY:
      for (Index_t k = 0; k < ghostDimZ; k++) {
        domain->ng_front[k*ghostDimY] += ghostData[k];
      }
      break;
    case EDGE_PXMY:
      for (Index_t k = 0; k < ghostDimZ; k++) {
        domain->ng_front[k*ghostDimY+(ghostDimY-1)] += ghostData[k];
      }
      break;
    case EDGE_MXPY:
      for (Index_t k = 0; k < ghostDimZ; k++) {
        domain->ng_back[k*ghostDimY] += ghostData[k];
      }
      break;
    case EDGE_MXMY:
      for (Index_t k = 0; k < ghostDimZ; k++) {
        domain->ng_back[k*ghostDimY+(ghostDimY-1)] += ghostData[k];
      }
      break;
  // y edges
    case EDGE_PXPZ:
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->ng_front[j] += ghostData[j];
      }
      break;
    case EDGE_PXMZ:
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->ng_front[((ghostDimZ-1)*ghostDimY+j)] += ghostData[j];
      }
      break;
    case EDGE_MXPZ:
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->ng_back[j] += ghostData[j];
      }
      break;
    case EDGE_MXMZ:
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->ng_back[((ghostDimZ-1)*ghostDimY+j)] += ghostData[j];
      }
      break;
  // x edges
    case EDGE_PYPZ:
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->ng_right[i] += ghostData[i];
      }
      break;
    case EDGE_PYMZ:
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->ng_right[((ghostDimZ-1)*ghostDimX+i)] += ghostData[i];
      }
      break;
    case EDGE_MYPZ:
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->ng_left[i] += ghostData[i];
      }
      break;
    case EDGE_MYMZ:
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->ng_left[((ghostDimZ-1)*ghostDimX+i)] += ghostData[i];
      }
      break;
    default:
      CkAbort("Error: Nodal Mass Ghost Direction\n");
  } 
#if CHARM_DEBUG
  CkPrintf("(%d,%d,%d): [%d] Nodes[%d] %d/%d\n",
           thisIndex.x, thisIndex.y, thisIndex.z, 
            thisIndex.z*chareDimY*chareDimX + thisIndex.y*chareDimX + thisIndex.x, dir,
           ghostNodeCount+1, numNodeGhosts);
#endif
  // Determine resume conditions
  if ((++ghostNodeCount) >= numNodeGhosts) {
    // Add in temporary force values
    for (Index_t k = 0; k < ghostDimZ; k++) {
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->m_nodalMass[k*ghostDimY*ghostDimX+j*ghostDimX] += domain->ng_front[(k*ghostDimY+j)];
        domain->m_nodalMass[k*ghostDimY*ghostDimX+j*ghostDimX+(ghostDimX-1)] += domain->ng_back[(k*ghostDimY+j)];
      }
    }
    for (Index_t k = 0; k < ghostDimZ; k++) {
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->m_nodalMass[k*ghostDimY*ghostDimX+i] += domain->ng_right[(k*ghostDimX+i)];
        domain->m_nodalMass[k*ghostDimY*ghostDimX+(ghostDimY-1)*ghostDimX+i] += domain->ng_left[(k*ghostDimX+i)];
      }
    }
    for (Index_t j = 0; j < ghostDimY; j++) {
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->m_nodalMass[j*ghostDimX+i] += domain->ng_up[(j*ghostDimX+i)];
        domain->m_nodalMass[(ghostDimZ-1)*ghostDimY*ghostDimX+j*ghostDimX+i] += domain->ng_down[(j*ghostDimX+i)];
      }
    }
    // Reset ghost count and node ghosts
    // for next iteration
    ghostNodeCount = 0;
    for (Index_t i = 0; i < ghostDimZ*ghostDimY; i++) {
      domain->ng_front[i] = Real_t(0.0);
      domain->ng_back[i] = Real_t(0.0);
    }
    for (Index_t i = 0; i < ghostDimZ*ghostDimX; i++) {
      domain->ng_right[i] = Real_t(0.0);
      domain->ng_left[i] = Real_t(0.0);
    }
    for (Index_t i = 0; i < ghostDimY*ghostDimX; i++) {
      domain->ng_up[i] = Real_t(0.0);
      domain->ng_down[i] = Real_t(0.0);
    }
    // Check in with main when done
    // Init timing parameters
    domain->deltaTime = DT_FIXED;
    CkCallback *cb = new CkCallback(CkIndex_Main::initCheckin(NULL),mainProxy);
    contribute(sizeof(Real_t), NULL, CkReduction::nop, *cb);
  }
}


#include "lulesh.def.h"
