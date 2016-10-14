#include <vector>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lulesh.h"

extern /*readonly*/ int elemDimX;
extern /*readonly*/ int elemDimY;
extern /*readonly*/ int elemDimZ;
extern /*readonly*/ int blockDimX;
extern /*readonly*/ int blockDimY;
extern /*readonly*/ int blockDimZ;
extern /*readonly*/ int ghostDimX;
extern /*readonly*/ int ghostDimY;
extern /*readonly*/ int ghostDimZ;
extern /*readonly*/ int chareDimX;
extern /*readonly*/ int chareDimY;
extern /*readonly*/ int chareDimZ;
extern /*readonly*/ double charesPerPE;
extern /*readonly*/ int numElems;
extern /*readonly*/ int numNodes;
extern /*readonly*/ int ghostElems;
extern /*readonly*/ int frontOffset;
extern /*readonly*/ int backOffset;
extern /*readonly*/ int rightOffset;
extern /*readonly*/ int leftOffset;
extern /*readonly*/ int upOffset;
extern /*readonly*/ int downOffset;
extern /*readonly*/ Real_t stopTime;
extern /*readonly*/ int lbFrequency;
extern /*readonly*/ int cpFrequency;

// Constant array needed in Hourglass calculations 
Real_t Gamma[4][8];


// -----------------------------------------
//
// Domain Initialization

void
Domain::InitializeLuleshData() {
  domain = new Dom;
  // Temporary values/indexes
  Real_t tx, ty, tz;
  Index_t nidx, zidx;
  // Determine behavior of ghost data
  // based on location of block
  numElemGhosts = 6;
  numNodeGhosts = 26;
  ghostNodeCount = 0;
  ghostElemCount = 0;
  if (thisIndex.x == 0 || thisIndex.x == chareDimX-1) {
    numElemGhosts--;
    numNodeGhosts -= 9;
    if (thisIndex.y == 0 || thisIndex.y == chareDimY-1) {
      numElemGhosts--;
      numNodeGhosts -= 6;
      if (thisIndex.z == 0 || thisIndex.z == chareDimZ-1) {
        numElemGhosts--;
        numNodeGhosts -= 4;
      }
    }
    else if (thisIndex.z == 0 || thisIndex.z == chareDimZ-1) {
      numElemGhosts--;
      numNodeGhosts -= 6;
    }
  }
  else if (thisIndex.y == 0 || thisIndex.y == chareDimY-1) {
    numElemGhosts--;
    numNodeGhosts -= 9;
    if (thisIndex.z == 0 || thisIndex.z == chareDimZ-1) {
      numElemGhosts--;
      numNodeGhosts -= 6;
    }
  }
  else if (thisIndex.z == 0 || thisIndex.z == chareDimZ-1) {
    numElemGhosts--;
    numNodeGhosts -= 9;
  }
  // Additional handling of ghost regions
  // if discretization spans entire domain
  if (thisIndex.x == 0 && thisIndex.x == chareDimX-1) {
    numElemGhosts--;
    numNodeGhosts -= 4;
    if (thisIndex.y == 0 && thisIndex.y == chareDimY-1) {
      numElemGhosts--;
      numNodeGhosts -= 2;
      if (thisIndex.z == 0 && thisIndex.z == chareDimZ-1) {
        numElemGhosts--;
        numNodeGhosts -= 1;
      }
    }
    else if (thisIndex.z == 0 && thisIndex.z == chareDimZ-1) {
      numElemGhosts--;
      numNodeGhosts -= 2;
    }
  }
  else if (thisIndex.y == 0 && thisIndex.y == chareDimY-1) {
    numElemGhosts--;
    numNodeGhosts -= 4;
    if (thisIndex.z == 0 && thisIndex.z == chareDimZ-1) {
      numElemGhosts--;
      numNodeGhosts -= 2;
    }
  }
  else if (thisIndex.z == 0 && thisIndex.z == chareDimZ-1) {
    numElemGhosts--;
    numNodeGhosts -= 4;
  }
  domain->ng_front.resize(3*ghostDimY*ghostDimZ,0.0);
  domain->ng_back.resize(3*ghostDimY*ghostDimZ,0.0);
  domain->ng_right.resize(3*ghostDimX*ghostDimZ,0.0);
  domain->ng_left.resize(3*ghostDimX*ghostDimZ,0.0);
  domain->ng_up.resize(3*ghostDimX*ghostDimY,0.0);
  domain->ng_down.resize(3*ghostDimX*ghostDimY,0.0);
  // Node Init (persistent)
  domain->m_x.resize(numNodes, 0.0);
  domain->m_y.resize(numNodes, 0.0);
  domain->m_z.resize(numNodes, 0.0);
  domain->m_xd.resize(numNodes, 0.0);
  domain->m_yd.resize(numNodes, 0.0);
  domain->m_zd.resize(numNodes, 0.0);
  domain->m_xdd.resize(numNodes, 0.0);
  domain->m_ydd.resize(numNodes, 0.0);
  domain->m_zdd.resize(numNodes, 0.0);
  domain->m_fx.resize(numNodes, 0.0);
  domain->m_fy.resize(numNodes, 0.0);
  domain->m_fz.resize(numNodes, 0.0);
  domain->m_nodalMass.resize(numNodes, 0.0);
  // Elem Init (persistent)
  domain->m_matElemlist.resize(numElems, 0);
  domain->m_nodelist.resize(NODES_PER_ELEM*numElems, 0);
  domain->m_lxim.resize(numElems, 0);
  domain->m_lxip.resize(numElems, 0);
  domain->m_letam.resize(numElems, 0);
  domain->m_letap.resize(numElems, 0);
  domain->m_lzetam.resize(numElems, 0);
  domain->m_lzetap.resize(numElems, 0);
  domain->m_elemBC.resize(numElems, 0);
  domain->m_e.resize(numElems,0.0);
  domain->m_p.resize(numElems,0.0);
  domain->m_q.resize(numElems,0.0);
  domain->m_ql.resize(numElems,0.0);
  domain->m_qq.resize(numElems,0.0);
  domain->m_v.resize(numElems, 1.0);
  domain->m_volo.resize(numElems, 0.0);
  domain->m_delv.resize(numElems, 0.0);
  domain->m_vdov.resize(numElems, 0.0);
  domain->m_arealg.resize(numElems, 0.0);
  domain->m_ss.resize(numElems, 0.0);
  domain->m_elemMass.resize(numElems, 0.0);
  // Elem Init (temporary)
  // Only delv_xi,eta,zeta require communication
  // So only those have additional elements
  domain->m_dxx.resize(numElems, 0.0);
  domain->m_dyy.resize(numElems, 0.0);
  domain->m_dzz.resize(numElems, 0.0);
  domain->m_delv_xi.resize(ghostElems, 0.0);
  domain->m_delv_eta.resize(ghostElems, 0.0);
  domain->m_delv_zeta.resize(ghostElems, 0.0);
  domain->m_delx_xi.resize(numElems, 0.0);
  domain->m_delx_eta.resize(numElems, 0.0);
  domain->m_delx_zeta.resize(numElems, 0.0);
  domain->m_vnew.resize(numElems, 0.0);
  domain->m_determ.resize(numElems, 0.0);

  // Initialize constant array needed in Hourglass calc

  Gamma[0][0] = Real_t( 1.);
  Gamma[0][1] = Real_t( 1.);
  Gamma[0][2] = Real_t(-1.);
  Gamma[0][3] = Real_t(-1.);
  Gamma[0][4] = Real_t(-1.);
  Gamma[0][5] = Real_t(-1.);
  Gamma[0][6] = Real_t( 1.);
  Gamma[0][7] = Real_t( 1.);
  Gamma[1][0] = Real_t( 1.);
  Gamma[1][1] = Real_t(-1.);
  Gamma[1][2] = Real_t(-1.);
  Gamma[1][3] = Real_t( 1.);
  Gamma[1][4] = Real_t(-1.);
  Gamma[1][5] = Real_t( 1.);
  Gamma[1][6] = Real_t( 1.);
  Gamma[1][7] = Real_t(-1.);
  Gamma[2][0] = Real_t( 1.);
  Gamma[2][1] = Real_t(-1.);
  Gamma[2][2] = Real_t( 1.);
  Gamma[2][3] = Real_t(-1.);
  Gamma[2][4] = Real_t( 1.);
  Gamma[2][5] = Real_t(-1.);
  Gamma[2][6] = Real_t( 1.);
  Gamma[2][7] = Real_t(-1.);
  Gamma[3][0] = Real_t(-1.);
  Gamma[3][1] = Real_t( 1.);
  Gamma[3][2] = Real_t(-1.);
  Gamma[3][3] = Real_t( 1.);
  Gamma[3][4] = Real_t( 1.);
  Gamma[3][5] = Real_t(-1.);
  Gamma[3][6] = Real_t( 1.);
  Gamma[3][7] = Real_t(-1.);
  
  // Initialize Nodal Coordinates
  nidx = 0;
  tz = MAX_NODE_POSITION/chareDimZ*thisIndex.z; 
  for (Index_t k = 0; k < ghostDimZ; k++) {
    ty = MAX_NODE_POSITION/chareDimY*thisIndex.y; 
    for (Index_t j = 0; j < ghostDimY; j++) {
      tx = MAX_NODE_POSITION/chareDimX*thisIndex.x; 
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->m_x[nidx] = tx;
        domain->m_y[nidx] = ty;
        domain->m_z[nidx] = tz;
#if INIT_DEBUG
        // Print out the node identifiers and its position
        CkPrintf("{%d,%d,%d} N(%d|%d): %f %f %f\n",
            thisIndex.x, thisIndex.y, thisIndex.z, nidx,
            thisIndex.z*(elemDimY+1)*(elemDimX+1)*blockDimZ
            + thisIndex.y*(elemDimX+1)*blockDimY + thisIndex.x*blockDimX
            + k*(elemDimX+1)*(elemDimX+1) + j*(elemDimX+1) + i,
            domain->m_x[nidx], domain->m_y[nidx], domain->m_z[nidx]);
#endif
        // Update entries
        nidx++;
        tx = MAX_NODE_POSITION*(thisIndex.x*blockDimX+i+1)/elemDimX;
      }
      ty = MAX_NODE_POSITION*(thisIndex.y*blockDimY+j+1)/elemDimY;
    }
    tz = MAX_NODE_POSITION*(thisIndex.z*blockDimZ+k+1)/elemDimZ;
  }

  // Embed hexahedral elements in nodal point lattice
  // (generate the node list for each element)
  nidx = 0;
  zidx = 0;
  for (Index_t k = 0; k < blockDimZ; k++) {
    for (Index_t j = 0; j < blockDimY; j++) {
      for (Index_t i = 0; i < blockDimX; i++) {
        Index_t *localNode = &domain->m_nodelist[NODES_PER_ELEM*zidx];
        localNode[0] = nidx                                      ;
        localNode[1] = nidx                                   + 1;
        localNode[2] = nidx                       + ghostDimX + 1;
        localNode[3] = nidx                       + ghostDimX    ;
        localNode[4] = nidx + ghostDimY*ghostDimX                ;
        localNode[5] = nidx + ghostDimY*ghostDimX             + 1;
        localNode[6] = nidx + ghostDimY*ghostDimX + ghostDimX + 1;
        localNode[7] = nidx + ghostDimY*ghostDimX + ghostDimX    ;
#if INIT_DEBUG
        // Print out the element identifier and its neighboring nodes
        CkPrintf("{%d,%d,%d} E[%d|%d]: %d %d %d %d %d %d %d %d\n",
            thisIndex.x, thisIndex.y, thisIndex.z, zidx,
            thisIndex.z*chareDimY*chareDimX*blockDimX*blockDimY*blockDimZ
            + thisIndex.y*chareDimX*blockDimX*blockDimY*blockDimZ
            + thisIndex.x*blockDimX*blockDimY*blockDimZ
            + k*blockDimY*blockDimX + j*blockDimX + i,
            localNode[0],localNode[1],localNode[2],localNode[3],
            localNode[4],localNode[5],localNode[6],localNode[7]);
#endif
        zidx++;
        nidx++;
      }
      nidx++;
    }
    nidx += ghostDimX;
  }

  // Create a Material IndexSet (entire domain same material for now)
  for (Index_t i = 0; i < numElems; i++) {
    domain->m_matElemlist[i] = i;
  }

  // Initialize Material parameters  
  domain->m_dtcourant = DT_COURANT;
  domain->m_dthydro = DT_HYDRO;

  domain->e_cut = E_CUT;
  domain->p_cut = P_CUT;
  domain->q_cut = Q_CUT;
  domain->u_cut = U_CUT;
  domain->v_cut = V_CUT;

  domain->hgcoef = HG_COEF;
  domain->ss4o3 = SS_4O3;

  domain->qstop = Q_STOP;
  domain->monoq_max_slope = MONOQ_MAX_SLOPE;
  domain->monoq_limiter_mult = MONOQ_LIM_MULT;
  domain->qlc_monoq = QLC_MONOQ;
  domain->qqc_monoq = QQC_MONOQ;
  domain->qqc = QQC;

  domain->pmin = P_MIN;
  domain->emin = E_MIN;
  domain->dvovmax = DVOV_MAX;
  domain->eosvmax = EOSV_MAX;
  domain->eosvmin = EOSV_MIN;
  domain->refdens = REF_DENS;

  // Initialize Field Data
  for (Index_t i = 0; i < numElems; i++) {
    Real_t Pos_local[8][3];
    Index_t *localNode = &domain->m_nodelist[NODES_PER_ELEM*i];
    for (Index_t lnode = 0; lnode < NODES_PER_ELEM; lnode++) {
      Index_t gnode = localNode[lnode];
      Pos_local[lnode][0] = domain->m_x[gnode];
      Pos_local[lnode][1] = domain->m_y[gnode];
      Pos_local[lnode][2] = domain->m_z[gnode];
    }
    // volume calculations
    Real_t volume = CalcElemVolume(Pos_local);
    domain->m_volo[i] = volume;
    domain->m_elemMass[i] = volume;
    for (Index_t j = 0; j < NODES_PER_ELEM; j++) {
      domain->m_nodalMass[localNode[j]] += volume/NODES_PER_ELEM;
    }
#if INIT_DEBUG
    CkPrintf("{%d,%d,%d} E[%d]: %e\n",
        thisIndex.x, thisIndex.y, thisIndex.z, i, domain->m_elemMass[i]);
#endif
  }

#if INIT_DEBUG
  nidx = 0;
  for (Index_t k = 0; k < ghostDimZ; k++) {
    for (Index_t j = 0; j < ghostDimY; j++) {
      for (Index_t i = 0; i < ghostDimX; i++) {
        CkPrintf("{%d,%d,%d} M(%d|%d): %e\n",
            thisIndex.x, thisIndex.y, thisIndex.z, nidx,
            thisIndex.z*(elemDimY+1)*(elemDimX+1)*blockDimZ
            + thisIndex.y*(elemDimX+1)*blockDimY + thisIndex.x*blockDimX
            + k*(elemDimX+1)*(elemDimX+1) + j*(elemDimX+1) + i, domain->m_nodalMass[nidx]);
        nidx++;
      }
    }
  }
#endif

  // Deposit Energy
  if(thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
    domain->m_e[0] = INITIAL_ENERGY;
  }

  // Setup Symmetry Nodesets
  // yz plane
  if (thisIndex.x == 0) {
    domain->m_symmX.resize(ghostDimY*ghostDimZ);
    nidx = 0;
    for (Index_t k = 0; k < ghostDimZ; k++) {
      for (Index_t j = 0; j < ghostDimY; j++) {
        domain->m_symmX[nidx] = k*ghostDimY*ghostDimX + j*ghostDimX;
        nidx++;
      }
    }
#if INIT_DEBUG
    std::string debugString;
    char buffer[30];
    sprintf(buffer, "{%d,%d,%d} Sx[%d]: ", thisIndex.x, thisIndex.y, thisIndex.z,
        thisIndex.z*chareDimY*chareDimX + thisIndex.y*chareDimX + thisIndex.x);
    debugString.append(buffer);
    for (Index_t i = 0; i < domain->m_symmX.size(); i++) {
      sprintf(buffer, "%d ", domain->m_symmX[i]);
      debugString.append(buffer);
    }
    CkPrintf("%s\n",debugString.c_str());
#endif
  }
  // xz plane
  if (thisIndex.y == 0) {
    domain->m_symmY.resize(ghostDimX*ghostDimZ);
    nidx = 0;
    for (Index_t k = 0; k < ghostDimZ; k++) {
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->m_symmY[nidx] = k*ghostDimY*ghostDimX + i;
        nidx++;
      }
    }
#if INIT_DEBUG
    std::string debugString;
    char buffer[30];
    sprintf(buffer, "{%d,%d,%d} Sy[%d]: ", thisIndex.x, thisIndex.y, thisIndex.z,
        thisIndex.z*chareDimY*chareDimX + thisIndex.y*chareDimX + thisIndex.x);
    debugString.append(buffer);
    for (Index_t j = 0; j < domain->m_symmY.size(); j++) {
      sprintf(buffer, "%d ", domain->m_symmY[j]);
      debugString.append(buffer);
    }
    CkPrintf("%s\n",debugString.c_str());
#endif
  }
  // xy plane
  if (thisIndex.z == 0) {
    domain->m_symmZ.resize(ghostDimX*ghostDimY);
    nidx = 0;
    for (Index_t j = 0; j < ghostDimY; j++) {
      for (Index_t i = 0; i < ghostDimX; i++) {
        domain->m_symmZ[nidx] = j*ghostDimX + i;
        nidx++;
      }
    }
#if INIT_DEBUG
    std::string debugString;
    char buffer[30];
    sprintf(buffer, "{%d,%d,%d} Sz[%d]: ", thisIndex.x, thisIndex.y, thisIndex.z,
        thisIndex.z*chareDimY*chareDimX + thisIndex.y*chareDimX + thisIndex.x);
    debugString.append(buffer);
    for (Index_t k = 0; k < domain->m_symmZ.size(); k++) {
      sprintf(buffer, "%d ", domain->m_symmZ[k]);
      debugString.append(buffer);
    }
    CkPrintf("%s\n",debugString.c_str());
#endif
  }

  // Set up element connectivity information
  // Also set up boundary condition information
  //   when updating the connectivity at the edges
  // Faces on "external" boundaries will be
  //   symmetry plane or free surface BCs
  // x -> xi
  domain->m_lxim[0] = 0;
  for (Index_t i = 1; i < numElems; i++) {
    domain->m_lxim[i] = i-1;
    domain->m_lxip[i-1] = i;
  }
  domain->m_lxip[numElems-1] = numElems-1;
  for (Index_t k = 0; k < blockDimZ; k++) {
    for (Index_t j = 0; j < blockDimY; j++) {
      domain->m_lxim[k*blockDimY*blockDimX+j*blockDimX] = backOffset + k*blockDimY+j;
      domain->m_lxip[k*blockDimY*blockDimX+j*blockDimX+(blockDimX-1)] = frontOffset + k*blockDimY+j;
      // xi BCs
      if (thisIndex.x == 0) {
        domain->m_elemBC[k*blockDimY*blockDimX+j*blockDimX] |= XI_M_SYMM;
      }
      if (thisIndex.x == chareDimX-1) {
        domain->m_elemBC[k*blockDimY*blockDimX+j*blockDimX+(blockDimX-1)] |= XI_P_FREE;
      }
    }
  }
  // y -> eta
  for (Index_t i = 0; i < blockDimX; i++) {
    domain->m_letam[i] = i;
    domain->m_letap[numElems-blockDimX+i] = numElems-blockDimX+i;
  }
  for (Index_t i = blockDimX; i < numElems; i++) {
    domain->m_letam[i] = i - blockDimX;
    domain->m_letap[i-blockDimX] = i;
  }
  for (Index_t k = 0; k < blockDimZ; k++) {
    for (Index_t i = 0; i < blockDimX; i++) {
      domain->m_letam[k*blockDimY*blockDimX+i] = leftOffset + k*blockDimX+i;
      domain->m_letap[k*blockDimY*blockDimX+(blockDimY-1)*blockDimX+i] = rightOffset + k*blockDimX+i;
      // eta BCs
      if (thisIndex.y == 0) {
        domain->m_elemBC[k*blockDimY*blockDimX+i] |= ETA_M_SYMM;
      }
      if (thisIndex.y == chareDimY-1) {
        domain->m_elemBC[k*blockDimY*blockDimX+(blockDimY-1)*blockDimX+i] |= ETA_P_FREE;
      }
    }
  }
  // z -> zeta
  for (Index_t i = 0; i < blockDimY*blockDimX; i++) {
    domain->m_lzetam[i] = i;
    domain->m_lzetap[numElems-blockDimY*blockDimX+i] = numElems-blockDimY*blockDimX+i;
  }
  for (Index_t i = blockDimY*blockDimX; i < numElems; i++) {
    domain->m_lzetam[i] = i - blockDimY*blockDimX;
    domain->m_lzetap[i-blockDimY*blockDimX] = i;
  }
  for (Index_t j = 0; j < blockDimY; j++) {
    for (Index_t i = 0; i < blockDimX; i++) {
      domain->m_lzetam[j*blockDimX+i] = downOffset + j*blockDimX+i;
      domain->m_lzetap[(blockDimZ-1)*blockDimY*blockDimX+j*blockDimX+i] = upOffset + j*blockDimX+i;
      // zeta BCs
      if (thisIndex.z == 0) {
        domain->m_elemBC[j*blockDimX+i] |= ZETA_M_SYMM;
      }
      if (thisIndex.z == chareDimZ-1) {
        domain->m_elemBC[(blockDimZ-1)*blockDimY*blockDimX+j*blockDimX+i] |= ZETA_P_FREE;
      }
    }
  }
#if INIT_DEBUG
  for (Index_t i = 0; i < numElems; i++) {
    CkPrintf("{%d,%d,%d} F[%d]: %d %d %d %d %d %d 0x%03x\n",
        thisIndex.x, thisIndex.y, thisIndex.z, i,
        domain->m_lxim[i], domain->m_lxip[i], domain->m_letam[i], domain->m_letap[i], domain->m_lzetam[i], domain->m_lzetap[i], domain->m_elemBC[i]);
  }
#endif
#if CHAPEL_MESH
  CkPrintf("%d %d\n", elemDimX*elemDimY*elemDimZ, ghostDimX*ghostDimY*ghostDimZ);
  for (Index_t i = 0; i < numNodes; i++) {
    CkPrintf("%f %f %f\n", domain->m_x[i], domain->m_y[i], domain->m_z[i]);
  }
  for (Index_t i = 0; i < numElems; i++) {
    Index_t *localNode = &domain->m_nodelist[NODES_PER_ELEM*i];
    CkPrintf("%d %d %d %d %d %d %d %d\n",
        localNode[0],localNode[1],localNode[2],localNode[3],
        localNode[4],localNode[5],localNode[6],localNode[7]);
  }
  for (Index_t i = 0; i < numElems; i++) {
    CkPrintf("%d %d %d %d %d %d\n",
        domain->m_lxim[i], domain->m_lxip[i], domain->m_letam[i], domain->m_letap[i], domain->m_lzetam[i], domain->m_lzetap[i]);
  }
  CkPrintf("%d\n", domain->m_symmX.size());
  for (Index_t i = 0; i < domain->m_symmX.size(); i++) {
    CkPrintf("%d\n", domain->m_symmX[i]);
  }
  CkPrintf("%d\n", domain->m_symmY.size());
  for (Index_t i = 0; i < domain->m_symmY.size(); i++) {
    CkPrintf("%d\n", domain->m_symmY[i]);
  }
  CkPrintf("%d\n", domain->m_symmZ.size());
  for (Index_t i = 0; i < domain->m_symmZ.size(); i++) {
    CkPrintf("%d\n", domain->m_symmZ[i]);
  }
  CkPrintf("%d\n", ghostDimX*ghostDimY + ghostDimX*(ghostDimZ-1) + (ghostDimY-1)*(ghostDimZ-1));
  for (Index_t i = 0; i < numNodes; i++) {
    if (domain->m_x[i] == MAX_NODE_POSITION ||
        domain->m_y[i] == MAX_NODE_POSITION ||
        domain->m_z[i] == MAX_NODE_POSITION) {
      CkPrintf("%d\n",i);
    }
  }
#endif
}

// Print Energy at end for corectness testing
void
Domain::printEnergy() {
  CkPrintf("Time: %.12e  Energy: %.12e\n", domain->totalTime, domain->m_e[0]);
  CkExit();
}

// -----------------------------------------
//
// Domain Computation Flow

void
Domain::LagrangeNodal_Forces() {
  Real_t hgcoef = domain->hgcoef;
  // Reset forces to zero
  for (Index_t i = 0; i < numNodes; i++) {
    domain->m_fx[i] = Real_t(0.0);
    domain->m_fy[i] = Real_t(0.0);
    domain->m_fz[i] = Real_t(0.0);
  }
  IntegrateStressForElems(); 
#if ITER_TIMING
  stressTime += elapsed(getticks(),tempTime);
  tempTime = getticks();
#endif
  if (hgcoef > Real_t(0.)) {
    CalcFBHourglassForceForElems(hgcoef);
  }
}

void
Domain::LagrangeNodal_Positions() {
  CalcAccelerationForNodes();
  ApplyAccelerationBoundaryConditionsForNodes();
  
  for (Index_t i = 0; i < numNodes; i++) {
    // Calculate new Velocity for the Node
    // velocity.x
    Real_t xdtmp = domain->m_xd[i] + domain->m_xdd[i]*domain->deltaTime;
    if (FABS(xdtmp) < domain->u_cut) { xdtmp = Real_t(0.0); }
    domain->m_xd[i] = xdtmp;
    // velocity.y
    Real_t ydtmp = domain->m_yd[i] + domain->m_ydd[i]*domain->deltaTime;
    if (FABS(ydtmp) < domain->u_cut) { ydtmp = Real_t(0.0); }
    domain->m_yd[i] = ydtmp;
    // velocity.z
    Real_t zdtmp = domain->m_zd[i] + domain->m_zdd[i]*domain->deltaTime;
    if (FABS(zdtmp) < domain->u_cut) { zdtmp = Real_t(0.0); }
    domain->m_zd[i] = zdtmp;
    
    // Calculate new Position for the Node
    domain->m_x[i] += domain->m_xd[i]*domain->deltaTime;
    domain->m_y[i] += domain->m_yd[i]*domain->deltaTime;
    domain->m_z[i] += domain->m_zd[i]*domain->deltaTime;
  }
}

void
Domain::LagrangeElements_KinVisc() {
  const Real_t dt = domain->deltaTime;

  // Element loop to do some stuff not included in the 'elemlib' function
  for (Index_t k = 0; k < numElems; k++) {
    Real_t Pos_local[8][3];
    Real_t Vel_local[8][3];
    Index_t *localNode = &domain->m_nodelist[k*NODES_PER_ELEM];
    // Collect local positions and velocities
    for (Index_t lnode = 0; lnode < NODES_PER_ELEM; lnode++) {
      Index_t gnode = localNode[lnode];
      Pos_local[lnode][0] = domain->m_x[gnode];
      Pos_local[lnode][1] = domain->m_y[gnode];
      Pos_local[lnode][2] = domain->m_z[gnode];
      Vel_local[lnode][0] = domain->m_xd[gnode];
      Vel_local[lnode][1] = domain->m_yd[gnode];
      Vel_local[lnode][2] = domain->m_zd[gnode];
    }

    CalcKinematicsForElems(k, Pos_local, Vel_local, dt);

    // Calc strain rate and apply as contraint (only done in FB element)
    Real_t vdov = domain->m_dxx[k] + domain->m_dyy[k] + domain->m_dzz[k];
    Real_t vdovthird = vdov/Real_t(3.0);
    // Make the rate of deformation tensor deviatoric
    domain->m_vdov[k] = vdov;
    domain->m_dxx[k] -= vdovthird;
    domain->m_dyy[k] -= vdovthird;
    domain->m_dzz[k] -= vdovthird;

    // Compute the Viscocity Gradients
    CalcMonotonicQGradientsForElems(k, Pos_local, Vel_local);
  }
}

void
Domain::LagrangeElements_QEOS() {
  const Real_t ptiny = P_TINY;
  const Real_t qstop = domain->qstop;
  const Real_t monoq_max_slope = domain->monoq_max_slope;
  const Real_t monoq_limiter_mult = domain->monoq_limiter_mult;
  const Real_t qlc_monoq = domain->qlc_monoq;
  const Real_t qqc_monoq = domain->qqc_monoq;
  const Real_t eosvmin = domain->eosvmin;
  const Real_t eosvmax = domain->eosvmax;

  const Real_t e_cut = domain->e_cut;
  const Real_t p_cut = domain->p_cut;
  const Real_t q_cut = domain->q_cut;
  const Real_t pmin  = domain->pmin;
  const Real_t emin  = domain->emin;
  const Real_t rho0  = domain->refdens;

  for (Index_t k = 0; k < numElems; k++) {
    CalcMonotonicQRegionForElems(k, qlc_monoq, qqc_monoq,
                        monoq_limiter_mult, monoq_max_slope,
                        ptiny, qstop);
    EvalEOSForElems(k, e_cut, p_cut, q_cut, eosvmin, eosvmax,
                    pmin, emin, rho0);
  }
#if ITER_TIMING
  viscTime += elapsed(getticks(),tempTime);
  tempTime = getticks();
#endif

#if ITER_TIMING
  matTime += elapsed(getticks(),tempTime);
  tempTime = getticks();
#endif
  for (Index_t k = 0; k < numElems; k++) {
    UpdateVolumesForElems(k);
  }
}

void
Domain::CalcTimeConstraintsForElems() {
  CalcCourantConstraintForElems();
  CalcHydroConstraintForElems();
  
  // Time Increment
  Real_t targetdt = stopTime - domain->totalTime;
  Real_t olddt = domain->deltaTime;

  Real_t newdt = DT_MAX;
  if (domain->m_dtcourant < newdt) {
    newdt = domain->m_dtcourant / Real_t(2.0);
  }
  if (domain->m_dthydro < newdt) {
    newdt = domain->m_dthydro * Real_t(2.0) / Real_t(3.0);
  }

  Real_t ratio = newdt / olddt;
  if (ratio >= Real_t(1.0)) {
    if (ratio < DT_MULT_LB) {
      newdt = olddt;
    }
    else if (ratio > DT_MULT_UB) {
      newdt = olddt * DT_MULT_UB;
    }
  }

  if (newdt > DT_MAX) {
    newdt = DT_MAX;
  }
  domain->deltaTime = newdt;

  // Try to prevent very small scaling on the next cycle
  if ((targetdt > domain->deltaTime) &&
      (targetdt < (Real_t(4.0) * domain->deltaTime / Real_t(3.0)))) {
    targetdt = Real_t(2.0) * domain->deltaTime / Real_t(3.0);
  }

  if (targetdt < domain->deltaTime) {
    domain->deltaTime = targetdt;
  }
}

// -----------------------------------------
//
// Nodal Functions
inline void
Domain::InitStressTermsForElems(Index_t k, Real_t *sigxx, Real_t *sigyy, Real_t *sigzz) {
  // Pull in the stresses appropriate to the hydro integration
  *sigxx = *sigyy = *sigzz =  - domain->m_p[k] - domain->m_q[k];
}

inline void
Domain::IntegrateStressForElems() {
  // loop over all elements
  for (Index_t k = 0; k < numElems; k++) {
    Real_t B[3][8]; // shape function derivatives
    Real_t x_local[8];
    Real_t y_local[8];
    Real_t z_local[8];

    const Index_t* const localNode = &domain->m_nodelist[NODES_PER_ELEM*k];
    
    // get nodal coordinates from global arrays and copy into local arrays
    for (Index_t lnode = 0; lnode < NODES_PER_ELEM; lnode++) {
      Index_t gnode = localNode[lnode];
      x_local[lnode] = domain->m_x[gnode];
      y_local[lnode] = domain->m_y[gnode];
      z_local[lnode] = domain->m_z[gnode];
    }

    // Volume calculation involves extra work for numerical consistency
    CalcElemShapeFunctionDerivatives(x_local, y_local, z_local,
        B, &domain->m_determ[k]);

    CalcElemNodeNormals(B[0], B[1], B[2], x_local, y_local, z_local);

    // Get Initial stress terms
    Real_t sigxx;
    Real_t sigyy;
    Real_t sigzz;
    InitStressTermsForElems(k, &sigxx, &sigyy, &sigzz);

    Real_t fx_local[8];
    Real_t fy_local[8];
    Real_t fz_local[8];
    SumElemStressesToNodeForces(B, sigxx, sigyy, sigzz,
                                fx_local, fy_local, fz_local);

    // copy nodal force contributions to global force arrray.
    for (Index_t lnode = 0; lnode < NODES_PER_ELEM; lnode++) {
      Index_t gnode = localNode[lnode];
      domain->m_fx[gnode] += fx_local[lnode];
      domain->m_fy[gnode] += fy_local[lnode];
      domain->m_fz[gnode] += fz_local[lnode];
    }
  }
}

inline void
Domain::CalcAccelerationForNodes() {
  for (Index_t i = 0; i < numNodes; i++) {
    domain->m_xdd[i] = domain->m_fx[i] / domain->m_nodalMass[i];
    domain->m_ydd[i] = domain->m_fy[i] / domain->m_nodalMass[i];
    domain->m_zdd[i] = domain->m_fz[i] / domain->m_nodalMass[i];
  }
}

inline void
Domain::ApplyAccelerationBoundaryConditionsForNodes() {
  // yz plane
  if (thisIndex.x == 0) {
    for (Index_t i = 0; i < ghostDimY*ghostDimZ; i++) {
      domain->m_xdd[domain->m_symmX[i]] = Real_t(0.0);
    }
  }
  // xz plane
  if (thisIndex.y == 0) {
    for (Index_t i = 0; i < ghostDimX*ghostDimZ; i++) {
      domain->m_ydd[domain->m_symmY[i]] = Real_t(0.0);
    }
  }
  // xy plane
  if (thisIndex.z == 0) {
    for (Index_t i = 0; i < ghostDimX*ghostDimY; i++) {
      domain->m_zdd[domain->m_symmZ[i]] = Real_t(0.0);
    }
  }
}


// -----------------------------------------
//
// Elemental Functions
inline void
Domain::CalcKinematicsForElems(const Index_t k, const Real_t Pos[8][3],
                               const Real_t Vel[8][3], const Real_t dt) {
  Real_t B[3][NODES_PER_ELEM] ; /** shape function derivatives */
  Real_t D[6] ;
  Real_t Pos_local[8][3];
  Real_t detJ = Real_t(0.0) ;

  Real_t volume;
  Real_t relativeVolume;

  // volume calculations
  volume = CalcElemVolume(Pos);
  relativeVolume = volume / domain->m_volo[k];
  domain->m_vnew[k] = relativeVolume;
  domain->m_delv[k] = relativeVolume - domain->m_v[k];

#if VOLUME_DEBUG
  if (relativeVolume <= Real_t(0.0)) {
    CkAbort(VOLUME_ERR);
  }
#endif

  // set characteristic length
  domain->m_arealg[k] = CalcElemCharacteristicLength(Pos, volume);

  Real_t dt2 = Real_t(0.5) * dt;
  for (Index_t j = 0; j < NODES_PER_ELEM; j++) {
    Pos_local[j][0] = Pos[j][0] - dt2 * Vel[j][0];
    Pos_local[j][1] = Pos[j][1] - dt2 * Vel[j][1];
    Pos_local[j][2] = Pos[j][2] - dt2 * Vel[j][2];
  }

  CalcElemShapeFunctionDerivatives(Pos_local, B, &detJ);
  CalcElemVelocityGradient(Vel, B, detJ, D);

  // put velocity gradient quantities into their global arrays.
  domain->m_dxx[k] = D[0];
  domain->m_dyy[k] = D[1];
  domain->m_dzz[k] = D[2];
}

inline void
Domain::UpdateVolumesForElems(const Index_t k) {
  Real_t vtmp = domain->m_vnew[k];
  if (FABS(vtmp - Real_t(1.0)) < domain->v_cut) {
    vtmp = Real_t(1.0);
  }
  domain->m_v[k] = vtmp;
}

// -----------------------------------------
//
// Timing Functions

inline void
Domain::CalcCourantConstraintForElems() {
  Real_t dtcourant = DT_COURANT;
  Index_t courant_elem = -1 ;
  Real_t qqc2 = Real_t(64.0) * domain->qqc * domain->qqc ;

  for (Index_t k = 0; k < numElems; k++) {
    Index_t indx = domain->m_matElemlist[k];

    Real_t dtf = domain->m_ss[indx] * domain->m_ss[indx] ;

    if (domain->m_vdov[indx] < Real_t(0.) ) {
      dtf = dtf
        + qqc2 * domain->m_arealg[indx] * domain->m_arealg[indx]
        * domain->m_vdov[indx] * domain->m_vdov[indx] ;
    }

    dtf = SQRT(dtf) ;

    dtf = domain->m_arealg[indx] / dtf ;

    /* determine minimum timestep with its corresponding elem */
    if (domain->m_vdov[indx] != Real_t(0.)) {
      if ( dtf < dtcourant ) {
        dtcourant = dtf ;
        courant_elem = indx ;
      }
    }
  }

  /* Don't try to register a time constraint if none of the elements
   * were active */
  if (courant_elem != -1) {
    domain->m_dtcourant = dtcourant ;
  }
}

inline void
Domain::CalcHydroConstraintForElems() {
  Real_t dthydro = DT_HYDRO;
  Index_t hydro_elem = -1 ;
  Real_t dvovmax = domain->dvovmax;

  for (Index_t k = 0; k < numElems; k++) {
    Index_t indx = domain->m_matElemlist[k];

    if (domain->m_vdov[indx] != Real_t(0.)) {
      Real_t dtdvov = dvovmax / (FABS(domain->m_vdov[indx])+DT_HYDRO_INV) ;
      if ( dthydro > dtdvov ) {
        dthydro = dtdvov ;
        hydro_elem = indx ;
      }
    }
  }

  if (hydro_elem != -1) {
    domain->m_dthydro = dthydro ;
  }
}

// -----------------------------------------
//
// Calculation Functions (nodal)
inline void
Domain::CalcElemNodeNormals(Real_t pfx[8], Real_t pfy[8], Real_t pfz[8],
    const Real_t x[8], const Real_t y[8], const Real_t z[8]) {
  for (Index_t i = 0 ; i < 8 ; ++i) {
    pfx[i] = Real_t(0.0);
    pfy[i] = Real_t(0.0);
    pfz[i] = Real_t(0.0);
  }
  /* evaluate face one: nodes 0, 1, 2, 3 */
  SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
      &pfx[1], &pfy[1], &pfz[1],
      &pfx[2], &pfy[2], &pfz[2],
      &pfx[3], &pfy[3], &pfz[3],
      x[0], y[0], z[0], x[1], y[1], z[1],
      x[2], y[2], z[2], x[3], y[3], z[3]);
  /* evaluate face two: nodes 0, 4, 5, 1 */
  SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
      &pfx[4], &pfy[4], &pfz[4],
      &pfx[5], &pfy[5], &pfz[5],
      &pfx[1], &pfy[1], &pfz[1],
      x[0], y[0], z[0], x[4], y[4], z[4],
      x[5], y[5], z[5], x[1], y[1], z[1]);
  /* evaluate face three: nodes 1, 5, 6, 2 */
  SumElemFaceNormal(&pfx[1], &pfy[1], &pfz[1],
      &pfx[5], &pfy[5], &pfz[5],
      &pfx[6], &pfy[6], &pfz[6],
      &pfx[2], &pfy[2], &pfz[2],
      x[1], y[1], z[1], x[5], y[5], z[5],
      x[6], y[6], z[6], x[2], y[2], z[2]);
  /* evaluate face four: nodes 2, 6, 7, 3 */
  SumElemFaceNormal(&pfx[2], &pfy[2], &pfz[2],
      &pfx[6], &pfy[6], &pfz[6],
      &pfx[7], &pfy[7], &pfz[7],
      &pfx[3], &pfy[3], &pfz[3],
      x[2], y[2], z[2], x[6], y[6], z[6],
      x[7], y[7], z[7], x[3], y[3], z[3]);
  /* evaluate face five: nodes 3, 7, 4, 0 */
  SumElemFaceNormal(&pfx[3], &pfy[3], &pfz[3],
      &pfx[7], &pfy[7], &pfz[7],
      &pfx[4], &pfy[4], &pfz[4],
      &pfx[0], &pfy[0], &pfz[0],
      x[3], y[3], z[3], x[7], y[7], z[7],
      x[4], y[4], z[4], x[0], y[0], z[0]);
  /* evaluate face six: nodes 4, 7, 6, 5 */
  SumElemFaceNormal(&pfx[4], &pfy[4], &pfz[4],
      &pfx[7], &pfy[7], &pfz[7],
      &pfx[6], &pfy[6], &pfz[6],
      &pfx[5], &pfy[5], &pfz[5],
      x[4], y[4], z[4], x[7], y[7], z[7],
      x[6], y[6], z[6], x[5], y[5], z[5]);
}
inline void
Domain::SumElemStressesToNodeForces(const Real_t B[3][8],
                                    const Real_t stress_xx,
                                    const Real_t stress_yy,
                                    const Real_t stress_zz,
                                    Real_t* const fx,
                                    Real_t* const fy,
                                    Real_t* const fz) {
  for (Index_t i = 0; i < 8; i++) {
    fx[i] = -( stress_xx * B[0][i] );
    fy[i] = -( stress_yy * B[1][i] );
    fz[i] = -( stress_zz * B[2][i] );
  }
}

inline void
Domain::CalcElemVolumeDerivative(Real_t dvdx[8], Real_t dvdy[8], Real_t dvdz[8],
    const Real_t x[8], const Real_t y[8],
    const Real_t z[8]) {
  VoluDer(x[1], x[2], x[3], x[4], x[5], x[7],
      y[1], y[2], y[3], y[4], y[5], y[7],
      z[1], z[2], z[3], z[4], z[5], z[7],
      &dvdx[0], &dvdy[0], &dvdz[0]);
  VoluDer(x[0], x[1], x[2], x[7], x[4], x[6],
      y[0], y[1], y[2], y[7], y[4], y[6],
      z[0], z[1], z[2], z[7], z[4], z[6],
      &dvdx[3], &dvdy[3], &dvdz[3]);
  VoluDer(x[3], x[0], x[1], x[6], x[7], x[5],
      y[3], y[0], y[1], y[6], y[7], y[5],
      z[3], z[0], z[1], z[6], z[7], z[5],
      &dvdx[2], &dvdy[2], &dvdz[2]);
  VoluDer(x[2], x[3], x[0], x[5], x[6], x[4],
      y[2], y[3], y[0], y[5], y[6], y[4],
      z[2], z[3], z[0], z[5], z[6], z[4],
      &dvdx[1], &dvdy[1], &dvdz[1]);
  VoluDer(x[7], x[6], x[5], x[0], x[3], x[1],
      y[7], y[6], y[5], y[0], y[3], y[1],
      z[7], z[6], z[5], z[0], z[3], z[1],
      &dvdx[4], &dvdy[4], &dvdz[4]);
  VoluDer(x[4], x[7], x[6], x[1], x[0], x[2],
      y[4], y[7], y[6], y[1], y[0], y[2],
      z[4], z[7], z[6], z[1], z[0], z[2],
      &dvdx[5], &dvdy[5], &dvdz[5]);
  VoluDer(x[5], x[4], x[7], x[2], x[1], x[3],
      y[5], y[4], y[7], y[2], y[1], y[3],
      z[5], z[4], z[7], z[2], z[1], z[3],
      &dvdx[6], &dvdy[6], &dvdz[6]);
  VoluDer(x[6], x[5], x[4], x[3], x[2], x[0],
      y[6], y[5], y[4], y[3], y[2], y[0],
      z[6], z[5], z[4], z[3], z[2], z[0],
      &dvdx[7], &dvdy[7], &dvdz[7]);
}

inline void
Domain::CalcFBHourglassForceForElems(Real_t hourg) {
  /*************************************************
   *
   *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
   *               force.
   *
   *************************************************/
  
  /*************************************************/
  /*    compute the hourglass modes */

  for(Index_t i2 = 0; i2 < numElems ; i2++) {

    Real_t   x1[8],   y1[8],   z1[8];
    Real_t dvdx[8], dvdy[8], dvdz[8];

    const Index_t *localNode = &domain->m_nodelist[NODES_PER_ELEM*i2];

    CollectDomainNodesToElemNodes(localNode, x1, y1, z1);
    CalcElemVolumeDerivative(dvdx, dvdy, dvdz, x1, y1, z1);

    Real_t determ = domain->m_volo[i2] * domain->m_v[i2];
  
    Real_t *fx_local, *fy_local, *fz_local;
    Real_t hgfx[8], hgfy[8], hgfz[8];

    Real_t coefficient;

    Real_t hourgam[8][4];
    Real_t xd1[8], yd1[8], zd1[8];

    Real_t volinv = Real_t(1.0)/determ;
    Real_t ss1, mass1, volume13;
    for(Index_t i1 = 0; i1 < 4; i1++) {

      Real_t hourmodx =
        x1[0] * Gamma[i1][0] + x1[1] * Gamma[i1][1] +
        x1[2] * Gamma[i1][2] + x1[3] * Gamma[i1][3] +
        x1[4] * Gamma[i1][4] + x1[5] * Gamma[i1][5] +
        x1[6] * Gamma[i1][6] + x1[7] * Gamma[i1][7];

      Real_t hourmody =
        y1[0] * Gamma[i1][0] + y1[1] * Gamma[i1][1] +
        y1[2] * Gamma[i1][2] + y1[3] * Gamma[i1][3] +
        y1[4] * Gamma[i1][4] + y1[5] * Gamma[i1][5] +
        y1[6] * Gamma[i1][6] + y1[7] * Gamma[i1][7];

      Real_t hourmodz =
        z1[0] * Gamma[i1][0] + z1[1] * Gamma[i1][1] +
        z1[2] * Gamma[i1][2] + z1[3] * Gamma[i1][3] +
        z1[4] * Gamma[i1][4] + z1[5] * Gamma[i1][5] +
        z1[6] * Gamma[i1][6] + z1[7] * Gamma[i1][7];

      hourgam[0][i1] = Gamma[i1][0] -  volinv*(dvdx[0] * hourmodx +
          dvdy[0] * hourmody +
          dvdz[0] * hourmodz );

      hourgam[1][i1] = Gamma[i1][1] -  volinv*(dvdx[1] * hourmodx +
          dvdy[1] * hourmody +
          dvdz[1] * hourmodz );

      hourgam[2][i1] = Gamma[i1][2] -  volinv*(dvdx[2] * hourmodx +
          dvdy[2] * hourmody +
          dvdz[2] * hourmodz );

      hourgam[3][i1] = Gamma[i1][3] -  volinv*(dvdx[3] * hourmodx +
          dvdy[3] * hourmody +
          dvdz[3] * hourmodz );

      hourgam[4][i1] = Gamma[i1][4] -  volinv*(dvdx[4] * hourmodx +
          dvdy[4] * hourmody +
          dvdz[4] * hourmodz );

      hourgam[5][i1] = Gamma[i1][5] -  volinv*(dvdx[5] * hourmodx +
          dvdy[5] * hourmody +
          dvdz[5] * hourmodz );

      hourgam[6][i1] = Gamma[i1][6] -  volinv*(dvdx[6] * hourmodx +
          dvdy[6] * hourmody +
          dvdz[6] * hourmodz );

      hourgam[7][i1] = Gamma[i1][7] -  volinv*(dvdx[7] * hourmodx +
          dvdy[7] * hourmody +
          dvdz[7] * hourmodz );
    }

    /* compute forces */
    /* store forces into h arrays (force arrays) */

    ss1 = domain->m_ss[i2];
    mass1 = domain->m_elemMass[i2];
    volume13 = CBRT(determ);

    Index_t n0si2 = localNode[0];
    Index_t n1si2 = localNode[1];
    Index_t n2si2 = localNode[2];
    Index_t n3si2 = localNode[3];
    Index_t n4si2 = localNode[4];
    Index_t n5si2 = localNode[5];
    Index_t n6si2 = localNode[6];
    Index_t n7si2 = localNode[7];

    xd1[0] = domain->m_xd[n0si2];
    xd1[1] = domain->m_xd[n1si2];
    xd1[2] = domain->m_xd[n2si2];
    xd1[3] = domain->m_xd[n3si2];
    xd1[4] = domain->m_xd[n4si2];
    xd1[5] = domain->m_xd[n5si2];
    xd1[6] = domain->m_xd[n6si2];
    xd1[7] = domain->m_xd[n7si2];

    yd1[0] = domain->m_yd[n0si2];
    yd1[1] = domain->m_yd[n1si2];
    yd1[2] = domain->m_yd[n2si2];
    yd1[3] = domain->m_yd[n3si2];
    yd1[4] = domain->m_yd[n4si2];
    yd1[5] = domain->m_yd[n5si2];
    yd1[6] = domain->m_yd[n6si2];
    yd1[7] = domain->m_yd[n7si2];

    zd1[0] = domain->m_zd[n0si2];
    zd1[1] = domain->m_zd[n1si2];
    zd1[2] = domain->m_zd[n2si2];
    zd1[3] = domain->m_zd[n3si2];
    zd1[4] = domain->m_zd[n4si2];
    zd1[5] = domain->m_zd[n5si2];
    zd1[6] = domain->m_zd[n6si2];
    zd1[7] = domain->m_zd[n7si2];

    coefficient = - hourg * Real_t(0.01) * ss1 * mass1 / volume13;

    CalcElemFBHourglassForce(xd1, yd1, zd1, hourgam,
        coefficient, hgfx, hgfy, hgfz);

    domain->m_fx[n0si2] += hgfx[0];
    domain->m_fy[n0si2] += hgfy[0];
    domain->m_fz[n0si2] += hgfz[0];

    domain->m_fx[n1si2] += hgfx[1];
    domain->m_fy[n1si2] += hgfy[1];
    domain->m_fz[n1si2] += hgfz[1];

    domain->m_fx[n2si2] += hgfx[2];
    domain->m_fy[n2si2] += hgfy[2];
    domain->m_fz[n2si2] += hgfz[2];

    domain->m_fx[n3si2] += hgfx[3];
    domain->m_fy[n3si2] += hgfy[3];
    domain->m_fz[n3si2] += hgfz[3];

    domain->m_fx[n4si2] += hgfx[4];
    domain->m_fy[n4si2] += hgfy[4];
    domain->m_fz[n4si2] += hgfz[4];

    domain->m_fx[n5si2] += hgfx[5];
    domain->m_fy[n5si2] += hgfy[5];
    domain->m_fz[n5si2] += hgfz[5];

    domain->m_fx[n6si2] += hgfx[6];
    domain->m_fy[n6si2] += hgfy[6];
    domain->m_fz[n6si2] += hgfz[6];

    domain->m_fx[n7si2] += hgfx[7];
    domain->m_fy[n7si2] += hgfy[7];
    domain->m_fz[n7si2] += hgfz[7];
  }
}

inline void
Domain::CalcElemFBHourglassForce(Real_t *xd, Real_t *yd, Real_t *zd,
    Real_t hourgam[8][4], Real_t coefficient,
    Real_t *hgfx, Real_t *hgfy, Real_t *hgfz) {
  Real_t hxx[4];
  for(Index_t i = 0; i < 4; i++) {

  hxx[i] =
    hourgam[0][i] * xd[0] + hourgam[1][i] * xd[1] +
    hourgam[2][i] * xd[2] + hourgam[3][i] * xd[3] +
    hourgam[4][i] * xd[4] + hourgam[5][i] * xd[5] +
    hourgam[6][i] * xd[6] + hourgam[7][i] * xd[7];
  }
  for(Index_t i = 0; i < 8; i++) {

  hgfx[i] = coefficient *
    (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
     hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
  }
  for(Index_t i = 0; i < 4; i++) {

  hxx[i] =
    hourgam[0][i] * yd[0] + hourgam[1][i] * yd[1] +
    hourgam[2][i] * yd[2] + hourgam[3][i] * yd[3] +
    hourgam[4][i] * yd[4] + hourgam[5][i] * yd[5] +
    hourgam[6][i] * yd[6] + hourgam[7][i] * yd[7];
  }

  for(Index_t i = 0; i < 8; i++) {

  hgfy[i] = coefficient *
    (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
     hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
  }

  for(Index_t i = 0; i < 4; i++) {

  hxx[i] =
    hourgam[0][i] * zd[0] + hourgam[1][i] * zd[1] +
    hourgam[2][i] * zd[2] + hourgam[3][i] * zd[3] +
    hourgam[4][i] * zd[4] + hourgam[5][i] * zd[5] +
    hourgam[6][i] * zd[6] + hourgam[7][i] * zd[7];
  }

  for(Index_t i = 0; i < 8; i++) {

  hgfz[i] = coefficient *
    (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
     hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
  }
}

// -----------------------------------------
//
// Calculation Functions (elemental)
inline Real_t
Domain::CalcElemCharacteristicLength(const Real_t Pos[8][3],
                                     const Real_t volume) {
  Real_t a, charLength = Real_t(0.0);

  a = AreaFace(Pos, 0, 1, 2, 3);
  charLength = std::max(a,charLength) ;

  a = AreaFace(Pos, 4, 5, 6, 7);
  charLength = std::max(a,charLength) ;

  a = AreaFace(Pos, 0, 1, 5, 4);
  charLength = std::max(a,charLength) ;

  a = AreaFace(Pos, 1, 2, 6, 5);
  charLength = std::max(a,charLength) ;

  a = AreaFace(Pos, 2, 3, 7, 6);
  charLength = std::max(a,charLength) ;

  a = AreaFace(Pos, 3, 0, 4, 7);
  charLength = std::max(a,charLength) ;

  charLength = Real_t(4.0) * volume / SQRT(charLength);

  return charLength;
}

inline void
Domain::CalcElemVelocityGradient(const Real_t vel[8][3], const Real_t b[][8],
                                 const Real_t detJ, Real_t* const d ) {
  const int x = 0;
  const int y = 1;
  const int z = 2;

  const Real_t inv_detJ = Real_t(1.0) / detJ ;
  Real_t dyddx, dxddy, dzddx, dxddz, dzddy, dyddz;
  const Real_t* const pfx = b[0];
  const Real_t* const pfy = b[1];
  const Real_t* const pfz = b[2];

  d[0] =   inv_detJ * (  pfx[0] * (vel[0][x]-vel[6][x])
               + pfx[1] * (vel[1][x]-vel[7][x])
               + pfx[2] * (vel[2][x]-vel[4][x])
               + pfx[3] * (vel[3][x]-vel[5][x]) );

  d[1] =   inv_detJ * (  pfy[0] * (vel[0][y]-vel[6][y])
               + pfy[1] * (vel[1][y]-vel[7][y])
               + pfy[2] * (vel[2][y]-vel[4][y])
               + pfy[3] * (vel[3][y]-vel[5][y]) );

  d[2]   = inv_detJ * (  pfz[0] * (vel[0][z]-vel[6][z])
               + pfz[1] * (vel[1][z]-vel[7][z])
               + pfz[2] * (vel[2][z]-vel[4][z])
               + pfz[3] * (vel[3][z]-vel[5][z]) );

  dyddx  = inv_detJ * (  pfx[0] * (vel[0][y]-vel[6][y])
             + pfx[1] * (vel[1][y]-vel[7][y])
             + pfx[2] * (vel[2][y]-vel[4][y])
             + pfx[3] * (vel[3][y]-vel[5][y]) );

  dxddy  = inv_detJ * (  pfy[0] * (vel[0][x]-vel[6][x])
             + pfy[1] * (vel[1][x]-vel[7][x])
             + pfy[2] * (vel[2][x]-vel[4][x])
             + pfy[3] * (vel[3][x]-vel[5][x]) );

  dzddx  = inv_detJ * (  pfx[0] * (vel[0][z]-vel[6][z])
             + pfx[1] * (vel[1][z]-vel[7][z])
             + pfx[2] * (vel[2][z]-vel[4][z])
             + pfx[3] * (vel[3][z]-vel[5][z]) );

  dxddz  = inv_detJ * (  pfz[0] * (vel[0][x]-vel[6][x])
             + pfz[1] * (vel[1][x]-vel[7][x])
             + pfz[2] * (vel[2][x]-vel[4][x])
             + pfz[3] * (vel[3][x]-vel[5][x]) );

  dzddy  = inv_detJ * (  pfy[0] * (vel[0][z]-vel[6][z])
             + pfy[1] * (vel[1][z]-vel[7][z])
             + pfy[2] * (vel[2][z]-vel[4][z])
             + pfy[3] * (vel[3][z]-vel[5][z]) );

  dyddz  = inv_detJ * (  pfz[0] * (vel[0][y]-vel[6][y])
             + pfz[1] * (vel[1][y]-vel[7][y])
             + pfz[2] * (vel[2][y]-vel[4][y])
             + pfz[3] * (vel[3][y]-vel[5][y]) );

  d[5]  = Real_t( .5) * ( dxddy + dyddx );
  d[4]  = Real_t( .5) * ( dxddz + dzddx );
  d[3]  = Real_t( .5) * ( dzddy + dyddz );
}

inline void
Domain::CalcMonotonicQGradientsForElems(Index_t k, const Real_t Pos[8][3],
                                        const Real_t Vel[8][3]) {
#define SUM4(a,b,c,d) (a + b + c + d)
  const Real_t ptiny = P_TINY;
  Real_t ax,ay,az ;
  Real_t dxv,dyv,dzv ;

  const int x = 0;
  const int y = 1;
  const int z = 2;

  Real_t vol = domain->m_volo[k] * domain->m_vnew[k];
  Real_t norm = Real_t(1.0) / ( vol + ptiny ) ;

  Real_t dxj = Real_t(-0.25)*
          ( SUM4( Pos[0][x], Pos[1][x], Pos[5][x], Pos[4][x] )
          - SUM4( Pos[3][x], Pos[2][x], Pos[6][x], Pos[7][x] )) ;
  Real_t dyj = Real_t(-0.25)*
          ( SUM4( Pos[0][y], Pos[1][y], Pos[5][y], Pos[4][y] )
          - SUM4( Pos[3][y], Pos[2][y], Pos[6][y], Pos[7][y] )) ;
  Real_t dzj = Real_t(-0.25)*
          ( SUM4( Pos[0][z], Pos[1][z], Pos[5][z], Pos[4][z] )
          - SUM4( Pos[3][z], Pos[2][z], Pos[6][z], Pos[7][z] )) ;

  Real_t dxi = Real_t(-0.25)*
          ( SUM4( Pos[1][x], Pos[2][x], Pos[6][x], Pos[5][x] )
          - SUM4( Pos[0][x], Pos[3][x], Pos[7][x], Pos[4][x] )) ;
  Real_t dyi = Real_t(-0.25)*
          ( SUM4( Pos[1][y], Pos[2][y], Pos[6][y], Pos[5][y] )
          - SUM4( Pos[0][y], Pos[3][y], Pos[7][y], Pos[4][y] )) ;
  Real_t dzi = Real_t(-0.25)*
          ( SUM4( Pos[1][z], Pos[2][z], Pos[6][z], Pos[5][z] )
          - SUM4( Pos[0][z], Pos[3][z], Pos[7][z], Pos[4][z] )) ;

  Real_t dxk = Real_t(-0.25)*
          ( SUM4( Pos[4][x], Pos[5][x], Pos[6][x], Pos[7][x] )
          - SUM4( Pos[0][x], Pos[1][x], Pos[2][x], Pos[3][x] )) ;
  Real_t dyk = Real_t(-0.25)*
          ( SUM4( Pos[4][y], Pos[5][y], Pos[6][y], Pos[7][y] )
          - SUM4( Pos[0][y], Pos[1][y], Pos[2][y], Pos[3][y] )) ;
  Real_t dzk = Real_t(-0.25)*
          ( SUM4( Pos[4][z], Pos[5][z], Pos[6][z], Pos[7][z] )
          - SUM4( Pos[0][z], Pos[1][z], Pos[2][z], Pos[3][z] )) ;

    /* find delvk and delxk ( i cross j ) */

    ax = dyi*dzj - dzi*dyj ;
    ay = dzi*dxj - dxi*dzj ;
    az = dxi*dyj - dyi*dxj ;

    domain->m_delx_zeta[k] = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

    ax *= norm ;
    ay *= norm ;
    az *= norm ;

  dxv = Real_t(-0.25)*
          ( SUM4( Vel[4][x], Vel[5][x], Vel[6][x], Vel[7][x] )
          - SUM4( Vel[0][x], Vel[1][x], Vel[2][x], Vel[3][x] )) ;
  dyv = Real_t(-0.25)*
          ( SUM4( Vel[4][y], Vel[5][y], Vel[6][y], Vel[7][y] )
          - SUM4( Vel[0][y], Vel[1][y], Vel[2][y], Vel[3][y] )) ;
  dzv = Real_t(-0.25)*
          ( SUM4( Vel[4][z], Vel[5][z], Vel[6][z], Vel[7][z] )
          - SUM4( Vel[0][z], Vel[1][z], Vel[2][z], Vel[3][z] )) ;

    domain->m_delv_zeta[k] = ax*dxv + ay*dyv + az*dzv ;

    /* find delxi and delvi ( j cross k ) */

    ax = dyj*dzk - dzj*dyk ;
    ay = dzj*dxk - dxj*dzk ;
    az = dxj*dyk - dyj*dxk ;

    domain->m_delx_xi[k] = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

    ax *= norm ;
    ay *= norm ;
    az *= norm ;

  dxv = Real_t(-0.25)*
        ( SUM4( Vel[1][x], Vel[2][x], Vel[6][x], Vel[5][x] )
        - SUM4( Vel[0][x], Vel[3][x], Vel[7][x], Vel[4][x] )) ;
  dyv = Real_t(-0.25)*
        ( SUM4( Vel[1][y], Vel[2][y], Vel[6][y], Vel[5][y] )
        - SUM4( Vel[0][y], Vel[3][y], Vel[7][y], Vel[4][y] )) ;
  dzv = Real_t(-0.25)*
        ( SUM4( Vel[1][z], Vel[2][z], Vel[6][z], Vel[5][z] )
        - SUM4( Vel[0][z], Vel[3][z], Vel[7][z], Vel[4][z] )) ;

    domain->m_delv_xi[k] = ax*dxv + ay*dyv + az*dzv ;

    /* find delxj and delvj ( k cross i ) */

    ax = dyk*dzi - dzk*dyi ;
    ay = dzk*dxi - dxk*dzi ;
    az = dxk*dyi - dyk*dxi ;

    domain->m_delx_eta[k] = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

    ax *= norm ;
    ay *= norm ;
    az *= norm ;

  dxv = Real_t(-0.25)*
        ( SUM4( Vel[0][x], Vel[1][x], Vel[5][x], Vel[4][x] )
        - SUM4( Vel[3][x], Vel[2][x], Vel[6][x], Vel[7][x] )) ;
  dyv = Real_t(-0.25)*
        ( SUM4( Vel[0][y], Vel[1][y], Vel[5][y], Vel[4][y] )
        - SUM4( Vel[3][y], Vel[2][y], Vel[6][y], Vel[7][y] )) ;
  dzv = Real_t(-0.25)*
        ( SUM4( Vel[0][z], Vel[1][z], Vel[5][z], Vel[4][z] )
        - SUM4( Vel[3][z], Vel[2][z], Vel[6][z], Vel[7][z] )) ;

    domain->m_delv_eta[k] = ax*dxv + ay*dyv + az*dzv ;
#undef SUM4
}

inline void
Domain::CalcMonotonicQRegionForElems(const Index_t k, const Real_t qlc_monoq,
                                     const Real_t qqc_monoq,
                                     const Real_t monoq_limiter_mult,
                                     const Real_t monoq_max_slope,
                                     const Real_t ptiny, const Real_t qstop) {
  // calculate the monotonic q for pure regions
  Real_t qlin, qquad ;
  Real_t phixi, phieta, phizeta ;
  Int_t bcMask = domain->m_elemBC[k];
  Real_t delvm, delvp ;

  /*  phixi     */
  Real_t norm = Real_t(1.) / ( domain->m_delv_xi[k] + ptiny ) ;

  switch (bcMask & XI_M) {
    case 0:         delvm = domain->m_delv_xi[domain->m_lxim[k]] ; break ;
    case XI_M_SYMM: delvm = domain->m_delv_xi[k] ;         break ;
    case XI_M_FREE: delvm = Real_t(0.0) ;          break ;
    default:        /* ERROR */ ;                  break ;
  }
  switch (bcMask & XI_P) {
    case 0:         delvp = domain->m_delv_xi[domain->m_lxip[k]]; break ;
    case XI_P_SYMM: delvp = domain->m_delv_xi[k] ;        break ;
    case XI_P_FREE: delvp = Real_t(0.0) ;         break ;
    default:        /* ERROR */ ;                 break ;
  }

  delvm = delvm * norm ;
  delvp = delvp * norm ;

  phixi = Real_t(.5) * ( delvm + delvp ) ;

  delvm *= monoq_limiter_mult ;
  delvp *= monoq_limiter_mult ;

  if ( delvm < phixi ) phixi = delvm ;
  if ( delvp < phixi ) phixi = delvp ;
  if ( phixi < Real_t(0.)) phixi = Real_t(0.) ;
  if ( phixi > monoq_max_slope) phixi = monoq_max_slope;

  /*  phieta     */
  norm = Real_t(1.) / ( domain->m_delv_eta[k] + ptiny ) ;

  switch (bcMask & ETA_M) {
    case 0:          delvm = domain->m_delv_eta[domain->m_letam[k]]; break ;
    case ETA_M_SYMM: delvm = domain->m_delv_eta[k] ;         break ;
    case ETA_M_FREE: delvm = Real_t(0.0) ;           break ;
    default:         /* ERROR */ ;                   break ;
  }
  switch (bcMask & ETA_P) {
    case 0:          delvp = domain->m_delv_eta[domain->m_letap[k]]; break ;
    case ETA_P_SYMM: delvp = domain->m_delv_eta[k] ;         break ;
    case ETA_P_FREE: delvp = Real_t(0.0) ;           break ;
    default:         /* ERROR */ ;                   break ;
  }

  delvm = delvm * norm ;
  delvp = delvp * norm ;

  phieta = Real_t(.5) * ( delvm + delvp ) ;

  delvm *= monoq_limiter_mult ;
  delvp *= monoq_limiter_mult ;

  if ( delvm  < phieta ) phieta = delvm ;
  if ( delvp  < phieta ) phieta = delvp ;
  if ( phieta < Real_t(0.)) phieta = Real_t(0.) ;
  if ( phieta > monoq_max_slope)  phieta = monoq_max_slope;

  /*  phizeta     */
  norm = Real_t(1.) / ( domain->m_delv_eta[k] + ptiny ) ;

  switch (bcMask & ETA_M) {
    case 0:          delvm = domain->m_delv_eta[domain->m_letam[k]]; break ;
    case ETA_M_SYMM: delvm = domain->m_delv_eta[k] ;         break ;
    case ETA_M_FREE: delvm = Real_t(0.0) ;           break ;
    default:         /* ERROR */ ;                   break ;
  }
  switch (bcMask & ETA_P) {
    case 0:          delvp = domain->m_delv_eta[domain->m_letap[k]]; break ;
    case ETA_P_SYMM: delvp = domain->m_delv_eta[k] ;         break ;
    case ETA_P_FREE: delvp = Real_t(0.0) ;           break ;
    default:         /* ERROR */ ;                   break ;
  }

  delvm = delvm * norm ;
  delvp = delvp * norm ;

  phieta = Real_t(.5) * ( delvm + delvp ) ;

  delvm *= monoq_limiter_mult ;
  delvp *= monoq_limiter_mult ;

  if ( delvm  < phieta ) phieta = delvm ;
  if ( delvp  < phieta ) phieta = delvp ;
  if ( phieta < Real_t(0.)) phieta = Real_t(0.) ;
  if ( phieta > monoq_max_slope)  phieta = monoq_max_slope;

  /*  phizeta     */
  norm = Real_t(1.) / ( domain->m_delv_zeta[k] + ptiny ) ;

  switch (bcMask & ZETA_M) {
    case 0:           delvm = domain->m_delv_zeta[domain->m_lzetam[k]]; break ;
    case ZETA_M_SYMM: delvm = domain->m_delv_zeta[k] ;          break ;
    case ZETA_M_FREE: delvm = Real_t(0.0) ;             break ;
    default:          /* ERROR */ ;                     break ;
  }
  switch (bcMask & ZETA_P) {
    case 0:           delvp = domain->m_delv_zeta[domain->m_lzetap[k]]; break ;
    case ZETA_P_SYMM: delvp = domain->m_delv_zeta[k] ;          break ;
    case ZETA_P_FREE: delvp = Real_t(0.0) ;             break ;
    default:          /* ERROR */ ;                     break ;
  }

  delvm = delvm * norm ;
  delvp = delvp * norm ;

  phizeta = Real_t(.5) * ( delvm + delvp ) ;

  delvm *= monoq_limiter_mult ;
  delvp *= monoq_limiter_mult ;

  if ( delvm   < phizeta ) phizeta = delvm ;
  if ( delvp   < phizeta ) phizeta = delvp ;
  if ( phizeta < Real_t(0.)) phizeta = Real_t(0.);
  if ( phizeta > monoq_max_slope  ) phizeta = monoq_max_slope;

  /* Remove length scale */

  if ( domain->m_vdov[k] > Real_t(0.) )  {
    qlin  = Real_t(0.) ;
    qquad = Real_t(0.) ;
  }
  else {
    Real_t delvxxi   = domain->m_delv_xi[k]   * domain->m_delx_xi[k]   ;
    Real_t delvxeta  = domain->m_delv_eta[k]  * domain->m_delx_eta[k]  ;
    Real_t delvxzeta = domain->m_delv_zeta[k] * domain->m_delx_zeta[k] ;

    if ( delvxxi   > Real_t(0.) ) delvxxi   = Real_t(0.) ;
    if ( delvxeta  > Real_t(0.) ) delvxeta  = Real_t(0.) ;
    if ( delvxzeta > Real_t(0.) ) delvxzeta = Real_t(0.) ;

    Real_t rho = domain->m_elemMass[k] / (domain->m_volo[k] * domain->m_vnew[k]) ;

    qlin = -qlc_monoq * rho *
      (  delvxxi   * (Real_t(1.) - phixi) +
         delvxeta  * (Real_t(1.) - phieta) +
         delvxzeta * (Real_t(1.) - phizeta)  ) ;

    qquad = domain->qqc_monoq * rho *
      (  delvxxi*delvxxi     * (Real_t(1.) - phixi*phixi) +
         delvxeta*delvxeta   * (Real_t(1.) - phieta*phieta) +
         delvxzeta*delvxzeta * (Real_t(1.) - phizeta*phizeta)  ) ;
  }

  domain->m_qq[k] = qquad ;
  domain->m_ql[k] = qlin  ;

  /* Don't allow excessive artificial viscocity */
  if (qlin > qstop) {
    CkAbort(QSTOP_ERR);
  }
}

inline void
Domain::EvalEOSForElems(const Index_t k, const Real_t e_cut,
                        const Real_t p_cut, const Real_t q_cut,
                        const Real_t eosvmin, const Real_t eosvmax,
                        const Real_t pmin, const Real_t emin,
                        const Real_t rho0) {
  Real_t vnewc = domain->m_vnew[k];
  Real_t e_old = domain->m_e[k];
  Real_t delvc = domain->m_delv[k];
  Real_t p_old = domain->m_p[k];
  Real_t q_old = domain->m_q[k];
  Real_t qq    = domain->m_qq[k];
  Real_t ql    = domain->m_ql[k];
  Real_t work  = Real_t(0.);

  Real_t compression = Real_t(1.) / vnewc - Real_t(1.);
  Real_t vchalf = vnewc - delvc * Real_t(.5);
  Real_t compHalfStep = Real_t(1.) / vchalf - Real_t(1.);

  if ((eosvmin != Real_t(0.)) && (vnewc < eosvmin)) {
    compHalfStep = compression;
  }
  if ((eosvmax != Real_t(0.)) && (vnewc > eosvmax)) {
    /* impossible due to calling func? */
    p_old = Real_t(0.);
    compression = Real_t(0.);
    compHalfStep = Real_t(0.);
  }

  Real_t p_new;
  Real_t e_new;
  Real_t q_new;
  Real_t bvc;
  Real_t pbvc;
  Real_t ss;

  CalcEnergyForElems(p_old, e_old,  q_old, compression, compHalfStep, vnewc,
                     work, delvc, pmin, p_cut, e_cut, q_cut, emin, qq, ql,
                     rho0, eosvmax, &p_new, &e_new, &q_new, &bvc, &pbvc);

  ss = CalcSoundSpeedForElem(vnewc, e_new, p_new, pbvc, bvc, rho0);

  domain->m_p[k] = p_new;
  domain->m_e[k] = e_new;
  domain->m_q[k] = q_new;
  domain->m_ss[k] = ss;
}

inline void
Domain::CalcEnergyForElems(const Real_t p_old, const Real_t e_old,
                           const Real_t q_old, const Real_t compression,
                           const Real_t compHalfStep, const Real_t vnewc,
                           const Real_t work, const Real_t delvc,
                           const Real_t pmin, const Real_t p_cut,
                           const Real_t e_cut, const Real_t q_cut,
                           const Real_t emin, const Real_t qq, const Real_t ql,
                           const Real_t rho0, const Real_t eosvmax,
                           Real_t *p_new, Real_t *e_new, Real_t *q_new,
                           Real_t *bvc, Real_t *pbvc) {
  Real_t sixth = Real_t(1.0) / Real_t(6.0) ;
  Real_t pHalfStep ;
  Real_t vhalf;
  Real_t ssc ;
  Real_t q_tilde ;

  *e_new = e_old - Real_t(0.5) * delvc * (p_old + q_old) + Real_t(0.5) * work;

  if (*e_new  < emin )  *e_new = emin ;

  CalcPressureForElems(*e_new, compHalfStep, vnewc, pmin, p_cut, eosvmax,
             &pHalfStep, bvc, pbvc);

  vhalf = Real_t(1.) / (Real_t(1.) + compHalfStep) ;

  ssc = CalcSoundSpeedForElem(vhalf, *e_new, pHalfStep, *pbvc, *bvc, rho0);

  *q_new = ( delvc > Real_t(0.0)) ? Real_t(0.0) : ssc * ql + qq ;

  *e_new = *e_new + Real_t(0.5) * delvc * (  Real_t(3.0) * (p_old + q_old)
                           - Real_t(4.0) * (pHalfStep  + *q_new))
            + Real_t(0.5) * work;

  if (FABS(*e_new) < e_cut)  *e_new = Real_t(0.)  ;
  if (     *e_new  < emin )  *e_new = emin ;

  CalcPressureForElems(*e_new, compression, vnewc, pmin, p_cut, eosvmax,
             p_new, bvc, pbvc );

  ssc = CalcSoundSpeedForElem(vnewc, *e_new, *p_new, *pbvc, *bvc, rho0);

  q_tilde = ( delvc > Real_t(0.0)) ? Real_t(0.0) : ssc * ql + qq ;

  *e_new = *e_new - (  Real_t(7.0) * (p_old     +  q_old)
                       - Real_t(8.0) * (pHalfStep + *q_new)
                       + (*p_new + q_tilde)) * delvc * sixth ;

  if (FABS(*e_new) < e_cut)   *e_new = Real_t(0.)  ;
  if (     *e_new  < emin )   *e_new = emin ;

  CalcPressureForElems(*e_new, compression, vnewc, pmin, p_cut, eosvmax,
             p_new, bvc, pbvc );

  if ( delvc <= Real_t(0.) ) {
    ssc = CalcSoundSpeedForElem(vnewc, *e_new, *p_new, *pbvc, *bvc, rho0);
    *q_new = (ssc * ql + qq) ;
    if (FABS(*q_new) < q_cut) *q_new = Real_t(0.) ;
  }
}

inline void
Domain::CalcPressureForElems(const Real_t e_old, const Real_t compression,
                             const Real_t vnewc, const Real_t pmin,
                             const Real_t p_cut, const Real_t eosvmax,
                             Real_t *p_new, Real_t *bvc, Real_t *pbvc) {
    Real_t c1s = Real_t(2.0)/Real_t(3.0) ;
        *bvc = c1s * (compression + Real_t(1.));
       *pbvc = c1s;
      *p_new = *bvc * e_old ;

    if (FABS(*p_new) <  p_cut)  *p_new = Real_t(0.0) ;
    if ( vnewc    >= eosvmax )  *p_new = Real_t(0.0) ;
    if (*p_new       <  pmin )  *p_new = pmin ;
}

inline Real_t
Domain::CalcSoundSpeedForElem(const Real_t new_vol, const Real_t enew,
                              const Real_t pnew, const Real_t pbvc,
                              const Real_t bvc, const Real_t rho0) {
  Real_t ss;
  ss = (   pbvc * enew + new_vol * new_vol * bvc * pnew) / rho0 ;
  ss = ( ss <= Real_t(SS_MIN) ) ? Real_t(SSC_MIN) : SQRT(ss)  ;
  return ss;
}

// -----------------------------------------
//
// Calculation Functions (helper/shared)
inline void
Domain::CalcElemShapeFunctionDerivatives(const Real_t* const x,
    const Real_t* const y,
    const Real_t* const z,
    Real_t b[][8], Real_t* const volume) {
  const Real_t x0 = x[0] ;   const Real_t x1 = x[1] ;
  const Real_t x2 = x[2] ;   const Real_t x3 = x[3] ;
  const Real_t x4 = x[4] ;   const Real_t x5 = x[5] ;
  const Real_t x6 = x[6] ;   const Real_t x7 = x[7] ;

  const Real_t y0 = y[0] ;   const Real_t y1 = y[1] ;
  const Real_t y2 = y[2] ;   const Real_t y3 = y[3] ;
  const Real_t y4 = y[4] ;   const Real_t y5 = y[5] ;
  const Real_t y6 = y[6] ;   const Real_t y7 = y[7] ;

  const Real_t z0 = z[0] ;   const Real_t z1 = z[1] ;
  const Real_t z2 = z[2] ;   const Real_t z3 = z[3] ;
  const Real_t z4 = z[4] ;   const Real_t z5 = z[5] ;
  const Real_t z6 = z[6] ;   const Real_t z7 = z[7] ;

  Real_t fjxxi, fjxet, fjxze;
  Real_t fjyxi, fjyet, fjyze;
  Real_t fjzxi, fjzet, fjzze;
  Real_t cjxxi, cjxet, cjxze;
  Real_t cjyxi, cjyet, cjyze;
  Real_t cjzxi, cjzet, cjzze;

  fjxxi = .125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
  fjxet = .125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
  fjxze = .125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

  fjyxi = .125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
  fjyet = .125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
  fjyze = .125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

  fjzxi = .125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
  fjzet = .125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
  fjzze = .125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );

  /* compute cofactors */
  cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
  cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
  cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

  cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
  cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
  cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

  cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
  cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
  cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

  /* calculate partials :
     this need only be done for l = 0,1,2,3   since , by symmetry ,
     (6,7,4,5) = - (0,1,2,3) .
   */
  b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
  b[0][1] =      cjxxi  -  cjxet  -  cjxze;
  b[0][2] =      cjxxi  +  cjxet  -  cjxze;
  b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
  b[0][4] = -b[0][2];
  b[0][5] = -b[0][3];
  b[0][6] = -b[0][0];
  b[0][7] = -b[0][1];

  b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
  b[1][1] =      cjyxi  -  cjyet  -  cjyze;
  b[1][2] =      cjyxi  +  cjyet  -  cjyze;
  b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
  b[1][4] = -b[1][2];
  b[1][5] = -b[1][3];
  b[1][6] = -b[1][0];
  b[1][7] = -b[1][1];

  b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
  b[2][1] =      cjzxi  -  cjzet  -  cjzze;
  b[2][2] =      cjzxi  +  cjzet  -  cjzze;
  b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
  b[2][4] = -b[2][2];
  b[2][5] = -b[2][3];
  b[2][6] = -b[2][0];
  b[2][7] = -b[2][1];

  /* calculate jacobian determinant (volume) */
  *volume = Real_t(8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
}

inline void
Domain::CalcElemShapeFunctionDerivatives(const Real_t  Pos[8][3], Real_t b[][8],
                                         Real_t* const volume) {
  const int x = 0;
  const int y = 1;
  const int z = 2;

  const Real_t x0 = Pos[0][x] ;   const Real_t x1 = Pos[1][x] ;
  const Real_t x2 = Pos[2][x] ;   const Real_t x3 = Pos[3][x] ;
  const Real_t x4 = Pos[4][x] ;   const Real_t x5 = Pos[5][x] ;
  const Real_t x6 = Pos[6][x] ;   const Real_t x7 = Pos[7][x] ;

  const Real_t y0 = Pos[0][y] ;   const Real_t y1 = Pos[1][y] ;
  const Real_t y2 = Pos[2][y] ;   const Real_t y3 = Pos[3][y] ;
  const Real_t y4 = Pos[4][y] ;   const Real_t y5 = Pos[5][y] ;
  const Real_t y6 = Pos[6][y] ;   const Real_t y7 = Pos[7][y] ;

  const Real_t z0 = Pos[0][z] ;   const Real_t z1 = Pos[1][z] ;
  const Real_t z2 = Pos[2][z] ;   const Real_t z3 = Pos[3][z] ;
  const Real_t z4 = Pos[4][z] ;   const Real_t z5 = Pos[5][z] ;
  const Real_t z6 = Pos[6][z] ;   const Real_t z7 = Pos[7][z] ;

  Real_t fjxxi, fjxet, fjxze;
  Real_t fjyxi, fjyet, fjyze;
  Real_t fjzxi, fjzet, fjzze;
  Real_t cjxxi, cjxet, cjxze;
  Real_t cjyxi, cjyet, cjyze;
  Real_t cjzxi, cjzet, cjzze;

  fjxxi = .125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
  fjxet = .125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
  fjxze = .125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

  fjyxi = .125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
  fjyet = .125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
  fjyze = .125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

  fjzxi = .125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
  fjzet = .125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
  fjzze = .125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );

  /* compute cofactors */
  cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
  cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
  cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

  cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
  cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
  cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

  cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
  cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
  cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

  /* calculate partials :
     this need only be done for l = 0,1,2,3   since , by symmetry ,
     (6,7,4,5) = - (0,1,2,3) .
   */
  b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
  b[0][1] =      cjxxi  -  cjxet  -  cjxze;
  b[0][2] =      cjxxi  +  cjxet  -  cjxze;
  b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
  b[0][4] = -b[0][2];
  b[0][5] = -b[0][3];
  b[0][6] = -b[0][0];
  b[0][7] = -b[0][1];

  b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
  b[1][1] =      cjyxi  -  cjyet  -  cjyze;
  b[1][2] =      cjyxi  +  cjyet  -  cjyze;
  b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
  b[1][4] = -b[1][2];
  b[1][5] = -b[1][3];
  b[1][6] = -b[1][0];
  b[1][7] = -b[1][1];

  b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
  b[2][1] =      cjzxi  -  cjzet  -  cjzze;
  b[2][2] =      cjzxi  +  cjzet  -  cjzze;
  b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
  b[2][4] = -b[2][2];
  b[2][5] = -b[2][3];
  b[2][6] = -b[2][0];
  b[2][7] = -b[2][1];

  /* calculate jacobian determinant (volume) */
  *volume = Real_t(8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
}

inline void
Domain::SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
    Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
    Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
    Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
    const Real_t x0, const Real_t y0, const Real_t z0,
    const Real_t x1, const Real_t y1, const Real_t z1,
    const Real_t x2, const Real_t y2, const Real_t z2,
    const Real_t x3, const Real_t y3, const Real_t z3) {
  Real_t bisectX0 = Real_t(0.5) * (x3 + x2 - x1 - x0);
  Real_t bisectY0 = Real_t(0.5) * (y3 + y2 - y1 - y0);
  Real_t bisectZ0 = Real_t(0.5) * (z3 + z2 - z1 - z0);
  Real_t bisectX1 = Real_t(0.5) * (x2 + x1 - x3 - x0);
  Real_t bisectY1 = Real_t(0.5) * (y2 + y1 - y3 - y0);
  Real_t bisectZ1 = Real_t(0.5) * (z2 + z1 - z3 - z0);
  Real_t areaX = Real_t(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
  Real_t areaY = Real_t(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
  Real_t areaZ = Real_t(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

  *normalX0 += areaX;
  *normalX1 += areaX;
  *normalX2 += areaX;
  *normalX3 += areaX;

  *normalY0 += areaY;
  *normalY1 += areaY;
  *normalY2 += areaY;
  *normalY3 += areaY;

  *normalZ0 += areaZ;
  *normalZ1 += areaZ;
  *normalZ2 += areaZ;
  *normalZ3 += areaZ;
}

inline void
Domain::CollectDomainNodesToElemNodes(const Index_t* elemToNode, Real_t elemX[8],
    Real_t elemY[8], Real_t elemZ[8]) {
  Index_t nd0i = elemToNode[0] ;
  Index_t nd1i = elemToNode[1] ;
  Index_t nd2i = elemToNode[2] ;
  Index_t nd3i = elemToNode[3] ;
  Index_t nd4i = elemToNode[4] ;
  Index_t nd5i = elemToNode[5] ;
  Index_t nd6i = elemToNode[6] ;
  Index_t nd7i = elemToNode[7] ;

  elemX[0] = domain->m_x[nd0i];
  elemX[1] = domain->m_x[nd1i];
  elemX[2] = domain->m_x[nd2i];
  elemX[3] = domain->m_x[nd3i];
  elemX[4] = domain->m_x[nd4i];
  elemX[5] = domain->m_x[nd5i];
  elemX[6] = domain->m_x[nd6i];
  elemX[7] = domain->m_x[nd7i];

  elemY[0] = domain->m_y[nd0i];
  elemY[1] = domain->m_y[nd1i];
  elemY[2] = domain->m_y[nd2i];
  elemY[3] = domain->m_y[nd3i];
  elemY[4] = domain->m_y[nd4i];
  elemY[5] = domain->m_y[nd5i];
  elemY[6] = domain->m_y[nd6i];
  elemY[7] = domain->m_y[nd7i];

  elemZ[0] = domain->m_z[nd0i];
  elemZ[1] = domain->m_z[nd1i];
  elemZ[2] = domain->m_z[nd2i];
  elemZ[3] = domain->m_z[nd3i];
  elemZ[4] = domain->m_z[nd4i];
  elemZ[5] = domain->m_z[nd5i];
  elemZ[6] = domain->m_z[nd6i];
  elemZ[7] = domain->m_z[nd7i];
}

inline Real_t
Domain::CalcElemVolume( const Real_t Pos[8][3] )
{
  const int x = 0;
  const int y = 1;
  const int z = 2;
  Real_t twelveth = Real_t(1.0)/Real_t(12.0);


  Real_t dx61 = Pos[6][x] - Pos[1][x];
  Real_t dy61 = Pos[6][y] - Pos[1][y];
  Real_t dz61 = Pos[6][z] - Pos[1][z];

  Real_t dx70 = Pos[7][x] - Pos[0][x];
  Real_t dy70 = Pos[7][y] - Pos[0][y];
  Real_t dz70 = Pos[7][z] - Pos[0][z];

  Real_t dx63 = Pos[6][x] - Pos[3][x];
  Real_t dy63 = Pos[6][y] - Pos[3][y];
  Real_t dz63 = Pos[6][z] - Pos[3][z];

  Real_t dx20 = Pos[2][x] - Pos[0][x];
  Real_t dy20 = Pos[2][y] - Pos[0][y];
  Real_t dz20 = Pos[2][z] - Pos[0][z];

  Real_t dx50 = Pos[5][x] - Pos[0][x];
  Real_t dy50 = Pos[5][y] - Pos[0][y];
  Real_t dz50 = Pos[5][z] - Pos[0][z];

  Real_t dx64 = Pos[6][x] - Pos[4][x];
  Real_t dy64 = Pos[6][y] - Pos[4][y];
  Real_t dz64 = Pos[6][z] - Pos[4][z];

  Real_t dx31 = Pos[3][x] - Pos[1][x];
  Real_t dy31 = Pos[3][y] - Pos[1][y];
  Real_t dz31 = Pos[3][z] - Pos[1][z];

  Real_t dx72 = Pos[7][x] - Pos[2][x];
  Real_t dy72 = Pos[7][y] - Pos[2][y];
  Real_t dz72 = Pos[7][z] - Pos[2][z];

  Real_t dx43 = Pos[4][x] - Pos[3][x];
  Real_t dy43 = Pos[4][y] - Pos[3][y];
  Real_t dz43 = Pos[4][z] - Pos[3][z];

  Real_t dx57 = Pos[5][x] - Pos[7][x];
  Real_t dy57 = Pos[5][y] - Pos[7][y];
  Real_t dz57 = Pos[5][z] - Pos[7][z];

  Real_t dx14 = Pos[1][x] - Pos[4][x];
  Real_t dy14 = Pos[1][y] - Pos[4][y];
  Real_t dz14 = Pos[1][z] - Pos[4][z];

  Real_t dx25 = Pos[2][x] - Pos[5][x];
  Real_t dy25 = Pos[2][y] - Pos[5][y];
  Real_t dz25 = Pos[2][z] - Pos[5][z];

#define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

  Real_t volume =
    TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
           dy31 + dy72, dy63, dy20,
           dz31 + dz72, dz63, dz20) +
    TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
           dy43 + dy57, dy64, dy70,
           dz43 + dz57, dz64, dz70) +
    TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
           dy14 + dy25, dy61, dy50,
           dz14 + dz25, dz61, dz50);

#undef TRIPLE_PRODUCT

  volume *= twelveth;

  return volume ;
}

inline void
Domain::VoluDer(const Real_t x0, const Real_t x1, const Real_t x2,
    const Real_t x3, const Real_t x4, const Real_t x5,
    const Real_t y0, const Real_t y1, const Real_t y2,
    const Real_t y3, const Real_t y4, const Real_t y5,
    const Real_t z0, const Real_t z1, const Real_t z2,
    const Real_t z3, const Real_t z4, const Real_t z5,
    Real_t* dvdx, Real_t* dvdy, Real_t* dvdz) {
  const Real_t twelfth = Real_t(1.0) / Real_t(12.0) ;

  *dvdx = (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
    (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
    (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5);

  *dvdy = - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
    (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
    (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5);

  *dvdz = - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
    (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
    (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5);
  *dvdx *= twelfth;
  *dvdy *= twelfth;
  *dvdz *= twelfth;
}

inline Real_t
Domain::AreaFace(const Real_t Pos[8][3], const int c0,
                 const int c1, const int c2, const int c3) {
  const int x = 0;
  const int y = 1;
  const int z = 2;

  Real_t fx = (Pos[c2][x] - Pos[c0][x]) - (Pos[c3][x] - Pos[c1][x]);
  Real_t fy = (Pos[c2][y] - Pos[c0][y]) - (Pos[c3][y] - Pos[c1][y]);
  Real_t fz = (Pos[c2][z] - Pos[c0][z]) - (Pos[c3][z] - Pos[c1][z]);

  Real_t gx = (Pos[c2][x] - Pos[c0][x]) + (Pos[c3][x] - Pos[c1][x]);
  Real_t gy = (Pos[c2][y] - Pos[c0][y]) + (Pos[c3][y] - Pos[c1][y]);
  Real_t gz = (Pos[c2][z] - Pos[c0][z]) + (Pos[c3][z] - Pos[c1][z]);

  Real_t area = (fx * fx + fy * fy + fz * fz) *
          (gx * gx + gy * gy + gz * gz) -
          (fx * gx + fy * gy + fz * gz) *
          (fx * gx + fy * gy + fz * gz);
  return area ;
}

