#ifndef __LULESH_CHARM_H__
#define __LULESH_CHARM_H__
#include "lulesh_defs.h"
#include "lulesh_types.h"
#include "cycle.h"
#include "lulesh.decl.h"

// Reduction Operations
void registerMinReal(void);
CkReductionMsg *minReal(int nMsg, CkReductionMsg **msgs);

void registerAveTime(void);
CkReductionMsg *aveTime(int nMsg, CkReductionMsg **msgs);

class Main : public CBase_Main {
  private:
    // Member Variables
    CProxy_Domain domains;
    int totalIterations;
    std::string fileOut;
    // Timing Variables
    double startSimTime;
    double endSimTime;
    double simTime;
    double iterTimers[11];
    // Member Functions
    void writeToFile();
  public:
    // Constructors & Destructors
    Main(CkArgMsg *msg);
    Main(CkMigrateMessage *msg);
    void pup(PUP::er &p);

    // Entry Methods
    void averageTimeCheckin(CkReductionMsg *msg);
    void initCheckin(CkReductionMsg *msg);
};

class Domain : public CBase_Domain {
  private:
    Domain_SDAG_CODE
    // Clocking Variables
    int iterations;
    double startIterTime;
    double endIterTime;
    double iterTime;
#if ITER_TIMING
    // Timing Variables
    ticks startIterTicks;
    double iterTicks;
    ticks tempTime;
    double stressTime;
    double hgTime;
    double posTime;
    double kinTime;
    double viscTime;
    double matTime;
    double volTime;
    double consTime;
#endif
    // Member Variables
    // Charm Communication
    int numElemGhosts;
    int numNodeGhosts;
    int ghostNodeCount;
    int ghostElemCount;
    // Domain Information
    struct Dom *domain;
    // Member Functions
    // Domain Control Flow
    void InitializeLuleshData();
    void LagrangeNodal_Forces();
    void LagrangeNodal_Positions();
    void LagrangeElements_KinVisc();
    void LagrangeElements_QEOS();
    void CalcTimeConstraintsForElems();
    // Nodal Functions
    void InitStressTermsForElems(Index_t k, Real_t *sigxx,
                                 Real_t *sigyy, Real_t *sigzz);
    void IntegrateStressForElems();
    void CalcAccelerationForNodes();
    void ApplyAccelerationBoundaryConditionsForNodes();
    // Elemental Functions
    void CalcKinematicsForElems(const Index_t k, const Real_t Pos[8][3],
                                const Real_t Vel[8][3], const Real_t dt);
    void CalcMonotonicQGradientsForElems(Index_t k, const Real_t Pos[8][3],
                                         const Real_t Vel[8][3]);
    void CalcMonotonicQRegionForElems(const Index_t k, const Real_t qlc_monoq,
                                      const Real_t qqc_monoq,
                                      const Real_t monoq_limiter_mult,
                                      const Real_t monoq_max_slope,
                                      const Real_t ptiny, const Real_t qstop);
    void EvalEOSForElems(const Index_t k, const Real_t e_cut,
                         const Real_t p_cut, const Real_t q_cut,
                         const Real_t eosvmin, const Real_t eosvmax,
                         const Real_t pmin, const Real_t emin,
                         const Real_t rho0);
    void UpdateVolumesForElems(const Index_t k);
    // Timing Functions
    void CalcCourantConstraintForElems();
    void CalcHydroConstraintForElems();
    // Calculation Functions (nodal)
    void CalcElemNodeNormals(Real_t pfx[8], Real_t pfy[8], Real_t pfz[8],
                             const Real_t x[8], const Real_t y[8], const Real_t z[8]);
    void SumElemStressesToNodeForces(const Real_t B[3][8],
                                     const Real_t stress_xx,
                                     const Real_t stress_yy,
                                     const Real_t stress_zz,
                                     Real_t* const fx,
                                     Real_t* const fy,
                                     Real_t* const fz);
    void CalcElemVolumeDerivative(Real_t dvdx[8], Real_t dvdy[8], Real_t dvdz[8],
                                  const Real_t x[8], const Real_t y[8],
                                  const Real_t z[8]);
    void CalcFBHourglassForceForElems(Real_t hourg);
    void CalcElemFBHourglassForce(Real_t *xd, Real_t *yd, Real_t *zd,
                                  Real_t hourgam[8][4], Real_t coefficient,
                                  Real_t *hgfx, Real_t *hgfy, Real_t *hgfz);
    // Calculation Functions (elemental)
    Real_t CalcElemCharacteristicLength(const Real_t Pos[8][3], const Real_t volume);
    void CalcElemVelocityGradient(const Real_t vel[8][3], const Real_t b[][8],
                                  const Real_t detJ, Real_t* const d );
    void CalcEnergyForElems(const Real_t p_old, const Real_t e_old,
                            const Real_t q_old, const Real_t compression,
                            const Real_t compHalfStep, const Real_t vnewc,
                            const Real_t work, const Real_t delvc,
                            const Real_t pmin, const Real_t p_cut,
                            const Real_t e_cut, const Real_t q_cut,
                            const Real_t emin, const Real_t qq, const Real_t ql,
                            const Real_t rho0, const Real_t eosvmax,
                            Real_t *p_new, Real_t *e_new, Real_t *q_new,
                            Real_t *bvc, Real_t *pbvc);
    void CalcPressureForElems(const Real_t e_old, const Real_t compression,
                              const Real_t vnewc, const Real_t pmin,
                              const Real_t p_cut, const Real_t eosvmax,
                              Real_t *p_new, Real_t *bvc, Real_t *pbvc);
    Real_t CalcSoundSpeedForElem(const Real_t new_vol, const Real_t enew,
                              const Real_t pnew, const Real_t pbvc,
                              const Real_t bvc, const Real_t rho0);
    // Calculation Functions (helper/shared)
    void CalcElemShapeFunctionDerivatives(const Real_t* const x,
                                          const Real_t* const y,
                                          const Real_t* const z,
                                          Real_t b[][8], Real_t* const volume);
    void CalcElemShapeFunctionDerivatives(const Real_t  Pos[8][3], Real_t b[][8],
                                          Real_t* const volume);
    void SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
                           Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
                           Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
                           Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
                           const Real_t x0, const Real_t y0, const Real_t z0,
                           const Real_t x1, const Real_t y1, const Real_t z1,
                           const Real_t x2, const Real_t y2, const Real_t z2,
                           const Real_t x3, const Real_t y3, const Real_t z3);
    void CollectDomainNodesToElemNodes(const Index_t* elemToNode, Real_t elemX[8],
                                       Real_t elemY[8], Real_t elemZ[8]);
    Real_t CalcElemVolume(const Real_t Pos[8][3]);
    void VoluDer(const Real_t x0, const Real_t x1, const Real_t x2,
                 const Real_t x3, const Real_t x4, const Real_t x5,
                 const Real_t y0, const Real_t y1, const Real_t y2,
                 const Real_t y3, const Real_t y4, const Real_t y5,
                 const Real_t z0, const Real_t z1, const Real_t z2,
                 const Real_t z3, const Real_t z4, const Real_t z5,
                 Real_t* dvdx, Real_t* dvdy, Real_t* dvdz);
    Real_t AreaFace(const Real_t Pos[8][3], const int c0,
                 const int c1, const int c2, const int c3);
    // Load Balancing
    void ResumeFromSync();
  public:
    // Constructors & Destructors
    Domain();
    Domain(CkMigrateMessage *msg);
    ~Domain();
    void pup(PUP::er &p);
    void startLB();
    
    // Entry Methods
    void beginIteration();
    void resumeNodeIteration();
    void resumeElemIteration();
    void sendNodeGhosts();
    void sendElemGhosts();
    void sendNodalMass();
    void processNodeGhosts(int dir, int width, int height,
                       Real_t ghostData[]);
    void updateForceGhosts();
    void processElemGhosts(int dir, int width, int height,
                       Real_t ghostData[]);
    void receiveNodalMass(int dir, int width, int height,
                       Real_t ghostData[]);
    void averageTime();
    void printEnergy();
};

#endif //__LULESH_CHARM_H__

