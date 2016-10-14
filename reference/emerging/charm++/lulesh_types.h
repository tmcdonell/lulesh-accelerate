#ifndef __LULESH_TYPES_H__
#define __LULESH_TYPES_H__
#include <math.h>

// Allow flexibility for arithmetic representations
// Could also support fixed point and interval arithmetic types
typedef float       real4;
typedef double      real8;
typedef long double real10; // 10 bytes on x86

typedef int   Index_t;  // array subscript and loop index
typedef real8 Real_t;   // floating point representation
typedef int   Int_t;    // integer representation

inline real4  SQRT(real4  arg) { return sqrtf(arg); }
inline real8  SQRT(real8  arg) { return sqrt(arg);  }
inline real10 SQRT(real10 arg) { return sqrtl(arg); }

inline real4  CBRT(real4  arg) { return cbrtf(arg); }
inline real8  CBRT(real8  arg) { return cbrt(arg);  }
inline real10 CBRT(real10 arg) { return cbrtl(arg); }

inline real4  FABS(real4  arg) { return fabsf(arg); }
inline real8  FABS(real8  arg) { return fabs(arg);  }
inline real10 FABS(real10 arg) { return fabsl(arg); }

struct Dom {
  // Simulation Time
  Real_t deltaTime;
  Real_t totalTime;
  // Ghosts
  std::vector<Real_t> ng_front;
  std::vector<Real_t> ng_back;
  std::vector<Real_t> ng_right;
  std::vector<Real_t> ng_left;
  std::vector<Real_t> ng_up;
  std::vector<Real_t> ng_down;
  // Node centered persistent
  std::vector<Real_t> m_x;
  std::vector<Real_t> m_y;
  std::vector<Real_t> m_z;
  std::vector<Real_t> m_xd;
  std::vector<Real_t> m_yd;
  std::vector<Real_t> m_zd;
  std::vector<Real_t> m_xdd;
  std::vector<Real_t> m_ydd;
  std::vector<Real_t> m_zdd;
  std::vector<Real_t> m_fx;
  std::vector<Real_t> m_fy;
  std::vector<Real_t> m_fz;
  std::vector<Real_t> m_nodalMass;
  // Node centered nodesets
  std::vector<Index_t> m_symmX;
  std::vector<Index_t> m_symmY;
  std::vector<Index_t> m_symmZ;
  // Elem centered persistent
  std::vector<Index_t> m_matElemlist;
  std::vector<Index_t> m_nodelist;
  std::vector<Index_t> m_lxim;
  std::vector<Index_t> m_lxip;
  std::vector<Index_t> m_letam;
  std::vector<Index_t> m_letap;
  std::vector<Index_t> m_lzetam;
  std::vector<Index_t> m_lzetap;
  std::vector<Int_t> m_elemBC;
  std::vector<Real_t> m_e;
  std::vector<Real_t> m_p;
  std::vector<Real_t> m_q;
  std::vector<Real_t> m_ql;
  std::vector<Real_t> m_qq;
  std::vector<Real_t> m_v;
  std::vector<Real_t> m_volo;
  std::vector<Real_t> m_delv;
  std::vector<Real_t> m_vdov;
  std::vector<Real_t> m_arealg;
  std::vector<Real_t> m_ss;
  std::vector<Real_t> m_elemMass;
  // Elem centered temporary
  std::vector<Real_t> m_dxx;
  std::vector<Real_t> m_dyy;
  std::vector<Real_t> m_dzz;
  std::vector<Real_t> m_delv_xi;
  std::vector<Real_t> m_delv_eta;
  std::vector<Real_t> m_delv_zeta;
  std::vector<Real_t> m_delx_xi;
  std::vector<Real_t> m_delx_eta;
  std::vector<Real_t> m_delx_zeta;
  std::vector<Real_t> m_vnew;
  std::vector<Real_t> m_determ;
  // Parameters
  Real_t u_cut;
  Real_t hgcoef;
  Real_t qstop;
  Real_t monoq_max_slope;
  Real_t monoq_limiter_mult;
  Real_t e_cut;
  Real_t p_cut;
  Real_t ss4o3;
  Real_t q_cut;
  Real_t v_cut;
  Real_t qlc_monoq;
  Real_t qqc_monoq;
  Real_t qqc;
  Real_t eosvmax;
  Real_t eosvmin;
  Real_t pmin;
  Real_t emin;
  Real_t dvovmax;
  Real_t refdens;
  Real_t m_dtcourant;
  Real_t m_dthydro;
};

#endif //__LULESH_TYPES_H__
