#ifndef __LULESH_DEFS_H__
#define __LULESH_DEFS_H__


//
// Application Configuration Options
//
#define CHARM_DEBUG false
#define INIT_DEBUG false
#define CHAPEL_MESH false
#define CHECKPOINT_DEBUG false

#define WRITE_TO_FILE false
#define LULESH_SHOW_PROGRESS false
#define ITER_TIMING false

#define APPLY_TIME_CONSTRAINTS true
#define USE_FIXED_DT true
#define MAX_ITERATIONS 500

#define PERFORM_LOAD_BALANCE false
#define LB_FREQUENCY 250

#define PERFORM_CHECKPOINT false
#define CP_FREQUENCY 350

//
// Application Specific Definitions
//
#define INITIAL_ENERGY 3.948746e+7
#define SSC_MIN        0.333333e-18
#define SS_MIN         1.111111e-36
#define P_TINY         1e-36
// Error
#define VOLUME_ERR "Volume Error"
#define QSTOP_ERR "Q Stop Error"
enum { VolumeError = -1, QStopError = -2 };
// Timing
#define STOP_TIME     1e-2
#define DT_FIXED      1e-7
#define DT_MAX        1e-2
#define DT_MULT_LB    1.1
#define DT_MULT_UB    1.2
#define DT_COURANT    1e+20
#define DT_HYDRO      1e+20
#define DT_HYDRO_INV  1e-20
// Mesh
#define MAX_NODE_POSITION 1.125
#define NODES_PER_ELEM    8
#define E_CUT             1e-7
#define P_CUT             1e-7
#define Q_CUT             1e-7
#define U_CUT             1e-7
#define V_CUT             1e-10
#define HG_COEF           3.0
#define SS_4O3            4.0/3.0
#define Q_STOP            1e+12
#define MONOQ_MAX_SLOPE   1.0
#define MONOQ_LIM_MULT    2.0
#define QLC_MONOQ         0.5
#define QQC_MONOQ         2.0/3.0
#define QQC               2.0
#define P_MIN             0.0
#define E_MIN            -1e+15
#define DVOV_MAX          0.1
#define EOSV_MAX          1e+9
#define EOSV_MIN          1e-9
#define REF_DENS          1.0
// Boundary Conditions
// 2 BCs on each of 6 hexahedral faces (12 bits)
#define XI_M        0x003
#define XI_M_SYMM   0x001
#define XI_M_FREE   0x002
#define XI_P        0x00c
#define XI_P_SYMM   0x004
#define XI_P_FREE   0x008

#define ETA_M       0x030
#define ETA_M_SYMM  0x010
#define ETA_M_FREE  0x020
#define ETA_P       0x0c0
#define ETA_P_SYMM  0x040
#define ETA_P_FREE  0x080

#define ZETA_M      0x300
#define ZETA_M_SYMM 0x100
#define ZETA_M_FREE 0x200
#define ZETA_P      0xc00
#define ZETA_P_SYMM 0x400
#define ZETA_P_FREE 0x800

#endif //__LULESH_DEFS_H__
