/*
 * Module       : visit
 * Copyright    : [2016] Trevor L. McDonell
 * License      : BSD3
 *
 * Maintainer   : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
 * Stability    : experimental
 * Portability  : non-portable (GHC extensions)
 *
 * The functions in this module can be used to write the result of each step of
 * the simulation to file. Data is written using the 'silo' library, which you
 * will need to install separately:
 *
 *   <https://wci.llnl.gov/simulation/computer-codes/silo>
 *
 * The resulting file(s) can then be used with the VisIt visualisation program:
 *
 *   <https://wci.llnl.gov/simulation/computer-codes/visit>
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include <silo.h>

/*
 * The type of mesh elements, as defined in src/Type.hs.
 */
#define DB_R DB_DOUBLE

typedef double    R;
typedef uint32_t  I;


static void writeDomainDB
(
    DBfile*     db,
    int         numElem,
    int         step,
    R           time,
    R*          x,            // nodal position
    R*          y,
    R*          z,
    R*          xd,           // nodal velocity
    R*          yd,
    R*          zd,
    R*          e,            // internal energy
    R*          p,            // pressure
    R*          v,            // relative volume
    R*          q             // viscosity
)
{
  int ok        = 0;
  int numNode   = numElem + 1;

  int numNode3  = numNode * numNode * numNode;
  int numElem3  = numElem * numElem * numElem;

  /*
   * Create an option list that will give some hints to VisIt for printing out
   * the cycle and time in the annotations
   */
  DBoptlist *optlist;

  /*
   * Write out the mesh connectivity in fully unstructured format.
   *
   * In the reference implementation the linear index of each of the 8
   * neighbouring nodes to an element were stored in the 'nodelist' array so
   * that they could be accessed conveniently. In Accelerate that is unnecessary
   * since we have multidimensional indexing built in. This function just
   * regenerates those indices into a flattened array.
   */
  int shapetype[1] = { DB_ZONETYPE_HEX };
  int shapesize[1] = { 8 };
  int shapecnt[1]  = { numElem3 };
  int *conn        = new int[ numElem3*8 ];

  int i = 0;
  int n = 0;

  for (int plane=0; plane < numElem; ++plane) {
    for (int row=0; row < numElem; ++row) {
      for (int col=0; col < numElem; ++col) {
        int *local = &conn[n];
        local[0] = i                                 ;
        local[1] = i                             + 1 ;
        local[2] = i                   + numNode + 1 ;
        local[3] = i                   + numNode     ;
        local[4] = i + numNode*numNode               ;
        local[5] = i + numNode*numNode           + 1 ;
        local[6] = i + numNode*numNode + numNode + 1 ;
        local[7] = i + numNode*numNode + numNode     ;

        n += 8;
        i += 1;
      }
      i += 1;
    }
    i += numNode ;
  }

  ok += DBPutZonelist2(db, "connectivity", numElem3, 3, conn, numElem3*8, 0,0,0, shapetype, shapesize, shapecnt, 1, NULL);
  delete [] conn ;

  /*
   * Write out the mesh coordinates associated with the mesh
   */
  const char  *coordnames[] = {"X", "Y", "Z"};
  R           *coords[3]    = { x, y, z };

  optlist = DBMakeOptlist(2);
  ok += DBAddOption(optlist, DBOPT_DTIME, &time);
  ok += DBAddOption(optlist, DBOPT_CYCLE, &step);
  ok += DBPutUcdmesh(db, "mesh", 3, coordnames, coords, numNode3, numElem3, "connectivity", 0, DB_R, optlist);
  ok += DBFreeOptlist(optlist);

  /*
   * Write out pressure, energy, relative volume, viscosity
   */
  ok += DBPutUcdvar1(db, "e", "mesh", e, numElem3, NULL, 0, DB_R, DB_ZONECENT, NULL);
  ok += DBPutUcdvar1(db, "p", "mesh", p, numElem3, NULL, 0, DB_R, DB_ZONECENT, NULL);
  ok += DBPutUcdvar1(db, "v", "mesh", v, numElem3, NULL, 0, DB_R, DB_ZONECENT, NULL);
  ok += DBPutUcdvar1(db, "q", "mesh", q, numElem3, NULL, 0, DB_R, DB_ZONECENT, NULL);

  /*
   * Write out nodal speed, velocities
   */
  R *speed = new R[numNode3];
  for(int i = 0; i < numNode3; ++i) {
    speed[i] = sqrt((xd[i]*xd[i])+(yd[i]*yd[i])+(zd[i]*zd[i]));
  }

  ok += DBPutUcdvar1(db, "xd",    "mesh", xd,    numNode3, NULL, 0, DB_R, DB_NODECENT, NULL);
  ok += DBPutUcdvar1(db, "yd",    "mesh", yd,    numNode3, NULL, 0, DB_R, DB_NODECENT, NULL);
  ok += DBPutUcdvar1(db, "zd",    "mesh", zd,    numNode3, NULL, 0, DB_R, DB_NODECENT, NULL);
  ok += DBPutUcdvar1(db, "speed", "mesh", speed, numNode3, NULL, 0, DB_R, DB_NODECENT, NULL);
  delete [] speed;

  if (ok != 0) {
    printf("Error writing output file\n");
  }
}

/*
 * Create a new silo output database and write all the domain data to it.
 */
extern "C" void writeDomain
(
    const char* basename,
    int         numElem,
    int         step,
    R           time,
    R*          x,            // nodal position
    R*          y,
    R*          z,
    R*          xd,           // nodal velocity
    R*          yd,
    R*          zd,
    R*          e,            // internal energy
    R*          p,            // pressure
    R*          v,            // relative volume
    R*          q             // viscosity
)
{
  DBfile *db = NULL;
  db         = (DBfile*) DBCreate(basename, DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5X);

  if (db) {
    writeDomainDB(db, numElem, step, time, x, y, z, xd, yd, zd, e, p, v, q);
    DBClose(db);
  }
  else {
    printf("Error creating output file\n");
  }
}

