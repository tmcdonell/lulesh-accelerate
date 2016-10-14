/*

                 Copyright (c) 2010.
      Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
                  LLNL-CODE-461231
                All rights reserved.

This file is part of LULESH, Version 1.0.
Please also read this link -- http://www.opensource.org/licenses/index.php

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.



Additional Language / Algorithm Notice

   This code is written in an array notation language called A++.  A++ is an array class
   library that is part of the ROSE Project.  ROSE is a tool for building source-to-source
   translators (sometimes confused with a preprocessor which is generally less sophisticated
   internally). ROSE is particularly useful in building custom tools that operate on source
   code for C, C99, C++, F77, F90, & F2003.  The primary user of A++ is the Overture Project
   and most of the documentation on A++ can be found there.  In particular, The Reference
   Manual and installation instructions for A++ can be found at:

   https://computation.llnl.gov/casc/Overture/

   in the Documents section of the web page.

   This version of the LULESH code avoids the use of global scope of accessing all arrays and
   variables.  Most functions receive a parameter list of all input-only, input-output, and
   output objects to help clarify data-dependencies and explain what each algorithm actually does.

*/


#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define BOUNDS_CHECK
#include "A++.h"



typedef float        real4 ;
typedef double       real8 ;
typedef long double  real10 ;  // 10 bytes on x86

typedef int    Index_t ;      // array subscript and loop index
typedef real8  Real_t ;       // floating point representation
typedef int    Int_t ;        // integer representation

enum { VolumeError   = -1,
        QStopError   = -2,
        FaceError    = -3,
        MatMultError = -4,
        PhiCalcError = -5   } ;


//  LULESH Constants

    //  Array Size Constants:        these constants will be accessed globally
    //                               NOT passed by parameter.   The first constant,
    //                               edgeElems, sets the problem size for the entire code.

    int edgeElems         = 45 ;
    int edgeNodes         = edgeElems + 1;
    int Node_count        = edgeNodes * edgeNodes * edgeNodes;
    int Elem_count        = edgeElems * edgeElems * edgeElems;
    int BCface_count      = edgeNodes * edgeNodes;
    int matZonelist_count = Elem_count;

    //  Array Index Constants:       these constants will be accessed globally
    //                               NOT passed by parameter.

    Index all;                                   //  empty index set matches whatever range needed
    Index    Node_index (1,Node_count);          //  Index of nodes
    Index    Elem_index (1,Elem_count);          //  Index of elements
    Index      BC_index (1,BCface_count);        //  Index of nodes on a boundary face
    Index   Coord_index (1,3);                   //  Index of X, Y, Z
    Index     Dir_index (1,3);                   //  Index of eta, xi, zeta
    Index  Strain_Modes (1,4);                   //  Index of hourglass modes
    Index      PM_index (1,2);                   //  Index of minus, plus
    Index    Face_crnrs (1,4);                   //
    Index    Face_index (1,6);                   //
    Index  Corner_index (1,8);                   //  corners of an element
    Index  MaterialList (1,matZonelist_count);   //  material index set
    Index   Dummy_index (1,1);                   //  Needed in A++ to change array dimensionality
    Index    Test_index (1,10);                  //  Debug tool

    Index    Elem3D_index (1,edgeElems);         //  Used for debug only
    Index    Node3D_index (1,edgeNodes);         //  Used for debug only


    //  Array Constants:      these constants will be accessed globally
    //                        NOT passed by parameter.


    doubleArray   Gamma (Corner_index, Strain_Modes);   //  coefficient matrix

    //  Material constants:  denoted by mc_ : these constants will be accessed globally
    //                                        NOT passed by parameter.

    Real_t  mc_dtfixed             = Real_t(-1.0e-7);           // fixed time increment
    Real_t  mc_dtmax               = Real_t(1.0e-2);            // maximum allowable time increment
    Real_t  mc_stoptime            = Real_t(1.0e-2);            // end time for simulation
    Real_t  mc_deltatimemultlb     = Real_t(1.1) ;
    Real_t  mc_deltatimemultub     = Real_t(1.2) ;
    Real_t  mc_e_cut               = Real_t(1.0e-7)  ;          // energy tolerance
    Real_t  mc_p_cut               = Real_t(1.0e-7)  ;          // pressure tolerance
    Real_t  mc_q_cut               = Real_t(1.0e-7)  ;          // q tolerance   (art. viscocity)
    Real_t  mc_u_cut               = Real_t(1.0e-7)  ;          // velocity tolerance
    Real_t  mc_v_cut               = Real_t(1.0e-10) ;          // relative volume tolerance
    Real_t  mc_ss4o3               = Real_t(4.0)/Real_t(3.0) ;
    Real_t  mc_hgcoef              = Real_t(3.0) ;              // hourglass control
    Real_t  mc_qstop               = Real_t(1.0e+12) ;          // excessive q indicator
    Real_t  mc_monoq_max_slope     = Real_t(1.0) ;
    Real_t  mc_monoq_limiter_mult  = Real_t(2.0) ;
    Real_t  mc_qlc_monoq           = Real_t(0.5) ;              // linear term coef for q
    Real_t  mc_qqc_monoq           = Real_t(2.0)/Real_t(3.0) ;  // quadratic term coef for q
    Real_t  mc_qqc                 = Real_t(2.0);
    Real_t  mc_eosvmax             = Real_t(1.0e+9);
    Real_t  mc_eosvmin             = Real_t(1.0e-9);
    Real_t  mc_pmin                = Real_t(0.) ;               // pressure floor
    Real_t  mc_emin                = Real_t(-1.0e+15) ;         // energy floor
    Real_t  mc_dvovmax             = Real_t(0.1) ;              // maximum allowable volume change
    Real_t  mc_refdens             = Real_t(1.0) ;              // reference density

   // Stuff needed for boundary conditions
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

//  LULESH Program Variables :   all of these quantities will be passed as parameters to all
//                               functions that use them, whether read-only, r/w, or write-only

    // Iteration-related Quantities

    Real_t  dtcourant              = Real_t(1.0e+20) ;          // courant constraint
    Real_t  dthydro                = Real_t(1.0e+20) ;          // volume change constraint
    Real_t  CurrentTime            = Real_t(0.);                // current time
    Int_t   cycle                  = 0 ;                        // iteration count for simulation
    Real_t  deltatime              = Real_t(1.0e-7);            // variable time increment

    //  Nodal Quantities

    doubleArray  coords (Node_index, Coord_index);               // X, Y, Z coordinates
    doubleArray     vel (Node_index, Coord_index);               // X, Y, Z velocities
    doubleArray   accel (Node_index, Coord_index);               // X, Y, Z accelerations
    doubleArray   force (Node_index, Coord_index);               // X, Y, Z Forces
    doubleArray  n_mass (Node_index);                            // nodal masses

    // Zonal Quantities

    intArray     ElemFaceToElem (Elem_index, Dir_index, PM_index);   //  map elem face to adj. elem
    intArray         ElemToNode (Elem_index, Corner_index);          //  map each elem to 8 nodes
    intArray            BCnodes (BC_index,    Coord_index);          //  map each elem to 8 nodes
    doubleArray         delta_v (Elem_index);                        //  old mesh.delv
    doubleArray          energy (Elem_index);                        //  old mesh.e
    doubleArray        pressure (Elem_index);                        //  old mesh.p
    doubleArray               q (Elem_index);                        //  old mesh.q
    doubleArray              ql (Elem_index);                        //  old mesh.ql
    doubleArray              qq (Elem_index);                        //  old mesh.qq
    doubleArray         rel_vol (Elem_index);                        //  old mesh.v  relative volume
    doubleArray         ref_vol (Elem_index);                        //  old mesh.volo    ref volume
    doubleArray         new_vol (Elem_index);                        //  old mesh.vnew    new volume
    doubleArray            vdov (Elem_index);                        //  old mesh.vdov
    doubleArray        char_len (Elem_index);                        //  old mesh.arealg
    doubleArray              ss (Elem_index);                        //  old mesh.ss - sound speed
    doubleArray          Z_mass (Elem_index);                        //  old mesh.zonalMass
    intArray             elemBC (Elem_index);                        //  bit map indicates elem BCs
    intArray        matElemlist (MaterialList);                      //  material index set



//  Utility functions  for common math ops and accessing element info (like coords. of an element)


static inline
doubleArray cross( const doubleArray &a,
                   const doubleArray &b  )
{
   doubleArray c (Coord_index);
   c(1) = a(2)*b(3) - b(2)*a(3);
   c(2) = a(3)*b(1) - b(3)*a(1);
   c(3) = a(1)*b(2) - b(1)*a(2);
   return c;
}


static inline
void CollectElemCoords(const int         &ElemId,
                       const intArray    &ElemToNode,
                       const doubleArray &coords,
                             doubleArray &ElemCoords)
{
  for( Index_t lnode=1 ; lnode<9 ; ++lnode )
    {
      Index_t          gnode = ElemToNode(ElemId,lnode  );
      ElemCoords(lnode, all) = coords(gnode, all);
   }
}


static inline
void CollectElemVelocities(const int         &ElemId,
                           const intArray    &ElemToNode,
                           const doubleArray &vel,
                                 doubleArray &ElemVels)
{
    for( Index_t lnode=1 ; lnode<9 ; ++lnode )
    {
      Index_t        gnode = ElemToNode(ElemId,lnode  );
      ElemVels(lnode, all) = vel(gnode, all);
    }
}


static inline
void Accum_Forces_To_Nodes(const int         &ElemId,
                           const doubleArray &elem_force,
                                 doubleArray &force)
{
   for( Index_t lnode=1 ; lnode<9 ; ++lnode )
    {
      Index_t        gnode = ElemToNode(ElemId,lnode  );
      force(gnode,all) +=  elem_force(lnode,all);
    }
}


static inline
doubleArray GetFace( const doubleArray  &elemCoords,
                     const Index_t      &FaceId      )
{
    doubleArray FaceCoords ( Face_crnrs, Coord_index);

    switch (FaceId) {
       case 1:  FaceCoords(1,all) = elemCoords(1, all);    // Face 1:  Corners 1,2,3,4
                FaceCoords(2,all) = elemCoords(2, all);
                FaceCoords(3,all) = elemCoords(3, all);
                FaceCoords(4,all) = elemCoords(4, all);
                break;
       case 2:  FaceCoords(1,all) = elemCoords(5, all);    // Face 2:  Corners 5,8,7,6
                FaceCoords(2,all) = elemCoords(8, all);
                FaceCoords(3,all) = elemCoords(7, all);
                FaceCoords(4,all) = elemCoords(6, all);
                break;
       case 3:  FaceCoords(1,all) = elemCoords(4, all);    // Face 3:  Corners 4,8,5,1
                FaceCoords(2,all) = elemCoords(8, all);
                FaceCoords(3,all) = elemCoords(5, all);
                FaceCoords(4,all) = elemCoords(1, all);
                break;
       case 4:  FaceCoords(1,all) = elemCoords(2, all);    // Face 4:  Corners 2,6,7,3
                FaceCoords(2,all) = elemCoords(6, all);
                FaceCoords(3,all) = elemCoords(7, all);
                FaceCoords(4,all) = elemCoords(3, all);
                break;
       case 5:  FaceCoords(1,all) = elemCoords(1, all);    // Face 5:  Corners 1,5,6,2
                FaceCoords(2,all) = elemCoords(5, all);
                FaceCoords(3,all) = elemCoords(6, all);
                FaceCoords(4,all) = elemCoords(2, all);
                break;
       case 6:  FaceCoords(1,all) = elemCoords(3, all);    // Face 6:  Corners 3,7,8,4
                FaceCoords(2,all) = elemCoords(7, all);
                FaceCoords(3,all) = elemCoords(8, all);
                FaceCoords(4,all) = elemCoords(4, all);
                break;
      default:  exit(FaceError);
    }
    return FaceCoords;
}


static inline
doubleArray AveFaceLocs(const doubleArray &A)
{
    doubleArray Ans ;
    Ans = (A(1, all) + A(2, all) + A(3, all) + A(4, all)) * Real_t(0.25);
    Ans.reshape(Coord_index);
    return Ans;
}


static inline
doubleArray MatMult(const doubleArray &A,
                    const doubleArray &B)
{
    Index IndexA ( 1, A.rows() );
    Index IndexB ( 1, B.cols() );
    doubleArray C (IndexA, IndexB);

    C = Real_t(0.0);
    if (A.cols() == B.rows() )
       { for (int i = 1; i <= A.rows(); ++ i) {
             for (int j = 1; j <= B.cols(); ++ j) {
                  C(i,j) = Real_t(0.0);
                  for (int k = 1; k<=A.cols(); k++) {
                      C(i,j) += A(i,k) * B(k,j);
                  }
             }
         }
        }
    else   {
        printf("MatMult error:  mismatched arrays \n");
        A.display("A array\n");
        B.display("B array\n");
        exit(MatMultError);
        }
    return C;
}


//  Start of LULESH specific code


static
void TimeIncrement(   const Real_t  &dtcourant,
                      const Real_t  &dthydro,
                            Real_t  &deltatime,
                            Real_t  &CurrentTime,
                            Int_t   &cycle)
{
   Real_t targetdt = mc_stoptime - CurrentTime ;

   if ((mc_dtfixed <= Real_t(0.0)) && (cycle != Int_t(0))) {
      Real_t ratio ;
      Real_t olddt = deltatime ;

      /* This will require a reduction in parallel */
      Real_t newdt = Real_t(1.0e+20) ;
      if (dtcourant < newdt) {
         newdt = dtcourant / Real_t(2.0) ;
      }
      if (dthydro < newdt) {
         newdt = dthydro * Real_t(2.0) / Real_t(3.0) ;
      }

      ratio = newdt / olddt ;
      if (ratio >= Real_t(1.0)) {
         if (ratio < mc_deltatimemultlb) {
            newdt = olddt ;
         }
         else if (ratio >   mc_deltatimemultub) {
            newdt = olddt * mc_deltatimemultub ;
         }
      }

      if (newdt > mc_dtmax) {
         newdt = mc_dtmax ;
      }
      deltatime = newdt ;
   }

   /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
   if ((targetdt > deltatime) &&
       (targetdt < (Real_t(4.0) * deltatime / Real_t(3.0))) ) {
        targetdt = (Real_t(2.0) * deltatime / Real_t(3.0)) ;
   }

   if (targetdt < deltatime) {
      deltatime = targetdt ;
   }

   CurrentTime += deltatime ;

   ++cycle ;
}


//  InitStressTermsForElems :   Moved content of procedure inline at the only call site.


static
void CalcElemShapeFunctionDerivatives( const  doubleArray  &ElemCoords,
                                              doubleArray  &b,
                                              Real_t       &volume )
{
  doubleArray diag71, diag64, diag82, diag53 ;
  doubleArray   fjxi,   fjet,   fjze         ;
  doubleArray   cjxi,   cjet,   cjze         ;
  doubleArray     p1,     p2,     p3,     p4 ;

  // compute diagonal differences

  diag71 = ElemCoords(7,all) - ElemCoords(1,all);
  diag64 = ElemCoords(6,all) - ElemCoords(4,all);
  diag82 = ElemCoords(8,all) - ElemCoords(2,all);
  diag53 = ElemCoords(5,all) - ElemCoords(3,all);

  // compute jacobians

  fjxi = .125 * (diag71 + diag64 - diag82 - diag53);
  fjet = .125 * (diag71 - diag64 + diag82 - diag53);
  fjze = .125 * (diag71 + diag64 + diag82 + diag53);

  fjxi.reshape(Coord_index);              //  These reshapes required by A++ because the
  fjet.reshape(Coord_index);              //   Diag and "f" arrays are really 2D (1:1 x 1:3)
  fjze.reshape(Coord_index);              //   and the next operations need   1D (1:3)

  // compute cofactors

  cjxi = cross( fjet, fjze);
  cjet = cross( fjze, fjxi);
  cjze = cross( fjxi, fjet);

  // calcualte partials
  //     this need only be done for l = 1,2,3,4   since , by symmetry ,
  //     (7,8,5,6) = - (1,2,3,4) .

  p1 =  -  cjxi  -  cjet  -  cjze;
  p2 =     cjxi  -  cjet  -  cjze;
  p3 =     cjxi  +  cjet  -  cjze;
  p4 =  -  cjxi  +  cjet  -  cjze;

  b(1,all) =                 p1.reshape(Dummy_index,Coord_index);  //  These reshapes required by A++
  b(2,all) =                 p2.reshape(Dummy_index,Coord_index);  //  to create a 2D array to match the
  b(3,all) =                 p3.reshape(Dummy_index,Coord_index);  //  target array shape.
  b(4,all) =                 p4.reshape(Dummy_index,Coord_index);
  b(5,all) =  Real_t(-1.0) * p3.reshape(Dummy_index,Coord_index);
  b(6,all) =  Real_t(-1.0) * p4.reshape(Dummy_index,Coord_index);
  b(7,all) =  Real_t(-1.0) * p1.reshape(Dummy_index,Coord_index);
  b(8,all) =  Real_t(-1.0) * p2.reshape(Dummy_index,Coord_index);

  // calculate jacobian determinant (volume)

  volume = Real_t(8.0) * sum(fjet*cjet);

}


static inline
doubleArray SumElemFaceNormal(  const doubleArray &FaceCoords)
{
   doubleArray bisect1;
   doubleArray bisect2;

   bisect1 = Real_t(0.5) * (   FaceCoords(4,all) + FaceCoords(3,all)
                             - FaceCoords(2,all) - FaceCoords(1,all)  );
   bisect2 = Real_t(0.5) * (   FaceCoords(3,all) + FaceCoords(2,all)
                             - FaceCoords(4,all) - FaceCoords(1,all)  );

   bisect1.reshape(Coord_index);
   bisect2.reshape(Coord_index);

   return  (Real_t(0.25) * cross(bisect1,bisect2));
}


static inline
void CalcElemNodeNormals(const  doubleArray &ElemCoords,
                                doubleArray &elemForce         )
{
   doubleArray N1 ;
   doubleArray N2 ;
   doubleArray N3 ;
   doubleArray N4 ;
   doubleArray N5 ;
   doubleArray N6 ;

   N1 = SumElemFaceNormal( GetFace(ElemCoords, 1) );  // Face 1:  Corners 1,2,3,4
   N2 = SumElemFaceNormal( GetFace(ElemCoords, 5) );  // Face 5:  Corners 1,5,6,2
   N3 = SumElemFaceNormal( GetFace(ElemCoords, 4) );  // Face 4:  Corners 2,6,7,3
   N4 = SumElemFaceNormal( GetFace(ElemCoords, 6) );  // Face 6:  Corners 3,7,8,4
   N5 = SumElemFaceNormal( GetFace(ElemCoords, 3) );  // Face 3:  Corners 4,8,5,1
   N6 = SumElemFaceNormal( GetFace(ElemCoords, 2) );  // Face 2:  Corners 5,8,7,6

   N1.reshape(Dummy_index, Coord_index);
   N2.reshape(Dummy_index, Coord_index);
   N3.reshape(Dummy_index, Coord_index);
   N4.reshape(Dummy_index, Coord_index);
   N5.reshape(Dummy_index, Coord_index);
   N6.reshape(Dummy_index, Coord_index);

   elemForce(1,all) = N1 + N2 + N5;
   elemForce(2,all) = N1 + N2 + N3;
   elemForce(3,all) = N1 + N3 + N4;
   elemForce(4,all) = N1 + N4 + N5;
   elemForce(5,all) = N2 + N5 + N6;
   elemForce(6,all) = N2 + N3 + N6;
   elemForce(7,all) = N3 + N4 + N6;
   elemForce(8,all) = N4 + N5 + N6;
}


static inline
void SumElemStressesToNodeForces( const  doubleArray &B,
                                  const  doubleArray &stress,
                                         doubleArray &elemForce      )
{
  const int xx = 1;
  const int yy = 2;
  const int zz = 3;

  for (Index_t Corner = 1; Corner <= 8; ++ Corner) {
      elemForce(Corner,xx) = -( stress(xx) * B( Corner,xx ) );
      elemForce(Corner,yy) = -( stress(yy) * B( Corner,yy ) );
      elemForce(Corner,zz) = -( stress(zz) * B( Corner,zz ) );
      };
}


static
void IntegrateStressForElems( const intArray     &ElemToNode,
                              const doubleArray  &coords,
                              const doubleArray  &Stress,
                                    doubleArray  &volume,
                                    doubleArray  &force     )
{
  doubleArray          B (Corner_index, Coord_index);   // shape function derivatives
  doubleArray ElemCoords (Corner_index, Coord_index);   // X, Y, Z, values for element corners
  doubleArray  elemForce (Corner_index, Coord_index);   // X, Y, Z, forces for element corners


  // loop over all elements -- this could be a parallel loop, except for update to global force array

  for( Index_t ElemId=1 ; ElemId<= ElemToNode.getLength(0) ; ++ElemId) {

    CollectElemCoords( ElemId, ElemToNode, coords, ElemCoords);

    CalcElemShapeFunctionDerivatives( ElemCoords, B, volume(ElemId));

    CalcElemNodeNormals( ElemCoords, B );

    doubleArray ElemStress;
    ElemStress = Stress(ElemId,all);
    ElemStress.reshape(Coord_index);

    SumElemStressesToNodeForces( B, ElemStress, elemForce );

    Accum_Forces_To_Nodes(ElemId, elemForce, force);
  }
}


static inline
doubleArray VoluDer( const doubleArray &N1,
                     const doubleArray &N2,
                     const doubleArray &N3,
                     const doubleArray &N4,
                     const doubleArray &N5,
                     const doubleArray &N7  )
{
   const Real_t twelfth = Real_t(1.0) / Real_t(12.0) ;
   doubleArray dvd (Coord_index);
   const int x = 1;
   const int y = 2;
   const int z = 3;

   dvd(x) =   (N2(y) + N3(y)) * (N1(z) + N2(z)) - (N1(y) + N2(y)) * (N2(z) + N3(z)) +
              (N1(y) + N5(y)) * (N4(z) + N5(z)) - (N4(y) + N5(y)) * (N1(z) + N5(z)) -
              (N3(y) + N7(y)) * (N4(z) + N7(z)) + (N4(y) + N7(y)) * (N3(z) + N7(z));

   dvd(y) = - (N2(x) + N3(x)) * (N1(z) + N2(z)) + (N1(x) + N2(x)) * (N2(z) + N3(z)) -
              (N1(x) + N5(x)) * (N4(z) + N5(z)) + (N4(x) + N5(x)) * (N1(z) + N5(z)) +
              (N3(x) + N7(x)) * (N4(z) + N7(z)) - (N4(x) + N7(x)) * (N3(z) + N7(z));

   dvd(z) = - (N2(y) + N3(y)) * (N1(x) + N2(x)) + (N1(y) + N2(y)) * (N2(x) + N3(x)) -
              (N1(y) + N5(y)) * (N4(x) + N5(x)) + (N4(y) + N5(y)) * (N1(x) + N5(x)) +
              (N3(y) + N7(y)) * (N4(x) + N7(x)) - (N4(y) + N7(y)) * (N3(x) + N7(x));

   dvd *= twelfth;

   dvd.reshape(Dummy_index,Coord_index);                // A++ fix to get correct array shape for
                                                        // CalcElemVolumeDerivative
   return dvd;
}


static
void CalcElemVolumeDerivative(const doubleArray  &Loc,
                                    doubleArray  &Elem_dvd    )
{
      doubleArray Node1, Node2, Node3, Node4, Node5, Node6, Node7, Node8;

      Node1 =Loc(1,all);   Node1.reshape(Coord_index);
      Node2 =Loc(2,all);   Node2.reshape(Coord_index);
      Node3 =Loc(3,all);   Node3.reshape(Coord_index);
      Node4 =Loc(4,all);   Node4.reshape(Coord_index);
      Node5 =Loc(5,all);   Node5.reshape(Coord_index);
      Node6 =Loc(6,all);   Node6.reshape(Coord_index);
      Node7 =Loc(7,all);   Node7.reshape(Coord_index);
      Node8 =Loc(8,all);   Node8.reshape(Coord_index);

      Elem_dvd(1, all) =  VoluDer(Node2, Node3,  Node4, Node5, Node6, Node8);
      Elem_dvd(2, all) =  VoluDer(Node3, Node4,  Node1, Node6, Node7, Node5);
      Elem_dvd(3, all) =  VoluDer(Node4, Node1,  Node2, Node7, Node8, Node6);
      Elem_dvd(4, all) =  VoluDer(Node1, Node2,  Node3, Node8, Node5, Node7);
      Elem_dvd(5, all) =  VoluDer(Node8, Node7,  Node6, Node1, Node4, Node2);
      Elem_dvd(6, all) =  VoluDer(Node5, Node8,  Node7, Node2, Node1, Node3);
      Elem_dvd(7, all) =  VoluDer(Node6, Node5,  Node8, Node3, Node2, Node4);
      Elem_dvd(8, all) =  VoluDer(Node7, Node6,  Node5, Node4, Node3, Node1);
}


static inline
void CalcElemFBHourglassForce( const Real_t       &coeff,
                               const doubleArray  &elem_vel,
                               const doubleArray  &HourGam,
                                     doubleArray  &elem_hgf     )
{
   doubleArray    Xpose1 = transpose(elem_vel);           Xpose1.setBase(1);
   doubleArray  MM_HG_EV = MatMult(Xpose1, HourGam );
   doubleArray    Xpose2 = transpose(MM_HG_EV);           Xpose2.setBase(1);
                elem_hgf = coeff * MatMult( HourGam, Xpose2 );
}


static
void FBHour(const    intArray  &ElemToNode,
            const doubleArray  &coords,
            const doubleArray  &vel,
            const doubleArray  &ss,
            const doubleArray  &z_mass,
            const doubleArray  &volume,
                  doubleArray  &force      )
{
//  Calculate the Flanagan-Belytschko anti-hourglass force.
//  Gamma is a globally-defined constant array of +1 and -1 coefficients

        doubleArray ElemCoords (Corner_index,  Coord_index);
        doubleArray   Elem_dvd (Corner_index,  Coord_index);
        doubleArray    HourMod ( Coord_index, Strain_Modes);
        doubleArray    HourGam (Corner_index, Strain_Modes);
        doubleArray   elem_vel (Corner_index,  Coord_index);
        doubleArray   elem_hgf (Corner_index,  Coord_index);
        doubleArray  Xpose;

   for( Index_t ElemId = 1; ElemId<=ElemToNode.getLength(0) ; ++ ElemId)  {

        CollectElemCoords( ElemId, ElemToNode, coords, ElemCoords);
        CalcElemVolumeDerivative( ElemCoords, Elem_dvd );

        Real_t volinv = Real_t(1.0)/volume(ElemId);

        Xpose = transpose(ElemCoords);
        Xpose.setBase(1);
        HourMod = MatMult( Xpose, Gamma );                       //  Coord_index   x  Strain_modes
        HourGam = Gamma - volinv * MatMult(Elem_dvd, HourMod );  //  Corner_Index  x  Strain_modes

        Real_t ss1       = ss(ElemId);
        Real_t mass1     = Z_mass(ElemId);
        Real_t volume13  = cbrt(volume(ElemId));
        Real_t coeff     = - mc_hgcoef * Real_t(0.01) * ss1 * mass1 / volume13;

        CollectElemVelocities( ElemId, ElemToNode, vel, elem_vel);
        CalcElemFBHourglassForce( coeff, elem_vel, HourGam, elem_hgf );
        Accum_Forces_To_Nodes( ElemId, elem_hgf, force);

    }
}



//  CalcHourglassControlForElems:  This procedure has been included in the next one below.
//  CalcVolumeForceForElems:       This procedure has been included in the next one below.


static
void CalcForceForNodes(const    intArray  &ElemToNode,
                       const doubleArray  &coords,
                       const doubleArray  &vel,
                       const doubleArray  &ss,
                       const doubleArray  &p,
                       const doubleArray  &q,
                       const doubleArray  &z_mass,
                       const doubleArray  &ref_vol,
                       const doubleArray  &rel_vol,
                             doubleArray  &force       )
{
   force = Real_t(0.0);
   if (ElemToNode.getLength(0) != 0) {

      doubleArray  stressvector;
      doubleArray  stress (Elem_index, Coord_index);
      doubleArray  volume (Elem_index);

      // pull in the stresses appropriate to the hydro integration
      // previously a call to InitStressTermsForElems

      stressvector  =  - p - q ;
      stress(all,1) =  stressvector.reshape(Elem_index, Dummy_index) ;
      stress(all,2) =  stressvector ;
      stress(all,3) =  stressvector ;

      // compute nodal volume and forces from material stresses.

      IntegrateStressForElems(ElemToNode, coords, stress, volume, force ) ;

      assert(min(volume) > Real_t(0.0));
      if ( min(volume) <= Real_t(0.0)) { exit(VolumeError) ; }    // check for negative elem volume

      // Compute and add hourglass forces to force results from material stress.
      // Following three lines replace the original call to LagrangeCalcHourglassControl

      volume = ref_vol * rel_vol;

      assert(min(rel_vol) > Real_t(0.0));

      if ( min(rel_vol) <= Real_t(0.0) ) {exit(VolumeError) ; }     // check for negative rel volume

      assert(mc_hgcoef > Real_t(0.0));

      if ( mc_hgcoef > Real_t(0.0) ) {

         FBHour( ElemToNode, coords, vel, ss, z_mass, volume, force ) ;
      }
   }
  //   Calculate Nodal Forces at domain boundaries
  //   problem->commSBN->Transfer(CommSBN::forces);
}


static inline
void CalcAccelerationForNodes( const doubleArray  &force,
                               const doubleArray  &n_mass,
                                     doubleArray  &accel  )
{
   for (int i = 1; i <= 3; ++i) {
        accel(all,i) = force(all,i) / n_mass;     //  loop over X, Y, Z coordinates
   }
}


static inline
void ApplyAccelerationBoundaryConditionsForNodes( const intArray     &BCnodes,
                                                        doubleArray  &accel   )
{
  Index_t numNodeBC = BCnodes.getLength(0); ;
  for(Index_t i=1 ; i<=numNodeBC ; ++i) {
     accel(BCnodes(i,1),1) = Real_t(0.0) ;
  };
  for(Index_t i=1 ; i<=numNodeBC ; ++i)  {
     accel(BCnodes(i,2),2) = Real_t(0.0) ;
  };
  for(Index_t i=1 ; i<=numNodeBC ; ++i)  {
     accel(BCnodes(i,3),3) = Real_t(0.0) ;
  };
}


static inline
void CalcVelocityForNodes(const      Real_t  &dt,
                          const doubleArray  &accel,
                                doubleArray  &vel  )
{
   vel += accel * dt;
   where  (fabs(vel) < mc_u_cut)
          { vel = Real_t(0.0) ; };
 }


static inline
void  CalcPositionForNodes(const     Real_t  &dt,
                          const doubleArray  &vel,
                                doubleArray  &coords)
{
   coords += vel * dt;
}


static
void LagrangeNodalLeapFrog( const    intArray  &ElemToNode,
                            const doubleArray  &ss,
                            const doubleArray  &p,
                            const doubleArray  &q,
                            const doubleArray  &Z_mass,
                            const doubleArray  &ref_vol,
                            const doubleArray  &rel_vol,
                            const doubleArray  &n_mass,
                            const intArray     &BCnodes,
                            const Real_t       &deltatime,
                                  doubleArray  &force,
                                  doubleArray  &accel,
                                  doubleArray  &vel,
                                  doubleArray  &coords)
{
  //  time of boundary condition evaluation is beginning of step for force and
  //  acceleration boundary conditions.

  CalcForceForNodes(ElemToNode, coords, vel, ss, p, q, Z_mass, ref_vol, rel_vol, force);
  CalcAccelerationForNodes(force, n_mass, accel);
  ApplyAccelerationBoundaryConditionsForNodes(BCnodes, accel);
  CalcVelocityForNodes(deltatime, accel, vel ) ;
  CalcPositionForNodes(deltatime, vel, coords );
  return;
}


 //  The following function replaces 2 Volume fctns in org. code

static
Real_t  CalcElemVolume( const doubleArray &ElemCoords )
{
  Real_t twelveth = Real_t(1.0)/Real_t(12.0);

  // compute diagonal differences
  doubleArray diag72 ;
  doubleArray diag81 ;
  doubleArray diag74 ;
  doubleArray diag31 ;
  doubleArray diag61 ;
  doubleArray diag75 ;
  doubleArray diag42 ;
  doubleArray diag83 ;
  doubleArray diag54 ;
  doubleArray diag68 ;
  doubleArray diag25 ;
  doubleArray diag36 ;

  diag72  = ElemCoords(7,all) - ElemCoords(2,all);
  diag81  = ElemCoords(8,all) - ElemCoords(1,all);
  diag74  = ElemCoords(7,all) - ElemCoords(4,all);
  diag31  = ElemCoords(3,all) - ElemCoords(1,all);
  diag61  = ElemCoords(6,all) - ElemCoords(1,all);
  diag75  = ElemCoords(7,all) - ElemCoords(5,all);
  diag42  = ElemCoords(4,all) - ElemCoords(2,all);
  diag83  = ElemCoords(8,all) - ElemCoords(3,all);
  diag54  = ElemCoords(5,all) - ElemCoords(4,all);
  diag68  = ElemCoords(6,all) - ElemCoords(8,all);
  diag25  = ElemCoords(2,all) - ElemCoords(5,all);
  diag36  = ElemCoords(3,all) - ElemCoords(6,all);

  diag72.reshape(Coord_index);
  diag81.reshape(Coord_index);
  diag74.reshape(Coord_index);
  diag31.reshape(Coord_index);
  diag61.reshape(Coord_index);
  diag75.reshape(Coord_index);
  diag42.reshape(Coord_index);
  diag83.reshape(Coord_index);
  diag54.reshape(Coord_index);
  diag68.reshape(Coord_index);
  diag25.reshape(Coord_index);
  diag36.reshape(Coord_index);

  Real_t volume =  sum( (diag42 + diag83) * cross(diag74, diag31)) +
                   sum( (diag54 + diag68) * cross(diag75, diag81)) +
                   sum( (diag25 + diag36) * cross(diag72, diag61)) ;

  volume *= twelveth;

  return volume ;
}


static inline
Real_t AreaFace( const doubleArray  &FaceCoords)
{
   doubleArray diff31 = (FaceCoords(3,all) - FaceCoords(1,all));
   doubleArray diff42 = (FaceCoords(4,all) - FaceCoords(2,all));
   doubleArray      f = diff31 - diff42;
   doubleArray      g = diff31 + diff42;

   Real_t        area = sum(f*f) * sum(g*g) - sum(f*g) * sum(f*g);
   return area ;
}


static inline
Real_t CharacteristicLength( const doubleArray  &ElemCoords,
                             const Real_t       &volume     )
{
   doubleArray area (Face_index);

   for (Index_t FaceId = 1; FaceId <= 6 ; ++FaceId )  {
       area(FaceId) = AreaFace(GetFace(ElemCoords, FaceId)) ;
   }

   Real_t charLength = Real_t(4.0) * volume / sqrt( max(area) );

   return charLength;
}


static
void  CalcElemVelocityGrandient( const doubleArray  &ElemVel,
                                 const doubleArray  &b,
                                 const Real_t       &volume,
                                       doubleArray  &d       )
{
   Index term (1,4);
   const doubleArray VelDif (term, Coord_index);

   VelDif(1,all) = ElemVel(1,all) - ElemVel(7,all);
   VelDif(2,all) = ElemVel(2,all) - ElemVel(8,all);
   VelDif(3,all) = ElemVel(3,all) - ElemVel(5,all);
   VelDif(4,all) = ElemVel(4,all) - ElemVel(6,all);

   doubleArray pf (term, Coord_index);
   doubleArray MMresult (Coord_index, Coord_index);
   Real_t inv_vol = Real_t(1.0) / volume ;

   pf       = b(term, all);
   doubleArray Xpose = transpose(pf);  Xpose.setBase(1);

   MMresult = inv_vol * MatMult( Xpose, VelDif );

   Index_t x = 1;
   Index_t y = 2;
   Index_t z = 3;

   d(1) =  MMresult(x,x) ;
   d(2) =  MMresult(y,y) ;
   d(3) =  MMresult(z,z) ;
   d(4) =  Real_t( 0.5) * ( MMresult(z,y) + MMresult(y,z) );       //  .5 * ( dzddy + dyddz )
   d(5) =  Real_t( 0.5) * ( MMresult(x,z) + MMresult(z,x) );       //  .5 * ( dxddz + dzddx )
   d(6) =  Real_t( 0.5) * ( MMresult(x,y) + MMresult(y,x) );       //  .5 * ( dxddy + dyddx )
}


static
void CalcKinematicsForElems( const intArray     &ElemToNode,
                             const doubleArray  &coords,
                             const doubleArray  &vel,
                             const doubleArray  &ref_vol,
                             const doubleArray  &rel_vol,
                             const      Real_t  &dt,
                                   doubleArray  &new_vol,
                                   doubleArray  &delta_v,
                                   doubleArray  &char_len,
                                   doubleArray  &Prin_Strains   )
{
  Index D_index (1,6);
  Index D_1To3  (1,3);
  doubleArray ElemCoords ( Corner_index, Coord_index );
  doubleArray   ElemVels ( Corner_index, Coord_index );
  doubleArray          B ( Corner_index, Coord_index );     // Shape function derivatives
  doubleArray          D ( D_index );

  for( Index_t ElemId = 1 ; ElemId<=ElemToNode.getLength(0) ; ++ ElemId )
  {
    Real_t volume ;
    Real_t relativeVolume ;
    CollectElemCoords( ElemId, ElemToNode, coords, ElemCoords);

                 volume  =  CalcElemVolume(ElemCoords);
         relativeVolume  =  volume / ref_vol(ElemId) ;
        new_vol(ElemId)  =  relativeVolume ;
        delta_v(ElemId)  =  relativeVolume - rel_vol(ElemId); ;
       char_len(ElemId)  =  CharacteristicLength( ElemCoords, volume);

    CollectElemVelocities( ElemId, ElemToNode, vel, ElemVels);
    ElemCoords  -= ElemVels * (dt * Real_t(0.5)) ;
    CalcElemShapeFunctionDerivatives( ElemCoords, B, volume );
    CalcElemVelocityGrandient( ElemVels, B, volume, D );

    // put velocity gradient quantities into global arrays.

    doubleArray Dtemp = D(D_1To3);
    Prin_Strains( ElemId, Coord_index) = Dtemp.reshape(Dummy_index,Coord_index);
  }
}


static
void CalcLagrangeElements( const intArray     &ElemToNode,
                           const doubleArray  &coords,
                           const doubleArray  &vel,
                           const doubleArray  &ref_vol,
                           const doubleArray  &rel_vol,
                           const      Real_t  &dt,
                                 doubleArray  &new_vol,
                                 doubleArray  &delta_v,
                                 doubleArray  &vdov,
                                 doubleArray  &char_len    )
{
  if (ElemToNode.getLength(0) > 0) {

      doubleArray Prin_Strains (Elem_index,  Coord_index);                 //  old dxx, dyy, dzz


      CalcKinematicsForElems(ElemToNode, coords, vel, ref_vol, rel_vol, dt,
                             new_vol, delta_v, char_len, Prin_Strains) ;    // results returned

      // element loop to do some stuff not included in the elemlib function.

      for( Index_t ElemId = 1 ; ElemId<=ElemToNode.getLength(0) ; ++ ElemId )
      {
        // calc strain rate and apply as constraint (only done in FB element)
        Real_t vdov_local = sum( Prin_Strains (ElemId, all) ) ;
        Real_t vdovthird  = vdov_local/Real_t(3.0) ;

        // make the rate of deformation tensor deviatoric
        vdov(ElemId) = vdov_local ;
        Prin_Strains (ElemId, all) -= vdovthird ;

        // See if any volumes are negative, and take appropriate action.
        if (new_vol(ElemId) <= Real_t(0.0)) { exit(VolumeError) ; }
      }
   }
}


static
void MonotonicQCalcVelocityGradient( const intArray     &ElemToNode,
                                     const doubleArray  &coords,
                                     const doubleArray  &ref_vol,
                                     const doubleArray  &new_vol,
                                     const doubleArray  &vel,
                                           doubleArray  &delta_v_grad,
                                           doubleArray  &delta_x)
{
   doubleArray  ElemC ( Corner_index, Coord_index );
   doubleArray  ElemV ( Corner_index, Coord_index );
   doubleArray  FaceDifsEta  (Coord_index);
   doubleArray  FaceDifsXi   (Coord_index);
   doubleArray  FaceDifsZeta (Coord_index);
   doubleArray  VelDifsEta   (Coord_index);
   doubleArray  VelDifsXi    (Coord_index);
   doubleArray  VelDifsZeta  (Coord_index);
   doubleArray             a (Coord_index);

   const Real_t ptiny = Real_t(1.e-36) ;

   Index_t  eta = 1;
   Index_t   xi = 2;
   Index_t zeta = 3;

   for( Index_t ElemId = 1 ; ElemId<=ElemToNode.getLength(0) ; ++ ElemId )
     {
      Real_t vol  = ref_vol(ElemId) * new_vol(ElemId) ;
      Real_t norm = Real_t(1.0) / ( vol + ptiny ) ;

      CollectElemCoords(     ElemId, ElemToNode, coords, ElemC );
      CollectElemVelocities( ElemId, ElemToNode, vel,    ElemV );

      //  Calc FaceDifs and VelDifs for the Element
      //  Original code had neg. coefficients for "eta".
      //  This version reverses the order of the faces instead.

      FaceDifsEta  =    AveFaceLocs ( GetFace( ElemC, 6))       // Face 6:  Corners 3,7,8,4
                      - AveFaceLocs ( GetFace( ElemC, 5)) ;     // Face 5:  Corners 1,5,6,2
      FaceDifsXi   =    AveFaceLocs ( GetFace( ElemC, 4))       // Face 4:  Corners 2,6,7,3
                      - AveFaceLocs ( GetFace( ElemC, 3)) ;     // Face 3:  Corners 4,8,5,1
      FaceDifsZeta =    AveFaceLocs ( GetFace( ElemC, 2))       // Face 2:  Corners 5,8,7,6
                      - AveFaceLocs ( GetFace( ElemC, 1)) ;     // Face 1:  Corners 1,2,3,4

       VelDifsEta  =    AveFaceLocs ( GetFace( ElemV, 6))
                      - AveFaceLocs ( GetFace( ElemV, 5)) ;
       VelDifsXi   =    AveFaceLocs ( GetFace( ElemV, 4))
                      - AveFaceLocs ( GetFace( ElemV, 3)) ;
       VelDifsZeta =    AveFaceLocs ( GetFace( ElemV, 2))
                      - AveFaceLocs ( GetFace( ElemV, 1)) ;

      a = cross(FaceDifsXi, FaceDifsEta );
      delta_x( ElemId, zeta)      = vol / sqrt( sum( a * a ) + ptiny) ;
      delta_v_grad( ElemId, zeta) = norm * sum( a * VelDifsZeta ) ;

      a = cross(FaceDifsEta, FaceDifsZeta );
      delta_x( ElemId,   xi)      = vol / sqrt( sum( a * a ) + ptiny) ;
      delta_v_grad( ElemId,   xi) = norm * sum( a * VelDifsXi ) ;

      a = cross(FaceDifsZeta, FaceDifsXi );
      delta_x( ElemId,  eta)      = vol / sqrt( sum( a * a ) + ptiny) ;
      delta_v_grad( ElemId,  eta) = norm * sum( a * VelDifsEta ) ;
   }
}


static
Real_t CalcPhi (  const Index_t      &ElemId,            // New routine used by CalcMonotonicQForElems
                  const intArray     &ElemFaceToElem,
                  const doubleArray  &delta_v_grad,
                  const Int_t        &dir,
                  const Int_t        &bc_M,
                  const Int_t        &bc_P )
{
    Real_t phi;
    Real_t delvm;
    Real_t delvp;
    Real_t ptiny   = Real_t(1.e-36) ;
    Real_t norm    = Real_t(1.0) / ( delta_v_grad(ElemId,dir) + ptiny ) ;

    switch (bc_M) {                                                         // minus calc
       case 0:   delvm = delta_v_grad(ElemFaceToElem(ElemId,dir,1),dir) ;  break ;  // internal Elem
       case 1:   delvm = delta_v_grad(ElemId,dir) ;                        break ;  // Symm boundary
       case 2:   delvm = Real_t(0.0) ;                                     break ;  // Free boundary
       default:          exit(PhiCalcError) ;                              break ;  // error
       }
    switch (bc_P) {                                                         // plus calc
       case 0:   delvp = delta_v_grad(ElemFaceToElem(ElemId,dir,2),dir) ;  break ;  // internal Elem
       case 4:   delvp = delta_v_grad(ElemId,dir) ;                        break ;  // Symm boundary
       case 8:   delvp = Real_t(0.0) ;                                     break ;  // Free boundary
       default:          exit(PhiCalcError) ;                              break ;  // error
       }

    delvm = delvm * norm ;
    delvp = delvp * norm ;

    phi = Real_t(0.5) * ( delvm + delvp ) ;

    delvm *= mc_monoq_limiter_mult ;
    delvp *= mc_monoq_limiter_mult ;

    if ( delvm < phi )             phi = delvm ;
    if ( delvp < phi )             phi = delvp ;
    if ( phi < Real_t(0.0))        phi = Real_t(0.0) ;
    if ( phi > mc_monoq_max_slope) phi = mc_monoq_max_slope;

    return phi;
}


static
void CalcMonotonicQForElems(  const intArray     &ElemFaceToElem,
                              const doubleArray  &delta_x,
                              const doubleArray  &delta_v_grad,
                              const doubleArray  &Z_mass,
                              const doubleArray  &vdov,
                              const intArray     &matElemlist,
                              const intArray     &elemBC,
                                    doubleArray  &qq,
                                    doubleArray  &ql                 )
{
   //  CalcMonotonicQForElems and CalcMonotonicQRegionForElems were combined into one routine
   //  because the former had essentially no significant content.

   Real_t              qlin, qquad, rho ;
   doubleArray         phi (Dir_index)  ;

   Int_t  eta = 1;     Int_t    bcXi_M;     Int_t    bcXi_P;
   Int_t   xi = 2;     Int_t   bcEta_M;     Int_t   bcEta_P;
   Int_t zeta = 3;     Int_t  bcZeta_M;     Int_t  bcZeta_P;

   for ( Index_t izone = 1 ; izone <=matElemlist.getLength(0) ; ++izone ) {

      Index_t  ElemId = matElemlist(izone);
      Int_t    bcMask = elemBC(ElemId) ;

        bcXi_M =  (bcMask &   XI_M)     ;       bcXi_P =  (bcMask &   XI_P)     ;
       bcEta_M =  (bcMask &  ETA_M) >> 4;      bcEta_P =  (bcMask &  ETA_P) >> 4;
      bcZeta_M =  (bcMask & ZETA_M) >> 8;     bcZeta_P =  (bcMask & ZETA_P) >> 8;

      phi(  xi) = CalcPhi ( ElemId, ElemFaceToElem, delta_v_grad,   xi,   bcXi_M,   bcXi_P);
      phi( eta) = CalcPhi ( ElemId, ElemFaceToElem, delta_v_grad,  eta,  bcEta_M,  bcEta_P);
      phi(zeta) = CalcPhi ( ElemId, ElemFaceToElem, delta_v_grad, zeta, bcZeta_M, bcZeta_P);

      // Remove length scale

      if ( vdov(ElemId) > Real_t(0.0) )  {

         qlin  = Real_t(0.0) ;
         qquad = Real_t(0.0) ;
      }
      else {

         doubleArray loc_delta_v   ;
         loc_delta_v = delta_v_grad(ElemId, all) * delta_x(ElemId, all);
         loc_delta_v.reshape(Dir_index);

         where (loc_delta_v > Real_t(0.0))  { loc_delta_v = Real_t(0.0); };

         rho    = Z_mass(ElemId) / (ref_vol(ElemId) * new_vol(ElemId)) ;
         qlin   = -mc_qlc_monoq * rho * sum( loc_delta_v * (Real_t(1.0) - phi));
         qquad  =  mc_qqc_monoq * rho * sum( loc_delta_v * loc_delta_v * (Real_t(1.0) - phi * phi));
      }

      qq(ElemId) = qquad ;
      ql(ElemId) = qlin  ;
   }
}


static inline
void CalcQForElems( const intArray     &ElemFaceToElem,
                    const intArray     &matElemlist,
                    const intArray     &elemBC,
                    const intArray     &ElemToNode,
                    const doubleArray  &coords,
                    const doubleArray  &ref_vol,
                    const doubleArray  &new_vol,
                    const doubleArray  &vel,
                    const doubleArray  &q,
                    const doubleArray  &Z_mass,
                    const doubleArray  &vdov,
                          doubleArray  &delta_v_grad,
                          doubleArray  &delta_x,
                          doubleArray  &qq,
                          doubleArray  &ql                 )
{

   MonotonicQCalcVelocityGradient( ElemToNode, coords, ref_vol, new_vol, vel,
                                   delta_v_grad, delta_x) ;                      // results returned

   // Transfer veloctiy gradients in the first order elements
   // problem->commElements->Transfer(CommElements::monoQ) ;

   CalcMonotonicQForElems(ElemFaceToElem, delta_x, delta_v_grad, Z_mass, vdov, matElemlist, elemBC,
                          qq, ql) ;                                        // results returned

   // Don't allow excessive artificial viscosity

   if (max(q) > mc_qstop) {exit(QStopError) ;}
}


static
void CalcElemPressure( const Real_t  &e_old_elem,            // Changed from original version
                       const Real_t  &compression_elem,      // This calculates one element's
                       const Real_t  &new_vol_elem,          // pressure.  Calling code must loop
                             Real_t  &p_new_elem,            // over all elements.
                             Real_t  &bvc_elem,
                             Real_t  &pbvc_elem        )
{
   Real_t      c1s = Real_t(2.0)/Real_t(3.0) ;
          bvc_elem = c1s * (compression_elem + Real_t(1.0));
         pbvc_elem = c1s;
        p_new_elem = bvc_elem * e_old_elem;

   if (  fabs(p_new_elem) <  mc_p_cut  )    { p_new_elem = Real_t(0.0) ; };
   if (      new_vol_elem >= mc_eosvmax)    { p_new_elem = Real_t(0.0) ; };
   if (        p_new_elem <  mc_pmin   )    { p_new_elem = mc_pmin      ; };
}


static inline
Real_t CalcSoundSpeedForElem( const Real_t  &new_vol_elem,
                              const Real_t  &enew_elem,
                              const Real_t  &pnew_elem,
                              const Real_t  &pbvc_elem,
                              const Real_t  &bvc_elem)
{
   Real_t  ss;
   ss = (   pbvc_elem * enew_elem
          + new_vol_elem * new_vol_elem * bvc_elem * pnew_elem) / mc_refdens ;
   ss =  ( ss <= Real_t(1.111111e-36) ) ? Real_t(.333333e-18) :  sqrt(ss)  ;
   return ss;
}


static
void CalcElemEnergy( const Real_t  &p_old_elem,                 // Changed from original version
                     const Real_t  &e_old_elem,                 // This calculates one element's
                     const Real_t  &q_old_elem,                 // energy.  Calling code must loop
                     const Real_t  &compression_elem,           // over all elements.
                     const Real_t  &compHalfStep_elem,
                     const Real_t  &new_vol_elem,
                     const Real_t  &work_elem,
                     const Real_t  &delta_v_elem,
                     const Real_t  &qq_elem,
                     const Real_t  &ql_elem,
                           Real_t  &p_new_elem,
                           Real_t  &e_new_elem,
                           Real_t  &q_new_elem,
                           Real_t  &bvc_elem,
                           Real_t  &pbvc_elem           )
{
   Real_t      sixth = Real_t(1.0) / Real_t(6.0) ;
   Real_t pHalfStep ;
   Real_t     vhalf ;
   Real_t       ssc ;
   Real_t   q_tilde ;

   e_new_elem = max( mc_emin,
                     e_old_elem
                     - Real_t(0.5) * delta_v_elem * (p_old_elem + q_old_elem)
                     + Real_t(0.5) * work_elem);

   CalcElemPressure( e_new_elem, compHalfStep_elem, new_vol_elem,
                     pHalfStep, bvc_elem, pbvc_elem);                //  result parameters

   vhalf = Real_t(1.0) / (Real_t(1.0) + compHalfStep_elem) ;

   ssc = CalcSoundSpeedForElem(vhalf, e_new_elem, pHalfStep, pbvc_elem, bvc_elem);

   q_new_elem = ( delta_v_elem > Real_t(0.0)) ? Real_t(0.0) : ssc * ql_elem + qq_elem ;

   e_new_elem = e_new_elem
                + Real_t(0.5) * delta_v_elem * (  Real_t(3.0) * (p_old_elem + q_old_elem)
                                                - Real_t(4.0) * (pHalfStep  + q_new_elem))
                + Real_t(0.5) * work_elem;

         if ( fabs(e_new_elem) < mc_e_cut )  { e_new_elem = Real_t(0.0)  ; };
         if (      e_new_elem  < mc_emin  )  { e_new_elem = mc_emin     ; };


   CalcElemPressure( e_new_elem, compression_elem, new_vol_elem,
                     p_new_elem, bvc_elem, pbvc_elem);

   ssc = CalcSoundSpeedForElem(new_vol_elem, e_new_elem, p_new_elem, pbvc_elem, bvc_elem);

   q_tilde =  (delta_v_elem > Real_t(0.0)) ? Real_t(0.0)  : ssc * ql_elem + qq_elem ;

   e_new_elem = e_new_elem - (   Real_t(7.0) * (p_old_elem + q_old_elem)
                               - Real_t(8.0) * (pHalfStep  + q_new_elem)
                               + (p_new_elem + q_tilde )  ) * delta_v_elem * sixth ;
          if ( fabs(e_new_elem) < mc_e_cut )  { e_new_elem = Real_t(0.0)  ; };
          if (      e_new_elem  < mc_emin  )  { e_new_elem = mc_emin     ; };

   CalcElemPressure( e_new_elem, compression_elem, new_vol_elem,
                     p_new_elem, bvc_elem, pbvc_elem);

   ssc = CalcSoundSpeedForElem(new_vol_elem, e_new_elem, p_new_elem, pbvc_elem, bvc_elem);

   q_new_elem = ssc * ql_elem + qq_elem ;
         if ( fabs(q_new_elem) < mc_q_cut )
                       {  q_new_elem = Real_t(0.0)  ; }
   return ;
}


// Function EvalEOSforElems merged into: ApplyMaterialPropertiesForElems


static
void ApplyMaterialPropertiesForElems(  const intArray     &matElemlist,
                                       const doubleArray  &new_vol,
                                       const doubleArray  &rel_vol,
                                       const doubleArray  &delta_v,
                                       const doubleArray  &qq,
                                       const doubleArray  &ql,
                                       const doubleArray  &vdov,
                                             doubleArray  &pressure,
                                             doubleArray  &energy,
                                             doubleArray  &q,
                                             doubleArray  &ss          )
{
   Real_t  new_volc;
   Real_t  rel_volc;
   Real_t  p_new;
   Real_t  p_old;
   Real_t  e_new;
   Real_t  e_old;
   Real_t  q_new;
   Real_t  q_old;
   Real_t  qq_old;
   Real_t  ql_old;
   Real_t  delta_v_old;
   Real_t  work;
   Real_t  compression;
   Real_t  compHalfStep;
   Real_t  bvc;
   Real_t  pbvc;
   Real_t  vc;
   Index_t zidx;

   for (Index_t i=1; i<=matElemlist.getLength(0); ++i) {    // this loop can be executed in parallel
      zidx         = matElemlist(i) ;
      new_volc     = (mc_eosvmin != Real_t(0.0))  ? max( new_vol(zidx), mc_eosvmin) : new_vol(zidx) ;
      new_volc     = (mc_eosvmax != Real_t(0.0))  ? min( new_volc,      mc_eosvmax) : new_volc      ;
      rel_volc     = (mc_eosvmin != Real_t(0.0))  ? max( rel_vol(zidx), mc_eosvmin) : rel_vol(zidx) ;
      rel_volc     = (mc_eosvmax != Real_t(0.0))  ? min( rel_volc,      mc_eosvmax) : rel_volc      ;

      if (rel_volc <= Real_t(0.0)) {exit(VolumeError) ;}

      p_old        = pressure(zidx) ;
      e_old        = energy(zidx);
      delta_v_old  = delta_v(zidx);
      q_old        = q(zidx);
      qq_old       = qq(zidx);
      ql_old       = ql(zidx);
      work         = Real_t(0.0) ;

      compression  = Real_t(1.0) /  new_volc - Real_t(1.0);
      compHalfStep = Real_t(1.0) / (new_volc - delta_v_old * Real_t(0.5) ) - Real_t(1.0);

      // Check for v > eosvmax or v < eosvmin

      if ( mc_eosvmin != Real_t(0.0) && (new_volc <= mc_eosvmin))  // impossible due to calling func?
               { compHalfStep = compression ;  }

      if ( mc_eosvmax != Real_t(0.0) && (new_volc >= mc_eosvmax))  // impossible due to calling func?
               { p_old        = Real_t(0.0) ;
                 compression  = Real_t(0.0) ;
                 compHalfStep = Real_t(0.0) ;
               }

      CalcElemEnergy( p_old, e_old,  q_old, compression, compHalfStep, new_volc,
                  work, delta_v_old, qq_old, ql_old,
                  p_new, e_new, q_new, bvc, pbvc       );                  //   result parameters

      pressure(zidx) = p_new ;
      energy(zidx)   = e_new ;
      q(zidx)        = q_new ;
      ss(zidx)       = CalcSoundSpeedForElem(new_volc, e_new, p_new, pbvc, bvc) ;

   }  // end loop over material zone list
}


static inline
void UpdateVolumesForElems( const doubleArray  &new_vol,
                                  doubleArray  &rel_vol )
{
   rel_vol  = new_vol;
   where  (fabs(rel_vol - Real_t(1.0)) < mc_v_cut )
          { rel_vol = Real_t(1.0); };
}


static
void LagrangeElements( const intArray     &ElemToNode,
                       const doubleArray  &coords,
                       const doubleArray  &vel,
                       const doubleArray  &ref_vol,
                       const doubleArray  &Z_mass,
                       const intArray     &ElemFaceToElem,
                       const intArray     &matElemlist,
                       const intArray     &elemBC,
                       const Real_t       &deltatime,
                             doubleArray  &new_vol,
                             doubleArray  &vdov,
                             doubleArray  &char_len,
                             doubleArray  &delta_v,
                             doubleArray  &qq,
                             doubleArray  &ql,
                             doubleArray  &pressure,
                             doubleArray  &energy,
                             doubleArray  &q,
                             doubleArray  &ss,
                             doubleArray  &rel_vol     )
{
  const Real_t dt = deltatime;
  doubleArray    delta_v_grad (Elem_index,    Dir_index);          //  old mesh.delv_xi ...
  doubleArray         delta_x (Elem_index,    Dir_index);          //  old mesh.delx_xi ...

  CalcLagrangeElements( ElemToNode, coords, vel, ref_vol, rel_vol, dt,
                        new_vol, delta_v, vdov, char_len) ;                  // results returned

  //   Calculate Q.  (Monotonic q option requires communication)

  CalcQForElems(  ElemFaceToElem, matElemlist, elemBC, ElemToNode, coords,
                  ref_vol, new_vol, vel, q, Z_mass, vdov,
                  delta_v_grad, delta_x, qq, ql) ;                           // results returned


  ApplyMaterialPropertiesForElems( matElemlist, new_vol, rel_vol,
                                   delta_v, qq, ql, vdov,
                                   pressure, energy, q, ss ) ;              // results returned

  UpdateVolumesForElems(  new_vol, rel_vol );
}


static
void CalcCourantConstraintForElems( const intArray     &matElemlist,
                                    const doubleArray  &char_len,
                                    const doubleArray  &ss,
                                    const doubleArray  &vdov,
                                          Real_t       &dtcourant)
{
   Real_t   dtcourant_tmp  = Real_t(1.0e+20) ;
   Real_t   dtf;
   Real_t   char_lenc;
   Index_t  courant_elem   = -1 ;
   Index_t  zidx;

   Real_t  qqc2 = Real_t(64.0) * mc_qqc * mc_qqc ;

   for (Index_t i = 1 ; i <= matElemlist.getLength(0) ; ++i) {
      zidx      = matElemlist(i) ;
      char_lenc = char_len(zidx);

      dtf = ( vdov(zidx) >= Real_t(0.0) )  ?
               ss(zidx) * ss(zidx)
            :  ss(zidx) * ss(zidx) + ( qqc2 * char_lenc * char_lenc * vdov(zidx) * vdov(zidx) );

      dtf = char_lenc / sqrt(dtf) ;

      /* determine minimum timestep with its corresponding elem */

      if (vdov(zidx) != Real_t(0.0)) {
         if ( dtf < dtcourant_tmp ) {
            dtcourant_tmp = dtf ;
            courant_elem = zidx ;
         }
      }
   }
   if (courant_elem != -1) { dtcourant = dtcourant_tmp ; }
   return ;
}


static inline
void CalcHydroConstraintForElems( const intArray     &matElemlist,
                                  const doubleArray  &vdov,
                                        Real_t       &dthydro     )
{
   Real_t dthydro_tmp = Real_t(1.0e+20) ;
   Index_t hydro_elem = -1 ;
   Int_t zidx;

   for (Index_t i = 1 ; i <= matElemlist.getLength(0) ; ++i) {
      zidx      = matElemlist(i) ;

      if (vdov(zidx) != Real_t(0.0)) {
         Real_t dtdvov = mc_dvovmax / (fabs(vdov(zidx))+Real_t(1.e-20)) ;
         if ( dthydro_tmp > dtdvov ) {
            dthydro_tmp = dtdvov ;
            hydro_elem = zidx ;
         }
      }
   }
   if (hydro_elem != -1) { dthydro = dthydro_tmp ; };

   return ;
}


static inline
void CalcTimeConstraintsForElems( const intArray     &matElemlist,
                                  const doubleArray  &char_len,
                                  const doubleArray  &ss,
                                  const doubleArray  &vdov,
                                        Real_t       &dtcourant,
                                        Real_t       &dthydro       )
{
   CalcCourantConstraintForElems( matElemlist, char_len, ss, vdov, dtcourant);
   CalcHydroConstraintForElems  ( matElemlist, vdov, dthydro);
}


static inline
void LagrangeLeapFrog(const intArray     &ElemToNode,
                      const intArray     &ElemFaceToElem,
                      const intArray     &matElemlist,
                      const intArray     &elemBC,
                      const intArray     &BCnodes,
                      const doubleArray  &Z_mass,
                      const doubleArray  &ref_vol,
                      const doubleArray  &n_mass,
                      const Real_t       &deltatime,
                            doubleArray  &force,
                            doubleArray  &accel,
                            doubleArray  &vel,
                            doubleArray  &coords,
                            doubleArray  &new_vol,
                            doubleArray  &vdov,
                            doubleArray  &char_len,
                            doubleArray  &delta_v,
                            doubleArray  &qq,
                            doubleArray  &ql,
                            doubleArray  &pressure,
                            doubleArray  &energy,
                            doubleArray  &q,
                            doubleArray  &ss,
                            doubleArray  &rel_vol,
                            Real_t       &dtcourant,
                            Real_t       &dthydro  )
{
   // calculate nodal forces, accelerations, velocities, positions, with applied
   // boundary conditions and slide surface considerations

   LagrangeNodalLeapFrog( ElemToNode, ss, pressure, q, Z_mass,
                          ref_vol, rel_vol, n_mass, BCnodes, deltatime,
                          force, accel, vel, coords);                              // return results

   // calculate zonal quantities (i.e. velocity gradient & q), and update material states

   LagrangeElements( ElemToNode, coords, vel, ref_vol, Z_mass,
                     ElemFaceToElem, matElemlist, elemBC, deltatime,
                     new_vol, vdov, char_len, delta_v,                             // return results
                     qq, ql, pressure, energy, q, ss, rel_vol);                    // more   results


   CalcTimeConstraintsForElems( matElemlist, char_len, ss, vdov,
                                dtcourant, dthydro );                              // return results

   // LagrangeRelease() ;  Creation/destruction of temps may be important to capture
}


int main(int argc, char *argv[])
{

   Index::setBoundsCheck (On);

   Real_t tx, ty, tz ;
   Index_t nidx, zidx ;

   Index_t faceElems = edgeElems * edgeElems ;
   Index_t edgeNodes = edgeElems+1 ;

   /****************************/
   /*   Initialize Sedov Mesh  */
   /****************************/

   /* construct a uniform box for this processor */

   Int_t   numElem = edgeElems * edgeElems * edgeElems ;
   Int_t   numNode = edgeNodes * edgeNodes * edgeNodes ;

   /* get run options to measure various metrics */

   //  Initialize constant array of coefficients:  Gamma

    Real_t plus1  = Real_t( 1.0);
    Real_t minus1 = Real_t(-1.0);

    Gamma(1,1) =  plus1;  Gamma(1,2)=   plus1;  Gamma(1,3)=   plus1;  Gamma(1,4)=  minus1;
    Gamma(2,1) =  plus1;  Gamma(2,2)=  minus1;  Gamma(2,3)=  minus1;  Gamma(2,4)=   plus1;
    Gamma(3,1) = minus1;  Gamma(3,2)=  minus1;  Gamma(3,3)=   plus1;  Gamma(3,4)=  minus1;
    Gamma(4,1) = minus1;  Gamma(4,2)=   plus1;  Gamma(4,3)=  minus1;  Gamma(4,4)=   plus1;
    Gamma(5,1) = minus1;  Gamma(5,2)=  minus1;  Gamma(5,3)=   plus1;  Gamma(5,4)=   plus1;
    Gamma(6,1) = minus1;  Gamma(6,2)=   plus1;  Gamma(6,3)=  minus1;  Gamma(6,4)=  minus1;
    Gamma(7,1) =  plus1;  Gamma(7,2)=   plus1;  Gamma(7,3)=   plus1;  Gamma(7,4)=   plus1;
    Gamma(8,1) =  plus1;  Gamma(8,2)=  minus1;  Gamma(8,3)=  minus1;  Gamma(8,4)=  minus1;



   // Basic Field Initialization

   energy   = Real_t(0.0);
   pressure = Real_t(0.0);
   rel_vol  = Real_t(1.0);

   vel   = Real_t(0.0);
   accel = Real_t(0.0);

   // initialize nodal coordinates

   zidx = 1;
   for (Index_t     plane=1; plane<= edgeNodes; ++plane) {
      for (Index_t    row=1; row  <= edgeNodes; ++row)   {
         for (Index_t col=1; col  <= edgeNodes; ++col)   {
            coords( zidx, 1) = Real_t(col  -1) ;                       // X coordinate
            coords( zidx, 2) = Real_t(row  -1) ;                       // Y coordinate
            coords( zidx, 3) = Real_t(plane-1) ;                       // Z coordinate
            zidx += 1;
         }
      }
   }
   coords *= Real_t(1.125)/Real_t(edgeElems) ;

   // embed hexehedral elements in nodal point lattice

   const Int_t      UP = edgeNodes*edgeNodes;
   const Int_t FORWARD = edgeNodes;
   const Int_t   RIGHT = 1;
   nidx = 1 ;
   zidx = 1 ;
   for (Index_t     plane=1; plane <= edgeElems; ++plane) {
      for (Index_t    row=1; row   <= edgeElems; ++row)   {
         for (Index_t col=1; col   <= edgeElems; ++col)   {
            ElemToNode( zidx, 1) = nidx                        ;
            ElemToNode( zidx, 2) = nidx                + RIGHT ;
            ElemToNode( zidx, 3) = nidx      + FORWARD + RIGHT ;
            ElemToNode( zidx, 4) = nidx      + FORWARD         ;
            ElemToNode( zidx, 5) = nidx + UP                   ;
            ElemToNode( zidx, 6) = nidx + UP           + RIGHT ;
            ElemToNode( zidx, 7) = nidx + UP + FORWARD + RIGHT ;
            ElemToNode( zidx, 8) = nidx + UP + FORWARD         ;
            ++nidx ;
            ++zidx ;
         }
         ++nidx ;
      }
      nidx += edgeNodes ;
   }

   // Create a material IndexSet (entire domain same material for now)

   for (Index_t i=1; i <= numElem; ++i) {
      matElemlist(i) = i ;
   }

   // initialize field data

   n_mass = 0.0;

   doubleArray ElemCoords(Corner_index, Coord_index);

   for (Index_t i=1; i<= numElem; ++i) {
       CollectElemCoords( i, ElemToNode, coords, ElemCoords);
       Real_t volume = CalcElemVolume(ElemCoords);                  // volume calculations
       ref_vol(i)    = volume ;
       Z_mass(i)     = volume ;
       for (Index_t j=1; j<=8; ++j) {
          n_mass(ElemToNode(i,j)) += volume / Real_t(8.0) ;        // distr. Elem mass to corners
      }
   }

   // deposit energy

   energy(1) = Real_t(3.948746e+7) ;

   // set up symmetry nodesets

   nidx = 1 ;
   for (Index_t i=0; i<edgeNodes; ++i) {
      Index_t planeInc = i*edgeNodes*edgeNodes ;
      Index_t rowInc   = i*edgeNodes ;
      for (Index_t j=0; j<edgeNodes; ++j) {
         BCnodes(nidx,1) = planeInc + j*edgeNodes + 1 ;        // x=0 Face
         BCnodes(nidx,2) = planeInc + j + 1 ;                  // y=0 Face
         BCnodes(nidx,3) = rowInc   + j + 1 ;                  // z=0 Face
         ++nidx ;
      }
   }

   // set up connectivity of each element face to a specific node
   // very tricky code to adapt when changing from 0-origin to 1-origin node numbering

   Int_t      eta = 1;
   Int_t       xi = 2;
   Int_t     zeta = 3;
   Int_t        M = 1;
   Int_t        P = 2;

   Int_t     x_offset = 1 ;
   Int_t     e_offset = edgeElems             ;
   Int_t     z_offset = edgeElems * edgeElems ;

   Int_t       xi_cut = numElem ;
   Int_t      eta_cut = numElem  - edgeElems ;
   Int_t     zeta_cut = numElem  - edgeElems * edgeElems  ;

   Index_t ElemId = 1;
   for (Index_t kk = 1; kk <= edgeElems; ++kk) {
       for (Index_t jj = 1; jj <= edgeElems; ++jj) {
          for (Index_t ii = 1; ii <= edgeElems; ++ii) {

             ElemFaceToElem( ElemId,   xi, M) =  (ElemId - x_offset) * ( ii > 1);
             ElemFaceToElem( ElemId,   xi, P) =  (ElemId + x_offset) * ( ii < edgeElems);
             ElemFaceToElem( ElemId,  eta, M) =  (ElemId - e_offset) * ( jj > 1);
             ElemFaceToElem( ElemId,  eta, P) =  (ElemId + e_offset) * ( jj < edgeElems);
             ElemFaceToElem( ElemId, zeta, M) =  (ElemId - z_offset) * ( kk > 1);
             ElemFaceToElem( ElemId, zeta, P) =  (ElemId + z_offset) * ( kk < edgeElems);
             ElemId += 1;
           };
       };
   };

   // set up boundary condition information
   // faces on "external" boundaries will be symmetry plane or free surface BCs

   elemBC = 0;

   for (Index_t i=0; i<edgeElems; ++i) {
      Index_t planeInc = i*edgeElems*edgeElems ;
      Index_t rowInc   = i*edgeElems ;
      for (Index_t j=0; j<edgeElems; ++j) {
         elemBC(planeInc+j*edgeElems+1)                     |= XI_M_SYMM   ;
         elemBC(planeInc+j*edgeElems+edgeElems)             |= XI_P_FREE   ;
         elemBC(planeInc+j+1)                               |= ETA_M_SYMM  ;
         elemBC(planeInc+j+edgeElems*edgeElems-edgeElems+1) |= ETA_P_FREE  ;
         elemBC(rowInc+j+1)                                 |= ZETA_M_SYMM ;
         elemBC(rowInc+j+numElem-edgeElems*edgeElems+1)     |= ZETA_P_FREE ;
      }
   }

   // timestep to solution

   while( (CurrentTime <  mc_stoptime) ) {
      TimeIncrement( dtcourant, dthydro, deltatime, CurrentTime, cycle) ;

      LagrangeLeapFrog( ElemToNode, ElemFaceToElem, matElemlist, elemBC, BCnodes,
                        Z_mass, ref_vol, n_mass, deltatime,
                        force, accel, vel, coords, new_vol, vdov, char_len,
                        delta_v, qq, ql, pressure, energy, q, ss,
                        rel_vol, dtcourant, dthydro  ) ;

      /* problem->commNodes->Transfer(CommNodes::syncposvel) ; */

   }

   return 0 ;
}

