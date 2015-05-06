{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ViewPatterns    #-}

module LULESH where

import Type

import Prelude                                          as P
import Data.Array.Accelerate                            as A
import Data.Array.Accelerate.Linear                     as L
import Data.Array.Accelerate.Control.Lens               as L hiding ( _1, _2, _3, _4, _5, _6, _7, _8, _9 )


-- | Get the indices of the nodes surrounding the element at a given index. This
-- follows the numbering convention:
--
-- <<../images/nodes.png>>
--
collectDomainNodesToElemNodes
    :: Acc (Field Point)
    -> Exp Ix
    -> Exp (Hexahedron Point)
collectDomainNodesToElemNodes mesh ix@(unlift -> Z :. z :. y :. x)
  = let n0 = mesh ! ix
        n1 = mesh ! index3 z     y     (x+1)
        n2 = mesh ! index3 z     (y+1) (x+1)
        n3 = mesh ! index3 z     (y+1) x
        n4 = mesh ! index3 (z+1) y      x
        n5 = mesh ! index3 (z+1) y     (x+1)
        n6 = mesh ! index3 (z+1) (y+1) (x+1)
        n7 = mesh ! index3 (z+1) (y+1) x
    in
    lift (n0, n1, n2, n3, n4, n5, n6, n7)


-- | Calculate nodal forces, accelerations, velocities, and positions, with
-- applied boundary conditions and slide surface considerations.
--
lagrangeNodal :: () -- Acc (Array DIM3 Node)
lagrangeNodal = ()

  -- Time of boundary condition evaluation is beginning of step for force and
  -- acceleration boundary conditions
  --
  -- calcForceForNodes

  -- calcAccelerationForNodes
  -- ApplyAccelerationBoundaryConditionsForNodes
  -- CalcVelocityForNodes
  -- CalcPositionForNodes



-- | Calculate the three-dimensional force vector F at each mesh node based on
-- the values of mesh variables at time t_n.
--
-- A volume force contribution is calculated within each mesh element. This is
-- then distributed to the surrounding nodes.
--
calcForceForNodes :: ()
calcForceForNodes = ()


-- | Calculate the volume force contribute for each hexahedral mesh element. The
-- main steps are:
--
--  1. Initialise stress terms for each element
--  2. Integrate the volumetric stress terms for each element
--  3. Calculate the hourglass control contribution for each element.
--
calcVolumeForceForElems :: ()
calcVolumeForceForElems =
  let
      -- Sum contributions to total stress tensor
      -- sigma     = initStressTermsForElems undefined undefined

      -- call elemlib stress integration loop to produce nodal forces from
      -- material stresses
      -- _         = integrateStressForElems sigma
  in
  ()


-- | Initialize stress terms for each element. Our assumption of an inviscid
-- isotropic stress tensor implies that the three principal stress components
-- are equal, and the shear stresses are zero. Thus, we initialize the diagonal
-- terms of the stress tensor sigma to −(p + q) in each element.
--
initStressTermsForElems
    :: Acc (Field Pressure)
    -> Acc (Field Viscosity)
    -> Acc (Field Sigma)
initStressTermsForElems =
  A.zipWith (\p (view _z -> q) -> let s = -p - q in lift (V3 s s s))


-- | Integrate the volumetric stress contributions for each element.
--
integrateStressForElems
    :: Exp (Hexahedron Point)
    -> Exp Sigma
    -> Exp (Hexahedron Force)
integrateStressForElems p sigma =
  let
      -- Volume calculation involves extra work for numerical consistency
      det    = calcElemShapeFunctionDerivatives p ^._1
      b      = calcElemNodeNormals p
      f      = sumElemStressesToNodeForces b sigma

      -- At this point, the forces at each of the corners of the hexahedron
      -- defining this element would be distributed to the nodal mesh. This
      -- corresponds to a global scatter operation.
      --
      -- Instead, we define the problem over the nodal mesh. At each node, visit
      -- the (up to) eight surrounding elements, and calculate the contribution
      -- from that element at this node.
  in
  f

{--
integrateStressForNode
    :: Acc (Field Pressure)
    -> Acc (Field Viscosity)
    -> 
--}


calcElemShapeFunctionDerivatives
    :: Exp (Hexahedron Point)                   -- node coordinates bounding this hexahedron
    -> Exp (Hexahedron Normal, Volume)          -- (shape function derivatives, jacobian determinant (volume))
calcElemShapeFunctionDerivatives p =
  let
      -- distance between diametrically opposed corners of the hexahedron
      d60       = p^._6 - p^._0
      d53       = p^._5 - p^._3
      d71       = p^._7 - p^._1
      d42       = p^._4 - p^._2

      -- what is this??
      fjxi      = 0.125 * ( d60 + d53 - d71 - d42 )
      fjet      = 0.125 * ( d60 - d53 + d71 - d42 )
      fjze      = 0.125 * ( d60 + d53 + d71 + d42 )

      -- calculate cofactors
      cjxi      =         cross fjet fjze
      cjet      = negate (cross fjxi fjze)
      cjze      =         cross fjxi fjet

      -- calculate partials
      -- By symmetry, [6,7,4,5] = - [0,1,2,3]
      b0        = - cjxi - cjet - cjze
      b1        =   cjxi - cjet - cjze
      b2        =   cjxi + cjet - cjze
      b3        = - cjxi + cjet - cjze
      b4        = -b2
      b5        = -b3
      b6        = -b0
      b7        = -b1

      -- calculate jacobian determinant (volume)
      volume    = 0.8 * dot fjet cjet
  in
  lift ((b0, b1, b2, b3, b4, b5, b6, b7), volume)


surfaceElemFaceNormal
    :: Exp (Quad Point)
    -> Exp Normal
surfaceElemFaceNormal p =
  let
      bisectx   = 0.5 * (p^._3 + p^._2 - p^._1 - p^._0)
      bisecty   = 0.5 * (p^._2 + p^._1 - p^._3 - p^._0)
  in
  0.25 * cross bisectx bisecty


calcElemNodeNormals
    :: Exp (Hexahedron Point)
    -> Exp (Hexahedron Normal)
calcElemNodeNormals p =
  let
      -- calculate the normal to each of the 6 faces of the hexahedron
      n0123     = surfaceElemFaceNormal (lift (p^._0, p^._1, p^._2, p^._3))
      n0451     = surfaceElemFaceNormal (lift (p^._0, p^._4, p^._5, p^._1))
      n1562     = surfaceElemFaceNormal (lift (p^._1, p^._5, p^._6, p^._2))
      n2673     = surfaceElemFaceNormal (lift (p^._2, p^._6, p^._7, p^._3))
      n3740     = surfaceElemFaceNormal (lift (p^._3, p^._7, p^._4, p^._0))
      n4765     = surfaceElemFaceNormal (lift (p^._4, p^._7, p^._6, p^._5))

      -- The normal at each node is the sum of the normals of the three faces
      -- that meet at that node
  in
  lift ( n0123 + n0451 + n3740
       , n0123 + n0451 + n1562
       , n0123 + n1562 + n2673
       , n0123 + n2673 + n3740
       , n0451 + n3740 + n4765
       , n0451 + n1562 + n4765
       , n1562 + n2673 + n4765
       , n2673 + n3740 + n4765
       )


sumElemStressesToNodeForces
    :: Exp (Hexahedron Normal)
    -> Exp Sigma
    -> Exp (Hexahedron Force)
sumElemStressesToNodeForces pf sigma =
  lift ( - sigma * pf^._0
       , - sigma * pf^._1
       , - sigma * pf^._2
       , - sigma * pf^._3
       , - sigma * pf^._4
       , - sigma * pf^._5
       , - sigma * pf^._6
       , - sigma * pf^._7
       )

