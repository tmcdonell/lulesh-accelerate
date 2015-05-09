{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ViewPatterns    #-}

module LULESH where

import Type

import Prelude                                          as P
import Data.Array.Accelerate                            as A
import Data.Array.Accelerate.Linear                     as L
import Data.Array.Accelerate.Control.Lens               as L hiding ( _1, _2, _3, _4, _5, _6, _7, _8, _9 )

import Options
import Domain
import Data.Array.Accelerate.Interpreter


-- -----------------------------------------------------------------------------
-- TESTING
-- -----------------------------------------------------------------------------

domain :: Domain
domain = initDomain defaultOpts

elemNodes :: Acc (Field (Hexahedron Point))
elemNodes =
  let m                 = A.use (mesh domain)
      Z :. z :. y :. x  = unlift (shape m)
      sh                = index3 (z-1) (y-1) (x-1)
  in
  A.generate sh (collectDomainNodesToElemNodes m)

stress :: Acc (Field Sigma)
stress =
  initStressTermsForElems (A.use $ pressure domain)
                          (A.use $ viscosity domain)

-- -----------------------------------------------------------------------------
-- END TESTING BLOCK
-- -----------------------------------------------------------------------------


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
--
--  2. Integrate the volumetric stress terms for each element
--
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


-- Calculate the shape function derivative for the element. This is used to
-- compute the velocity gradient of the element.
--
calcElemShapeFunctionDerivatives
    :: Exp (Hexahedron Point)                   -- node coordinates bounding this hexahedron
    -> Exp (Hexahedron Normal, Volume)          -- (shape function derivatives, jacobian determinant (volume))
calcElemShapeFunctionDerivatives p =
  let
      -- (distance between diametrically opposed corners of the hexahedron??)
      d60       = p^._6 - p^._0
      d53       = p^._5 - p^._3
      d71       = p^._7 - p^._1
      d42       = p^._4 - p^._2

      -- what is this??
      -- 0.125 = 1/6
      fj_xi     = 0.125 * ( d60 + d53 - d71 - d42 )
      fj_eta    = 0.125 * ( d60 - d53 + d71 - d42 )
      fj_zeta   = 0.125 * ( d60 + d53 + d71 + d42 )

      -- calculate cofactors (= determinant??)
      cj_xi     =         cross fj_eta fj_zeta
      cj_eta    = negate (cross fj_xi  fj_zeta)
      cj_zeta   =         cross fj_xi  fj_eta

      -- calculate partials
      -- By symmetry, [6,7,4,5] = - [0,1,2,3]
      b0        = - cj_xi - cj_eta - cj_zeta
      b1        =   cj_xi - cj_eta - cj_zeta
      b2        =   cj_xi + cj_eta - cj_zeta
      b3        = - cj_xi + cj_eta - cj_zeta
      b4        = -b2
      b5        = -b3
      b6        = -b0
      b7        = -b1

      -- calculate jacobian determinant (volume)
      volume    = 0.8 * dot fj_eta cj_eta
  in
  lift ((b0, b1, b2, b3, b4, b5, b6, b7), volume)



-- | Calculate normal vectors at element nodes, as an interpolation of element
-- face normals.
--
--  1. The normal at each node of the element is initially zero
--
--  2. Enumerate all six faces of the element. For each face, calculate a normal
--     vector, scale the magnitude by one quarter, and sum the scaled vector
--     into each of the four nodes of the element corresponding to a face.
--
calcElemNodeNormals
    :: Exp (Hexahedron Point)
    -> Exp (Hexahedron Normal)
calcElemNodeNormals p =
  let
      -- Calculate a face normal
      --
      surfaceElemFaceNormal :: Exp (Quad Point) -> Exp Normal
      surfaceElemFaceNormal p =
        let
            bisectx   = 0.5 * (p^._3 + p^._2 - p^._1 - p^._0)
            bisecty   = 0.5 * (p^._2 + p^._1 - p^._3 - p^._0)
        in
        0.25 * cross bisectx bisecty

      -- The normals at each of the six faces of the hexahedron.
      --
      -- The direction that we trace out the coordinates forming a face is such
      -- that it points towards the inside the hexahedron (RH-rule)
      --
      n0123     = surfaceElemFaceNormal (lift (p^._0, p^._1, p^._2, p^._3))
      n0451     = surfaceElemFaceNormal (lift (p^._0, p^._4, p^._5, p^._1))
      n1562     = surfaceElemFaceNormal (lift (p^._1, p^._5, p^._6, p^._2))
      n2673     = surfaceElemFaceNormal (lift (p^._2, p^._6, p^._7, p^._3))
      n3740     = surfaceElemFaceNormal (lift (p^._3, p^._7, p^._4, p^._0))
      n4765     = surfaceElemFaceNormal (lift (p^._4, p^._7, p^._6, p^._5))

      -- The normal at each node is then the sum of the normals of the three
      -- faces that meet at that node.
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


-- | Sum force contribution in element to local vector for each node around
-- element.
--
sumElemStressesToNodeForces
    :: Exp (Hexahedron Normal)
    -> Exp Sigma
    -> Exp (Hexahedron Force)
sumElemStressesToNodeForces pf sigma =
  over each (\x -> -sigma * x) pf


-- | Calculate the volume derivatives for an element. Starting with a formula
-- for the volume of a hexahedron, take the derivative of that volume formula
-- with respect to the coordinates at one of the nodes. By symmetry, the formula
-- for one node can be applied to each of the other seven nodes
--
calcElemVolumeDerivative
    :: Exp (Hexahedron Point)
    -> Exp (Hexahedron (V3 R))
calcElemVolumeDerivative p =
  let
      volumeDerivative :: Exp (V3 R, V3 R, V3 R, V3 R, V3 R, V3 R) -> Exp (V3 R)
      volumeDerivative p =
        let p01 = p^._0 + p^._1
            p12 = p^._1 + p^._2
            p04 = p^._0 + p^._4
            p34 = p^._3 + p^._4
            p25 = p^._2 + p^._5
            p35 = p^._3 + p^._5
        in
        (1/12) * cross p12 p01 + cross p04 p34 + cross p35 p25
  in
  lift ( volumeDerivative (lift (p^._1, p^._2, p^._3, p^._4, p^._5, p^._7))
       , volumeDerivative (lift (p^._0, p^._1, p^._2, p^._7, p^._4, p^._6))
       , volumeDerivative (lift (p^._3, p^._0, p^._1, p^._6, p^._7, p^._5))
       , volumeDerivative (lift (p^._2, p^._3, p^._0, p^._5, p^._6, p^._4))
       , volumeDerivative (lift (p^._7, p^._6, p^._5, p^._0, p^._3, p^._1))
       , volumeDerivative (lift (p^._4, p^._7, p^._6, p^._1, p^._0, p^._2))
       , volumeDerivative (lift (p^._5, p^._4, p^._7, p^._2, p^._1, p^._3))
       , volumeDerivative (lift (p^._6, p^._5, p^._4, p^._3, p^._2, p^._0))
       )


-- Calculate the hourglass control contribution for each element. For each
-- element:
--
--  1. Gather the node coordinates for that element.
--
--  2. Calculate the element volume derivative.
--
--  3. Perform a diagnosis check for any element volumes <= zero
--
--  4. Compute the Flanagan-Belytschko hourglass control force for each element.
--  This is described in the paper:
--
--    [1] "A uniform strain hexahedron and quadrilateral with orthogonal
--        hourglass control", Flanagan, D. P. and Belytschko, T. International
--        Journal for Numerical Methods in Engineering, (17) 5, May 1981.
--
calcHourglassControlForElems
    :: ()
calcHourglassControlForElems = ()


calcFBHourglassForceForElems
    :: Exp (Hexahedron Point)   -- from collectDomainNodesToElemNodes
    -> Exp (Hexahedron Velocity)
    -> Exp Volume               -- from calcElemShapeFunctionDerivatives
    -> Exp (Hexahedron (V3 R))  -- from calcElemVolumeDerivative
    -> Exp (Hexahedron Force)
calcFBHourglassForceForElems pos vel vol dvol =
  let
      -- Hourglass base vectors, from [1] table 1. This defines the hourglass
      -- patterns for a unit cube.
      --
      gamma :: Exp (Hexahedron (V4 R))
      gamma = constant
        ( V4 ( 1) ( 1) ( 1) (-1)
        , V4 ( 1) (-1) (-1) ( 1)
        , V4 (-1) (-1) ( 1) (-1)
        , V4 (-1) ( 1) (-1) ( 1)
        , V4 (-1) (-1) ( 1) ( 1)
        , V4 (-1) ( 1) (-1) (-1)
        , V4 ( 1) ( 1) ( 1) ( 1)
        , V4 ( 1) (-1) (-1) (-1)
        )

      -- Compute hourglass modes
      --
      hourgam :: Exp (Hexahedron (V4 R))
      hourgam =
        let hg :: Exp (V4 R) -> Exp (Point) -> Exp (V3 R) -> Exp (V4 R)
            hg g p dv   = (1 - volinv * dot dv p) *^ g

            volinv      = 1 / vol
        in
        lift ( hg (gamma^._0) (pos^._0) (dvol^._0)
             , hg (gamma^._1) (pos^._1) (dvol^._1)
             , hg (gamma^._2) (pos^._2) (dvol^._2)
             , hg (gamma^._3) (pos^._3) (dvol^._3)
             , hg (gamma^._4) (pos^._4) (dvol^._4)
             , hg (gamma^._5) (pos^._5) (dvol^._5)
             , hg (gamma^._6) (pos^._6) (dvol^._6)
             , hg (gamma^._7) (pos^._7) (dvol^._7)
             )

      -- Compute forces
      coefficient = undefined -- hourg * 0.01 * ss1 * mass1 / volume13

      f           = undefined -- calcElemFBHourglassForce vel hourgam coefficient

{--
      hourmodx =
            x8n[i3]   * gamma[i1][0] +
            x8n[i3+1] * gamma[i1][1] +
            x8n[i3+2] * gamma[i1][2] +
            x8n[i3+3] * gamma[i1][3] +
            x8n[i3+4] * gamma[i1][4] +
            x8n[i3+5] * gamma[i1][5] +
            x8n[i3+6] * gamma[i1][6] +
            x8n[i3+7] * gamma[i1][7];

      hourmody =
            y8n[i3]   * gamma[i1][0] +
            y8n[i3+1] * gamma[i1][1] +
            y8n[i3+2] * gamma[i1][2] +
            y8n[i3+3] * gamma[i1][3] +
            y8n[i3+4] * gamma[i1][4] +
            y8n[i3+5] * gamma[i1][5] +
            y8n[i3+6] * gamma[i1][6] +
            y8n[i3+7] * gamma[i1][7];

      hourmodz =
            z8n[i3]   * gamma[i1][0] +
            z8n[i3+1] * gamma[i1][1] +
            z8n[i3+2] * gamma[i1][2] +
            z8n[i3+3] * gamma[i1][3] +
            z8n[i3+4] * gamma[i1][4] +
            z8n[i3+5] * gamma[i1][5] +
            z8n[i3+6] * gamma[i1][6] +
            z8n[i3+7] * gamma[i1][7];

      hourgam0[i1] = gamma[i1][0] -  volinv * ( dvdx[i3  ] * hourmodx + dvdy[i3  ] * hourmody + dvdz[i3  ] * hourmodz );
      hourgam1[i1] = gamma[i1][1] -  volinv * ( dvdx[i3+1] * hourmodx + dvdy[i3+1] * hourmody + dvdz[i3+1] * hourmodz );
      hourgam2[i1] = gamma[i1][2] -  volinv * ( dvdx[i3+2] * hourmodx + dvdy[i3+2] * hourmody + dvdz[i3+2] * hourmodz );
      hourgam3[i1] = gamma[i1][3] -  volinv * ( dvdx[i3+3] * hourmodx + dvdy[i3+3] * hourmody + dvdz[i3+3] * hourmodz );
      hourgam4[i1] = gamma[i1][4] -  volinv * ( dvdx[i3+4] * hourmodx + dvdy[i3+4] * hourmody + dvdz[i3+4] * hourmodz );
      hourgam5[i1] = gamma[i1][5] -  volinv * ( dvdx[i3+5] * hourmodx + dvdy[i3+5] * hourmody + dvdz[i3+5] * hourmodz );
      hourgam6[i1] = gamma[i1][6] -  volinv * ( dvdx[i3+6] * hourmodx + dvdy[i3+6] * hourmody + dvdz[i3+6] * hourmodz );
      hourgam7[i1] = gamma[i1][7] -  volinv * ( dvdx[i3+7] * hourmodx + dvdy[i3+7] * hourmody + dvdz[i3+7] * hourmodz );
--}
  in
  f


calcElemFBHourglassForce
    :: Exp R
    -> Exp (Hexahedron Velocity)
    -> Exp (Hexahedron (V4 R))
    -> Exp (Hexahedron Force)
calcElemFBHourglassForce coefficient vel hourgam =
  let
      h00, h01, h02, h03 :: Exp (V3 R)
      h00 = P.sum $ P.zipWith (*^) (hourgam ^.. (each._x)) (vel ^.. each)
      h01 = P.sum $ P.zipWith (*^) (hourgam ^.. (each._y)) (vel ^.. each)
      h02 = P.sum $ P.zipWith (*^) (hourgam ^.. (each._z)) (vel ^.. each)
      h03 = P.sum $ P.zipWith (*^) (hourgam ^.. (each._w)) (vel ^.. each)

      hh :: Exp (V4 (V3 R))
      hh  = lift (V4 h00 h01 h02 h03)

      hg :: Exp (V4 R) -> Exp Force
      hg h = coefficient *^ (P.sum $ P.zipWith (*^) (h^..each) (hh^..each))
  in
  over each hg hourgam

{--

   h00 = hourgam0[i00] * xd[0] + hourgam1[i00] * xd[1] + hourgam2[i00] * xd[2] + hourgam3[i00] * xd[3] + hourgam4[i00] * xd[4] + hourgam5[i00] * xd[5] + hourgam6[i00] * xd[6] + hourgam7[i00] * xd[7];
   h00 = hourgam0[i00] * yd[0] + hourgam1[i00] * yd[1] + hourgam2[i00] * yd[2] + hourgam3[i00] * yd[3] + hourgam4[i00] * yd[4] + hourgam5[i00] * yd[5] + hourgam6[i00] * yd[6] + hourgam7[i00] * yd[7];
   h00 = hourgam0[i00] * zd[0] + hourgam1[i00] * zd[1] + hourgam2[i00] * zd[2] + hourgam3[i00] * zd[3] + hourgam4[i00] * zd[4] + hourgam5[i00] * zd[5] + hourgam6[i00] * zd[6] + hourgam7[i00] * zd[7];

   h01 = hourgam0[i01] * xd[0] + hourgam1[i01] * xd[1] + hourgam2[i01] * xd[2] + hourgam3[i01] * xd[3] + hourgam4[i01] * xd[4] + hourgam5[i01] * xd[5] + hourgam6[i01] * xd[6] + hourgam7[i01] * xd[7];
   h01 = hourgam0[i01] * yd[0] + hourgam1[i01] * yd[1] + hourgam2[i01] * yd[2] + hourgam3[i01] * yd[3] + hourgam4[i01] * yd[4] + hourgam5[i01] * yd[5] + hourgam6[i01] * yd[6] + hourgam7[i01] * yd[7];
   h01 = hourgam0[i01] * zd[0] + hourgam1[i01] * zd[1] + hourgam2[i01] * zd[2] + hourgam3[i01] * zd[3] + hourgam4[i01] * zd[4] + hourgam5[i01] * zd[5] + hourgam6[i01] * zd[6] + hourgam7[i01] * zd[7];

   h02 = hourgam0[i02] * xd[0] + hourgam1[i02] * xd[1] + hourgam2[i02] * xd[2] + hourgam3[i02] * xd[3] + hourgam4[i02] * xd[4] + hourgam5[i02] * xd[5] + hourgam6[i02] * xd[6] + hourgam7[i02] * xd[7];
   h02 = hourgam0[i02] * yd[0] + hourgam1[i02] * yd[1] + hourgam2[i02] * yd[2] + hourgam3[i02] * yd[3] + hourgam4[i02] * yd[4] + hourgam5[i02] * yd[5] + hourgam6[i02] * yd[6] + hourgam7[i02] * yd[7];
   h02 = hourgam0[i02] * zd[0] + hourgam1[i02] * zd[1] + hourgam2[i02] * zd[2] + hourgam3[i02] * zd[3] + hourgam4[i02] * zd[4] + hourgam5[i02] * zd[5] + hourgam6[i02] * zd[6] + hourgam7[i02] * zd[7];

   h03 = hourgam0[i03] * xd[0] + hourgam1[i03] * xd[1] + hourgam2[i03] * xd[2] + hourgam3[i03] * xd[3] + hourgam4[i03] * xd[4] + hourgam5[i03] * xd[5] + hourgam6[i03] * xd[6] + hourgam7[i03] * xd[7];
   h03 = hourgam0[i03] * yd[0] + hourgam1[i03] * yd[1] + hourgam2[i03] * yd[2] + hourgam3[i03] * yd[3] + hourgam4[i03] * yd[4] + hourgam5[i03] * yd[5] + hourgam6[i03] * yd[6] + hourgam7[i03] * yd[7];
   h03 = hourgam0[i03] * zd[0] + hourgam1[i03] * zd[1] + hourgam2[i03] * zd[2] + hourgam3[i03] * zd[3] + hourgam4[i03] * zd[4] + hourgam5[i03] * zd[5] + hourgam6[i03] * zd[6] + hourgam7[i03] * zd[7];

   hgfx[0] = coefficient * (hourgam0[i00] * h00 + hourgam0[i01] * h01 + hourgam0[i02] * h02 + hourgam0[i03] * h03);
   hgfy[0] = coefficient * (hourgam0[i00] * h00 + hourgam0[i01] * h01 + hourgam0[i02] * h02 + hourgam0[i03] * h03);
   hgfz[0] = coefficient * (hourgam0[i00] * h00 + hourgam0[i01] * h01 + hourgam0[i02] * h02 + hourgam0[i03] * h03);




   hgfx[0] = coefficient * (hourgam0[i00] * h00 + hourgam0[i01] * h01 + hourgam0[i02] * h02 + hourgam0[i03] * h03);
   hgfx[1] = coefficient * (hourgam1[i00] * h00 + hourgam1[i01] * h01 + hourgam1[i02] * h02 + hourgam1[i03] * h03);
   hgfx[2] = coefficient * (hourgam2[i00] * h00 + hourgam2[i01] * h01 + hourgam2[i02] * h02 + hourgam2[i03] * h03);
   hgfx[3] = coefficient * (hourgam3[i00] * h00 + hourgam3[i01] * h01 + hourgam3[i02] * h02 + hourgam3[i03] * h03);
   hgfx[4] = coefficient * (hourgam4[i00] * h00 + hourgam4[i01] * h01 + hourgam4[i02] * h02 + hourgam4[i03] * h03);
   hgfx[5] = coefficient * (hourgam5[i00] * h00 + hourgam5[i01] * h01 + hourgam5[i02] * h02 + hourgam5[i03] * h03);
   hgfx[6] = coefficient * (hourgam6[i00] * h00 + hourgam6[i01] * h01 + hourgam6[i02] * h02 + hourgam6[i03] * h03);
   hgfx[7] = coefficient * (hourgam7[i00] * h00 + hourgam7[i01] * h01 + hourgam7[i02] * h02 + hourgam7[i03] * h03);

   hgfy[0] = coefficient * (hourgam0[i00] * h00 + hourgam0[i01] * h01 + hourgam0[i02] * h02 + hourgam0[i03] * h03);
   hgfy[1] = coefficient * (hourgam1[i00] * h00 + hourgam1[i01] * h01 + hourgam1[i02] * h02 + hourgam1[i03] * h03);
   hgfy[2] = coefficient * (hourgam2[i00] * h00 + hourgam2[i01] * h01 + hourgam2[i02] * h02 + hourgam2[i03] * h03);
   hgfy[3] = coefficient * (hourgam3[i00] * h00 + hourgam3[i01] * h01 + hourgam3[i02] * h02 + hourgam3[i03] * h03);
   hgfy[4] = coefficient * (hourgam4[i00] * h00 + hourgam4[i01] * h01 + hourgam4[i02] * h02 + hourgam4[i03] * h03);
   hgfy[5] = coefficient * (hourgam5[i00] * h00 + hourgam5[i01] * h01 + hourgam5[i02] * h02 + hourgam5[i03] * h03);
   hgfy[6] = coefficient * (hourgam6[i00] * h00 + hourgam6[i01] * h01 + hourgam6[i02] * h02 + hourgam6[i03] * h03);
   hgfy[7] = coefficient * (hourgam7[i00] * h00 + hourgam7[i01] * h01 + hourgam7[i02] * h02 + hourgam7[i03] * h03);

   hgfz[0] = coefficient * (hourgam0[i00] * h00 + hourgam0[i01] * h01 + hourgam0[i02] * h02 + hourgam0[i03] * h03);
   hgfz[1] = coefficient * (hourgam1[i00] * h00 + hourgam1[i01] * h01 + hourgam1[i02] * h02 + hourgam1[i03] * h03);
   hgfz[2] = coefficient * (hourgam2[i00] * h00 + hourgam2[i01] * h01 + hourgam2[i02] * h02 + hourgam2[i03] * h03);
   hgfz[3] = coefficient * (hourgam3[i00] * h00 + hourgam3[i01] * h01 + hourgam3[i02] * h02 + hourgam3[i03] * h03);
   hgfz[4] = coefficient * (hourgam4[i00] * h00 + hourgam4[i01] * h01 + hourgam4[i02] * h02 + hourgam4[i03] * h03);
   hgfz[5] = coefficient * (hourgam5[i00] * h00 + hourgam5[i01] * h01 + hourgam5[i02] * h02 + hourgam5[i03] * h03);
   hgfz[6] = coefficient * (hourgam6[i00] * h00 + hourgam6[i01] * h01 + hourgam6[i02] * h02 + hourgam6[i03] * h03);
   hgfz[7] = coefficient * (hourgam7[i00] * h00 + hourgam7[i01] * h01 + hourgam7[i02] * h02 + hourgam7[i03] * h03);
--}

