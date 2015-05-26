{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE RebindableSyntax    #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ViewPatterns        #-}

module LULESH where

import Domain
import Options
import Type

import Prelude                                          as P hiding ( (<*) )
import Data.Array.Accelerate                            as A
import Data.Array.Accelerate.Linear                     as L
import Data.Array.Accelerate.Control.Lens               as L hiding ( _1, _2, _3, _4, _5, _6, _7, _8, _9, at, ix, use )


-- -----------------------------------------------------------------------------
-- TESTING
-- -----------------------------------------------------------------------------

import Data.Array.Accelerate.Interpreter                as I

domain :: Domain
domain = initDomain defaultOpts

elemNodes :: Acc (Field (Hexahedron Position))
elemNodes =
  let m                 = A.use (mesh domain)
      Z :. z :. y :. x  = unlift (shape m)
      sh                = index3 (z-1) (y-1) (x-1)
  in
  A.generate sh (collectToElem m)

-- stress :: Acc (Field Sigma)
-- stress =
--   initStressTermsForElems (A.use $ pressure domain)
--                           (A.use $ viscosity domain)


step1 :: Domain -> Acc (Field Force)
step1 Domain{..}
  = calcForceForNodes (constant hgcoef)
                      (use mesh)
                      (use velocity)
                      (use pressure)
                      (use viscosity)
                      (use volume)
                      (use volume_ref)
                      (use ss)
                      (use elemMass)


-- -----------------------------------------------------------------------------
-- END TESTING BLOCK
-- -----------------------------------------------------------------------------

-- Utilities
-- ---------


-- Use rebindable syntax to improve readability of branches in scalar code.
--
ifThenElse :: Elt t => Exp Bool -> Exp t -> Exp t -> Exp t
ifThenElse x y z
  = x ? (y, z)


-- | Get the values at the nodes surrounding the element of a given index. This
-- follows the numbering convention:
--
-- <<../images/nodes.png>>
--
collectToElem
    :: Elt a
    => Acc (Field a)
    -> Exp Ix
    -> Exp (Hexahedron a)
collectToElem mesh ix@(unlift -> Z :. z :. y :. x)
  = let
        n0 = mesh ! ix
        n1 = mesh ! index3 z     y     (x+1)
        n2 = mesh ! index3 z     (y+1) (x+1)
        n3 = mesh ! index3 z     (y+1) x
        n4 = mesh ! index3 (z+1) y      x
        n5 = mesh ! index3 (z+1) y     (x+1)
        n6 = mesh ! index3 (z+1) (y+1) (x+1)
        n7 = mesh ! index3 (z+1) (y+1) x
    in
    lift (n0, n1, n2, n3, n4, n5, n6, n7)


-- Given a set of values at the nodes surrounding each element, form a nodal
-- grid by combining all those elements.
--
distributeToNode
    :: forall a. Elt a
    => (Exp a -> Exp a -> Exp a)
    -> Exp a
    -> Acc (Field (Hexahedron a))
    -> Acc (Field a)
distributeToNode f zero arr =
  let
      numElem           = unindex3 (shape arr) ^. _2
      numNode           = numElem + 1
      sh'               = index3 numNode numNode numNode

      at :: Exp Int -> Exp Int -> Exp Int
         -> Lens' (Exp (Hexahedron a)) (Exp a)
         -> Exp a
      at z y x node =
        if 0 <=* z &&* z <* numElem &&* 0 <=* y &&* y <* numElem &&* 0 <=* x &&* x <* numElem
           then (arr ! index3 z y x) ^. node
           else zero

      mesh :: Exp Ix -> Exp a
      mesh (unlift -> Z :. z :. y :. x) =
        at z     y     x     _0 `f`
        at z     y     (x-1) _1 `f`
        at z     (y-1) (x-1) _2 `f`
        at z     (y-1) x     _3 `f`
        at (z-1) y     x     _4 `f`
        at (z-1) y     (x-1) _5 `f`
        at (z-1) (y-1) (x-1) _6 `f`
        at (z-1) (y-1) x     _7
  in
  generate sh' mesh


-- Lagrange Leapfrog Algorithm
-- ===========================

-- | 'lagrangeLeapFrog' advances the solution from t_n to t_{n+1} over the time
-- increment delta_t. The process of advance the solution is comprised of two
-- major parts:
--
--   1. Advance variables on the nodal mesh; and
--   2. Advance the element variables
--
lagrangeLeapFrog :: ()
lagrangeLeapFrog = ()
  -- lagrangeNodal
  -- lagrangeElements
  -- calcTimeConstraintsForElems


-- Advance Node Quantities
-- -----------------------

-- | Advance the nodal mesh variables, primarily the velocity and position. The
-- main steps are:
--
--   1. Calculate the nodal forces: 'calcForceForNodes'
--   2. Calculate nodal accelerations: 'calcAccelerationForNodes'
--   3. Apply acceleration boundary conditions ('applyAccelerationBoundaryConditionsForNodes, but called from (2))
--   4. Integrate nodal accelerations to obtain updated velocities: 'calcVelocityForNodes'
--   5. Integrate nodal velocities to obtain updated positions: 'calcPositionForNodes'
--
lagrangeNodal
    :: Exp R                    -- hourglass coefficient
    -> Exp R                    -- acceleration cutoff
    -> Exp R                    -- timestep
    -> Acc (Field Position)
    -> Acc (Field Velocity)
    -> Acc (Field Pressure)
    -> Acc (Field Viscosity)
    -> Acc (Field Volume)       -- relative volume
    -> Acc (Field Volume)       -- reference volume
    -> Acc (Field R)            -- speed of sound
    -> Acc (Field Mass)         -- element mass
    -> Acc (Field Mass)         -- nodal mass
    -> ( Acc (Field Position), Acc (Field Velocity) )
lagrangeNodal hgcoef ucut dt position velocity pressure viscosity volume volumeRef soundSpeed elemMass nodalMass =
  let
      -- Time of boundary condition evaluation is beginning of step for force
      -- and acceleration boundary conditions
      force             = calcForceForNodes hgcoef position velocity pressure viscosity volume volumeRef soundSpeed elemMass
      acceleration      = calcAccelerationForNodes force nodalMass
      velocity'         = calcVelocityForNodes dt ucut velocity acceleration
      position'         = calcPositionForNodes dt position velocity'
  in
  (position', velocity')


-- | Calculate the three-dimensional force vector F at each mesh node based on
-- the values of mesh variables at time t_n.
--
-- A volume force contribution is calculated within each mesh element. This is
-- then distributed to the surrounding nodes.
--
calcForceForNodes
    :: Exp R
    -> Acc (Field Position)
    -> Acc (Field Velocity)
    -> Acc (Field Pressure)
    -> Acc (Field Viscosity)
    -> Acc (Field Volume)       -- volume
    -> Acc (Field Volume)       -- reference valume
    -> Acc (Field R)
    -> Acc (Field Mass)
    -> Acc (Field Force)
calcForceForNodes hgcoef position velocity pressure viscosity volume volumeRef soundSpeed elemMass
  = distributeToNode (+) 0
  $ calcVolumeForceForElems hgcoef position velocity pressure viscosity volume volumeRef soundSpeed elemMass


-- | Calculate the volume force contribute for each hexahedral mesh element. The
-- main steps are:
--
--  1. Initialise stress terms for each element
--  2. Integrate the volumetric stress terms for each element
--  3. Calculate the hourglass control contribution for each element.
--
calcVolumeForceForElems
    :: Exp R
    -> Acc (Field Position)
    -> Acc (Field Velocity)
    -> Acc (Field Pressure)
    -> Acc (Field Viscosity)
    -> Acc (Field Volume)
    -> Acc (Field Volume)
    -> Acc (Field R)
    -> Acc (Field Mass)
    -> Acc (Field (Hexahedron Force))
calcVolumeForceForElems hgcoef position velocity pressure viscosity volume volumeRef soundSpeed elemMass =
  let
      numNode           = unindex3 (shape position) ^. _2
      numElem           = numNode - 1
      sh                = index3 numElem numElem numElem

      -- sum contributions to total stress tensor
      sigma             = A.zipWith initStressTermsForElems pressure viscosity

      -- calculate nodal forces from element stresses
      (stress, _determ) = A.unzip
                        $ A.zipWith integrateStressForElems
                                    (generate sh (collectToElem position))
                                    sigma

      -- TODO: check for negative element volume
      -- A.any (<=* 0) determ --> error

      -- Calculate the hourglass control contribution for each element
      hourglass         = generate sh $ \ix ->
        let pos         = collectToElem position ix
            vel         = collectToElem velocity ix

            v           = volume     ! ix
            volo        = volumeRef  ! ix
            ss          = soundSpeed ! ix
            mass        = elemMass   ! ix
        in
        calcHourglassControlForElems pos vel volo v ss mass hgcoef

      -- Add the nodal forces
      combine :: Exp (Hexahedron Force) -> Exp (Hexahedron Force) -> Exp (Hexahedron Force)
      combine x y       = lift ( x^._0 + y^._0
                               , x^._1 + y^._1
                               , x^._2 + y^._2
                               , x^._3 + y^._3
                               , x^._4 + y^._4
                               , x^._5 + y^._5
                               , x^._6 + y^._6
                               , x^._7 + y^._7 )
  in
  A.zipWith combine stress hourglass


-- | Initialize stress terms for each element. Our assumption of an inviscid
-- isotropic stress tensor implies that the three principal stress components
-- are equal, and the shear stresses are zero. Thus, we initialize the diagonal
-- terms of the stress tensor sigma to −(p + q) in each element.
--
initStressTermsForElems
    :: Exp Pressure
    -> Exp Viscosity
    -> Exp Sigma
initStressTermsForElems p (view _z -> q) =
  let s = -p - q
  in  lift (V3 s s s)


-- | Integrate the volumetric stress contributions for each element.
--
-- In the reference LULESH code, the forces at each of the corners of the
-- hexahedron defining this element would be distributed to the nodal mesh. This
-- corresponds to a global scatter operation.
--
-- Instead, we just return all the values directly, and the individual
-- contributions to the nodes will be combined in a different step.
--
integrateStressForElems
    :: Exp (Hexahedron Position)
    -> Exp Sigma
    -> Exp (Hexahedron Force, Volume)
integrateStressForElems p sigma =
  let
      -- Volume calculation involves extra work for numerical consistency
      det    = calcElemShapeFunctionDerivatives p ^._1
      b      = calcElemNodeNormals p
      f      = sumElemStressesToNodeForces b sigma
  in
  lift (f, det)


-- Calculate the shape function derivative for the element. This is used to
-- compute the velocity gradient of the element.
--
calcElemShapeFunctionDerivatives
    :: Exp (Hexahedron Position)                -- node coordinates bounding this hexahedron
    -> Exp (Hexahedron Force, Volume)           -- (shape function derivatives, jacobian determinant (volume))
calcElemShapeFunctionDerivatives p =
  let
      -- compute diagonal differences
      d60       = p^._6 - p^._0
      d53       = p^._5 - p^._3
      d71       = p^._7 - p^._1
      d42       = p^._4 - p^._2

      -- compute jacobians
      fj_xi     = 0.125 * ( d60 + d53 - d71 - d42 )
      fj_eta    = 0.125 * ( d60 - d53 + d71 - d42 )
      fj_zeta   = 0.125 * ( d60 + d53 + d71 + d42 )

      -- calculate cofactors (= determinant??)
      cj_xi     = cross fj_eta  fj_zeta
      cj_eta    = cross fj_zeta fj_xi
      cj_zeta   = cross fj_xi   fj_eta

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
    :: Exp (Hexahedron Position)
    -> Exp (Hexahedron Normal)
calcElemNodeNormals p =
  let
      -- Calculate a face normal
      --
      surfaceElemFaceNormal :: Exp (Quad Position) -> Exp Normal
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
  over each (\x -> -sigma * x) pf       -- interesting shorthand to map over a tuple


-- | Calculate the volume derivatives for an element. Starting with a formula
-- for the volume of a hexahedron, take the derivative of that volume formula
-- with respect to the coordinates at one of the nodes. By symmetry, the formula
-- for one node can be applied to each of the other seven nodes
--
calcElemVolumeDerivative
    :: Exp (Hexahedron Position)
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


-- Calculate the hourglass control contribution for each element.
--
-- For each element:
--
--  1. Gather the node coordinates for that element.
--  2. Calculate the element volume derivative.
--  3. Perform a diagnosis check for any element volumes <= zero
--  4. Compute the Flanagan-Belytschko hourglass control force for each element.
--     This is described in the paper:
--
--     [1] "A uniform strain hexahedron and quadrilateral with orthogonal
--         hourglass control", Flanagan, D. P. and Belytschko, T. International
--         Journal for Numerical Methods in Engineering, (17) 5, May 1981.
--
calcHourglassControlForElems
    :: Exp (Hexahedron Position)
    -> Exp (Hexahedron Velocity)
    -> Exp Volume                       -- relative volume
    -> Exp Volume                       -- reference volume
    -> Exp R                            -- speed of sound
    -> Exp Mass                         -- mass
    -> Exp R
    -> Exp (Hexahedron Force)
calcHourglassControlForElems pos vel volo v ss mass hourg =
  let dvol      = calcElemVolumeDerivative pos
      determ    = volo * v
  in
  if hourg >* 0
     then calcFBHourglassForceForElems pos vel determ dvol ss mass hourg
     else constant (0,0,0,0,0,0,0,0)


calcFBHourglassForceForElems
    :: Exp (Hexahedron Position)
    -> Exp (Hexahedron Velocity)
    -> Exp Volume
    -> Exp (Hexahedron (V3 R))          -- from calcElemVolumeDerivative
    -> Exp R                            -- speed of sound
    -> Exp Mass                         -- mass
    -> Exp R
    -> Exp (Hexahedron Force)
calcFBHourglassForceForElems pos vel determ dvol ss mass hourg =
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
        let hg :: Exp (V4 R) -> Exp (Position) -> Exp (V3 R) -> Exp (V4 R)
            hg g p dv   = (1 - volinv * dot dv p) *^ g

            volinv      = 1 / determ
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
      cbrt x      = x ** (1/3)          -- cube root
      coefficient = - hourg * 0.01 * ss * mass / cbrt determ
  in
  calcElemFBHourglassForce coefficient vel hourgam


calcElemFBHourglassForce
    :: Exp R
    -> Exp (Hexahedron Velocity)
    -> Exp (Hexahedron (V4 R))
    -> Exp (Hexahedron Force)
calcElemFBHourglassForce coefficient vel hourgam =
  let
      -- TLM: this looks like a small matrix multiplication?

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


-- | Calculate the three-dimensional acceleration vector at each mesh node, and
-- apply the symmetry boundary conditions.
--
calcAccelerationForNodes
    :: Acc (Field Force)
    -> Acc (Field Mass)
    -> Acc (Field Acceleration)
calcAccelerationForNodes force mass
  = applyAccelerationBoundaryConditionsForNodes
  $ A.zipWith (^/) force mass


-- | Applies symmetry boundary conditions at nodes on the boundaries of the
-- mesh. This sets the normal component of the acceleration vector at the
-- boundary to zero. This implies that the normal component of the velocity
-- vector will remain constant in time.
--
-- Recall that the benchmark Sedov problem is spherically-symmetric and that we
-- simulate it in a cubic domain containing a single octant of the sphere. To
-- maintain spherical symmetry of the domain, we apply symmetry boundary
-- conditions along the faces of the cubic domain that contact the planes
-- separating the octants of the sphere. This forces the normal component of the
-- velocity vector to be zero along these boundary faces for all time, since
-- they were initialised to zero.
--
applyAccelerationBoundaryConditionsForNodes
    :: Acc (Field Acceleration)
    -> Acc (Field Acceleration)
applyAccelerationBoundaryConditionsForNodes acc =
  generate (shape acc) $ \ix ->
    let
        Z :. z :. y :. x        = unlift ix
        V3 xd yd zd             = unlift $ acc ! ix
    in
    lift $ V3 (x ==* 0 ? (0, xd))
              (y ==* 0 ? (0, yd))
              (z ==* 0 ? (0, zd))


-- | Integrate the acceleration at each node to advance the velocity at the
-- node.
--
-- Note that the routine applies a cutoff to each velocity vector value.
-- Specifically, if a value is below some prescribed threshold the term is set
-- to zero. The reason for this cutoff is to prevent spurious mesh motion which
-- may arise due to floating point roundoff error when the velocity is near
-- zero.
--
calcVelocityForNodes
    :: Exp R
    -> Exp R
    -> Acc (Field Velocity)
    -> Acc (Field Acceleration)
    -> Acc (Field Velocity)
calcVelocityForNodes dt ucut u ud
  = A.map (over each (\x -> abs x <* ucut ? (0,x)))
  $ integrate dt u ud

-- | Integrate the velocity at each node to advance the position of the node
--
calcPositionForNodes
    :: Exp R
    -> Acc (Field Position)
    -> Acc (Field Velocity)
    -> Acc (Field Position)
calcPositionForNodes = integrate


-- | Euler integration
--
integrate
    :: Exp R
    -> Acc (Field (V3 R))
    -> Acc (Field (V3 R))
    -> Acc (Field (V3 R))
integrate dt
  = A.zipWith (\x xd -> x + xd ^* dt)


-- Advance Element Quantities
-- --------------------------

-- | Advance element quantities, primarily pressure, internal energy, and
-- relative volume. The artificial viscosity in each element is also calculated
-- here. The main steps are:
--
--   1. Calculate element quantities based on nodal kinematic quantities
--   2. Calculate element artificial viscosity terms
--   3. Apply material properties in each element needed to calculate updated
--      pressure and internal energy.
--   4. Compute updated element volume
--
lagrangeElements
    :: Exp R                    -- timestep
    -> Acc (Field Position)
    -> Acc (Field Velocity)
    -> Acc (Field Volume)
    -> Acc (Field Mass)
    -> ()
lagrangeElements dt position velocity volume _ =
  ()


-- | Calculate various element quantities that are based on the new kinematic
-- node quantities position and velocity.
--
calcLagrangeElements
    :: Exp R                    -- timestep
    -> Acc (Field Position)
    -> Acc (Field Velocity)
    -> Acc (Field Volume)       -- relative volume
    -> Acc (Field Volume)       -- reference volume
    -> ()
calcLagrangeElements dt pos vel relVol refVol =
  ()


-- | Calculate terms in the total strain rate tensor epsilon_tot that are used
-- to compute the terms in the deviatoric strain rate tensor epsilon.
--
calcKinematicsForElems
    :: Exp R                    -- timestep
    -> Acc (Field Position)
    -> Acc (Field Velocity)
    -> Acc (Field Volume)       -- relative volume
    -> Acc (Field Volume)       -- reference volume
    -> ()
calcKinematicsForElems dt pos vel relVol refVol =
  ()


-- | Calculate the volume of an element given the nodal coordinates
--
calcElemVolume
    :: Exp (Hexahedron Position)
    -> Exp Volume
calcElemVolume p =
  let
      -- compute diagonal differences
      d61       = p^._6 - p^._1
      d70       = p^._7 - p^._0
      d63       = p^._6 - p^._3
      d20       = p^._2 - p^._0
      d50       = p^._5 - p^._0
      d64       = p^._6 - p^._4
      d31       = p^._3 - p^._1
      d72       = p^._7 - p^._2
      d43       = p^._4 - p^._3
      d57       = p^._5 - p^._7
      d14       = p^._1 - p^._4
      d25       = p^._2 - p^._5
  in
  (1/12) * ( triple (d31 + d72) d63 d20
           + triple (d43 + d57) d64 d70
           + triple (d14 + d25) d61 d50
           )

