{-# LANGUAGE RecordWildCards #-}

module Domain where

import Type
import Options

import Prelude                                  as P
import Data.Array.Accelerate                    as A
import Data.Array.Accelerate.Linear             hiding ( Epsilon )

-- The size of the simulation box
--
_WIDTH, _HEIGHT, _DEPTH :: Fractional a => a
_WIDTH  = 1.125
_HEIGHT = 1.125
_DEPTH  = 1.125


-- TODO: Pack the node-centred / elem-centred fields into a data structure
--
-- type Vec3 a = (a,a,a)
-- type Node a = (Vec3 a, Vec3 a, Vec3 a, Vec3 a, a)
-- type Elem a = (a, a, Vec3 a, a, a, a, a, a, a)


-- The domain is the primary data structure for LULESH.
--
-- A mesh is partitioned into a logically rectangular collection called the
-- Domain. All data for variables represented by this portion of the mesh are
-- contained within the Domain data structure. Note that nodes on domain
-- boundaries are shared by neighbouring domains.
--
-- Thus, we have a three level hierarchy for computation: the problem level;
-- which is the union of all domains (i.e., the whole mesh); the domain level;
-- and the element level. The problem may contain one or more domains, and a
-- domain may contain one or more elements. This allows for flexible scalability
-- from a single element problem to a very large problem with an arbitrary
-- number of elements, where the computation is defined over the same physical
-- volume of space.
--
data Domain = Domain
  {
    -- Node centred
    -- ------------
    mesh                :: Field Position               -- (x, y, z)
  , velocity            :: Field Velocity               -- (xd, yd, zd)
  , acceleration        :: Field Acceleration           -- (xdd, ydd, zdd)
  , force               :: Field Force                  -- (fx, fy, fz)
  , nodeMass            :: Field Mass

--  , symmetry            :: Vector (I,I,I)             -- symmetry plane nodesets

    -- Element centred
    -- ---------------
--  , numReg              :: ?? I
--  , cost                :: ?? I                    -- imbalance cost
--  , regElemSize         :: ?? [I]                  -- size of region sets
--  , regNumList          :: ?? [I]                  -- region number per domain element
--  , regElemList         :: ?? [Acc (Vector I)]     -- region index set

-- , nodeList           :: ??                     -- element to node connectivity
--
--   ... element connectivity across each face
--   ... symmetric/free surface flags for each element face

  , energy              :: Field Energy                 -- e
  , pressure            :: Field Pressure               -- p
  , viscosity           :: Field Viscosity              -- (qq, ql, q)
  , strains             :: Field Epsilon                -- principal strains (dxx, dyy, dzz)

--  , viscosity           :: Array DIM3 R         -- q
--  , viscosity_linear    :: Vector R             -- ql
--  , viscosity_quadratic :: Vector R             -- qq
    -- TLM: viscosity = qq^2 + ql + q ??

  , volume              :: Field Volume         -- v (relative)
  , volume_ref          :: Field Volume         -- volo (reference volume)
  , volume_dov          :: Array DIM3 R         -- volume derivative over volume

  , arealg              :: Array DIM3 R         -- characteristic length of an element

    -- Courant-Friedrichs-Lewy (CFL) condition determines the maximum size of
    -- each time increment based on the shortest distance across any mesh
    -- element, and the speed of sound in the material of that element.
  , ss                  :: Array DIM3 R         -- "speed of sound"

  , elemMass            :: Array DIM3 R         -- mass

    -- Simulation variables
    -- --------------------
  , iteration           :: Int                  -- simulation iteration
  , time                :: R                    -- current simulation time
  , dt                  :: R                    -- variable time increment
  }
  deriving Show


data Parameters = Parameters
  {
    -- Simulation constants
    -- --------------------
    e_cut               :: Exp R                -- energy tolerance
  , p_cut               :: Exp R                -- pressure tolerance
  , q_cut               :: Exp R                -- viscosity tolerance
  , v_cut               :: Exp R                -- relative volume tolerance
  , u_cut               :: Exp R                -- velocity tolerance
  , p_min               :: Exp R                -- pressure floor
  , e_min               :: Exp R                -- energy floor

  , dt_scaling          :: (Exp R, Exp R)
  , dt_max              :: Exp R                -- maximum allowable time increment
  , t_end               :: Exp R                -- end time of the simulation

    -- Other constants
    -- ---------------
  , hgcoef              :: Exp R                -- hourglass coefficient
  , qstop               :: Exp R                -- excessive viscosity indicator
  , monoq_max_slope     :: Exp R
  , monoq_limiter       :: Exp R
  , qlc_monoq           :: Exp R                -- linear term coefficient for viscosity
  , qqc_monoq           :: Exp R                -- quadratic term for coefficient for viscosity
  , qqc                 :: Exp R
  , eosvmax             :: Exp R
  , eosvmin             :: Exp R
  , dvovmax             :: Exp R                -- maximum allowable volume change
  , ref_dens            :: Exp R                -- reference density
  }
  deriving Show


-- Initialise a domain with the starting configuration
--
initDomain :: Options -> Domain
initDomain Options{..} =
  let
      numElem   = optSize
      numNode   = numElem + 1

      ev        = initElemVolume numElem
      efill v   = fromFunction (Z :. numElem :. numElem :. numElem) (const v)
      nfill v   = fromFunction (Z :. numNode :. numNode :. numNode) (const v)

      n000      = nfill (V3 0 0 0)
      e0        = efill 0
  in
  Domain
  {
    -- node centred
    mesh                = initMesh numElem
  , velocity            = n000
  , acceleration        = n000
  , force               = n000
  , nodeMass            = initNodeMass numElem

    -- element centred
  , energy              = initEnergy numElem
  , pressure            = e0
  , viscosity           = e0
  , volume              = efill 1
  , volume_ref          = ev
  , volume_dov          = e0
  , arealg              = e0
  , ss                  = e0
  , elemMass            = ev

    -- simulation parameters
  , iteration           = 0
  , time                = 0
  , dt                  = 1.0e-7
  }

parameters :: Parameters
parameters = Parameters
  {
    e_cut               = 1.0e-7
  , p_cut               = 1.0e-7
  , q_cut               = 1.0e-7
  , u_cut               = 1.0e-7
  , v_cut               = 1.0e-7
  , p_min               = 0
  , e_min               = -1.0e15

  , dt_scaling          = (1.1, 1.2)
  , dt_max              = 1.0e-2
  , t_end               = 1.0e-2

  , hgcoef              = 3.0
  , qstop               = 1.0e12
  , monoq_max_slope     = 1.0
  , monoq_limiter       = 2.0
  , qlc_monoq           = 0.5
  , qqc_monoq           = 2.0/3.0
  , qqc                 = 2.0
  , eosvmax             = 1.0e+9
  , eosvmin             = 1.0e-9
  , dvovmax             = 0.1
  , ref_dens            = 1.0
  }



-- Deposit some energy at the origin. The simulation is symmetric so we only
-- simulate one quadrant, being sure to maintain the boundary conditions.
--
initEnergy :: (Elt a, Fractional a) => Int -> Array DIM3 a
initEnergy numElem
  = fromFunction (Z :. numElem :. numElem :. numElem)
  $ \ix -> case ix of
             Z :. 0 :. 0 :. 0 -> 3.948746e+7
             _                -> 0

-- Initialise the nodal coordinates to a regular hexahedron mesh. The
-- coordinates of the mesh nodes will change as the simulation progresses.
--
-- We don't need the nodal point lattice to record the indices of our neighbours
-- because we have native multidimensional arrays.
--
initMesh :: (Elt a, Fractional a) => Int -> Array DIM3 (V3 a)
initMesh numElem
  = let numNode                 = numElem + 1
        n                       = P.fromIntegral numElem
        f (Z :. k :. j :. i)    =
          let x = _WIDTH  * P.fromIntegral i / n
              y = _HEIGHT * P.fromIntegral j / n
              z = _DEPTH  * P.fromIntegral k / n
          in
          V3 x y z
    in
    fromFunction (Z :. numNode :. numNode :. numNode) f

-- Initialise the volume of each element.
--
-- Since we begin with a regular hexahedral mesh we just compute the volume
-- directly and initialise all elements to that value.
--
initElemVolume :: (Elt a, Fractional a) => Int -> Array DIM3 a
initElemVolume numElem
  = let w = _WIDTH  / P.fromIntegral numElem
        h = _HEIGHT / P.fromIntegral numElem
        d = _DEPTH  / P.fromIntegral numElem
        v = w * h * d
    in
    fromFunction (Z :. numElem :. numElem :. numElem) (const v)

-- Initialise the mass at each node. This is the average of the contribution of
-- each of the surrounding elements.
--
-- Again, since we begin with a regular mesh, we just compute this value directly.
--
initNodeMass :: (Elt a, Fractional a) => Int -> Array DIM3 a
initNodeMass numElem
  = let numNode = numElem + 1

        w = _WIDTH  / P.fromIntegral numElem
        h = _HEIGHT / P.fromIntegral numElem
        d = _DEPTH  / P.fromIntegral numElem
        v = w * h * d

        at (Z :. z :. y :. x) =
          if 0 <= z && z < numElem && 0 <= y && y < numElem && 0 <= x && x < numElem
             then v
             else 0

        -- This corresponds to the node -> surrounding elements index mapping
        neighbours (Z :. z :. y :. x)
          = ( at (Z :. z   :. y   :. x)
            + at (Z :. z   :. y   :. x-1)
            + at (Z :. z   :. y-1 :. x-1)
            + at (Z :. z   :. y-1 :. x)
            + at (Z :. z-1 :. y   :. x)
            + at (Z :. z-1 :. y   :. x-1)
            + at (Z :. z-1 :. y-1 :. x-1)
            + at (Z :. z-1 :. y-1 :. x)
            ) / 8
    in
    fromFunction (Z :. numNode :. numNode :. numNode) neighbours

