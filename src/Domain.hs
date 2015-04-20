{-# LANGUAGE RecordWildCards #-}

module Domain where

import Options

import Prelude                                  as P
import Data.Array.Accelerate                    as A

type R = Double         -- Floating point representation
type I = Int32          -- Indices and subscripts


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
    mesh                :: Array DIM3 (R,R,R)   -- (x, y, z)
  , velocity            :: Vector (R,R,R)       -- (xd, yd, zd)
  , acceleration        :: Vector (R,R,R)       -- (xdd, ydd, zdd)
  , force               :: Vector (R,R,R)       -- (fx, fy, fz)
  , nodalMass           :: Vector R
  , symmetry            :: Vector (I,I,I)       -- symmetry plane nodesets

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

  , energy              :: Array DIM3 R         -- e            -- Array DIM3 R
  , pressure            :: Vector R             -- p
  , viscosity           :: Vector R             -- q
  , viscosity_linear    :: Vector R             -- ql
  , viscosity_quadratic :: Vector R             -- qq
    -- TLM: viscosity = qq^2 + ql + q ??

  , volume              :: Vector R             -- v (relative)
  , volume_ref          :: Vector R             -- volo
  , volume_dov          :: Vector R             -- volume derivative over volume

  , arealg              :: Vector R             -- characteristic length of an element

    -- Courant-Friedrichs-Lewy (CFL) condition determines the maximum size of
    -- each time increment based on the shortest distance across any mesh
    -- element, and the speed of sound in the material of that element.
  , ss                  :: Vector R             -- "speed of sound"

  , elemMass            :: Vector R             -- mass

    -- Cutoffs (constants)
    -- -------------------
  , e_cut               :: R                    -- energy tolerance
  , p_cut               :: R                    -- pressure tolerance
  , q_cut               :: R                    -- viscosity tolerance
  , v_cut               :: R                    -- relative volume tolerance
  , u_cut               :: R                    -- velocity tolerance

    -- Other constants
    -- ---------------
  , hgcoef              :: R                    -- hourglass coefficient
  , qstop               :: R                    -- excessive viscosity indicator
  , monoq_max_slope     :: R
  , monoq_limiter       :: R
  , qlc_monoq           :: R                    -- linear term coefficient for viscosity
  , qqc_monoq           :: R                    -- quadratic term for coefficient for viscosity
  , qqc                 :: R
  , eosvmax             :: R
  , eosvmin             :: R
  , pmin                :: R                    -- pressure floor
  , emin                :: R                    -- energy floor
  , dvovmax             :: R                    -- maximum allowable volume change
  , refdens             :: R                    -- reference density

    -- Simulation variables
    -- --------------------
  , iteration           :: Int                  -- simulation iteration
  , time                :: R                    -- current simulation time
  , time_end            :: R                    -- simulation stop time
  , dt                  :: R                    -- variable time increment
  , dt_scaling          :: (R,R)                -- lower and upper bounds to scale dt by
  , dt_max              :: R                    -- maximum allowable time increment
  , dt_courant          :: R                    -- courant constraint
  , dt_hydro            :: R                    -- volume change constraint
  }
  deriving Show


-- Initialise a domain with the starting configuration
--
initDomain :: Options -> Domain
initDomain Options{..} =
  let
--      edgeElems         = optSize
--      edgeNodes         = optSize + 1
--
--      numElems          = edgeElems * edgeElems * edgeElems
--      numNodes          = edgeNodes * edgeNodes * edgeNodes
  in
  Domain
  {
    -- node centred
    mesh                = initMesh optSize

    -- element centred
  , energy              = initEnergy optSize

    -- constants
  , e_cut               = 1.0e-7
  , p_cut               = 1.0e-7
  , q_cut               = 1.0e-7
  , u_cut               = 1.0e-7
  , v_cut               = 1.0e-7

  , hgcoef              = 3.0
  , qstop               = 1.0e12
  , monoq_max_slope     = 1.0
  , monoq_limiter       = 2.0
  , qlc_monoq           = 0.5
  , qqc_monoq           = 2.0/3.0
  , qqc                 = 2.0
  , eosvmax             = 1.0e+9
  , eosvmin             = 1.0e-9
  , pmin                = 0
  , emin                = -1.0e15
  , dvovmax             = 0.1
  , refdens             = 1.0

    -- simulation parameters
  , iteration           = 0
  , time                = 0
  , time_end            = 1.0e-2
  , dt                  = 1.0e-7
  , dt_scaling          = (1.1, 1.2)
  , dt_max              = 1.0e-2
  , dt_courant          = 1.0e20
  , dt_hydro            = 1.0e20
  }


-- Deposit some energy at the origin. The simulation is symmetric and we only
-- simulate one quadrant, being sure to maintain the boundary conditions.
--
initEnergy :: Int -> Array DIM3 R
initEnergy numEdges
  = fromFunction (Z :. numEdges :. numEdges :. numEdges)
  $ \ix -> case ix of
             Z :. 0 :. 0 :. 0 -> 3.948746e+7
             _                -> 0

-- Initialise the nodal coordinates to a regular hexahedron mesh. The
-- coordinates of the mesh nodes will change as the simulation progresses.
--
-- We don't need the nodal point lattice to record the indices of our neighbours
-- because we have native multidimensional arrays.
--
initMesh :: Int -> Array DIM3 (R,R,R)
initMesh numEdges
  = let numNodes                = numEdges + 1
        n                       = P.fromIntegral numEdges
        f (Z :. k :. j :. i)    =
          let x = 1.125 * P.fromIntegral i / n
              y = 1.125 * P.fromIntegral j / n
              z = 1.125 * P.fromIntegral k / n
          in
          (x, y, z)
    in
    fromFunction (Z :. numNodes :. numNodes :. numNodes) f


