
module Domain where

import Options

import Prelude                          as P
import Data.Array.Accelerate            as A

type R = Double                         -- Floating point representation
type I = Int32                          -- array subscript and loop index

data Error
  = Volume
  | QStop
  deriving (Eq, Enum, Show)


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
    coordinate          :: Vector (R,R,R)       -- (x, y, z)
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

  , energy              :: Vector R             -- e
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
  , ss4o3               :: R
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
  , time_stop           :: R                    -- simulation stop time
  , dt                  :: R                    -- variable time increment
  , dt_mult_bounds      :: (R,R)                -- lower and upper bounds to scale dt by
  , dt_max              :: R                    -- maximum allowable time increment
  , dt_courant          :: R                    -- courant constraint
  , dt_hydro            :: R                    -- volume change constraint
  }
  deriving Show


-- Initialise a domain with the starting configuration
--
initDomain :: Options -> Domain
initDomain = error "TODO: initDomain"

