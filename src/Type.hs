
module Type where

import Data.Array.Accelerate
import Data.Array.Accelerate.Linear

import Control.Lens                             ( Lens, Field1, Field2, Field3, Field4, Field5, Field6, Field7, Field8 )
import qualified Control.Lens.Tuple             as L


type R                  = Double        -- Floating point representation
type I                  = Int32         -- Indices and subscripts
type Ix                 = DIM3

type Field a            = Array Ix a
type Quad a             = (a, a, a, a)
type Hexahedron a       = (a, a, a, a, a, a, a, a)

-- Nodal quantities
type Position           = V3 R          -- position vector (old: x,y,z)
type Velocity           = V3 R          -- velocity vector (old: xd, yd, zd)
type Acceleration       = V3 R          -- acceleration vector (old: xdd, ydd, zdd)
type Force              = V3 R          -- force vector (old: fx, fy, fz)
type Mass               = R             -- nodal mass (old: nodalMass)

-- Element quantities
type Pressure           = R             -- pressure (old: p)
type Energy             = R             -- internal energy (old: e)
type Volume             = R             -- relative volume (old: v)
type Viscosity          = V3 R          -- artificial viscosity (old: qq, ql, q)
type Epsilon            = V3 R          -- diagonal terms of deviatoric strain (old: dxx, dyy, dzz)

-- Other useful type synonyms
type Normal             = V3 R          -- normal vector
type Sigma              = V3 R          -- stress term
type Timestep           = R


-- Node numbering in LULESH is zero-based. Rename the one-based tuple accessors
-- to match the diagrams.

_0 :: Field1 s t a b => Lens s t a b
_0 = L._1

_1 :: Field2 s t a b => Lens s t a b
_1 = L._2

_2 :: Field3 s t a b => Lens s t a b
_2 = L._3

_3 :: Field4 s t a b => Lens s t a b
_3 = L._4

_4 :: Field5 s t a b => Lens s t a b
_4 = L._5

_5 :: Field6 s t a b => Lens s t a b
_5 = L._6

_6 :: Field7 s t a b => Lens s t a b
_6 = L._7

_7 :: Field8 s t a b => Lens s t a b
_7 = L._8

