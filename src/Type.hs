
module Type where

import Data.Array.Accelerate
import Data.Array.Accelerate.Linear

import Control.Lens                             ( Lens, Field1, Field2, Field3, Field4, Field5, Field6, Field7, Field8 )
import qualified Control.Lens.Tuple             as L


type Quad a             = (a, a, a, a)
type Hexahedron a       = (a, a, a, a, a, a, a, a)

type R = Double                         -- Floating point representation
type I = Int32                          -- Indices and subscripts

type Ix                 = DIM3
type Field a            = Array Ix a
type Point              = V3 R
type Normal             = V3 R

type Pressure           = R             -- The quantities that we are interested in
type Volume             = R
type Viscosity          = V3 R
type Sigma              = V3 R
type Force              = V3 R


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

