{-# LANGUAGE RebindableSyntax #-}
{-# LANGUAGE ViewPatterns     #-}

module Init where

import Options
import Type

import Data.Array.Accelerate                    as A hiding ( fromInteger, V3 )
import Data.Array.Accelerate.Linear             as A

import Prelude                                  ( fromInteger )
import qualified Prelude                        as P


-- | Deposit some energy at the origin. The simulation is symmetric so we only
-- simulate one quadrant, being sure to maintain the boundary conditions.
--
initEnergy :: Int -> Acc (Field Energy)
initEnergy numElem =
  let
      sh = constant (Z :. numElem :. numElem :. numElem)

      f :: Exp Ix -> Exp Energy
      f (unlift -> Z :. z :. y :. x) =
        if z == 0 && y == 0 && x == 0
           then 3.948746e+7
           else 0
  in
  A.generate sh f


-- | Initialise the nodal coordinates to a regular hexahedron mesh. The
-- coordinates of the mesh nodes will change as the simulation progresses.
--
-- We don't need the nodal point lattice to record the indices of our neighbours
-- because we have native multidimensional arrays.
--
initMesh :: Int -> Acc (Field Position)
initMesh numElem =
  let
      numNode   = numElem + 1
      sh        = constant (Z :. numNode :. numNode :. numNode)
      n         = P.fromIntegral numElem

      f :: Exp Ix -> Exp Position
      f (unlift -> Z :. k :. j :. i) =
        let x = _WIDTH  * A.fromIntegral i / n
            y = _HEIGHT * A.fromIntegral j / n
            z = _DEPTH  * A.fromIntegral k / n
        in
        lift (V3 x y z)
  in
  A.generate sh f


-- | Initialise the volume of each element.
--
-- Since we begin with a regular hexahedral mesh we just compute the volume
-- directly and initialise all elements to that value.
--
initElemVolume :: Int -> Acc (Field Volume)
initElemVolume numElem =
  let
      sh = constant (Z :. numElem :. numElem :. numElem)
      w  = _WIDTH  / P.fromIntegral numElem
      h  = _HEIGHT / P.fromIntegral numElem
      d  = _DEPTH  / P.fromIntegral numElem
      v  = w * h * d
  in
  A.fill sh v


-- | Initialise the mass at each node. This is the average of the contribution
-- of each of the surrounding elements.
--
-- Again, since we begin with a regular mesh, we just compute this value
-- directly, but we could equivalently read from the array of element volumes.
--
initNodeMass :: Int -> Acc (Field Mass)
initNodeMass numElem =
  let
      numNode   = numElem + 1
      sh        = constant (Z :. numNode :. numNode :. numNode)

      w = _WIDTH  / P.fromIntegral numElem
      h = _HEIGHT / P.fromIntegral numElem
      d = _DEPTH  / P.fromIntegral numElem
      v = w * h * d

      at z y x =
        if 0 <= z && z < constant numElem &&
           0 <= y && y < constant numElem &&
           0 <= x && x < constant numElem
           then v
           else 0

      -- This corresponds to the node -> surrounding elements index mapping
      neighbours :: Exp Ix -> Exp Mass
      neighbours (unlift -> Z :. z :. y :. x)
        = ( at z      y     x
          + at z      y    (x-1)
          + at z     (y-1) (x-1)
          + at z     (y-1)  x
          + at (z-1)  y     x
          + at (z-1)  y    (x-1)
          + at (z-1) (y-1) (x-1)
          + at (z-1) (y-1)  x
          ) / 8
  in
  generate sh neighbours

