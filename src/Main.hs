{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ViewPatterns #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing #-}
-- |
-- Module       : Main
-- Copyright    : [2015] Trevor L. McDonell
-- License      : BSD3
--
-- Maintainer   : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability    : experimental
-- Portability  : non-portable (GHC extensions)
--
-- This applications is an implementation of the Livermore Unstructured
-- Lagrangian Explicit Shock Hydrodynamics
-- (<https://codesign.llnl.gov/lulesh.php LULESH>) mini-app in
-- <https://github.com/AccelerateHS/accelerate Accelerate>.
--
-- This shock hydrodynamics challenge problem was originally defined and
-- implemented by LLNL as one of five challenge problems in the DARPA UHPC
-- program and has since become a widely studied proxy application in DOE
-- co-design efforts for exascale.
--
-- <<https://codesign.llnl.gov/images/sedov-3d-LLNL.png>>
--

module Main where

import Domain
import Init
import LULESH
import Options
import Time
import Timing                                           ( time )
import Type
import qualified Backend                                as B

import Data.Array.Accelerate                            as A
import Data.Array.Accelerate.Array.Sugar                as S
import Data.Array.Accelerate.Linear                     as A
import Data.Array.Accelerate.Control.Lens               as L hiding ( _1, _2, _3, _4, _5, _6, _7, _8, _9 )

import Prelude                                          as P hiding ( (<*) )
import Control.Exception
import System.IO
import Text.Printf


main :: IO ()
main = do
  (opts,_)      <- parseArgs

  let
      run       = B.run  (opts ^. optBackend)
      run3      = B.run3 (opts ^. optBackend)
      numElem   = opts ^. optSize
      numNode   = numElem + 1
      maxSteps  = constant (view optMaxSteps opts)

      -- Initialise the primary data structures
      x0        = initMesh numElem
      dx0       = A.fill (constant (Z:.numNode:.numNode:.numNode)) (constant (V3 0 0 0))
      e0        = initEnergy numElem
      p0        = zeros
      q0        = zeros
      ss0       = zeros
      mN0       = initNodeMass numElem
      v0        = initElemVolume numElem
      vrel0     = A.fill (constant (Z:.numElem:.numElem:.numElem)) 1
      zeros     = A.fill (constant (Z:.numElem:.numElem:.numElem)) 0
      dt0       = unit 1.0e-7
      t0        = unit 0
      n0        = unit 0

      initial :: Acc Domain
      initial = lift (x0, dx0, e0, p0, q0, vrel0, ss0, (t0, dt0, n0))

      lulesh :: Acc (Field Volume)
             -> Acc (Field Mass)
             -> Acc Domain
             -> Acc Domain
      lulesh v0 mN0 dom0 =
          awhile
            -- loop condition
            (\domain -> A.zipWith (&&*) (A.map (<* t_end parameters) (domain^._7._0))
                                        (A.map (<* maxSteps)         (domain^._7._2)))
            -- loop body
            (\domain ->
                let
                    (x, dx, e, p, q, v, ss, r) = unlift domain
                    (t, dt, n)                 = unlift (r :: Acc (Scalar Time, Scalar Time, Scalar Int))

                    (x', dx', e', p', q', v', ss', dtc, dth)
                        = lagrangeLeapFrog parameters (the dt) x dx e p q v v0 ss v0 mN0

                    (t', dt')
                        = timeIncrement parameters t dt dtc dth

                    n'  = A.map (+1) n
                in
                lift (x', dx', e', p', q', v', ss', (t', dt', n')))
            dom0

  printf "Running problem size     : %d^3\n" numElem
  printf "Total number of elements : %d\n\n" (numElem * numElem * numElem)

  printf "Initialising accelerate...            " >> hFlush stdout
  (compute, t1)         <- time (evaluate $ run3 lulesh)
  print t1

  printf "Initialising data...                  " >> hFlush stdout
  ((v0,m0,dom0), t2)    <- time (evaluate $ run (lift (v0, mN0, initial)))
  print t2

  printf "Running simulation...                 " >> hFlush stdout
  (result, t3)          <- time (evaluate $ compute v0 m0 dom0)
  printf "%s\n\n" (show t3)

  let energy     = result ^._2
      iterations = result ^._7._2
      sh         = arrayShape energy

  printf "Run completed\n"
  printf "   Iteration count     : %d\n"     (iterations `indexArray` Z)
  printf "   Final origin energy : %.6e\n\n" (energy     `indexArray` (Z:.0:.0:.0))

  let go !j !k !maxRelDiff !maxAbsDiff !totalAbsDiff
        | j >= numElem  = (maxRelDiff, maxAbsDiff, totalAbsDiff)
        | k >= numElem  = go (j+1) (j+2) maxRelDiff maxAbsDiff totalAbsDiff
        | otherwise     =
            let x       = energy `indexArray` S.fromIndex sh (j * numElem + k)
                y       = energy `indexArray` S.fromIndex sh (k * numElem + j)

                diff    = abs (x - y)
                rel     = diff / y
            in
            go j (k+1) (maxRelDiff `max` rel) (maxAbsDiff `max` diff) (totalAbsDiff + diff)

      (relDiff, absDiff, totalDiff) = go 0 1 0 0 0

  printf "Testing Plane 0 of Energy Array\n"
  printf "   Maximum relative difference : %.6e\n" relDiff
  printf "   Maximum absolute difference : %.6e\n" absDiff
  printf "   Total absolute difference   : %.6e\n" totalDiff

