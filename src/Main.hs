{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE CPP          #-}
{-# LANGUAGE ViewPatterns #-}
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

import Backend
import Domain
import Init
import LULESH
import Options
import Timing                                           ( time )
import Type
import VisIt

import Data.Array.Accelerate                            as A
import Data.Array.Accelerate.Array.Sugar                as S
import Data.Array.Accelerate.Linear                     as A
import Data.Array.Accelerate.Control.Lens               as L hiding ( _1, _2, _3, _4, _5, _6, _7, _8, _9 )
import Data.Array.Accelerate.Debug                      ( beginMonitoring )

import Prelude                                          as P hiding ( (<*) )
import System.IO
import Text.Printf


main :: IO ()
main = do
  beginMonitoring
  (opts,_)      <- parseArgs

  let backend   = opts ^. optBackend
      numElem   = opts ^. optSize
      numNode   = numElem + 1
      maxSteps  = view optMaxSteps opts

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

      write     = case opts ^. optOutputPath of
                    Nothing -> \_ _   -> return ()
                    Just p  -> \i dom -> writeDomain p i dom

      initial :: Acc Domain
      initial = lift (x0, dx0, e0, p0, q0, vrel0, ss0, t0, dt0)

      step :: Acc Domain -> Acc Domain
      step domain =
        let
            (x, dx, e, p, q, v, ss, t, dt) = unlift domain

            (x', dx', e', p', q', v', ss', dtc, dth)
                = lagrangeLeapFrog parameters (the dt) x dx e p q v v0 ss v0 mN0

            (t', dt')
                = timeIncrement parameters t dt dtc dth
        in
        lift (x', dx', e', p', q', v', ss', t', dt')

      simulate :: Int -> Domain -> IO (Domain, Int)
      simulate !i !dom@(_, _, _, _, _, _, _, t, _)
        | i                >= maxSteps         = return (dom, i)
        | t `indexArray` Z >= t_end parameters = return (dom, i)
        | otherwise                            = do
            let !dom' = go dom
            write i dom
            simulate (i+1) dom'
        where
          go = run1 backend step

{--
      lulesh :: Acc (Field Volume)
             -> Acc (Field Mass)
             -> Acc Domain
             -> Acc Domain
      lulesh v0 mN0 dom0 =
          awhile
            -- loop condition
            undefined
            -- (\domain -> A.zipWith (&&*) (A.map (<* t_end parameters)  (elapsed domain))
            --                             (A.map (<* constant maxSteps) (iteration domain)))
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
--}

  -- Problem description
  --
  printf "Running problem size     : %d^3\n" numElem
  printf "Total number of elements : %d\n\n" (numElem * numElem * numElem)

  -- Initialise the accelerate computation by performing a single step
  --
  printf "Initialising accelerate...            " >> hFlush stdout
  (dom0, t1) <- time $ do
    let dom0 = run backend initial
    r <- simulate (maxSteps-1) dom0
    r `seq` return dom0
  print t1

  -- Run the simulation proper...
  --
  printf "Running simulation...                 " >> hFlush stdout
  ((result, iterations), t2) <- time $ simulate 0 dom0
  printf "%s\n\n" (show t2)

  -- Final results summary
  --
  let energy            = result ^._2
      sh                = arrayShape energy

      go !j !k !maxRelDiff !maxAbsDiff !totalAbsDiff
        | j >= numElem  = (maxRelDiff, maxAbsDiff, totalAbsDiff)
        | k >= numElem  = go (j+1) (j+2) maxRelDiff maxAbsDiff totalAbsDiff
        | otherwise     =
            let x       = energy `indexArray` S.fromIndex sh (j * numElem + k)
                y       = energy `indexArray` S.fromIndex sh (k * numElem + j)
                --
                diff    = abs (x - y)
                rel     = diff / y
            in
            go j (k+1) (maxRelDiff `P.max` rel) (maxAbsDiff `P.max` diff) (totalAbsDiff + diff)

      (relDiff, absDiff, totalDiff) = go 0 1 0 0 0

  printf "Run completed\n"
  printf "   Iteration count     : %d\n"     iterations
  printf "   Final origin energy : %.6e\n\n" (energy `indexArray` (Z:.0:.0:.0))

  printf "Testing Plane 0 of Energy Array\n"
  printf "   Maximum relative difference : %.6e\n" relDiff
  printf "   Maximum absolute difference : %.6e\n" absDiff
  printf "   Total absolute difference   : %.6e\n" totalDiff

