{-# LANGUAGE BangPatterns    #-}
{-# LANGUAGE CPP             #-}
{-# LANGUAGE GADTs           #-}
{-# LANGUAGE TemplateHaskell #-}
#if defined(ACCELERATE_LLVM_NATIVE_BACKEND) && __GLASGOW_HASKELL__ < 806
{-# OPTIONS_GHC -fplugin=Data.Array.Accelerate.LLVM.Native.Plugin #-}
#endif
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

import Data.Array.Accelerate                            as A hiding ( Ord(..), V3 )
import Data.Array.Accelerate.Sugar.Shape                as S
import Data.Array.Accelerate.Linear                     as A
import Data.Array.Accelerate.Control.Lens               as L hiding ( _1, _2, _3, _4, _5, _6, _7, _8, _9 )

import Prelude                                          as P
import System.IO
import Text.Printf

#if __GLASGOW_HASKELL__ >= 800
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
import qualified Data.Array.Accelerate.LLVM.Native      as CPU
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
import qualified Data.Array.Accelerate.LLVM.PTX         as PTX
#endif
#endif


main :: IO ()
main = do
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
      mN0       = run backend $ initNodeMass numElem
      v0        = run backend $ initElemVolume numElem
      vrel0     = A.fill (constant (Z:.numElem:.numElem:.numElem)) 1
      zeros     = A.fill (constant (Z:.numElem:.numElem:.numElem)) 0
      dt0       = unit 1.0e-7
      t0        = unit 0

      write     = case opts ^. optOutputPath of
                    Nothing -> \_ _   -> return ()
                    Just p  -> \i dom -> writeDomain p i dom

      initial :: Acc Domain
      initial = lift (x0, dx0, e0, p0, q0, vrel0, ss0, t0, dt0)

      simulate :: Int -> Domain -> IO (Domain, Int)
      simulate !i !dom@(_, _, _, _, _, _, _, t, _)
        | i                >= maxSteps         = return (dom, i)
        | t `indexArray` Z >= t_end parameters = return (dom, i)
        | otherwise                            = do
            let !dom' = step dom
            write i dom
            simulate (i+1) dom'
        where
          !step = case backend of
#if __GLASGOW_HASKELL__ >= 800
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
            CPU -> $( CPU.runQ lulesh ) mN0 v0
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
            PTX -> $( PTX.runQ lulesh ) mN0 v0
#endif
#endif
            _   -> run1 backend (lulesh (A.use mN0) (A.use v0))

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

