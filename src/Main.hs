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
import Type
-- import qualified Backend                                as B

import Prelude                                          as P hiding ( (<*) )
import Data.Array.Accelerate                            as A
import Data.Array.Accelerate.Linear                     as A
import Data.Array.Accelerate.Control.Lens               as L hiding ( _1, _2, _3, _4, _5, _6, _7, _8, _9, at, ix, use )


main :: IO ()
main = do
  (opts,_)      <- parseArgs

  let
--      run       = B.run  (opts ^. optBackend)
      numElem   = opts ^. optSize
      numNode   = numElem + 1

      -- We don't have loop-invariant code motion. This forces 'x' to be
      -- evaluated before applying it in 'f'
      licm x f  = (id >-> f) x

      -- Initialise the primary data structures
      x0        = initMesh numElem
      dx0       = A.fill (constant (Z:.numNode:.numNode:.numNode)) (constant (V3 0 0 0))
      e0        = initEnergy numElem
      p0        = zeros
      q0        = zeros
      ss0       = zeros
      vrel0     = A.fill (constant (Z:.numElem:.numElem:.numElem)) 1
      zeros     = A.fill (constant (Z:.numElem:.numElem:.numElem)) 0

      -- Timestep to solution
      dt0       = unit 1.0e-7
      t0        = unit 0

      lulesh    =
        licm (initElemVolume numElem) $ \v0  ->
        licm (initNodeMass numElem)   $ \mN0 ->
          awhile
            (\domain -> A.map (<* t_end parameters) (domain ^._7))
            (\domain ->
                let
                    (x, dx, e, p, q, v, ss, t, dt) = unlift domain

                    (x', dx', e', p', q', v', ss', dtc, dth)
                        = lagrangeLeapFrog parameters (the dt) x dx e p q v v0 ss v0 mN0

                    (t', dt')
                        = timeIncrement parameters t dt dtc dth
                in
                lift (x', dx', e', p', q', v', ss', t', dt'))
            (lift (x0, dx0, e0, p0, q0, vrel0, ss0, t0, dt0))

  print lulesh

