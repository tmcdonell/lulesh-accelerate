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


