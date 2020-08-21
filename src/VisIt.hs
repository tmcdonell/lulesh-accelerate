{-# LANGUAGE CPP                      #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# OPTIONS_GHC -fno-warn-unused-imports #-}

module VisIt where

import Type
import Domain

import Data.Array.Accelerate.Array.Data
import Data.Array.Accelerate.Control.Lens               as L hiding ( _1, _2, _3, _4, _5, _6, _7, _8, _9 )

import Text.Printf
import Foreign.C
import Foreign.Ptr


-- | Write the simulation state to file, suitable for use with the 'VisIt'
-- visualisation program. The input file path should have space (via %d) to
-- write the simulation step as part of the name.
--
writeDomain :: FilePath -> Int -> Domain -> IO ()
#ifndef ACCELERATE_VISIT
writeDomain _ _ _ = error "visualisation is not enabled: recompile with -fvisit"
#else
writeDomain basename step domain =
  withCString name $ \c_name -> do
    c_writeDomain c_name
                  (fromIntegral numElem)
                  (fromIntegral step)
                  elapsed
                  p_x
                  p_y
                  p_z
                  p_xd
                  p_yd
                  p_zd
                  p_e
                  p_p
                  p_v
                  p_q

    touchArrayData ad_pos
    touchArrayData ad_vel
    touchArrayData ad_e
    touchArrayData ad_p
    touchArrayData ad_q
    touchArrayData ad_v
  where
    name            = printf basename step

    elapsed         = (domain^._7) ! Z

    (_,numElem)     = s

    Array _ ad_pos  = domain ^._0
    Array _ ad_vel  = domain ^._1
    Array s ad_e    = domain ^._2
    Array _ ad_p    = domain ^._3
    Array _ ad_q    = domain ^._4
    Array _ ad_v    = domain ^._5

    ((((), p_x),  p_y),  p_z)   = ptrsOfArrayData ad_pos
    ((((), p_xd), p_yd), p_zd)  = ptrsOfArrayData ad_vel
    p_e                         = ptrsOfArrayData ad_e
    p_p                         = ptrsOfArrayData ad_p
    p_q                         = ptrsOfArrayData ad_q
    p_v                         = ptrsOfArrayData ad_v


foreign import ccall "writeDomain"
  c_writeDomain :: CString
                -> CInt
                -> CInt
                -> R
                -> Ptr R
                -> Ptr R
                -> Ptr R
                -> Ptr R
                -> Ptr R
                -> Ptr R
                -> Ptr R
                -> Ptr R
                -> Ptr R
                -> Ptr R
                -> IO ()
#endif

