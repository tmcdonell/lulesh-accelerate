{-# LANGUAGE RecordWildCards #-}

module Domain where

import Type

import Prelude                                  as P
import Data.Array.Accelerate                    as A


type Domain =
    ( Field Position
    , Field Velocity
    , Field Energy
    , Field Pressure
    , Field Viscosity
    , Field Volume
    , Field SoundSpeed
    , Scalar Time                 -- current simulation time
    , Scalar Time                 -- delta t
    )

data Parameters = Parameters
  {
    -- Simulation constants
    -- --------------------
    e_cut               :: Exp R                -- energy tolerance
  , p_cut               :: Exp R                -- pressure tolerance
  , q_cut               :: Exp R                -- viscosity tolerance
  , v_cut               :: Exp R                -- relative volume tolerance
  , u_cut               :: Exp R                -- velocity tolerance
  , p_min               :: Exp R                -- pressure floor
  , e_min               :: Exp R                -- energy floor

  , dt_scaling          :: (Exp R, Exp R)
  , dt_max              :: Exp R                -- maximum allowable time increment
  , t_end               :: R                    -- end time of the simulation

    -- Other constants
    -- ---------------
  , hgcoef              :: Exp R                -- hourglass coefficient
  , qstop               :: Exp R                -- excessive viscosity indicator
  , monoq_max_slope     :: Exp R
  , monoq_limiter       :: Exp R
  , qlc_monoq           :: Exp R                -- linear term coefficient for viscosity
  , qqc_monoq           :: Exp R                -- quadratic term for coefficient for viscosity
  , qqc                 :: Exp R
  , eosvmax             :: Exp R
  , eosvmin             :: Exp R
  , dvovmax             :: Exp R                -- maximum allowable volume change
  , ref_dens            :: Exp R                -- reference density
  }
  deriving Show


parameters :: Parameters
parameters = Parameters
  {
    e_cut               = 1.0e-7
  , p_cut               = 1.0e-7
  , q_cut               = 1.0e-7
  , u_cut               = 1.0e-7
  , v_cut               = 1.0e-7
  , p_min               = 0
  , e_min               = -1.0e15

  , dt_scaling          = (1.1, 1.2)
  , dt_max              = 1.0e-2
  , t_end               = 1.0e-2

  , hgcoef              = 3.0
  , qstop               = 1.0e12
  , monoq_max_slope     = 1.0
  , monoq_limiter       = 2.0
  , qlc_monoq           = 0.5
  , qqc_monoq           = 2.0/3.0
  , qqc                 = 2.0
  , eosvmax             = 1.0e+9
  , eosvmin             = 1.0e-9
  , dvovmax             = 0.1
  , ref_dens            = 1.0
  }

