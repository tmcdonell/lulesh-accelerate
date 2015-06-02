{-# LANGUAGE RebindableSyntax #-}
{-# LANGUAGE RecordWildCards  #-}

module Time where

import Domain
import Type
import Util

import Prelude                          as P hiding ( (<*) )
import Data.Array.Accelerate            as A


-- | Compute the time increment for the current loop iteration. We aim for a
-- "target" timestep value which completes the simulation in the next step, but
-- is only allowed to change from the previous value by a certain amount,
-- subject to courant and hydro constraints.
--
timeIncrement
    :: Parameters
    -> Exp R                    -- current simulation time
    -> Exp Timestep             -- old timestep
    -> Exp Timestep
    -> Exp Timestep
    -> Exp Timestep
timeIncrement Parameters{..} t_now dt_old dt_courant dt_hydro =
  let
      dt_end    = t_end - t_now
      (lb,ub)   = dt_scaling

      -- try to prevent very small scaling on the next cycle
      target    = if dt_end >* step &&* dt_end <* 4 * step / 3
                     then 2 * step / 3
                     else dt_end

      -- increment the previous timestep by a small amount
      step      = min dt_new dt_max

      c1        = 1.0e20
      c2        = if dt_courant <* c1 then dt_courant / 2   else c1
      c3        = if dt_hydro   <* c2 then dt_hydro * 2 / 3 else c2

      ratio     = c3 / dt_old
      dt_new    = caseof ratio
                [ (\r -> r >=* 1 &&* r <* lb, dt_old)
                , (\r -> r >=* 1 &&* r >* ub, dt_old * ub)
                ]
                c3
  in
  min step target

