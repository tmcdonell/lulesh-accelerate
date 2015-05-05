{-# LANGUAGE RecordWildCards #-}

module Time where

import Domain


-- Compute the time increment for the current loop iteration. We aim for a
-- "target" timestep value which completes the simulation in the next step, but
-- is only allowed to change from the previous value by a certain amount, is is
-- subject to courant and hydro constraints.
--
timeIncrement :: Domain -> R
timeIncrement Domain{..}
  = step `min` target
  where
    -- try to prevent very small scaling on the next cycle
    target
      | dt_end > step && dt_end < 4 * step / 3 = 2 * step / 3
      | otherwise                              = dt_end
      where
        dt_end = time_end - time

    -- time increment can from the previous value by a small amount
    step
      | iteration == 0 = dt
      | otherwise      = dt' `min` dt_max
      where
        -- magic numbers from lulesh.cc:TimeIncrement()
        c1 = 1.0e20
        c2 = if dt_courant < c1
                then dt_courant / 2
                else c1
        c3 = if dt_hydro < c2
                then dt_hydro * 2 / 3
                else c2

        (lb,ub)       = dt_scaling
        r             = c3 / dt
        dt'
          | r >= 1 && r < lb  = dt      -- TLM: * lb ??
          | r >= 1 && r > ub  = dt * ub
          | otherwise         = c3

