{-# LANGUAGE BangPatterns #-}
--
-- Timing utilities
--
module Timing where

import Numeric
import Data.Time
import System.CPUTime
import Control.Applicative
import Prelude


-- Time -----------------------------------------------------------------------
-- | Abstract representation of process time.
--
data Time = Time
  { cpuTime     :: !Integer
  , wallTime    :: !UTCTime
  }
  deriving Show

data ElapsedTime = ElapsedTime !Integer !NominalDiffTime

instance Show ElapsedTime where
  show (ElapsedTime cpu wall) =
    "wall: " ++ showFFloatSIBase (Just 3) 1000 (fromRational (toRational wall) :: Double) "s, " ++
    "cpu: "  ++ showFFloatSIBase (Just 3) 1000 (fromIntegral cpu * 1.0e-12 :: Double) "s"

-- | Get the current time
--
getTime :: IO Time
getTime =
  Time <$> getCPUTime <*> getCurrentTime


-- | Difference of two times
--
elapsedTime :: Time -> Time -> ElapsedTime
elapsedTime (Time c1 w1) (Time c2 w2) = ElapsedTime (c1-c2) (diffUTCTime w1 w2)


-- Timing operations ----------------------------------------------------------

-- | Time some IO action.
--
--   Make sure to deepseq the result before returning it from the action. If you
--   don't do this then there's a good chance that you'll just pass a suspension
--   out of the action, and the computation time will be zero.
--
time :: IO a -> IO (a, ElapsedTime)
{-# NOINLINE time #-}
time it = do
  start <- getTime
  !r    <- it
  end   <- getTime
  return (r, elapsedTime end start)


-- | Show a signed 'RealFloat' value using SI unit prefixes. In the call to:
--
-- > showFFloatSIBase prec base val
--
-- If @prec@ is @'Nothing'@ the value is shown to full precision, and if @prec@
-- is @'Just' d@, then at most @d@ digits are shown after the decimal place.
-- Here @base@ represents the increment size between multiples of the original
-- unit. For measures in base-10 this will be 1000 and for values in base-2 this
-- is usually 1024, for example when measuring seconds versus bytes,
-- respectively.
--
showFFloatSIBase :: RealFloat a => Maybe Int -> a -> a -> ShowS
showFFloatSIBase prec !base !k
  = showString
  $ case pow of
      4   -> with "T"
      3   -> with "G"
      2   -> with "M"
      1   -> with "k"
      -1  -> with "m"
      -2  -> with "µ"
      -3  -> with "n"
      -4  -> with "p"
      _   -> showGFloat prec k " "      -- no unit or unhandled SI prefix
  where
    !k'         = k / (base ^^ pow)
    !pow        = floor (logBase base k) :: Int
    with unit   = showFFloat prec k' (' ':unit)

