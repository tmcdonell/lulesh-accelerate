
module Options where

data Options = Options
  {
    maxSteps    :: Int                  -- maximum number of cycles to run
  , size        :: Int                  -- length of cube mesh along each dimension
  , numRegions  :: Int                  -- number of distinct regions
  , balance     :: Int                  -- load balance between regions of a domain
  , cost        :: Int                  -- extra cost of more expensive regions
  , outputFile  :: Maybe FilePath       -- output viz files
  , progress    :: Bool                 -- print out progress
  , help        :: Bool                 -- print help message
  }
  deriving (Eq, Show)

defaultOpts :: Options
defaultOpts = Options
  {
    maxSteps    = 9999999
  , size        = 30
  , numRegions  = 11
  , balance     = 1
  , cost        = 1
  , outputFile  = Nothing
  , progress    = False
  , help        = False
  }

parseOptions :: [String] -> Options
parseOptions _ = defaultOpts

