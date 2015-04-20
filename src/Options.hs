
module Options where

data Options = Options
  {
    optMaxSteps         :: Int                  -- maximum number of cycles to run
  , optSize             :: Int                  -- length of cube mesh along each dimension
  , optNumRegions       :: Int                  -- number of distinct regions
  , optBalance          :: Int                  -- load balance between regions of a domain
  , optCost             :: Int                  -- extra cost of more expensive regions
  , optOutputFile       :: Maybe FilePath       -- output viz files
  , optProgress         :: Bool                 -- print out progress
  , optHelp             :: Bool                 -- print help message
  }
  deriving (Eq, Show)

defaultOpts :: Options
defaultOpts = Options
  {
    optMaxSteps         = 9999999
  , optSize             = 30
  , optNumRegions       = 11
  , optBalance          = 1
  , optCost             = 1
  , optOutputFile       = Nothing
  , optProgress         = False
  , optHelp             = False
  }

parseOptions :: [String] -> Options
parseOptions _ = defaultOpts

