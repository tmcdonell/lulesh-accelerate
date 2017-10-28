{-# LANGUAGE CPP             #-}
{-# LANGUAGE TemplateHaskell #-}

module Options
  where

import Backend

import Data.List                                ( intercalate )
import Control.Lens
import System.Exit
import System.Environment
import System.Console.GetOpt
import Text.PrettyPrint.ANSI.Leijen

#if !MIN_VERSION_accelerate(1,2,0)
import Data.Array.Accelerate.Debug
endif


data Options = Options
  {
    _optBackend         :: Backend              -- backend to execute with
  , _optMaxSteps        :: Int                  -- maximum number of cycles to run
  , _optSize            :: Int                  -- length of cube mesh along each dimension
--  , _optNumRegions      :: Int                  -- number of distinct regions
--  , _optBalance         :: Int                  -- load balance between regions of a domain
--  , _optCost            :: Int                  -- extra cost of more expensive regions
  , _optOutputPath      :: Maybe FilePath       -- output simulation state suitable for visualisation
  , _optProgress        :: Bool                 -- print out progress
  , _optHelp            :: Bool                 -- print help message
  }
  deriving (Eq, Show)

$(makeLenses ''Options)


-- The size of the simulation box
--
_WIDTH, _HEIGHT, _DEPTH :: Fractional a => a
_WIDTH  = 1.125
_HEIGHT = 1.125
_DEPTH  = 1.125


-- Default simulation options
--
defaultOptions :: Options
defaultOptions = Options
  {
    _optBackend         = defaultBackend
  , _optMaxSteps        = 9999999
  , _optSize            = 30
--  , _optNumRegions      = 11
--  , _optBalance         = 1
--  , _optCost            = 1
  , _optOutputPath      = Nothing
  , _optProgress        = False
  , _optHelp            = False
  }


options :: [OptDescr (Options -> Options)]
options =
  availableBackends optBackend ++
  [
    Option  "n" ["size"]
            (ReqArg (set optSize . read) "INT")
            "problem size in each dimension"

  , Option  []  ["max-steps"]
            (ReqArg (set optMaxSteps . read) "INT")
            "maximum number of iterations to run"

#ifdef ACCELERATE_VISIT
  , Option  []  ["output"]
            (ReqArg (set optOutputPath . Just) "FILE")
            "base name to output VisIt files"
#endif

  , Option  []  ["progress"]
            (NoArg (set optProgress True))
            "print simulation progress"

  , Option  "h?" ["help"]
            (NoArg (set optHelp True))
            "show help message"
  ]


-- | Format a (console) string as bold text. Assume the user has configured
-- their terminal colours to something that looks good (and avoids the light vs.
-- dark background debate).
--
title :: String -> String
title = show . bold . text


-- | Generate the list of available (and the selected) Accelerate backends.
--
fancyHeader :: Options -> [String] -> [String] -> String
fancyHeader opts hdr ftr = intercalate "\n" (hdr ++ body ++ ftr)
  where
    active this         = if this == show (view optBackend opts) then "*" else ""
    (ss,bs,ds)          = unzip3 $ map (\(b,d) -> (active b, b, d)) $ concatMap extract (availableBackends optBackend)
    table               = zipWith3 paste (sameLen ss) (sameLen bs) ds
    paste x y z         = "  " ++ x ++ "  " ++ y ++ "  " ++ z
    sameLen xs          = flushLeft ((maximum . map length) xs) xs
    flushLeft n xs      = [ take n (x ++ repeat ' ') | x <- xs ]
    --
    extract (Option _ los _ descr) =
      let losFmt  = intercalate ", " los
      in  case lines descr of
            []          -> [(losFmt, "")]
            (x:xs)      -> (losFmt, x) : [ ("",x') | x' <- xs ]
    --
    body                = title "Available backends:"
                        : table


-- | Parse the command line arguments and return the resulting options
-- structure, together with any unrecognised or unprocessed options.
--
-- Anything following a "--" is not processed.
--
parseArgs :: IO (Options, [String])
parseArgs = do
#if !MIN_VERSION_accelerate(1,2,0)
  accInit
#endif
  args <- getArgs

  let
      (argv, rest) =
        let (x,y) = span (/= "--") args
        in  (x, dropWhile (== "--") y)

      helpMsg []  = helpMsg'
      helpMsg err = unlines [ concat err, helpMsg' ]

      helpMsg'    = unlines
        [ fancyHeader defaultOptions header []
        , ""
        , usageInfo (title "Options:") options
        ]

      (opts, nonopts) = case getOpt' Permute options argv of
          (o,n,u,[])  -> (foldr id defaultOptions o, n ++ u)
          (_,_,_,err) -> error (helpMsg err)

  if view optHelp opts
     then putStr   (helpMsg []) >> exitSuccess
     else putStrLn (fancyHeader opts header footer)

  return (opts, nonopts ++ rest)


header, footer :: [String]
header =
  [ "accelerate-lulesh (c) [2015] The Accelerate Team"
  , ""
  , "Usage: accelerate-lulesh [OPTIONS]"
  , ""
  ]

footer =
  [ "" ]

