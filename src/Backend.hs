{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE CPP                 #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ViewPatterns        #-}

module Backend where

import Prelude                                          as P
import Control.Lens
import System.Console.GetOpt

import Data.Array.Accelerate                            as A
import qualified Data.Array.Accelerate.Interpreter      as Interp
#ifdef ACCELERATE_CUDA_BACKEND
import qualified Data.Array.Accelerate.CUDA             as CUDA
#endif
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
import qualified Data.Array.Accelerate.LLVM.Native      as CPU
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
import qualified Data.Array.Accelerate.LLVM.PTX         as PTX
#endif


-- | The set of backends available to execute the program.
--
data Backend = Interpreter
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
             | CPU
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
             | PTX
#endif
#ifdef ACCELERATE_CUDA_BACKEND
             | CUDA
#endif
  deriving (P.Eq, P.Bounded)

-- | The default backend to use
--
defaultBackend :: Backend
defaultBackend = maxBound


-- The choice of show instance is important because this will be used to
-- generate the command line flag.
--
instance Show Backend where
  show Interpreter      = "interpreter"
#ifdef ACCELERATE_CUDA_BACKEND
  show CUDA             = "cuda"
#endif
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
  show CPU              = "llvm-cpu"
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
  show PTX              = "llvm-gpu"
#endif

-- | Execute Accelerate expressions
--
run :: Arrays a => Backend -> Acc a -> a
run Interpreter = Interp.run
#ifdef ACCELERATE_CUDA_BACKEND
run CUDA        = CUDA.run
#endif
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
run CPU         = CPU.run
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
run PTX         = PTX.run
#endif


run1 :: (Arrays a, Arrays b) => Backend -> (Acc a -> Acc b) -> a -> b
run1 Interpreter = Interp.run1
#ifdef ACCELERATE_CUDA_BACKEND
run1 CUDA        = CUDA.run1
#endif
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
run1 CPU         = CPU.run1
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
run1 PTX         = PTX.run1
#endif


run2 :: (Arrays a, Arrays b, Arrays c) => Backend -> (Acc a -> Acc b -> Acc c) -> a -> b -> c
run2 backend f =
  let !go = run1 backend (A.uncurry f)
  in \x y -> go (x, y)

run3 :: forall a b c d. (Arrays a, Arrays b, Arrays c, Arrays d)
     => Backend
     -> (Acc a -> Acc b -> Acc c -> Acc d)
     -> a -> b -> c -> d
run3 backend f =
  let !go = run1 backend (uncurry3 f)
  in  \x y z -> go (x, y, z)

uncurry3
    :: Unlift f (f a, f b, f c)
    => (f a -> f b -> f c -> f d)
    -> f (Plain (f a), Plain (f b), Plain (f c))
    -> f d
uncurry3 f argv = let (a, b, c) = unlift argv in f a b c


-- | The set of available backnds. This will be used for both the command line
-- options as well as the fancy header generation.
--
availableBackends :: (ASetter' options Backend) -> [OptDescr (options -> options)]
availableBackends optBackend =
  [ Option  [] [show Interpreter]
            (NoArg (set optBackend Interpreter))
            "reference implementation (sequential)"

#ifdef ACCELERATE_CUDA_BACKEND
  , Option  [] [show CUDA]
            (NoArg (set optBackend CUDA))
            "implementation for NVIDIA GPUs (parallel)"
#endif
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
  , Option  [] [show CPU]
            (NoArg (set optBackend CPU))
            "LLVM based implementation for multicore CPUs (parallel)"
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
  , Option  [] [show PTX]
            (NoArg (set optBackend PTX))
            "LLVM based implementation for NVIDIA GPUs (parallel)"
#endif
  ]

