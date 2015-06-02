{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE RebindableSyntax    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ViewPatterns        #-}

module Util where

import Type

import Prelude                                          as P hiding ( (<*) )
import Data.Array.Accelerate                            as A hiding ( transpose )
import Data.Array.Accelerate.Control.Lens               as L hiding ( _1, _2, _3, _4, _5, _6, _7, _8, _9, at, ix, use )


-- Use rebindable syntax to improve readability of branches in scalar code.
--
ifThenElse :: Elt t => Exp Bool -> Exp t -> Exp t -> Exp t
ifThenElse x y z
  = x ? (y, z)


-- | Get the values at the nodes surrounding the element of a given index. This
-- follows the numbering convention:
--
-- <<../images/nodes.png>>
--
collectToElem
    :: Elt a
    => Acc (Field a)
    -> Exp Ix
    -> Exp (Hexahedron a)
collectToElem mesh ix@(unlift -> Z :. z :. y :. x)
  = let
        n0 = mesh ! ix
        n1 = mesh ! index3 z     y     (x+1)
        n2 = mesh ! index3 z     (y+1) (x+1)
        n3 = mesh ! index3 z     (y+1) x
        n4 = mesh ! index3 (z+1) y      x
        n5 = mesh ! index3 (z+1) y     (x+1)
        n6 = mesh ! index3 (z+1) (y+1) (x+1)
        n7 = mesh ! index3 (z+1) (y+1) x
    in
    lift (n0, n1, n2, n3, n4, n5, n6, n7)


-- Given a set of values at the nodes surrounding each element, form a nodal
-- grid by combining all those elements.
--
-- This is used to transform a local scatter operation (from elements to the
-- surrounding nodes) into a global gather operation.
--
distributeToNode
    :: forall a. Elt a
    => (Exp a -> Exp a -> Exp a)
    -> Exp a
    -> Acc (Field (Hexahedron a))
    -> Acc (Field a)
distributeToNode f zero arr =
  let
      numElem           = indexHead (shape arr)
      numNode           = numElem + 1
      sh'               = index3 numNode numNode numNode

      at :: Exp Int -> Exp Int -> Exp Int
         -> Lens' (Exp (Hexahedron a)) (Exp a)
         -> Exp a
      at z y x node =
        if 0 <=* z &&* z <* numElem &&* 0 <=* y &&* y <* numElem &&* 0 <=* x &&* x <* numElem
           then (arr ! index3 z y x) ^. node
           else zero

      mesh :: Exp Ix -> Exp a
      mesh (unlift -> Z :. z :. y :. x) =
        at z     y     x     _0 `f`
        at z     y     (x-1) _1 `f`
        at z     (y-1) (x-1) _2 `f`
        at z     (y-1) x     _3 `f`
        at (z-1) y     x     _4 `f`
        at (z-1) y     (x-1) _5 `f`
        at (z-1) (y-1) (x-1) _6 `f`
        at (z-1) (y-1) x     _7
  in
  generate sh' mesh


-- | Extract one of the six "faces" of a hexahedron.
--
-- The numbering of the faces is arbitrary, but the direction the corners trace
-- out is such that the face normal is towards the inside of the hexahedron.
--
collectFace :: Elt a => Int -> Exp (Hexahedron a) -> Exp (Quad a)
collectFace 0 p = lift (p^._0, p^._1, p^._2, p^._3)
collectFace 1 p = lift (p^._0, p^._4, p^._5, p^._1)
collectFace 2 p = lift (p^._1, p^._5, p^._6, p^._2)
collectFace 3 p = lift (p^._2, p^._6, p^._7, p^._3)
collectFace 4 p = lift (p^._3, p^._7, p^._4, p^._0)
collectFace 5 p = lift (p^._4, p^._7, p^._6, p^._5)
collectFace _ _ = error "collectFace: there are only six faces on a hexahedron"

