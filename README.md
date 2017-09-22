LULESH-Accelerate
=================

Implementation of the Livermore Unstructured Lagrangian Explicit Shock
Hydrodynamics (LULESH) mini-app using [Accelerate][accelerate].

[LULESH][lulesh] represents a typical hydrodynamics code such as [ALE3D][ale3d],
but is a highly simplified application, hard-coded to solve the Sedov blast
problem on an unstructured hexahedron mesh.

![What LULESH models](images/ale3d.gif)


Benchmarks
----------

Macbook Pro 10,1 (Mid 2012)

  - Intel Core i7-3720QM CPU @ 2.60GHz
  - NVIDIA GeForce GT 650M GPU @ 900MHz
  - macOS 10.11.6
  - 30 x 30 x 30 grid
  - 8 threads

| Implementation | Compiler | Time (s) | SLOC |
-----------------|----------|---------:|-----:|
| [LULESH-OMP.cc](reference/C/LULESH-OMP.cc) | Clang 5.0 | 15.1 | 2400 |
| [LULESH-OMP.cc](reference/C/LULESH-OMP.cc) | ICC 17.0.1 | 14.4 | 2400 |
| [CUDA](reference/CUDA/lulesh-kepler-singlegpu) | CUDA 8.0.46 | 8.82 | 3000 |
| [accelerate-llvm-native-1.1.0.0](http://hackage.haskell.org/package/accelerate-llvm-native-1.1.0.0) | GHC 8.2.1, LLVM 5.0 | 6.28 | 1200 |
| [accelerate-llvm-ptx-1.1.0.0](http://hackage.haskell.org/package/accelerate-llvm-ptx-1.1.0.0) | GHC 8.2.1, LLVM 5.0 | 8.35 | 1200 |


  [accelerate]:         https://github.com/AccelerateHS/accelerate
  [ale3d]:              https://wci.llnl.gov/simulation/computer-codes/ale3d
  [lulesh]:             https://codesign.llnl.gov/lulesh.php

