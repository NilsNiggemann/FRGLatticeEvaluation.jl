# FRGLatticeEvaluation.jl: Fourier Transforms on lattices and helper functions for evaluating the FRG flow

## Usage
```
using FRGLatticeEvaluation
using SpinFRGLattices

L1 = LatticeInfo(Pyrochlore.getPyrochlore(7),Pyrochlore)
coupl = L1.System.couplings
Fourier = getFourier(coupl,L1)
FFT = interpolatedFT(coupl,L1,64)


```