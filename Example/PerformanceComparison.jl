using FRGLatticeEvaluation

let
    L1 = LatticeInfo(Pyrochlore.getPyrochlore(16,[1,2,3]),Pyrochlore)

    coupl = - L1.System.couplings

    Fourier = getFourier(coupl,L1)
    tFFT = @elapsed FFT = getLatticeFFT(coupl,L1,64)
    tscanFFT = @elapsed getkMax(FFT,res = 64,ext = 4pi)
    tFFT_total = tFFT + tscanFFT

    @info "Performance for LatticeFFT" tFFT tscanFFT tFFT_total
    
    
    tConstruct = @elapsed FT = getnaiveLatticeFT(coupl,L1)
    
    tscan = @elapsed getkMax(FT,res = 64,ext = 4pi)

    tNaive = tConstruct + tscan

    @info "Speedup for FFT" tNaive speedup = tNaive/ tFFT_total

end