using FRGLatticeEvaluation

let
    L1 = LatticeInfo(Pyrochlore.getPyrochlore(12,[1,2,3]),Pyrochlore)

    coupl = - L1.System.couplings

    Fourier = getFourier(coupl,L1)
    tFFT = @elapsed FFT = getLatticeFFT(coupl,L1,64)
    tscan = @elapsed getkMax(FFT,res = 64,ext = 4pi)

    @info "Performance for LatticeFFT" tFFT tscan total = tFFT + tscan
    
    
    tConstruct = @elapsed FT = getnaiveLatticeFT(coupl,L1)
    
    tscan = @elapsed getkMax(FT,res = 64,ext = 4pi)

    @info "Performance for conventional FT" tConstruct tscan total = tConstruct + tscan

end