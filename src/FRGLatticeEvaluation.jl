module FRGLatticeEvaluation
    
    using SpinFRGLattices,StaticArrays,LaTeXStrings

    using SLEEFPirates:cos_fast

    include("Fourier.jl")
    export AbstractLattice, LatticeInfo,FourierTransform, Fourier2D,Fourier3D,  getFlow, getMaxFlow, pointPath, fetchKPath,getkMax,CorrelationLength
        
    include("helpers.jl")
    export hhlplane,normalTohhlplane,xyplane,zzerocut,sphereplane,hhllabels,xylabels,zzerolabels,normalTohhllabels,strd,strd3,getPairs,DimerResponse,DimerResponse_Keles,getFourier


end # module FRGLatticeEvaluation
