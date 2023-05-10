module FRGLatticeEvaluation
    
    using SpinFRGLattices,StaticArrays,LaTeXStrings,Optim,LinearAlgebra,LoopVectorization

    using SLEEFPirates
    include("Fourier.jl")
    export AbstractLattice, LatticeInfo,FourierTransform, FourierStruct2D, FourierStruct3D, FourierStruct, Fourier2D,Fourier3D,  getFlow, getMaxFlow, pointPath, fetchKPath,getkMax,CorrelationLength,maxCorrelationLength,getkMaxOptim
        
    include("helpers.jl")
    export hhlplane,normalTohhlplane,xyplane,zzerocut,sphereplane,hhllabels,xylabels,zzerolabels,normalTohhllabels,strd,strd3,getPairs,DimerResponse,DimerResponse_Keles,getFourier


end # module FRGLatticeEvaluation
