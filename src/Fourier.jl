struct FourierInfo{Dim}
    pairs::Vector{Int}
    Rij_vec::Vector{SVector{Dim,Float64}}
end
abstract type AbstractLattice end
Base.@kwdef struct LatticeInfo{BasisType,RvecType,FunctionType,Dim} <: AbstractLattice
    System::Geometry
    Basis::BasisType
    NLen::Int = System.NLen
    Npairs::Int = System.Npairs
    NUnique::Int = System.NUnique
    PairList::Vector{RvecType} = System.PairList
    PairTypes::Vector{sitePair} = System.PairTypes
    SiteList::Vector{RvecType} = unique(SpinFRGLattices.sortedPairList(NLen,Basis)[1])
    UnitCell::Vector{RvecType} = [SpinFRGLattices.getRvec(b,Basis) for b in Basis.b]
    pairToInequiv::FunctionType
    FourierInfos::FourierInfo{Dim} = PrecomputeFourier(UnitCell,SiteList,PairList,PairTypes,pairToInequiv,Basis)
end
LatticeInfo(System::Geometry,Basis::Basis_Struct,pairToInequiv::Function ) = LatticeInfo(;System,Basis,pairToInequiv)
LatticeInfo(System::Geometry,Mod::Module,PTI =Mod.pairToInequiv ) = LatticeInfo(System=System,Basis = Mod.Basis,pairToInequiv = PTI)

getDim(::Basis_Struct_2D) = 2
getDim(::Basis_Struct_3D) = 3
getDim(L::LatticeInfo) = getDim(L.Basis)

function PrecomputeFourier(UnitCell::AbstractVector{T},SiteList::AbstractVector{T},PairList::AbstractVector{T},PairTypes,pairToInequiv,Basis) where T <: Rvec
    Rij_vec = SVector{getDim(Basis),Float64}[]
    pairs = Int[]
    for i_site in UnitCell
        Ri = getCartesian(i_site,Basis)
        # println(Ri)
        for j_site in SiteList # site summation
            Rj = getCartesian(j_site,Basis)
            Rij = Ri - Rj
            R_Ref,ij = pairToInequiv(i_site,j_site) #Map j to correct pair so that we may use Chi_0,j'
            xi = getSiteType(R_Ref,Basis)
            pair = MapToPair(xi,ij,PairList,PairTypes)
            if pair !== 0
                push!(Rij_vec,Rij)
                push!(pairs,pair)
            end
        end
    end
    return FourierInfo(pairs,Rij_vec)
end

@inline function dot(v1::AbstractVector{T},v2::AbstractVector{T}) where T
    res = zero(T)
    @inbounds @simd for i in eachindex(v1,v2)
        res+= v1[i]*v2[i]
    end
    return res
end

@inline function FourierTransform_prec(k::AbstractVector,Chi_R, NCell,pairs,Rij_vec)
    Chi_k = 0.
    # @inbounds @simd for (p,Rij) in zip(pairs,Rij_vec)
    @inbounds @simd for i in eachindex(pairs,Rij_vec)
        p = pairs[i]
        Rij = Rij_vec[i]
        Chi_k += cos(dot(k,Rij)) * Chi_R[p]
    end
    return 1/NCell * Chi_k
end

@inline function FourierTransform(k::AbstractVector,Chi_R, NCell,pairs,Rij_vec)
    Chi_k = 0.
    # @inbounds @simd for (p,Rij) in zip(pairs,Rij_vec)
    @inbounds @simd for i in eachindex(pairs,Rij_vec)
        p = pairs[i]
        Rij = Rij_vec[i]
        
        Chi_k += cos_fast(dot(k,Rij)) * Chi_R[p]
    end
    return 1/NCell * Chi_k
end


@inline FourierTransform(k,Chi_R, Lattice::AbstractLattice) = FourierTransform(k,Chi_R, Lattice.Basis.NCell,Lattice.FourierInfos.pairs,Lattice.FourierInfos.Rij_vec) 

@inline FourierTransform_prec(k,Chi_R, Lattice::AbstractLattice) = FourierTransform_prec(k,Chi_R, Lattice.Basis.NCell,Lattice.FourierInfos.pairs,Lattice.FourierInfos.Rij_vec) 

"""Constructs a function χ(q) for given real space chi and the lattice structure and returns it"""
function getFourier(Chi_R::AbstractArray,Lattice::AbstractLattice)
    @inline FT(q::StaticVector) = FourierTransform(Float64.(q),Chi_R,Lattice.Basis.NCell,Lattice.FourierInfos.pairs,Lattice.FourierInfos.Rij_vec)
    @inline FT(args::Tuple) = FT(StaticVector(args))
    @inline FT(args...) = FT(SA[args...])
    @inline FT(q) = throw(ArgumentError("Please use a static vector or a tuple"))
    return FT
end

"""Returns 2D Fourier trafo in plane as specified by the "regionfunc" function. Eg for a plot in the xy plane we can use plane = (ki,kj) -> SA[ki,kj] """
function Fourier2D(Chi_R::AbstractArray,regionfunc::Function,Lattice::AbstractLattice;res=100,ext = pi,minext = -ext,kwargs...)
    karray = range(minext,stop = ext,length = res)
    Chi_k = Fourier2D(Chi_R,karray,karray ,regionfunc,Lattice;kwargs...)
    return karray,Chi_k
end

function Fourier2D(Chi_R::AbstractArray,x::AbstractVector,y::AbstractVector, regionfunc::Function,Lattice::AbstractLattice;fast = true)
    Chi_k = zeros(length(x),length(y))
    FT = fast ? FourierTransform : FourierTransform_prec
    Threads.@threads for j in eachindex(y)
        kj = y[j]
        for (i,ki) in enumerate(x)
            Chi_k[i,j] = FT(regionfunc(ki,kj),Chi_R,Lattice)
            # isnan(Chi_k[i,j]) && println((i,j,ki,kj))
        end
    end
    return Chi_k
end

"""Returns 3D Fourier trafo"""
function Fourier3D(Chi_R::AbstractArray,Lattice::AbstractLattice;res=50,ext = pi,minext = -ext,kwargs...)
    k = range(minext,stop = ext,length = res)
    k, Fourier3D(Chi_R,Lattice,k,k,k;kwargs...)
end

function Fourier3D(Chi_R::AbstractArray,Lattice::AbstractLattice,kx_vec::AbstractVector,ky_vec::AbstractVector,kz_vec::AbstractVector;fast = true)
    Chi_k = zeros(length(kx_vec),length(ky_vec),length(kz_vec))
    FT = fast ? FourierTransform : FourierTransform_prec
    Threads.@threads for iz in eachindex(kz_vec)
        kz = kz_vec[iz]
        for (iy,ky) in enumerate(ky_vec),(ix,kx) in enumerate(kx_vec)
            Chi_k[ix,iy,iz] = FT(SA[kx,ky,kz],Chi_R,Lattice)
        end
    end
    return Chi_k
end

##


##
function getFlow(k::StaticArray,Chi_LR,Lambdas,Lattice)
    flow = similar(Lambdas)
    @inline FT(k,Chi) = FourierTransform(k,Chi,Lattice)
    for i in eachindex(Lambdas,flow)
        flow[i] =  @views FT(k,Chi_LR[i,:])
    end
    return flow
end

function getFullFourier(Lattice::LatticeInfo{Basis_Struct_2D,Rvec_2D,F,2};kwargs...) where {F} 
    @inline FT(Chi) = Fourier2D(Chi,xyplane,Lattice;kwargs...)
end
function getFullFourier(Lattice::LatticeInfo{Basis_Struct_3D,Rvec_3D,F,3};kwargs...) where {F} 
    @inline FT(Chi) = Fourier3D(Chi,Lattice;kwargs...)
end

function getMaxFlow(Chi_LR::AbstractMatrix,Lambdas,Lattice::LatticeInfo;  res = 50,ext = 4pi,max=maximum(Lambdas))
    index = findfirst(<=(max),Lambdas) #Todo: make this work when lambda is inversely sorted
    Lamnew = Lambdas[index:end]
    flow = similar(Lamnew)
    Chislice = Chi_LR[index:end,:]'
    FT(Chi) = maximum(getFullFourier(Lattice,res = res,ext=ext)(Chi)[2])
    for i in eachindex(Lamnew)
        flow[i] = @views FT(Chislice[:,i])
    end
    return Lamnew,flow
end

function getkMax(Chi_R::AbstractVector,Lattice::LatticeInfo,ext = 4pi,res = 50;kwargs...) 
    FT = getFullFourier(Lattice;ext,res,kwargs...)
    k, Chik = FT(Chi_R)
    maxpos =  Tuple(argmax(Chik))
    kmax = SA[k[[maxpos...]]...]
end


##
function pointPath(p1::StaticArray,p2::StaticArray,res)
    Path = Vector{typeof(p1)}(undef,res)
    for i in eachindex(Path)
        Path[i] = p1 + i/res*(p2 -p1)
    end
    return Path
end
"""res contains the number of points along -pi,pi"""
function fetchKPath(points,res = 100)
    Path = Vector{typeof(points[begin])}(undef,0)
    # Path = []
    PointIndices = [1]
    for i in eachindex(points[begin:end-1])
        p1 = points[i]
        p2 = points[i+1]
        append!(Path,pointPath(p1,p2,round(Int,norm(p1-p2)/2pi * res)))
        append!(PointIndices,length(Path)) # get indices corresponding to points
    end
    return PointIndices,Path
end
##

##
"""Compute correlation length according to Sandvik's definition
A. W. Sandvik, AIP Conf. Proc. 1297, 135 (2010).
"""
CorrelationLength(Chi::Function,Q::AbstractVector,qa::AbstractVector) = 1/norm(Q-qa) * sqrt(Chi(Q)/Chi(qa)-1)

function CorrelationLength(Chi::AbstractArray,Q::AbstractVector,direction::AbstractVector,Lattice::LatticeInfo) 
    ChiFunc = getFourier(Chi,Lattice)
    qa = Q+direction/norm(direction) * 2pi/Lattice.System.NLen
    CorrelationLength(ChiFunc,Q,qa)
end

function CorrelationLength(Chi::AbstractArray,direction::AbstractVector,Lattice::LatticeInfo;kwargs...) 
    Q = getkMax(Chi,Lattice;kwargs...)
    CorrelationLength(Chi,Q,direction,Lattice)
end
