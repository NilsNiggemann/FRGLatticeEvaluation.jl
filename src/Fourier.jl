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
getDim(L::AbstractLattice) = getDim(L.Basis)

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

struct FourierStruct3D{T<:Real}
    FourierInfos::FourierInfo{3}
    Sij_vec::Vector{T}
    Rx_vec::Vector{T}
    Ry_vec::Vector{T}
    Rz_vec::Vector{T}
    NCell::Int
end
struct FourierStruct2D{T<:Real}
    FourierInfos::FourierInfo{2}
    Sij_vec::Vector{T}
    Rx_vec::Vector{T}
    Ry_vec::Vector{T}
    NCell::Int
end

function splitRij(Rij::Vector{SVector{Dim,T}}) where {T<:Real,Dim}
    return [[r[i] for r in Rij] for i in 1:Dim]
end

FourierStruct(F::FourierInfo{2},Sij,Rx_vec,Ry_vec,NCell) = FourierStruct2D(F,Sij,Rx_vec,Ry_vec,NCell)
FourierStruct(F::FourierInfo{3},Sij,Rx_vec,Ry_vec,Rz_vec,NCell) = FourierStruct3D(F,Sij,Rx_vec,Ry_vec,Rz_vec,NCell)

function FourierStruct(S_pairs::AbstractVector,Lattice::AbstractLattice)
    Sij = S_pairs[Lattice.FourierInfos.pairs]
    Rij_vec = Lattice.FourierInfos.Rij_vec

    R = splitRij(Rij_vec)
    return FourierStruct(Lattice.FourierInfos,Sij,R...,Lattice.Basis.NCell)
end

@inline function (F::FourierStruct3D)(kx::T,ky::T,kz::T) where T <: Real
    Chi_k = zero(T)
    (;Rx_vec,Ry_vec,Rz_vec,Sij_vec,NCell) = F
    @turbo for i in eachindex(Rx_vec,Ry_vec,Rz_vec,Sij_vec)
        Rx = Rx_vec[i]
        Ry = Ry_vec[i]
        Rz = Rz_vec[i]
        Chi_k += cos(kx*Rx+ky*Ry+kz*Rz) * Sij_vec[i]
    end
    return Chi_k/NCell
end

(F::FourierStruct3D)(k::AbstractVector{T}) where T <: Real = F(k...)

(F::FourierStruct2D)(k::AbstractVector{T}) where T <: Real = F(k...)

@inline function (F::FourierStruct3D)(kx::T,ky::T) where T <: Real
    Chi_k = zero(T)
    (;Rx_vec,Ry_vec,Sij_vec,NCell) = F
    @turbo for i in eachindex(Rx_vec,Ry_vec,Sij_vec)
        Rx = Rx_vec[i]
        Ry = Ry_vec[i]
        Chi_k += cos(kx*Rx+ky*Ry) * Sij_vec[i]
    end
    return Chi_k/NCell
end

@inline FourierTransform(k,Chi_R, Lattice::AbstractLattice) = FourierTransform(k,Chi_R, Lattice.Basis.NCell,Lattice.FourierInfos.pairs,Lattice.FourierInfos.Rij_vec) 

@inline FourierTransform_prec(k,Chi_R, Lattice::AbstractLattice) = FourierTransform_prec(k,Chi_R, Lattice.Basis.NCell,Lattice.FourierInfos.pairs,Lattice.FourierInfos.Rij_vec) 

"""Constructs a function χ(q) for given real space chi and the lattice structure and returns it"""
function getFourier(Chi_R::AbstractArray,Lattice::AbstractLattice)
    FTStr = FourierStruct(Chi_R,Lattice)
    @inline FT(q::SVector) = FTStr(q)
    @inline FT(qx::Real,qy::Real) = FT(SA[qx,qy])
    @inline FT(qx::Real,qy::Real,qz::Real) = FT(SA[qx,qy,qz])

    @inline FT(q) = throw(ArgumentError("Please use a static vector or function arguments"))
    return FT
end

"""Returns 2D Fourier trafo in plane as specified by the "regionfunc" function. Eg for a plot in the xy plane we can use plane = (ki,kj) -> SA[ki,kj] """
function Fourier2D(Chi_R::AbstractArray,regionfunc::Function,Lattice::AbstractLattice;res=100,ext = pi,minext = -ext,kwargs...)
    karray = range(minext,stop = ext,length = res)
    Chi_k = Fourier2D(Chi_R,karray,karray ,regionfunc,Lattice;kwargs...)
    return karray,Chi_k
end

function Fourier2D(Chi_R::AbstractArray,x::AbstractVector,y::AbstractVector, regionfunc::Function,Lattice::AbstractLattice)
    Chi_k = zeros(length(x),length(y))
    FT = FourierStruct(Chi_R,Lattice)
    Threads.@threads for j in eachindex(y)
        kj = y[j]
        for (i,ki) in enumerate(x)
            Chi_k[i,j] = FT(regionfunc(ki,kj))
        end
    end
    return Chi_k
end

"""Returns 3D Fourier trafo"""
function Fourier3D(Chi_R::AbstractArray,Lattice::AbstractLattice;res=50,ext = pi,minext = -ext,kwargs...)
    k = range(minext,stop = ext,length = res)
    k, Fourier3D(Chi_R,Lattice,k,k,k;kwargs...)
end

function Fourier3D(Chi_R::AbstractArray,Lattice::AbstractLattice,kx_vec::AbstractVector,ky_vec::AbstractVector,kz_vec::AbstractVector)
    Chi_k = zeros(length(kx_vec),length(ky_vec),length(kz_vec))
    FT = FourierStruct(Chi_R,Lattice)
    Threads.@threads for iz in eachindex(kz_vec)
        kz = kz_vec[iz]
        Threads.@threads for iy in eachindex(ky_vec)
            ky = ky_vec[iy]
            for ix in eachindex(kx_vec)
                kx = kx_vec[ix]
                Chi_k[ix,iy,iz] = FT(kx,ky,kz)
            end
        end
    end
    return Chi_k
end

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

function getMaxFlow(Chi_LR::AbstractMatrix,Lambdas,Lattice::AbstractLattice;  res = 50,ext = 4pi,max=maximum(Lambdas))
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


function getkMax(Chi_R::AbstractVector,Lattice::AbstractLattice,ext = 4pi,res = 50;kwargs...) 
    FT = getFullFourier(Lattice;ext,res,kwargs...)
    k, Chik = FT(Chi_R)
    maxpos =  Tuple(argmax(Chik))
    kmax = SA[k[[maxpos...]]...]
end


function getkMaxOptim(Chiq,kguess::T;kwargs...) where T <: SVector
    kguess = convert(Vector,kguess)
    MinusChi(k) = -Chiq(convert(T,k))
    res = optimize(MinusChi,kguess;kwargs...)
    if !Optim.converged(res)
        @warn "optimization did not converge" res res.x_converged res.f_converged res.iterations res.stopped_by
    end
    k = res.minimizer

    return convert(T,k)
end

"""Given a vector of real space susceptibilities and a lattice, compute the wavevector of the maximum susceptibility by first optimizing over a coarse k-space grid of resolution res and then refining the result with a gradient descent"""
function getkMaxOptim(Chi_R::AbstractVector,Lattice::AbstractLattice;ext = 4pi,res = 50,kwargs...)
    kguess = getkMax(Chi_R,Lattice,ext,res)
    Chiq = getFourier(Chi_R,Lattice)
    return getkMaxOptim(Chiq,kguess;kwargs...)
end

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
CorrelationLength(Chi::Function, Q::AbstractVector, qa::AbstractVector) = 1 / norm(Q - qa) * sqrt(abs(Chi(Q) / Chi(qa) - 1))

function CorrelationLength(Chi::AbstractArray,Q::AbstractVector,direction::AbstractVector,Lattice::AbstractLattice) 
    ChiFunc = getFourier(Chi,Lattice)
    qa = Q+direction/norm(direction) * 2pi/Lattice.System.NLen
    CorrelationLength(ChiFunc,Q,qa)
end

function CorrelationLength(Chi::AbstractArray,direction::AbstractVector,Lattice::AbstractLattice;kwargs...) 
    Q = getkMax(Chi,Lattice;kwargs...)
    CorrelationLength(Chi,Q,direction,Lattice)
end

"""Returns the maximum correlation length by averaging over a sphere in reciprocal space.
Example: maxCorrelationLength(SA[π,π,π],Chi,NLen,res=50) averages over all angles θ, and ϕ such that the interval [0,2π] contains 50 points.
"""
function maxCorrelationLength(Q::SVector, chi::Function, NLen; kwargs...)
    dirs = unitsphereGenerator(Q;kwargs...)
    qa(d) = Q + d / norm(d) * 2π / NLen
    return maximum(CorrelationLength(chi, Q, qa(d)) for d in dirs)
end

sphere(r,θ,ϕ) = SA[r*sin(θ)*cos(ϕ), r*sin(θ)*sin(ϕ), r*cos(θ)]
sphere(r,ϕ) = SA[r*cos(ϕ), r*sin(θ)]

function unitsphereGenerator(::SVector{3,<:Real}; res=50)
    theta = range(0, π, res ÷ 2)
    phi = range(0, 2π, res)
    return (sphere(1,t,p) for t in theta for p in phi)
end

function unitsphereGenerator(::SVector{2,<:Real}; res=50)
    phi = range(0, 2π, res)
    return (sphere(1,p) for p in phi)
end

maxCorrelationLength(Q::SVector, chi::Function; kwargs...) = maxCorrelationLength(Q, chi, chi.Lattice.NLen;kwargs...)
