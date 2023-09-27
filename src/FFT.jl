"""Given a lattice geometry, prepare a list of all possible pairs of sites and the corresponding inequivalent pair"""
function getCorrelationPairs(UnitCell::AbstractVector{R},SiteList::AbstractVector{R},PairList::AbstractVector{R},PairTypes,pairToInequiv,Basis) where {R <: Rvec}
    Ri_vec = empty(PairList)
    Rj_vec = empty(PairList)

    pairs = Int[]
    for i_site in UnitCell
        # println(Ri)
        for j_site in SiteList # site summation
            R_Ref,ij = pairToInequiv(i_site,j_site) #Map j to correct pair so that we may use Chi_0,j'
            xi = getSiteType(R_Ref,Basis)
            pair = MapToPair(xi,ij,PairList,PairTypes)
            if pair !== 0
                push!(Ri_vec,i_site)
                push!(Rj_vec,j_site)
                push!(pairs,pair)
            end
        end
    end
    return (;Ri_vec,Rj_vec,pairs)
end

isInUnitCell(R::Rvec_2D) = R.n1 == 0 && R.n2 == 0
isInUnitCell(R::Rvec_3D) = R.n1 == 0 && R.n2 == 0 && R.n3 == 0

function getCorrelationPairs(pairNumberDict::AbstractDict{Tuple{R,R},I}) where {I<:Integer,R <: Rvec}
    UCPair(R1,R2,pair) = isInUnitCell(R1)
    # UCPair(x) = UCPair(first(first(x)),second(first(x)),second(x))
    Ri_vec = R[]
    Rj_vec = R[]
    pairs = Int[]
    for ((R1,R2),pair) in pairNumberDict
        if UCPair(R1,R2,pair)
            push!(Ri_vec,R1)
            push!(Rj_vec,R2)
            push!(pairs,pair)
        end
    end

    return (;Ri_vec,Rj_vec,pairs)
end

"""Given lists of paired sites Rk and Rj, and their susceptibility Chi_ij, separate the susceptibility into a Matrix containing the sublattices to be used for FFT.
"""
function separateSublattices(Ri_vec::AbstractVector{Rvec_2D},Rj_vec::AbstractVector{Rvec_2D},Chi_ij::AbstractVector{T}) where {T}
    NCell = length(unique(x.b for x in Rj_vec))
    
    maxn(R::Rvec_2D) = maximum(abs,(R.n1,R.n2))

    L = maximum(maxn,Rj_vec)
    Lbox = 2L+1

    offset = L+1
    coords(R::Rvec_2D) = CartesianIndex(offset .+ (R.n1, R.n2))
    function chi_ij(α,β) 

        chi = zeros(T,Lbox,Lbox)
        for (i,χ) in enumerate(Chi_ij)
            if Ri_vec[i].b == α && Rj_vec[i].b == β
                chi[coords(Rj_vec[i])] = χ
            end
        end
        return chi
    end
    return [chi_ij(α,β) for α in 1:NCell, β in 1:NCell]
end

function separateSublattices(Ri_vec::AbstractVector{Rvec_3D},Rj_vec::AbstractVector{Rvec_3D},Chi_ij::AbstractArray{T}) where {T}
    NCell = length(unique(x.b for x in Rj_vec))
    
    maxn(R::Rvec_3D) = maximum(abs,(R.n1,R.n2,R.n3))

    L = maximum(maxn,Rj_vec)
    Lbox = 2L+1

    offset = L+1
    coords(R::Rvec_3D) = CartesianIndex(offset .+ (R.n1, R.n2, R.n3))

    function chi_ij(α,β) 

        chi = zeros(T,Lbox,Lbox,Lbox)
        for (i,χ) in enumerate(Chi_ij)
            if Ri_vec[i].b == α && Rj_vec[i].b == β
                chi[coords(Rj_vec[i])] = χ
            end
        end
        return chi
    end
    return [chi_ij(α,β) for α in 1:NCell, β in 1:NCell]
end


getCorrelationPairs(Lat::AbstractLattice) = getCorrelationPairs(Lat.UnitCell,Lat.SiteList,Lat.PairList,Lat.PairTypes,Lat.pairToInequiv,Lat.Basis)

getCorrelationPairs(Lattice::NamedTuple) = getCorrelationPairs(Lattice.pairNumberDict)
getCorrelationPairs(Lattice::ReducedLattice) = getCorrelationPairs(Lattice.pairNumberDict)

import LatticeFFTs.getLatticeFFT

LatticeFFTs.getLatticeFFT(S_ab::AbstractMatrix,Basis::Basis_Struct_3D,args...;kwargs...)  = getLatticeFFT(S_ab,[Basis.a1 Basis.a2 Basis.a3],Basis.b,args...;kwargs...)
LatticeFFTs.getLatticeFFT(S_ab::AbstractMatrix,Basis::Basis_Struct_2D,args...;kwargs...) = getLatticeFFT(S_ab,[Basis.a1 Basis.a2],Basis.b,args...;kwargs...)

getNaiveLatticeFT(S_ab::AbstractMatrix,Basis::Basis_Struct_3D,args...;kwargs...)  = naiveLatticeFT(S_ab,[Basis.a1 Basis.a2 Basis.a3],Basis.b,args...;kwargs...)
getNaiveLatticeFT(S_ab::AbstractMatrix,Basis::Basis_Struct_2D,args...;kwargs...) = naiveLatticeFT(S_ab,[Basis.a1 Basis.a2],Basis.b,args...;kwargs...)

getNPairs(redLat::Union{ReducedLattice,NamedTuple}) = length(unique(values(redLat.pairNumberDict)))

function setCouplings!(S_ab::AbstractMatrix{<:AbstractArray{T}},couplings,coupMap) where T
    for (Sab,coup) in zip(S_ab,coupMap)
        for (i,m) in enumerate(coup)
            if m !== 0
                Sab[i] = couplings[m]
            else
                Sab[i] = zero(T)
            end
        end
    end
    return S_ab
end

function getLatticeFFT(ChiR::AbstractVector,Lattice,args...;kwargs...)
    CorrelationPairs = getCorrelationPairs(Lattice)
    (;Ri_vec,Rj_vec,pairs) = CorrelationPairs
    S_ab = separateSublattices(Ri_vec,Rj_vec,ChiR[pairs])

    return getLatticeFFT(S_ab,Lattice.Basis,args...;kwargs...)
end

struct LatticeFTplan{T<:LatticeFT}
    FT::T
end

function planNaiveLatticeFT(Lattice,args...;kwargs...)
    CorrelationPairs = getCorrelationPairs(Lattice)
    (;Ri_vec,Rj_vec,pairs) = CorrelationPairs
    NPairs = getNPairs(Lattice)
    RealSpaceInds = 1:NPairs
    S_ab_inds = separateSublattices(Ri_vec,Rj_vec,RealSpaceInds[pairs])
    FT = getNaiveLatticeFT(S_ab_inds,Lattice.Basis,args...;kwargs...)
    plan = LatticeFTplan(FT)

    reducePlan!(plan)
    return plan
end


function reducePlan!(plan::LatticeFTplan)
    S_ab = plan.FT.S
    for S in S_ab
        delinds = findall(iszero,S.Sij)
        deleteat!(S.Sij,delinds)
        for rij_x in S.Rij
            deleteat!(rij_x,delinds)
        end
    end
    return plan
end

function convertTypes(T::Type,plan::LatticeFTplan)
    LatticeFTplan(convertTypes(T,plan.FT))
end

function convertTypes(T::Type,FT::LatticeFFTs.LatticeFT)
    LatticeFFTs.LatticeFT(convertTypes.(T,FT.S))
end

function convertTypes(T::Type,FT::LatticeFFTs.naiveSubLatticeFT)
    Sij = convertTypes(T,FT.Sij)
    Rij = convertTypes.(T,FT.Rij)
    LatticeFFTs.naiveSubLatticeFT(Sij,Rij)
end

convertTypes(T::Type,x::AbstractArray) = convert.(T,x)
convertTypes(T::Type,x) = convert(T,x)

function getNaiveLatticeFT(ChiR::AbstractVector,plan::LatticeFTplan,args...;kwargs...)
    Sq = convertTypes(eltype(ChiR),plan).FT
    setCouplings!(Sq,ChiR,plan)
    return Sq
end

getNaiveLatticeFT(ChiR::AbstractVector,Lattice,args...;kwargs...) = getNaiveLatticeFT(ChiR,planNaiveLatticeFT(Lattice,args...;kwargs...))

function setCouplings!(Sq::LatticeFFTs.LatticeFT,couplings::AbstractVector{T},plan) where T
    for (Sab,coup) in zip(Sq.S,plan.FT.S)
        for (i,m) in enumerate(coup.Sij)
            if m !== 0
                Sab.Sij[i] = couplings[m]
            else
                Sab.Sij[i] = zero(T)
            end
        end
    end
    return Sq
end

getDim(ChikFunction::PhaseShiftedFFT) = length(ChikFunction.PhaseVector)
getDim(::naiveSubLatticeFT{N,V,T}) where {N,V,T} = N

getDim(ChikFunction::LatticeFT) = getDim(ChikFunction[1,1])

function getkMax(ChikFunction::AbstractLatticeFT;res = 120,ext = 4pi,kwargs...) 
    dim = getDim(ChikFunction)
    kRange = Iterators.product((range(start = -ext,stop = ext,length = res) for _ in 1:dim)...)
    
    return argmax(ChikFunction,SVector(k) for k in kRange)
end

# function getkMax(ChikFunction::AbstractLatticeFT;res = 60,ext = 4pi,kwargs...) 
#     dim = getDim(ChikFunction)
#     kRange_1 = range(start = -ext,stop = ext,length = res)
#     kRange = collect(Iterators.product((kRange_1 for _ in 1:dim)...))

#     Chik = zeros((res for _ in 1:dim)...)

#     # for (ik,k) in zip(CartesianIndices(Chik),kRange)
#     Threads.@threads for (ik) in CartesianIndices(Chik)
#         k = kRange[ik]
#         Chik[ik] = ChikFunction(k...)
#     end
#     maxpos =  Tuple(argmax(Chik))
#     kmax = SVector{dim}(kRange_1[[maxpos...]])
# end