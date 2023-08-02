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
function separateSublattices(Ri_vec::AbstractVector{RType},Rj_vec::AbstractVector{RType},Chi_ij) where RType <: Rvec
    NCell = length(unique(x.b for x in Rj_vec))
    
    maxn(R::Rvec_3D) = maximum(abs,(R.n1,R.n2,R.n3))
    maxn(R::Rvec_2D) = maximum(abs,(R.n1,R.n2))

    L = maximum(maxn,Rj_vec)
    Lbox = 2L+1

    ChiBox(::Type{Rvec_3D}) = zeros(Lbox,Lbox,Lbox)
    ChiBox(::Type{Rvec_2D}) = zeros(Lbox,Lbox)

    offset = L+1
    coords(R::Rvec_3D) = CartesianIndex(offset .+ (R.n1, R.n2, R.n3))
    coords(R::Rvec_2D) = CartesianIndex(offset .+ (R.n1, R.n2))

    function chi_ij(α,β) 

        chi = ChiBox(RType)
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

getCorrelationPairs(Lattice::NamedTuple{(:pairNumberDict, :Basis)}) = getCorrelationPairs(Lattice.pairNumberDict)

import LatticeFFTs.getLatticeFFT
getLatticeFFT(S_ab,Basis::Basis_Struct,args...;kwargs...)  = getLatticeFFT(S_ab,[Basis.a1 Basis.a2 Basis.a3],Basis.b,args...;kwargs...)
getLatticeFFT(S_ab,Basis::Basis_Struct_2D,args...;kwargs...) = getLatticeFFT(S_ab,[Basis.a1 Basis.a2],Basis.b,args...;kwargs...)

getnaiveLatticeFT(S_ab,Basis::Basis_Struct,args...;kwargs...)  = naiveLatticeFT(S_ab,[Basis.a1 Basis.a2 Basis.a3],Basis.b,args...;kwargs...)
getnaiveLatticeFT(S_ab,Basis::Basis_Struct_2D,args...;kwargs...) = naiveLatticeFT(S_ab,[Basis.a1 Basis.a2],Basis.b,args...;kwargs...)

function getLatticeFFT(ChiR::AbstractVector,Lattice,args...;kwargs...)
    CorrelationPairs = getCorrelationPairs(Lattice)
    (;Ri_vec,Rj_vec,pairs) = CorrelationPairs
    S_ab = separateSublattices(Ri_vec,Rj_vec,ChiR[pairs])

    return getLatticeFFT(S_ab,Lattice.Basis,args...;kwargs...)
end

function getnaiveLatticeFT(ChiR::AbstractVector,Lattice,args...;kwargs...)
    CorrelationPairs = getCorrelationPairs(Lattice)
    (;Ri_vec,Rj_vec,pairs) = CorrelationPairs
    S_ab = separateSublattices(Ri_vec,Rj_vec,ChiR[pairs])

    return getnaiveLatticeFT(S_ab,Lattice.Basis,args...;kwargs...)
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