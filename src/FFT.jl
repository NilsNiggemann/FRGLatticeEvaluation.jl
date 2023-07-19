function separateSublattices(Ri_vec::AbstractVector{RType},Rj_vec::AbstractVector{RType},Chi_ij) where RType <: Rvec
    NCell = length(unique(x.b for x in Rj_vec))
    Lbox = 2L+1+2*padding

    ChiBox(::Type{Rvec_3D}) = zeros(Lbox,Lbox,Lbox)
    ChiBox(::Type{Rvec_2D}) = zeros(Lbox,Lbox)

    maxn(R::Rvec_3D) = maximum(abs,(R.n1,R.n2,R.n3))
    maxn(R::Rvec_2D) = maximum(abs,(R.n1,R.n2))

    L = maximum(maxn,Rj_vec)
    offset = L+1+padding
    coords(R::Rvec_3D) = offset .+ CartesianIndex(R.n1, R.n2, R.n3)
    coords(R::Rvec_2D) = offset .+ CartesianIndex(R.n1, R.n2)

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

getCorrelationPairs(Lat::AbstractLattice) = getCorrelationPairs(Lat.UnitCell,Lat.SiteList,Lat.PairList,Lat.PairTypes,Lat.pairToInequiv,Lat.Basis)