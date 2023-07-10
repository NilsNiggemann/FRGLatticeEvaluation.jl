using FFTViews, Interpolations
using FFTViews:fftfreq,fftshift,ifftshift,fft,rfft
abstract type AbstractPadding end
struct AutomaticPadding <: AbstractPadding end

"""
Given a function ChiR, return the Fourier transform via FFT.
Assumptions:
    - ChiR is given in lattice coordinates, i.e. Chi[n1,n2,n3] corresponds to chi(0,n1*a1+n2*a2+n3*a3) 
    - The central index corresponds to chi(0,0). 
    - each dimension is of odd length
"""
function getFFT(ChiR)
    # @assert all(isodd.(size(ChiR))) "works only for odd length of dimensions, but dims(ChiR) = $(size(ChiR))"
    chik = FFTView(fft(ifftshift(ChiR)))
    return chik
end

function padSusc(ChiR::AbstractArray{T},::AutomaticPadding) where T <: Number
    dims = size(ChiR)
    Nearest2Power = Tuple(ceil(Int,2^log2(d)) for d in dims)
    minSize = 64
    newDims =  Tuple(max(n,minSize) for n in Nearest2Power)
    return padSusc(ChiR,newDims)
end

function padSusc(ChiR::AbstractArray{T},newDims::Tuple) where T <: Number
    dims = size(ChiR)
    newDims =  max.(dims,newDims)
    shift = CartesianIndex(( (d ≠ nd) && iseven(nd)  for (d,nd) in zip(dims,newDims))...)

    Origin = CartesianIndex((newDims .- dims).÷ 2 )

    # @info "" dims newDims shift Origin
    PaddedChiR = zeros(T,newDims)
    for I in CartesianIndices(ChiR)
        newI = I + shift + Origin
        PaddedChiR[newI] = ChiR[I]
    end
    return PaddedChiR
end

function padSusc(ChiR::AbstractArray{T},val::Integer) where T <: Number
    val <= 0 && return ChiR
    padSusc(ChiR,Tuple(val for s in size(ChiR)))
end

function getInterpolatedFFT(Chi_ij::AbstractArray{<:Real},padding = 0)
    Chi_ij = padSusc(Chi_ij,padding)
    nk = Tuple(0:N for N in size(Chi_ij))
    FFT = real(getFFT(Chi_ij)[nk...])
    k = Tuple(2π/N .* (0:N) for N in size(Chi_ij))
    chik = Interpolations.interpolate(k,FFT, Gridded(Linear()))
    chik = extrapolate(chik,Periodic(OnGrid()))
    return chik
end

function separateSublattices(Ri_vec::AbstractVector{Rvec_3D},Rj_vec::AbstractVector{Rvec_3D},Chi_ij)
    NCell = length(unique(x.b for x in Rj_vec))
    maxn(R) = maximum(abs,(R.n1,R.n2,R.n3))
    L = maximum(maxn,Rj_vec)
    padding = 10
    # chi_ij(α,β) = [χ for (i,χ) in enumerate(Chi_ij) if Ri_vec[i].b == α && Rj_vec[i].b == β ]

    function chi_ij(α,β) 
        Lbox = 2L+1+2*padding

        chi = zeros(Lbox,Lbox,Lbox)
        offset = L+1+padding
        for (i,χ) in enumerate(Chi_ij)
            if Ri_vec[i].b == α && Rj_vec[i].b == β
                chi[ offset+ Rj_vec[i].n1, offset+ Rj_vec[i].n2, offset+ Rj_vec[i].n3] = χ
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


function interpolatedChi(Chi_ab,Basis)
    chiab = [linear_interpolation(getk(chi),getFFT(chi)) for chi in Chi_ab]

    chi(α,β,k) = exp(1im*k'*(Basis.b[α]-Basis.b[β]))* chiab[α,β]((Basis.T*k)...)

    return [k->chi(α,β,k) for α in 1:Basis.NCell, β in 1:Basis.NCell]
end