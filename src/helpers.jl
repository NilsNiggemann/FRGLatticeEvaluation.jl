strd(x,n=2) = string(round(x,digits = n))
strd3(x) = string(round(x,digits = 3))

function getNumberFromName(Name,subName)
    res_string = split(Name,subName*"=")[end]
    for i in length(res_string):-1:1
        N = tryparse(Int,res_string[1:i])
        if N !== nothing
            return N
        end
    end
    error("Could not get ", subName, "from string ",Name)
end

# planes for Fourier space
@inline hhlplane(x,z) = SA[x,x,z]
@inline normalTohhlplane(x,y) = SA[x,-x,y]
@inline xyplane(x,y) = SA[x,y]
@inline zzerocut(x,y) = SA[x,y,0]
@inline function sphereplane(origin,radius)
    sphere(θ,ϕ) = radius*SA[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ),cos(θ)] +origin
end

xylabels() = Dict([:xlabel => L"qx",:ylabel => L"qy"])
zzerolabels() = Dict([:xlabel => L"[h00]",:ylabel => L"[0l0]"])
hhllabels() = Dict([:xlabel => L"[hh0]",:ylabel => L"[00l]"])
normalTohhllabels() = (xlabel = L"[h\bar{h}0]",ylabel = L"[00l]")

getPair(R1::Rvec,R2::Rvec,Lattice) = SpinFRGLattices.pairNumber(R1,R2,Lattice.System.PairList,Lattice.System.PairTypes,Lattice.Basis,Lattice.pairToInequiv)

function getPairs(R1,R2,UnitCellTranslation::StaticVector,Lattice::LatticeInfo)
    (;System,Basis,pairToInequiv) = Lattice
    Pair1 = getPair(R1,R2,Lattice)
    R1_pr = translation(R1,UnitCellTranslation,Basis)
    R2_pr = translation(R2,UnitCellTranslation,Basis)
    Pair2 = getPair(R1_pr,R2_pr,Lattice)
    return Pair1,Pair2
end

function DimerResponse(Chi1,Chi2,delta)
    return @. (Chi1 - Chi2) / (Chi1 + Chi2) / delta
end


function DimerResponse(R1::Rvec,R2::Rvec,R1prime::Rvec,R2prime::Rvec,Chi_LR,delta,Lattice::LatticeInfo)
    Pair1,Pair2 = getPair(R1,R2,Lattice),getPair(R1prime,R2prime,Lattice)
    @views Chi1 = Chi_LR[:,Pair1]
    @views Chi2 = Chi_LR[:,Pair2]
    DimerResponse(Chi1,Chi2,delta)
end
function DimerResponse(R1::Rvec,R2::Rvec,UnitCellTranslation::StaticVector,Chi_LR,delta,Lattice::LatticeInfo)
    Pair1,Pair2 = getPairs(R1,R2,UnitCellTranslation,Lattice)
    @views Chi1 = Chi_LR[:,Pair1]
    @views Chi2 = Chi_LR[:,Pair2]
    DimerResponse(Chi1,Chi2,delta)
end


function DimerResponse_Keles(Chi1,Chi2,delta)
    return @. (Chi1 - Chi2) / (Chi2) / delta
end

function DimerResponse_Keles(R1::Rvec,R2::Rvec,Chi_LR::AbstractMatrix,Chi_LRNew::AbstractMatrix,delta,Lattice::LatticeInfo)
    Chiold = DimerFlow(R1,R2,Chi_LR,Lattice)
    ChiNew = DimerFlow(R1,R2,Chi_LRNew,Lattice)
    DimerResponse_Keles(ChiNew,Chiold,delta)
end

function DimerFlow(R1::Rvec,R2::Rvec,Chi_LR::AbstractMatrix,Lattice::LatticeInfo)
    Pair = getPair(R1,R2,Lattice)
    @views Chi1 = Chi_LR[:,Pair]
end