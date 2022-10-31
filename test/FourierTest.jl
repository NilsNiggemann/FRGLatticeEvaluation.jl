function test_fourier_onsite(pairs::AbstractVector{I},Rijvec::AbstractVector{RT},OnsitePairs::AbstractVector{I}) where {I<:Integer,RT <: StaticVector}
    Onsite1 = sort(findall(x-> x in OnsitePairs,pairs))
    Onsite2 = sort(findall(x-> norm(x) ≈ 0.,Rijvec))
    @testset "testing for number of onsite bonds" begin
        @test length(Onsite1) == length(Onsite2)
    end
    @testset "testing onsite bonds in Fourier infos" begin
        for (o1,o2) in zip(Onsite1,Onsite2)
            @test o1 == o2
        end
    end
end

test_fourier_onsite(Lat::LatticeInfo) = test_fourier_onsite(Lat.FourierInfos.pairs,Lat.FourierInfos.Rij_vec,Lat.System.OnsitePairs)

function test_fourier_pairs(pairs::AbstractVector{I},PairList::AbstractVector{RT}) where {I<:Integer,RT <: Rvec}
    @testset "testing that all Pairs are in Fourier Info" begin
        Inds = sort(unique(pairs))
        @test Inds == collect(eachindex(PairList))
    end
end

test_fourier_pairs(Lat::LatticeInfo) = test_fourier_pairs(Lat.FourierInfos.pairs,Lat.System.PairList)
