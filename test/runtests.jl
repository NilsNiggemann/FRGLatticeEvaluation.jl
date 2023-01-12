using SpinFRGLattices
using FRGLatticeEvaluation
using Test
using SpinFRGLattices.StaticArrays

include("FourierTest.jl")

@testset "SquareLattice" begin
    L1 = LatticeInfo(SquareLattice.getSquareLattice(7),SquareLattice)
    test_fourier_onsite(L1)
    test_fourier_pairs(L1)
end

@testset "Pyrochlore" begin
    L1 = LatticeInfo(Pyrochlore.getPyrochlore(7),Pyrochlore)
    test_fourier_onsite(L1)
    test_fourier_pairs(L1)
end

@testset "SquareKagome" begin
    L1 = LatticeInfo(SquareKagome.getSquareKagome(7),SquareKagome)
    test_fourier_onsite(L1)
    test_fourier_pairs(L1)
end

