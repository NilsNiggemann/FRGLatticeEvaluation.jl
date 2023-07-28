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
##
@testset "LatticeFFTs" begin

    Points2D = Dict([
        "Γ" => (0,0),
        "K" => (4π/3,0),
        "M" => (π, -π/√3),
        "(3π,3π)" => (3π,3π),
        "(√(22),100/3" => (√(22),100/3)
    ]
    )
    Points3D = Dict([
        "Γ" => (0,0,0),
        "X" => (0,2π,0),
        "L" => (π, π,π),
        "K" => (3π/2, 3π/2,0),
        "(4π,4π,0)" => (4π,4π,0),
        "(√(22),100/3,5.1)" => (√(22),100/3,5.1)
    ]
    )

    @testset "2D" begin
        L1 = LatticeInfo(SquareKagome.getSquareKagome(7,(1,2,3,4)),SquareKagome)
        coupl = L1.System.couplings
        Fourier = getFourier(coupl,L1)
        FFT = interpolatedFT(coupl,L1,64)
        @testset "k = $key" for (key,val) in Points2D
            k = SVector(val)
            @test Fourier(k) ≈ FFT(k) atol = 1e-6
        end
    end
    @testset "3D" begin
        L1 = LatticeInfo(Pyrochlore.getPyrochlore(7),Pyrochlore)
        coupl = L1.System.couplings
        Fourier = getFourier(coupl,L1)
        FFT = interpolatedFT(coupl,L1,64)
        @testset "k = $key" for (key,val) in Points3D
            k = SVector(val)
            @test Fourier(k) ≈ FFT(k) atol = 1e-4
            @test FFT[1,2](k) ≈ FFT[2,1](k)' atol = 1e-16
            @test FFT[3,4](k) ≈ FFT[4,3](k)' atol = 1e-16
        end
    end
end