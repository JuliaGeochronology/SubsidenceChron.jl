using SubsidenceChron
using Test, Statistics, StatsBase

@testset "Decompaction" begin
    yₚ = zeros(3)
    y = [1,2,3]
    SubsidenceChron.decompact!(yₚ, y, 0.32, 1, 1, 1)
    SubsidenceChron.decompact!(yₚ, y, 0.32, 1, 2, 1)
    @test yₚ ≈ [0.0, 1.143612863313294, 2.1821135232335362]
end
@testset "Subsidence" begin include("testSubsidence.jl") end
@testset "Utilities" begin include("testUtilities.jl") end
@testset "Eruption / deposition age distributions" begin include("testDist.jl") end
@testset "Strat only" begin include("testStratOnly.jl") end
@testset "Radiocarbon" begin include("testRadiocarbon.jl") end
@testset "Coupled model" begin include("testCoupled.jl") end
