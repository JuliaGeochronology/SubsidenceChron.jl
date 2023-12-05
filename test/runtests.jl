using SubsidenceChron
using Test, Statistics
const make_plots = get(ENV, "MAKE_PLOTS", false) == "true"

@testset "Decompaction" begin
    yₚ = zeros(3)
    y = [1,2,3]
    SubsidenceChron.decompact!(yₚ, y, 0.32, 1, 1, 1)
    SubsidenceChron.decompact!(yₚ, y, 0.32, 1, 2, 1)
    @test yₚ ≈ [0.0, 1.143612863313294, 2.1821135232335362]
end
@testset "Subsidence" begin include("testSubsidence.jl") end
@testset "Strat only" begin include("testStratOnly.jl") end
@testset "Coupled model" begin include("testCoupled.jl") end
