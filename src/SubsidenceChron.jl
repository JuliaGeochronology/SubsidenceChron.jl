module SubsidenceChron

    using Reexport
    @reexport using NaNStatistics
    @reexport using StatGeochemBase
    @reexport using Chron

    # Basic statistics and UI resources
    using ProgressMeter: @showprogress, Progress, update!
    using Distributions
    using LoopVectorization: @turbo
    using IfElse: ifelse

    # Custom objects for holding SubsidenceChron age data
    include("Objects.jl")

    # Decompaction and backstripping
    include("DecompactBackstrip.jl")

    # Age-depth modelling
    include("SubsidenceStratMetropolis.jl")

    # Structs
    export StratData, WaterDepth, ThermalSubsidenceParameters, SubsidenceStratAgeModel

    # High-level functions
    export DecompactBackstrip, SubsidenceStratMetropolis

    # Additional utility functions
    function find_formation_depths(formation, thickness)
        depth = cumsum([0; thickness])
        unique_formations = unique(formation)
        unique_formation_tops = [depth[findfirst(x->x==n, formation)] for n in unique_formations]
        unique_formation_bottoms = [depth[findlast(x->x==n, formation)+1] for n in unique_formations]
        return unique_formations, unique_formation_tops, unique_formation_bottoms
    end
    export find_formation_depths

end # module
