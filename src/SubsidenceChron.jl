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

    # Utility functions
    include("Utilities.jl")

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

    # Other exported functions
    export find_formation_depths, find_formation_heights

end # module
