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

    #Functions for subsidence modeling
    include("SubsidenceStratMetropolis.jl")

    # Structs
    export StratData, ThermalSubsidenceParameters, SubsidenceStratAgeModel

    # High-level functions
    export DecompactBackstrip, SubsidenceStratMetropolis

end # module
