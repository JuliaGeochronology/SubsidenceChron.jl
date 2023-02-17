## --- Custom objects for holding SubsidenceChron data

    # A type of object to hold data about the stratigraphy
    mutable struct StratData
        Lithology::Array{String}
        Thickness::Array{Float64,1}
    end

    function NewStratData(nLayers)
        strat = StratData(
            fill("", nLayers),
            fill(NaN,nLayers),
        )
        return strat
    end
    export NewStratData

    # A type of object to hold data about the thermal subsidence parameters
    mutable struct ThermalSubsidenceParameters
        Param::Array{Float64,1}
        Sigma::Array{Float64,1}
    end

    function NewThermalSubsidenceParameters()
        therm = ThermalSubsidenceParameters(
            fill(NaN,2),
            fill(NaN,2),
        )
        return therm
    end
    export NewThermalSubsidenceParameters

    struct SubsidenceStratAgeModel
        Height::Array{Float64,1}
        Age::Array{Float64,1}
        Age_sigma::Array{Float64,1}
        Age_Median::Array{Float64,1}
        Age_025CI::Array{Float64,1}
        Age_975CI::Array{Float64,1}
        Beta::Array{Float64,1}
        Beta_sigma::Array{Float64,1}
        Beta_Median::Array{Float64,1}
        Beta_025CI::Array{Float64,1}
        Beta_975CI::Array{Float64,1}
        T0::Array{Float64,1}
        T0_sigma::Array{Float64,1}
        T0_Median::Array{Float64,1}
        T0_025CI::Array{Float64,1}
        T0_975CI::Array{Float64,1}
    end
