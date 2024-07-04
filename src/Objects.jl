## --- Custom objects for holding SubsidenceChron data

    # A type of object to hold data about the stratigraphy
    mutable struct StratData
        Lithology::Array{String}
        Thickness::Array{Float64,1}
    end
    function StratData(nLayers)
        StratData(
            fill("", nLayers),
            fill(NaN,nLayers),
        )
    end

    mutable struct WaterDepth
        DepthID::Array{String}
        Thickness::Array{Float64,1}
    end
    function WaterDepth(wd_nLayers)
        WaterDepth(
            fill("", wd_nLayers),
            fill(NaN, wd_nLayers),
        )
    end

    # A type of object to hold data about the thermal subsidence parameters
    mutable struct ThermalSubsidenceParameters
        Param::Array{Float64,1}
        Sigma::Array{Float64,1}
    end
    function ThermalSubsidenceParameters()
        ThermalSubsidenceParameters(
            fill(NaN,2),
            fill(NaN,2),
        )
    end

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
        lithosphere::Array{Float64,1}
        lithosphere_sigma::Array{Float64,1}
        lithosphere_Median::Array{Float64,1}
        lithosphere_025CI::Array{Float64,1}
        lithosphere_975CI::Array{Float64,1}
    end

    struct SubsidenceStratAgeModel_Height
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
        TSHeight::Array{Float64,1}
        TSHeight_sigma::Array{Float64,1}
        TSHeight_Median::Array{Float64,1}
        TSHeight_025CI::Array{Float64,1}
        TSHeight_975CI::Array{Float64,1}
        lithosphere::Array{Float64,1}
        lithosphere_sigma::Array{Float64,1}
        lithosphere_Median::Array{Float64,1}
        lithosphere_025CI::Array{Float64,1}
        lithosphere_975CI::Array{Float64,1}
    end

# For backwards compatibility
const NewStratData = StratData
const NewWaterDepth = WaterDepth
const NewThermalSubsidenceParameters = ThermalSubsidenceParameters
export NewStratData, NewWaterDepth, NewThermalSubsidenceParameters
