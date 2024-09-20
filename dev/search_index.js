var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SubsidenceChron","category":"page"},{"location":"#SubsidenceChron","page":"Home","title":"SubsidenceChron","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SubsidenceChron.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SubsidenceChron]","category":"page"},{"location":"#SubsidenceChron.SubsidenceStratMetropolis-Tuple{ChronAgeData, StratAgeModelConfiguration, ThermalSubsidenceParameters, Vararg{Any, 5}}","page":"Home","title":"SubsidenceChron.SubsidenceStratMetropolis","text":"SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_depths, Sμ, Sσ, beta_ip, t0_ip;\n    subsidencebotom=minimum(smpl.Height),\n    subsidencetop=maximum(smpl.Height),\n    lithosphere = Normal(125000, 100),\n)\n\nRuns the main SubsidenceChron.jl age-depth model routine given a set of Gaussian age constraints specified in the smpl struct, an age-depth model configuration specified in the config struct, thermal subsidence parameters defined in the therm struct, decompation and backstripping outputs in the form of subsidence_strat_depths, Sμ, and Sσ, and prior estimates for stretching factor Beta (beta_ip) and time of thermal subsedence onset (t0_ip).\n\nExamples:\n\n(subsmdl, agedist, lldist, beta_t0dist, lldist_burnin) = SubsidenceStratMetropolis(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ_corr, Beta_sigma/10, T0_sigma/10)\n\n\n\n\n\n","category":"method"},{"location":"#SubsidenceChron.SubsidenceStratMetropolis_Height-Tuple{ChronAgeData, StratAgeModelConfiguration, ThermalSubsidenceParameters, Vararg{Any, 5}}","page":"Home","title":"SubsidenceChron.SubsidenceStratMetropolis_Height","text":"SubsidenceStratMetropolis_Height(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_depths, Sμ, Sσ, beta_ip, t0_ip;\n    subsidencetop=maximum(smpl.Height),\n    lithosphere = Normal(125000, 100),\n)\n\nRuns the main SubsidenceChron.jl age-depth model routine given a set of Gaussian age constraints specified in the smpl struct, an age-depth model configuration specified in the config struct, thermal subsidence parameters defined in the therm struct, decompation and backstripping outputs in the form of subsidence_strat_depths, Sμ, and Sσ, and prior estimates for stretching factor Beta (beta_ip) and time of thermal subsedence onset (t0_ip).\n\nExamples:\n\n(subsmdl, agedist, lldist, beta_tsdist, lldist_burnin, zdist) = SubsidenceStratMetropolis_Height(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ_corr, Beta_sigma/10, T0_sigma/10)\n\n\n\n\n\n","category":"method"},{"location":"#SubsidenceChron.find_formation_depths-Tuple{Any, Any}","page":"Home","title":"SubsidenceChron.find_formation_depths","text":"find_formation_depths(formation, thickness)\n\nFind the top and bottom depth of each unique formation given a list of stratigraphic units and thicknesses\n\nExamples:\n\nds = importdataset(joinpath(\"..\", \"examples\", \"Svalbard_highres.csv\"), importas=:Tuple)\nformations, tops bottoms = find_formation_depths(ds.Formation, ds.Thickness)\n\n\n\n\n\n","category":"method"},{"location":"#SubsidenceChron.find_formation_heights-Tuple{Any, Any}","page":"Home","title":"SubsidenceChron.find_formation_heights","text":"find_formation_heights(formation, thickness)\n\nFind the top and bottom height of each unique formation given a list of stratigraphic units and thicknesses\n\nExamples:\n\nds = importdataset(joinpath(\"..\", \"examples\", \"Svalbard_highres.csv\"), importas=:Tuple)\nformations, tops bottoms = find_formation_heights(ds.Formation, ds.Thickness)\n\n\n\n\n\n","category":"method"},{"location":"#SubsidenceChron.subsidenceparams-Tuple{AbstractString}","page":"Home","title":"SubsidenceChron.subsidenceparams","text":"subsidenceparams(lithology_string::AbstractString)\n\nReturn a tuple of the porosity-depth coefficient [1/m], surface porosity, and solid density of the given lithology. Supported lithologies include:\n\nshale\nblack shale\nsiltstone\nshaly sandstone\nsandstone\nchalk\nlimestone\ndolostone\nanhydrite\nquartzite\ndiamictite\nconglomerate\nbreccia\ndiabase\nbasalt\nandesite\nrhyolite\ntuff\n\n\n\n\n\n","category":"method"}]
}
