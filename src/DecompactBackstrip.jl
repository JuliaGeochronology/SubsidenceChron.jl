## --- Subsidence parameters
"""
```julia
subsidenceparams(lithology_string::AbstractString)
```

Return a tuple of the porosity-depth coefficient [1/m], surface porosity, and solid density of the given lithology.
Supported lithologies include:
```
shale
black shale
siltstone
shaly sandstone
sandstone
chalk
limestone
dolostone
anhydrite
quartzite
diamictite
conglomerate
breccia
diabase
basalt
andesite
rhyolite
tuff
```
"""
function subsidenceparams(lithology_string::AbstractString)
    # Turn string into lowercase symbol
    lithology = Symbol(lowercase(lithology_string))

    # Look up properties:
    #   Porosity-depth coefficient c [km^-1] has a lower bound of 0, since porosity should not decrease with depth
    #   Initial porosity ϕ₀ [unitless] has a lower bound of 0 and an upper bound of 1, definitionally
    #   Density ρg [kg/m^3] reflects that of the solid only (i.e., density at zero porosity)
    if lithology === :shale
        # From Allen & Allen 2013, Table 9.1
        c_dist = truncated(Normal(0.51, 0.15), 0, Inf) #
        ϕ₀_dist = truncated(Normal(0.63, 0.15), 0, 1)
        ρg = 2720
    elseif lithology === Symbol("black shale")
        # As normal shale, but more compressible and lower solid density
        c_dist = truncated(Normal(0.60, 0.15), 0, Inf) #
        ϕ₀_dist = truncated(Normal(0.63, 0.15), 0, 1)
        ρg = 2400
    elseif lithology === :siltstone
        c_dist = truncated(Normal(0.39, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.56, 0.1), 0, 1)
        ρg = 2650
    elseif lithology === Symbol("shaly sandstone")
        # From Allen & Allen 2013, Table 9.1
        c_dist = truncated(Normal(0.40, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.56, 0.1), 0, 1)
        ρg = 2650
    elseif lithology === :sandstone
        # From Allen & Allen 2013, Table 9.1
        c_dist = truncated(Normal(0.27, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.49, 0.1), 0, 1)
        ρg = 2650
    elseif lithology === :chalk
        # From Allen & Allen 2013, Table 9.1
        c_dist = truncated(Normal(0.71, 0.15), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.7, 0.15), 0, 1)
        ρg = 2710
    elseif lithology === :limestone
        c_dist = truncated(Normal(0.6, 0.2), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.4, 0.17), 0, 1)
        ρg = 2710
    elseif lithology === :dolostone
        c_dist = truncated(Normal(0.6, 0.2), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.2, 0.1), 0, 1)
        ρg = 2870
    elseif lithology === :anhydrite
        c_dist = truncated(Normal(0.2, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.05, 0.05), 0, 1)
        ρg = 2960
    elseif lithology === :quartzite
        # As sandstone, but lower porosity
        c_dist = truncated(Normal(0.27, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.15, 0.1), 0, 1)
        ρg = 2650
    elseif lithology === :diamictite
        c_dist = truncated(Normal(0.51, 0.15), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.63, 0.15), 0, 1)
        ρg = 2720
    elseif lithology === :conglomerate
        c_dist = truncated(Normal(0.51, 0.15), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.63, 0.15), 0, 1)
        ρg = 2720
    elseif lithology === :breccia
        c_dist = truncated(Normal(0.51, 0.15), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.63, 0.15), 0, 1)
        ρg = 2720
    elseif lithology === :diabase
        # As basalt, but lower initial porosity
        c_dist = truncated(Normal(0.5, 0.2), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.03, 0.05), 0, 1)
        ρg = 2960
    elseif lithology === :basalt
        c_dist = truncated(Normal(0.5, 0.2), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.095, 0.1), 0, 1)
        ρg = 2960
    elseif lithology === :andesite
        c_dist = truncated(Normal(0.5, 0.2), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.095, 0.1), 0, 1)
        ρg = 2650
    elseif lithology === :rhyolite
        c_dist = truncated(Normal(0.5, 0.2), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.275, 0.14), 0, 1)
        ρg = 2510
    elseif lithology === :tuff
        # As Rhyolite, but more porous and compressible
        c_dist = truncated(Normal(0.65, 0.3), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.35, 0.2), 0, 1)
        ρg = 2550
    else # fallback, if unknown
        @warn "lithology $lithology_string not recognized, using default porosity parameters"
        c_dist = truncated(Normal(0.5, 0.3), 0, Inf)
        ϕ₀_dist = Uniform(0, 1)
        ρg = 2700
    end
    # Divide by 1000 to convert c_dist from 1/km to 1/m
    return c_dist/1000, ϕ₀_dist, ρg
end
function subsidenceparams(lithology::Vector{<:AbstractString})
    # porosity depth coefficient(c)
    c_dist = similar(lithology, Distribution)
    # surface porosity (ϕ₀)
    ϕ₀_dist = similar(lithology, Distribution)
    # sediment grain density (ρg)
    ρg = similar(lithology, Float64)

    # Find the correct c, ϕ₀, and ρg for each layer based on lithology
    for i in eachindex(lithology)
        c_dist[i], ϕ₀_dist[i], ρg[i] = subsidenceparams(lithology[i])
    end
    return c_dist, ϕ₀_dist, ρg
end

## --- Stratigraphic MCMC model without hiata # # # # # # # # # # # # # # # # #

# Part 1: Decompaction and Backstripping

# Define the decompaction function
function decompact!(yₚ, y, ϕ₀, c, n, m; niterations = 10)
    # y[n+1] = y₂ = present-day depth for base of layer n
    # y[n] = y₁ = present-day depth for top of layer n
    δy = y[n+1] - y[n]
    # yₚ[n-m+2] = y₂' = depth for base of layer n during decompaction (at t=tₘ)
    # yₚ[n-m+1] = y₁' = depth for top of layer n during decompaction (at t=tₘ)
    yₚ[n-m+2] = yₚ[n-m+1] + δy
    # In this for loop, the initial approximatation for the thickness of a given layer is the present-day thickness of that layer (calculated above)
    # The true thickness is then reached by running the decompaction equation multiple times (niterations)
    e₁ = exp(-c*yₚ[n-m+1])
    e₃₄ = exp(-c*y[n]) - exp(-c*y[n+1])
    @inbounds for idx=1:niterations
        # yₚ[n-m+2] = yₚ[n-m+1] + δy + ϕ₀/c * ((exp(-c*yₚ[n-m+1]) - exp(-c*yₚ[n-m+2])) - (exp(-c*y[n]) - exp(-c*y[n+1])))
        yₚ[n-m+2] = yₚ[n-m+1] + δy + ϕ₀/c * ((e₁ - exp(-c*yₚ[n-m+2])) - e₃₄)
    end
end

# Decompaction and backstripping (Method 1: with water depth inputs)
function DecompactBackstrip(strat::StratData, wd::WaterDepth, nsims, res)

    # Import data from csv and assign parameters for each lithology
        lithology_inputs = strat.Lithology
        height_inputs = cumsum([0; strat.Thickness])
        nlayer_input = length(strat.Thickness)
        subsidence_strat_depths = [0:res:maximum(height_inputs);]
        model_nlayer = length(subsidence_strat_depths)-1

    # Allocate parameters as distributions; each element/distribution represents a layer
        # c: porosity depth coefficient, ϕ₀: surface porosity, ρg: sediment grain density
        c_dist, ϕ₀_dist, ρg = subsidenceparams(lithology_inputs)

    # Prep for decompaction and backstripping MC
        # Define parameters for decompaction (number of simulations; water and mantle densities)
        ρw = 1000
        ρm = 3330

        # Allocate depth matricies (rows are strat horizons, columns are timesteps)
        # Decompacted depth
        Y = fill(1E-18, (model_nlayer+1, model_nlayer+1))

        # Allocate decompaction parameter matricies (for all modeled layers with a resolution = res)
        c_highres = Array{Float64,1}(undef, model_nlayer)
        ϕ₀_highres = Array{Float64,1}(undef, model_nlayer)
        ρg_highres = Array{Float64,1}(undef, model_nlayer)

        # Allocate paleo water depth matrix (for all modeled layers with a resolution = res)
        paleo_wd_dist = Array{Float64,1}(undef, nsims)

        # Allocate porosity, density and tectonic subsidence matricies
        # Porosity of a strat unit at any depth
        ϕ_avg = fill(1E-18, (model_nlayer, model_nlayer+1))
        # Bulk density of a single layer
        ρ_bulk = fill(1E-18, (model_nlayer, model_nlayer+1))
        # Intermediate step - bulk density*thickness of a single layer
        m_bulk = fill(1E-18, (model_nlayer, model_nlayer+1))
        # Bulk density of the entire column
        ρ_bulk_column = Array{Float64,1}(undef, model_nlayer+1)
        # Tectonic subsidence # this is the only one that will need to propagate outside of the loop
        Sₜ = Array{Float64,2}(undef, model_nlayer+1, nsims) # meters

    # MC for decompaction and backstripping
        # Visualize the progress in terminal
        print("Decompaction and Backstripping: ", nsims, " steps\n")
        pgrs = Progress(nsims, desc="Decompaction and Backstripping...")
        pgrs_interval = ceil(Int,10)
        c = rand.(c_dist)
        ϕ₀ = rand.(ϕ₀_dist)

        for sim = 1:nsims
            # Randomly select c and ϕ₀ vectors from the distributions for each input layer
            @. c = rand(c_dist)
            @. ϕ₀ = rand(ϕ₀_dist)

            # Propagate these selections to every model layers; all model layers from the same input layer get the same c and ϕ₀ values
            @inbounds for i = 1:nlayer_input
                for j = 1:model_nlayer
                    if subsidence_strat_depths[j+1]>height_inputs[i]
                        c_highres[j]=c[i]
                        ϕ₀_highres[j]=ϕ₀[i]
                        ρg_highres[j]=ρg[i]
                    end
                end
            end

            # Fill the first column with modern observed values (present-day depths)
            Y[:,1] .= subsidence_strat_depths
            # Fill the first row with zeros
            Y[1,:] .= 0

            # Decompact
            # i = time steps during decompaction, which runs from 2 to layer_count b/c column 1 is present day
            # j = layer number, which runs from i to layer_count b/c the ith column begins with y₁' of layer i
            for i = 2:model_nlayer
                for j = i:model_nlayer
                    decompact!(view(Y,:,i), subsidence_strat_depths, ϕ₀_highres[j], c_highres[j], j, i)
                end
            end

            # Import water depth correction data (water depth categories based on sedimentological evidences)
            wd_id_inputs = wd.DepthID
            wd_height_inputs = cumsum([0; wd.Thickness])
            wd_nlayer_input = length(wd.Thickness)

            # Allocate matrix for paleo-water depth distributions (for all water depth input layers)
            paleo_wd_dist = Array{Distribution,1}(undef, wd_nlayer_input)
            # Define the water depth distributions for all input layers based on their water depth classifications
            for i = 1:wd_nlayer_input
                if wd_id_inputs[i] == "Exposure"
                    paleo_wd_dist[i] = Uniform(0, 0.01)
                elseif wd_id_inputs[i] == "FWWB"
                    paleo_wd_dist[i] = Uniform(5,15)
                elseif wd_id_inputs[i] == "SWB"
                    paleo_wd_dist[i] = Uniform(15,40)
                elseif wd_id_inputs[i] == "BWB"
                    paleo_wd_dist[i] = Uniform(40,100)
                end
            end

            # Randomly select paleo water depth value (from the distributions) for each input layer
            paleo_wd = rand.(paleo_wd_dist)
            # Allocate paleo_wd matrix for all modeled layers with a resolution = res
            paleo_wd_highres = Array{Float64,1}(undef, model_nlayer+1)

            # Setup this matrix - paleo water depth for the first (topmost) modeled layer should be the same as the paleo water depth for the first input layer
            paleo_wd_highres[1] = paleo_wd[1]
            # Propagate these selections to every model layers; all model layers from the same input layer get the same paleo water depth
            for i = 1:wd_nlayer_input
                for j = 1:model_nlayer+1
                    if subsidence_strat_depths[j]>wd_height_inputs[i]
                        paleo_wd_highres[j]=paleo_wd[i]
                    end
                end
            end

            # Remove the effect of sediment load - same indexing logic as the decompaction loop
            # i = time steps (columns) = 1:layer_count b/c need to calculate these parameters for the present day column/layers
            # j = layer number = i:layer_count b/c the ith column begins with y₁' of layer i
            @inbounds for i = 1:model_nlayer+1
                for j = i:model_nlayer
                    # Average porosity of all modeled layers
                    ϕ_avg[j-i+1,i] = (ϕ₀_highres[j]/c_highres[j])*(exp(-c_highres[j]*Y[j-i+1,i])-exp(-c_highres[j]*Y[j-i+2,i]))/(Y[j-i+2,i]-Y[j-i+1,i])
                    # Bulk density of all modeled layers
                    ρ_bulk[j-i+1,i] = ρw*ϕ_avg[j-i+1,i]+(1-ϕ_avg[j-i+1,i])*ρg_highres[j]
                    # Mass of all modeled layers
                    m_bulk[j-i+1,i] = (Y[j-i+2,i]-Y[j-i+1,i])*ρ_bulk[j-i+1,i]
                end
            end
            # (Mass = 0 for the last column, which represents the moment right before the deposition of the bottommost layer)
            m_bulk[:,end] .= 0

            # Backstripping calculations
            @inbounds for i = 1:model_nlayer+1
                # Bulk density of all the columns (aka bulk densities of the whole sediment column for all timesteps)
                # maximum([:,i]) = total depth of column at time step i
                ρ_bulk_column[i] = sum(view(m_bulk,:,i))/maximum(view(Y,:,i))
                # Tectonic subsidence for all timesteps; record results from all simulations
                Sₜ[i,sim] = Y[model_nlayer+2-i,i]*((ρm-ρ_bulk_column[i])/(ρm-ρw)).+paleo_wd_highres[i]
            end

            mod(sim,pgrs_interval)==0 && update!(pgrs, sim)
        end
        update!(pgrs,nsims)

        # Calculate summary statistics (mean and standard deviation)
        Sμ = nanmean(Sₜ, dim=2)
        Sσ = nanstd(Sₜ, dim=2)

    return Sₜ, Sμ, Sσ, subsidence_strat_depths
end

# Decompaction and backstripping (Method 2: without water depth inputs)
function DecompactBackstrip(strat::StratData, nsims, res)

    # Import data from csv and assign parameters for each lithology
        lithology_inputs = strat.Lithology
        height_inputs = cumsum([0; strat.Thickness])
        nlayer_input = length(strat.Thickness)
        subsidence_strat_depths = [0:res:maximum(height_inputs);]
        model_nlayer = length(subsidence_strat_depths)-1

    # Allocate parameters as distributions; each element/distribution represents a layer
        # c: porosity depth coefficient, ϕ₀: surface porosity, ρg: sediment grain density
        c_dist, ϕ₀_dist, ρg = subsidenceparams(lithology_inputs)

    # Prep for decompaction and backstripping MC
        # Define parameters for decompaction (number of simulations; water and mantle densities)
        ρw = 1000
        ρm = 3330

        # Allocate depth matricies (rows are strat horizons, columns are timesteps)
        # Decompacted depth
        Y = fill(1E-18, (model_nlayer+1, model_nlayer+1))

        # Allocate decompaction parameter matricies (for all modeled layers with a resolution = res)
        c_highres = Array{Float64,1}(undef, model_nlayer)
        ϕ₀_highres = Array{Float64,1}(undef, model_nlayer)
        ρg_highres = Array{Float64,1}(undef, model_nlayer)

        #Allocate porosity, density and tectonic subsidence matricies
        # porosity of a strat unit at any depth
        ϕ_avg = fill(1E-18, (model_nlayer, model_nlayer+1))
        # bulk density of a single layer
        ρ_bulk = fill(1E-18, (model_nlayer, model_nlayer+1))
        # intermediate step - bulk density*thickness of a single layer
        m_bulk = fill(1E-18, (model_nlayer, model_nlayer+1))
        # bulk density of the entire column
        ρ_bulk_column = Array{Float64,1}(undef, model_nlayer+1)
        # tectonic subsidence # this is the only one that will need to propagate outside of the loop
        Sₜ = Array{Float64,2}(undef, model_nlayer+1, nsims) # meters

    # MC for decompaction and backstripping
        # Visualize the progress in terminal
        print("Decompaction and Backstripping: ", nsims, " steps\n")
        pgrs = Progress(nsims, desc="Decompaction and Backstripping...")
        pgrs_interval = ceil(Int,10)
        c = rand.(c_dist)
        ϕ₀ = rand.(ϕ₀_dist)

        for sim = 1:nsims
            # Randomly select c and ϕ₀ vectors from the distributions for each input layer
            @. c = rand(c_dist)
            @. ϕ₀ = rand(ϕ₀_dist)

            # Propagate these selections to every model layers; all model layers from the same input layer get the same c and ϕ₀ values
            @inbounds for i = 1:nlayer_input
                for j = 1:model_nlayer
                    if subsidence_strat_depths[j+1]>height_inputs[i]
                        c_highres[j]=c[i]
                        ϕ₀_highres[j]=ϕ₀[i]
                        ρg_highres[j]=ρg[i]
                    end
                end
            end

            # Fill the first column with modern observed values (present-day depths)
            Y[:,1] .= subsidence_strat_depths
            # Fill the first row with zeros
            Y[1,:] .= 0

            # Decompact
            # i = time steps during decompaction, which runs from 2 to layer_count b/c column 1 is present day
            # j = layer number, which runs from i to layer_count b/c the ith column begins with y₁' of layer i
            for i = 2:model_nlayer
                for j = i:model_nlayer
                    decompact!(view(Y,:,i), subsidence_strat_depths, ϕ₀_highres[j], c_highres[j], j, i)
                end
            end

            # Remove the effect of sediment load - same indexing logic as the decompaction loop
            # i = time steps (columns) = 1:layer_count b/c need to calculate these parameters for the present day column/layers
            # j = layer number = i:layer_count b/c the ith column begins with y₁' of layer i
            @inbounds for i = 1:model_nlayer+1
                for j = i:model_nlayer
                    # Average porosity of all modeled layers
                    ϕ_avg[j-i+1,i] = (ϕ₀_highres[j]/c_highres[j])*(exp(-c_highres[j]*Y[j-i+1,i])-exp(-c_highres[j]*Y[j-i+2,i]))/(Y[j-i+2,i]-Y[j-i+1,i])
                    # Bulk density of all modeled layers
                    ρ_bulk[j-i+1,i] = ρw*ϕ_avg[j-i+1,i]+(1-ϕ_avg[j-i+1,i])*ρg_highres[j]
                    # Mass of all modeled layers
                    m_bulk[j-i+1,i] = (Y[j-i+2,i]-Y[j-i+1,i])*ρ_bulk[j-i+1,i]
                end
            end
            # (Mass = 0 for the last column, which represents the moment right before the deposition of the bottommost layer)
            m_bulk[:,end] .= 0

            # Backstripping calculations
            for i = 1:model_nlayer+1
                # Bulk density of all the columns (aka bulk densities of the whole sediment column for all timesteps)
                # maximum(Y_corr[:,i]) = total depth of column at time step i
                ρ_bulk_column[i] = sum(view(m_bulk,:,i))/maximum(view(Y,:,i))
                # Tectonic subsidence for all timesteps; record results from all simulations
                Sₜ[i,sim] = Y[model_nlayer+2-i,i]*((ρm-ρ_bulk_column[i])/(ρm-ρw))
            end

            mod(sim,pgrs_interval)==0 && update!(pgrs, sim)
        end
        update!(pgrs,nsims)

        # Calculate summary statistics (mean and standard deviation)
        Sμ = nanmean(Sₜ, dim=2)
        Sσ = nanstd(Sₜ, dim=2)

    return Sₜ, Sμ, Sσ, subsidence_strat_depths
end
