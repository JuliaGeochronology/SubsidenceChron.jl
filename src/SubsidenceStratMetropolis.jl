## --- Subsidence parameters

function subsidenceparams(lithology::AbstractString)
    if lithology == "Shale"
        # c has a lower bound of 0 b/c porosity at depth should not be greater than porosity at surface
        c_dist = truncated(Normal(0.51, 0.15), 0, Inf)
        # ϕ₀ has a lower bound of 0 and an upper bound of 1 b/c of the definition of porosity
        ϕ₀_dist = truncated(Normal(0.63, 0.15), 0, 1)
        ρg = 2720
    elseif lithology == "Siltstone"
        c_dist = truncated(Normal(0.39, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.56, 0.1), 0, 1)
        ρg = 2640 # Might need to update this number
    elseif lithology == "Sandstone"
        c_dist = truncated(Normal(0.27, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.49, 0.1), 0, 1)
        ρg = 2650
    elseif lithology == "Chalk"
        c_dist = truncated(Normal(0.71, 0.15), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.7, 0.15), 0, 1)
        ρg = 2710
    elseif lithology == "Limestone"
        c_dist = truncated(Normal(0.6, 0.2), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.4, 0.17), 0, 1)
        ρg = 2710
    elseif lithology == "Dolostone"
        c_dist = truncated(Normal(0.6, 0.2), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.2, 0.1), 0, 1)
        ρg = 2870
    elseif lithology == "Anhydrite"
        c_dist = truncated(Normal(0.2, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.05, 0.05), 0, 1)
        ρg = 2960
    elseif lithology == "Quartzite"
        c_dist = truncated(Normal(0.3, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.2, 0.1), 0, 1)
        ρg = 2650
     elseif lithology == "Diabase"
        c_dist = truncated(Normal(0.65, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.285, 0.1), 0, 1)
        ρg = 2960
    elseif lithology == "Rhyolite"
        c_dist = truncated(Normal(0.65, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.275, 0.1), 0, 1)
        ρg = 2510
    elseif lithology == "Diamictite"
        c_dist = truncated(Normal(0.51, 0.15), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.63, 0.15), 0, 1)
        ρg = 2720
    elseif lithology == "Conglomerate"
        c_dist = truncated(Normal(0.51, 0.15), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.63, 0.15), 0, 1)
        ρg = 2720
    elseif lithology == "Breccia"
        c_dist = truncated(Normal(0.51, 0.15), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.63, 0.15), 0, 1)
        ρg = 2720
    elseif lithology == "Basalt"
        c_dist = truncated(Normal(0.65, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.05, 0.1), 0, 1)
        ρg = 2980
    elseif lithology == "Andesite"
        c_dist = truncated(Normal(0.65, 0.1), 0, Inf)
        ϕ₀_dist = truncated(Normal(0.095, 0.1), 0, 1)
        ρg = 2650
    else # fallback, if unknown
        @warn "lithology $lithology not recognized, using default porosity parameters"
        c_dist = truncated(Normal(0.5, 0.3), 0, Inf)
        ϕ₀_dist = Uniform(0, 1)
        ρg = 2700
    end
    return c_dist, ϕ₀_dist, ρg
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
        model_strat_heights = [0:res:maximum(height_inputs);]
        model_nlayer = length(model_strat_heights)-1

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
        Sₜ_km = Array{Float64,2}(undef, model_nlayer+1, nsims)

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
                    if model_strat_heights[j+1]>height_inputs[i]
                        c_highres[j]=c[i]
                        ϕ₀_highres[j]=ϕ₀[i]
                        ρg_highres[j]=ρg[i]
                    end
                end
            end

            # Fill the first column with modern observed values (present-day depths)
            Y[:,1] .= model_strat_heights
            # Fill the first row with zeros
            Y[1,:] .= 0

            # Decompact
            # i = time steps during decompaction, which runs from 2 to layer_count b/c column 1 is present day
            # j = layer number, which runs from i to layer_count b/c the ith column begins with y₁' of layer i
            for i = 2:model_nlayer
                for j = i:model_nlayer
                    decompact!(view(Y,:,i), model_strat_heights, ϕ₀_highres[j], c_highres[j], j, i)
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
                    paleo_wd_dist[i] = Uniform(0, 0.00001)
                elseif wd_id_inputs[i] == "FWWB"
                    paleo_wd_dist[i] = Uniform(0.005,0.015)
                elseif wd_id_inputs[i] == "SWB"
                    paleo_wd_dist[i] = Uniform(0.015,0.04)
                elseif wd_id_inputs[i] == "BWB"
                    paleo_wd_dist[i] = Uniform(0.04,0.1)
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
                    if model_strat_heights[j]>wd_height_inputs[i]
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
                Sₜ_km[i,sim] = Y[model_nlayer+2-i,i]*((ρm-ρ_bulk_column[i])/(ρm-ρw)).+paleo_wd_highres[i]
            end

            mod(sim,pgrs_interval)==0 && update!(pgrs, sim)
        end
        update!(pgrs,nsims)

        # Calculate summary statistics (mean and standard deviation)
        Sₜ = Sₜ_km .* 1000
        Sμ = nanmean(Sₜ, dim=2)
        Sσ = nanstd(Sₜ, dim=2)

    return Sₜ, Sμ, Sσ, model_strat_heights
end

# Decompaction and backstripping (Method 2: without water depth inputs)
function DecompactBackstrip(strat::StratData, nsims, res)

    # Import data from csv and assign parameters for each lithology
        lithology_inputs = strat.Lithology
        height_inputs = cumsum([0; strat.Thickness])
        nlayer_input = length(strat.Thickness)
        model_strat_heights = [0:res:maximum(height_inputs);]
        model_nlayer = length(model_strat_heights)-1

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
        Sₜ_km = Array{Float64,2}(undef, model_nlayer+1, nsims)

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
                    if model_strat_heights[j+1]>height_inputs[i]
                        c_highres[j]=c[i]
                        ϕ₀_highres[j]=ϕ₀[i]
                        ρg_highres[j]=ρg[i]
                    end
                end
            end

            # Fill the first column with modern observed values (present-day depths)
            Y[:,1] .= model_strat_heights
            # Fill the first row with zeros
            Y[1,:] .= 0

            # Decompact
            # i = time steps during decompaction, which runs from 2 to layer_count b/c column 1 is present day
            # j = layer number, which runs from i to layer_count b/c the ith column begins with y₁' of layer i
            for i = 2:model_nlayer
                for j = i:model_nlayer
                    decompact!(view(Y,:,i), model_strat_heights, ϕ₀_highres[j], c_highres[j], j, i)
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
                Sₜ_km[i,sim] = Y[model_nlayer+2-i,i]*((ρm-ρ_bulk_column[i])/(ρm-ρw))
            end

            mod(sim,pgrs_interval)==0 && update!(pgrs, sim)
        end
        update!(pgrs,nsims)

        # Calculate summary statistics (mean and standard deviation)
        Sₜ = Sₜ_km .* 1000
        # reverse!(Sₜ, dims=1)
        # Sₜ .- Sₜ[findfirst(active), :]'
        Sμ = nanmean(Sₜ, dim=2)
        Sσ = nanstd(Sₜ, dim=2)

    return Sₜ, Sμ, Sσ, model_strat_heights
end


## --- Part 2: Age-depth modelling

# Define the function that calculates the log likelihood of the thermal subsidence fit
function subsidence_ll(E₀, τ, model_St, model_St_sigma, model_t, beta_t0)
    @assert eachindex(model_St) == eachindex(model_St_sigma) == eachindex(model_t)
    β, T₀ = beta_t0
    ll = zero(float(eltype(model_St_sigma)))
    @turbo for i ∈ eachindex(model_t)
        # Calculate subsidence_model_heights given this unique smooth thermal subsidence model at each age in model_t
        x = (E₀*β/pi)*sin(pi/β)*(1 -exp(-(T₀ -model_t[i])/τ))
        mu = model_St[i]
        sigma = model_St_sigma[i]
        # Turn that into a log likelihood using some age uncertainty of the curve
        ll -= (x-mu)*(x-mu) / (2*sigma*sigma)
    end
    return ll

    return log_likelihood
end

# Part 2a: Modified StratMetropolis for extensional basins - without hiatus

function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, model_strat_heights, Sμ, Sσ, beta_ip, t0_ip;
        subsidencebottom=minimum(smpl.Height),
        subsidencetop=maximum(smpl.Height)
    )

    # Run stratigraphic MCMC model
    print("Generating stratigraphic age-depth model...\n")

    # Define thermal subsidence model parameters
        y_litho= 125000
        ρ_mantle = 3330
        ρ_water = 1000
        αᵥ = 3.28*10^(-5)
        T_mantle = 1333
        τ = 50 #Myr
        E₀ = (4*y_litho*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water))

    # Stratigraphic age constraints
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)

    # Model configuration -- read from struct
        resolution = config.resolution
        bounding = config.bounding
        nsteps = config.nsteps
        burnin = config.burnin
        sieve = config.sieve
        model_heights = bottom:resolution:top

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = round((top-bottom)*bounding/resolution)*resolution
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [nanmean(Age_sigma)/10; Age_sigma; nanmean(Age_sigma)/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
        end

        active_height_t = bottom .<= model_heights .<= top
        npoints = length(model_heights)

    # STEP 1: calculate log likelihood of the modeled ages (and heights) in the initial proposal
        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]

        #Age_twosided = Array{Float64}(undef, count(==(0), Age))
        #Age_sigma_twosided = Array{Float64}(undef, count(==(0), Age))

        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
            #if Age_sidedness[i] == 0
            #    for j = 1:length(Age_twosided)
            #        Age_twosided[j] = Age[i]
            #        Age_sigma_twosided[j] = Age_sigma[i]
            #    end
            #end
        end

        ll = normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

    # STEP 2: calculate log likelihood of the subsidence model parameters in the initial proposal
        # Define thermal subsidence model parameters
        ideal_subs_parameters = therm.Param
        ideal_subs_parameters_sigma = therm.Sigma

        # Initial proposal for subsidence parameters - randomly pick a set of values from the distribution
        subs_parameters = [ideal_subs_parameters[1]+beta_ip, ideal_subs_parameters[2]+t0_ip]
        # Calculate log likelihood of this initial proposal for subsidence parameters
        ll += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)

    # STEP 3: calculate log likelihood for the fit of the thermal subsidence curve in the initial proposal
        subsidence_height_t = subsidencebottom .<= model_heights .<= subsidencetop
        ts_model_ages = reverse(model_ages[subsidence_height_t])
        input_subsidence_t = subsidencebottom .<= (-model_strat_heights*1000) .<= subsidencetop
        ts_Sμ = Sμ[input_subsidence_t]
        ts_Sσ = Sσ[input_subsidence_t]
        ll += subsidence_ll(E₀, τ, ts_Sμ, ts_Sσ, ts_model_ages, subs_parameters)/length(ts_Sμ)

    # Preallocate variables for MCMC proposals
        llₚ = ll
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)
        subs_parametersₚ = copy(subs_parameters)

    # Run burnin
        # acceptancedist = fill(false,burnin)
        lldist_burnin = Array{Float64}(undef,burnin÷1000)

        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ .+= randn.() .* ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points that are still stratigraphically below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points that are still stratigraphically above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            ts_model_agesₚ = reverse(model_agesₚ[subsidence_height_t])
            llₚ += subsidence_ll(E₀, τ, ts_Sμ, ts_Sσ, ts_model_agesₚ, subs_parametersₚ)/length(ts_Sμ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
                # acceptancedist[i] = true
            end
            if mod(n,1000) == 0
                lldist_burnin[n÷1000] = ll
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

    # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        beta_t0dist = Array{Float64}(undef,2,nsteps)

    # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            ts_model_agesₚ = reverse(model_agesₚ[subsidence_height_t])
            llₚ += subsidence_ll(E₀, τ, ts_Sμ, ts_Sσ, ts_model_agesₚ, subs_parametersₚ)/length(ts_Sμ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
                beta_t0dist[:,n÷sieve] .= subs_parameters
                #predicted_ages[n÷sieve] = beta_t0dist[2]+τ*log(1-(Sμ[?]*pi)/(E₀*beta_t0dist[1]*sin(pi/beta_t0dist[2])))
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

    # Crop the result
        agedist = agedist[active_height_t,:]
        subsmdl = SubsidenceStratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            nanpctile(agedist,2.5,dim=2), # 2.5th percentile
            nanpctile(agedist,97.5,dim=2), # 97.5th percentile
            nanmean(beta_t0dist[1,:],dim=1), # Mean beta
            nanstd(beta_t0dist[1,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[1,:],dim=1), # Median beta
            nanpctile(beta_t0dist[1,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[1,:],97.5,dim=1), # 97.5th percentile
            nanmean(beta_t0dist[2,:],dim=1), # Mean T0
            nanstd(beta_t0dist[2,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[2,:],dim=1), # Median T0
            nanpctile(beta_t0dist[2,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[2,:],97.5,dim=1) # 97.5th percentile
        )
    return subsmdl, agedist, lldist, beta_t0dist, lldist_burnin
end

#= Haven't checked the hiatus part yet - will review this for the next round of changes
## --- Stratigraphic MCMC model with hiatus # # # # # # # # # # # # # # # # #


# Part 2b: Modified StratMetropolis for extensional basins with hiatus (when we only know the strat height and do not have any previous knowledge about the duration)

function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, model_strat_heights, Sμ, Sσ, hiatus_height, beta_ip, t0_ip)

    # Run stratigraphic MCMC model
    print("Generating stratigraphic age-depth model...\n")

    # Thermal subsidence model parameters
        y_litho= 125000
        ρ_mantle = 3330
        ρ_water = 1000
        αᵥ = 3.28*10^(-5)
        T_mantle = 1333
        τ = 50 #Myr
        E₀ = (4*y_litho*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water))

    # Stratigraphic age constraints
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)
        model_heights = copy(-model_strat_heights[2:end]).*1000

    # Model configuration -- read from struct
        resolution = config.resolution
        bounding = config.bounding
        nsteps = config.nsteps
        burnin = config.burnin
        sieve = config.sieve

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = round((top-bottom)*bounding)
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [nanmean(Age_sigma)/10; Age_sigma; nanmean(Age_sigma)/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
        end

        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

    # STEP 1: calculate log likelihood of the modeled ages (and heights) in the initial proposal
        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end

        ll = normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

    # STEP 2: calculate log likelihood of the subsidence model parameters in the initial proposal
        # Define thermal subsidence model parameters
        ideal_subs_parameters = therm.Param
        ideal_subs_parameters_sigma = therm.Sigma

        # Initial proposal for subsidence parameters - randomly pick a set of values from the distribution
        subs_parameters = [ideal_subs_parameters[1]+beta_ip, ideal_subs_parameters[2]+t0_ip]
        # Calculate log likelihood of this initial proposal for subsidence parameters
        ll += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)

    # STEP 3: calculate log likelihood for the fit of the thermal subsidence curve in the initial proposal
        active_height_t_hiatus = bottom .<= model_heights .< hiatus_height
        ts_model_ages_hiatus = reverse(model_ages[active_height_t_hiatus])
        model_heights_all = copy(-model_strat_heights).*1000
        ts_Sμ_hiatus = Sμ[bottom .<= model_heights_all .< hiatus_height]
        ts_Sσ_hiatus = Sσ[bottom .<= model_heights_all .< hiatus_height]
        ll += subsidence_ll(E₀, τ, ts_Sμ_hiatus, ts_Sσ_hiatus, ts_model_ages_hiatus, subs_parameters)/(length(ts_Sμ))

    # Preallocate variables for MCMC proposals
        llₚ = ll
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)
        subs_parametersₚ = copy(subs_parameters)

    # Run burnin
        # acceptancedist = fill(false,burnin)
        lldist_burnin = Array{Float64}(undef,burnin÷10)

        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points that are still stratigraphically below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points that are still stratigraphically above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            ts_model_ages_hiatusₚ = reverse(model_agesₚ[active_height_t_hiatus])
            llₚ += subsidence_ll(E₀, τ, ts_Sμ_hiatus, ts_Sσ_hiatus, ts_model_ages_hiatusₚ, subs_parametersₚ)/(length(ts_Sμ))

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
                # acceptancedist[i] = true
            end
            if mod(n,100) == 0
                lldist_burnin[[n÷100]] = ll
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

    # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        beta_t0dist = Array{Float64}(undef,2,nsteps)

    # Run the model
        h = plot(xlabel="Age", ylabel="Height", framestyle=:box)

        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            ts_model_ages_hiatusₚ = reverse(model_agesₚ[active_height_t_hiatus])
            llₚ += subsidence_ll(E₀, τ, ts_Sμ_hiatus, ts_Sσ_hiatus, ts_model_ages_hiatusₚ, subs_parametersₚ)/(length(ts_Sμ))

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
                beta_t0dist[:,n÷sieve] .= subs_parameters
                #predicted_ages[n÷sieve] = beta_t0dist[2]+τ*log(1-(Sμ[?]*pi)/(E₀*beta_t0dist[1]*sin(pi/beta_t0dist[2])))
                if mod(n, sieve*4000) ==0
                    ts_model_heights = Sμ[6:86]
                    plot!(h, model_ages[active_height_t], reverse(ts_model_heights), label="", alpha=0.5, yflip = true, xflip = true)
                    yₛ = (3165.647578*subs_parameters[1]/pi)*sin(pi/subs_parameters[1]).*(1 .-exp.(-(subs_parameters[2] .-model_ages[active_height_t_hiatus])./50))
                    # y should be a smooth curve determined by the subsidence model with the current t0, beta, etc.
                    plot!(h, model_ages[active_height_t_hiatus], yₛ, label="", alpha=0.5)
                    scatter!(h, smpl.Age, [223.6125, 580.5430, 840.8018, 975.1254, 1030.9710, 1072.7511], label="data",seriestype=:scatter,color=:black)
                end
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

        savefig(h, "SubsidenceModelAgeComparison_hiatus.pdf")

    # Crop the result
        agedist = agedist[active_height_t,:]

        subsmdl = SubsidenceStratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            nanpctile(agedist,2.5,dim=2), # 2.5th percentile
            nanpctile(agedist,97.5,dim=2), # 97.5th percentile
            nanmean(beta_t0dist[1,:],dim=1), # Mean beta
            nanstd(beta_t0dist[1,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[1,:],dim=1), # Median beta
            nanpctile(beta_t0dist[1,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[1,:],97.5,dim=1), # 97.5th percentile
            nanmean(beta_t0dist[2,:],dim=1), # Mean T0
            nanstd(beta_t0dist[2,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[2,:],dim=1), # Median T0
            nanpctile(beta_t0dist[2,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[2,:],97.5,dim=1) # 97.5th percentile
        )

    return subsmdl, agedist, lldist, beta_t0dist, lldist_burnin
end


# Part 2c: Modified StratMetropolis for extensional basins with hiata (when we know the duration of hiata relatively precisely)

function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, model_strat_heights, Sμ, Sσ, hiatus::HiatusData)

    # Run stratigraphic MCMC model, with hiata
    print("Generating stratigraphic age-depth model...\n")

    # Thermal subsidence model parameters
        y_litho= 125000
        ρ_mantle = 3330
        ρ_water = 1000
        αᵥ = 3.28*10^(-5)
        T_mantle = 1333
        τ = 50 #Myr
        E₀ = (4*y_litho*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water))

    # Stratigraphic age constraints. Type assertions for stability
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)
        model_heights = copy(-model_strat_heights[2:end]).*1000

    # Model configuration -- read from struct
        resolution = config.resolution
        bounding = config.bounding
        nsteps = config.nsteps
        burnin = config.burnin
        sieve = config.sieve

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = round((top-bottom)*bounding)
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [nanmean(Age_sigma)/10; Age_sigma; nanmean(Age_sigma)/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

    # STEP 1: calculate log likelihood of the modeled ages (and heights) in the initial proposal
        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end
        ll = normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

    # STEP 2: calculate log likelihood of the subsidence model parameters in the initial proposal
        # Define thermal subsidence model parameters - read from struct
        ideal_subs_parameters = therm.Param
        ideal_subs_parameters_sigma = therm.Sigma
        # Initial proposal for subsidence parameters - randomly pick a set of values from the distribution
        subs_parameters = [ideal_subs_parameters[1]+0.01, ideal_subs_parameters[2]-1]
        # Calculate log likelihood of this initial proposal
        ll += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)

    # STEP 3: calculate log likelihood for the fit of the thermal subsidence curve in the initial proposal
        ll += subsidence_ll(E₀, τ, Sμ, Sσ, model_ages, subs_parameters)/(length(model_strat_heights)-1)

        # Ensure there is only one effective hiatus at most for each height node
        closest_hiatus = findclosestabove((hiatus.Height::Array{Float64,1}),model_heights)
        closest_hiatus_unique = unique(closest_hiatus)
        Hiatus_height = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration_sigma = Array{Float64}(undef,size(closest_hiatus_unique))
        for i=1:length(closest_hiatus_unique)
            t = closest_hiatus.==closest_hiatus_unique[i]
            Hiatus_height[i] = mean((hiatus.Height::Array{Float64,1})[t])
            Hiatus_duration[i] = sum((hiatus.Duration::Array{Float64,1})[t])
            Hiatus_duration_sigma[i] = sqrt(sum((hiatus.Duration_sigma::Array{Float64,1})[t].^2))
        end

    # STEP 4: add log likelihood for hiatus duration
        duration = @. min(model_ages[closest_hiatus_unique - 1] - model_ages[closest_hiatus_unique], Hiatus_duration)
        ll += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration)

    # Preallocate variables for MCMC proposals
        llₚ=ll
        chosen_point=0
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        durationₚ = copy(duration)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)
        subs_parametersₚ = copy(subs_parameters)

    # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else

                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r

                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points that are still stratigraphically below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points that are still stratigraphically above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if Hiatus_height_uncert>0
                    #     closest_hiatus = findclosestabove(h.Height+randn(size(h.Height)).*Hiatus_height_uncert,heights)
                    # end
                    if any(closest_hiatus_unique.==chosen_point)
                        chosen_point -= 1
                        model_agesₚ[chosen_point] = model_ages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if model_agesₚ[i] < model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if model_agesₚ[i] > model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += subsidence_ll(E₀, τ, Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(model_strat_heights)-1)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = min(model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique], Hiatus_duration)
            llₚ += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(duration, durationₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
                # acceptancedist[i] = true
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

    # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        beta_t0dist = Array{Float64}(undef,2,nsteps)
        hiatusdist = Array{Float64}(undef,length(duration),nsteps)

    # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else

                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r

                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points that are still stratigraphically below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points that are still stratigraphically above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if Hiatus_height_uncert>0
                    #     closest_hiatus = findclosestabove(h.Height+randn(size(h.Height)).*Hiatus_height_uncert,heights)
                    # end
                    if any(closest_hiatus_unique.==chosen_point)
                        chosen_point -= 1
                        model_agesₚ[chosen_point] = model_ages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if model_agesₚ[i] < model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if model_agesₚ[i] > model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += subsidence_ll(E₀, τ, Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(model_strat_heights)-1)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = min(model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique], Hiatus_duration)
            llₚ += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(duration, durationₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
                beta_t0dist[:,n÷sieve] .= subs_parameters
                hiatusdist[n÷sieve] = duration[1]
                # predicted_ages[n÷sieve] = beta_t0dist[2]+τ*log(1-(Sμ[?]*pi)/(E₀*beta_t0dist[1]*sin(pi/beta_t0dist[2])))
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

    # Crop the result
        agedist = agedist[active_height_t,:]

        subsmdl = SubsidenceStratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            nanpctile(agedist,2.5,dim=2), # 2.5th percentile
            nanpctile(agedist,97.5,dim=2), # 97.5th percentile
            nanmean(beta_t0dist[1,:],dim=1), # Mean beta
            nanstd(beta_t0dist[1,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[1,:],dim=1), # Median beta
            nanpctile(beta_t0dist[1,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[1,:],97.5,dim=1), # 97.5th percentile
            nanmean(beta_t0dist[2,:],dim=1), # Mean T0
            nanstd(beta_t0dist[2,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[2,:],dim=1), # Median T0
            nanpctile(beta_t0dist[2,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[2,:],97.5,dim=1) # 97.5th percentile
        )
    return subsmdl, agedist, lldist, hiatusdist, beta_t0dist
end
=#

## --- End of File
