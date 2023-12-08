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

# To avoid allocations when indexing by a boolean vector
function copyat!(dest, t, src)
    @assert eachindex(t) == eachindex(src)
    iₙ = firstindex(dest)
    iₗ = lastindex(dest)
    @inbounds for iₛ in eachindex(src)
        if t[iₛ]
            dest[iₙ] = src[iₛ]
            iₙ += 1
            iₙ > iₗ && break
        end
    end
    return dest
end
function reversecopyat!(dest, t, src)
    @assert eachindex(t) == eachindex(src)
    i₀ = firstindex(dest)
    iₙ = lastindex(dest)
    @inbounds for iₛ in eachindex(src)
        if t[iₛ]
            dest[iₙ] = src[iₛ]
            iₙ -= 1
            iₙ < i₀ && break
        end
    end
    return dest
end

function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_heights, Sμ, Sσ, beta_ip, t0_ip;
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
        # Only include two-sided age constraints in fitting
        t = Age_Sidedness .== 0
        (a,b) = hcat(fill!(similar(Height[t]), 1), Height[t]) \ Age[t]
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

        # Initial proposal ll
        ll = normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

    # STEP 2: calculate log likelihood of the subsidence model parameters in the initial proposal
        # Define thermal subsidence model parameters
        ideal_subs_parameters = therm.Param
        ideal_subs_parameters_sigma = therm.Sigma

        # Initial proposal for subsidence parameters - randomly pick a set of values from the distribution
        subs_parameters = [ideal_subs_parameters[1]+beta_ip, ideal_subs_parameters[2]+t0_ip]
        # Add subsidence parameter likelihood to the initial proposal ll
        ll += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)

    # STEP 3: calculate log likelihood for the fit of the thermal subsidence curve in the initial proposal
        subsidence_height_t = subsidencebottom .<= model_heights .<= subsidencetop
        ts_model_ages = model_ages[subsidence_height_t] # N.B., subsidence stuff goes top down
        @info "Subsidence active for $(count(subsidence_height_t)) model horizons"

        heightconversion = if smpl.Height_Unit == "km"
            1
        elseif smpl.Height_Unit == "m"
            1000
        elseif smpl.Height_Unit == "cm"
            100_000
        elseif smpl.Height_Unit == "mm"
            1_000_000
        else
            1
        end
        equivalent_strat_height = if all(x->!(x>0), smpl.Height)
            -subsidence_strat_heights*heightconversion
        else
            sectionthickness = maximum(subsidence_strat_heights)*heightconversion
            sectionthickness .- subsidence_strat_heights*heightconversion
        end
        # Find the subsidence model horizons that match each age model horizon within the range where subsidence modelling is active
        closest_subsidence = findclosest(model_heights[subsidence_height_t], equivalent_strat_height)
        @info "Found $(length(closest_subsidence)) closest model horizons"
        # N.B. this will reverse ts_Sμ and ts_Sσ to match the (bottom-up) order in model_heights 
        ts_Sμ = Sμ[closest_subsidence]
        ts_Sσ = Sσ[closest_subsidence]

        # Add subsidence likelihood to initial proposal ll
        ll += subsidence_ll(E₀, τ, ts_Sμ, ts_Sσ, ts_model_ages, subs_parameters)/length(ts_Sμ)

    # Preallocate variables for MCMC proposals
        llₚ = ll
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)
        subs_parametersₚ = copy(subs_parameters)
        ts_model_agesₚ = copy(ts_model_ages)

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
            for i in eachindex(subs_parametersₚ, ideal_subs_parameters_sigma)
                subs_parametersₚ[i] += randn() * ideal_subs_parameters_sigma[i]
            end

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

            copyat!(ts_model_agesₚ, subsidence_height_t, model_agesₚ)
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
            for i in eachindex(subs_parametersₚ, ideal_subs_parameters_sigma)
                subs_parametersₚ[i] += randn() * ideal_subs_parameters_sigma[i]
            end

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

            copyat!(ts_model_agesₚ, subsidence_height_t, model_agesₚ)
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

function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_heights, Sμ, Sσ, hiatus_height, beta_ip, t0_ip)

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
        model_heights = copy(-subsidence_strat_heights[2:end]).*1000

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
        model_heights_all = copy(-subsidence_strat_heights).*1000
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

function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_heights, Sμ, Sσ, hiatus::HiatusData)

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
        model_heights = copy(-subsidence_strat_heights[2:end]).*1000

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
        ll += subsidence_ll(E₀, τ, Sμ, Sσ, model_ages, subs_parameters)/(length(subsidence_strat_heights)-1)

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
            llₚ += subsidence_ll(E₀, τ, Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(subsidence_strat_heights)-1)
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
            llₚ += subsidence_ll(E₀, τ, Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(subsidence_strat_heights)-1)
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
