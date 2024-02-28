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

function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_depths, Sμ, Sσ, beta_ip, t0_ip;
        subsidencebottom=minimum(smpl.Height),
        subsidencetop=maximum(smpl.Height),
        y_lithosphere= 125000, # Meters!
        τ = 50 #Myr
    )

    # Run stratigraphic MCMC model
    print("Generating stratigraphic age-depth model...\n")

    # Define thermal subsidence model parameters
        ρ_mantle = 3330
        ρ_water = 1000
        αᵥ = 3.28*10^(-5)
        T_mantle = 1333
        E₀ = (4*y_lithosphere*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water)) # Also meters!

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
            0.001
        elseif smpl.Height_Unit == "m"
            1.
        elseif smpl.Height_Unit == "cm"
            100.
        elseif smpl.Height_Unit == "mm"
            1_000.
        else
            1.
        end
        equivalent_strat_height = if all(x->!(x>0), smpl.Height)
            -subsidence_strat_depths*heightconversion
        else
            sectionthickness = maximum(subsidence_strat_depths)*heightconversion
            sectionthickness .- subsidence_strat_depths*heightconversion
        end
        # Find the subsidence model horizons that match each age model horizon within the range where subsidence modelling is active
        closest_subsidence = findclosest(model_heights[subsidence_height_t], equivalent_strat_height)
        closest_subsidencebottom = findclosest(subsidencebottom, equivalent_strat_height)

        @info "Found $(length(closest_subsidence)) closest model horizons"
        # N.B. this will reverse ts_Sμ and ts_Sσ to match the (bottom-up) order in model_heights 
        ts_Sμ = Sμ[closest_subsidence] .- Sμ[closest_subsidencebottom] # We only consider subsidence that happens after `subsidencebottom`
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
        betadist = beta_t0dist[1,:]
        t0dist = beta_t0dist[2,:]
        subsmdl = SubsidenceStratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            nanpctile(agedist,2.5,dim=2), # 2.5th percentile
            nanpctile(agedist,97.5,dim=2), # 97.5th percentile
            nanmean(betadist,dim=1), # Mean beta
            nanstd(betadist,dim=1), # Standard deviation
            nanmedian(betadist,dim=1), # Median beta
            nanpctile(betadist,2.5,dim=1), # 2.5th percentile
            nanpctile(betadist,97.5,dim=1), # 97.5th percentile
            nanmean(t0dist,dim=1), # Mean T0
            nanstd(t0dist,dim=1), # Standard deviation
            nanmedian(t0dist,dim=1), # Median T0
            nanpctile(t0dist,2.5,dim=1), # 2.5th percentile
            nanpctile(t0dist,97.5,dim=1) # 97.5th percentile
        )
    return subsmdl, agedist, lldist, beta_t0dist, lldist_burnin
end

## --- Stratigraphic MCMC model with hiatus # # # # # # # # # # # # # # # # #

