## --- Part 2: Age-depth modeling

# Define the function that calculates the log likelihood of the thermal subsidence fit
# Method 1: Strat position for rift-drift transition is known (i.e., prior = age)
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
end
# Method 2: Strat position for rift-drift transition is unknown (i.e., prior = strat height)
#=
# Attemp 1: Only modified the subsidence_ll function slightly; majority of additional calculations are placed outside of this function
# This didn't work - takes too long to run
function subsidence_strat_ll(E₀, τ, model_St, model_St_sigma, model_t, beta)
    @assert eachindex(model_St) == eachindex(model_St_sigma) == eachindex(model_t)
    β = beta
    ll = zero(float(eltype(model_St_sigma)))
    @turbo for i ∈ eachindex(model_t)
        # Calculate subsidence_model_heights given this unique smooth thermal subsidence model at each age in model_t
        x = (E₀*β/pi)*sin(pi/β)*(1 -exp(-(model_t[1] -model_t[i])/τ))
        mu = model_St[i]
        sigma = model_St_sigma[i]
        # Turn that into a log likelihood using some age uncertainty of the curve
        ll -= (x-mu)*(x-mu) / (2*sigma*sigma)
    end
    return ll
end

# Attempt 2: Put all additional calculations in this new function subsidence_ll_strat (to avoid allocating a set of variables that changes length during each simulation)
# This didn't work - burnin would take 2.5 days...
function subsidence_ll_strat(E₀, τ, Sμ, Sσ, model_ages, model_h, s_top, equivalent_strat_height, beta_ts)
    ll = zero(float(eltype(Sσ)))
    β, h = beta_ts
    
    subsidence_height_t = h .<= model_h .<= s_top
    ts_model_ages = model_ages[subsidence_height_t]
    closest_subsidence = findclosest(model_h[subsidence_height_t], equivalent_strat_height)
    closest_subsidencebottom = findclosest(h, equivalent_strat_height)
    # N.B. this will reverse ts_Sμ and ts_Sσ to match the (bottom-up) order in model_heights 
    ts_Sμ = Sμ[closest_subsidence] .- Sμ[closest_subsidencebottom] # We only consider subsidence that happens after `subsidencebottom`
    ts_Sσ = Sσ[closest_subsidence]
    @assert eachindex(ts_Sμ) == eachindex(ts_Sσ) == eachindex(ts_model_ages)

    @turbo for i ∈ eachindex(ts_model_ages)
        # Calculate subsidence_model_heights given this unique smooth thermal subsidence model at each age in model_t
        x = (E₀*β/pi)*sin(pi/β)*(1 -exp(-(ts_model_ages[1] -ts_model_ages[i])/τ))
        mu = ts_Sμ[i]
        sigma = ts_Sσ[i]
        # Turn that into a log likelihood using some age uncertainty of the curve
        ll -= (x-mu)*(x-mu) / (2*sigma*sigma)
    end
    return ll/length(ts_Sμ)
    return log_likelihood
end
=#

# Attemp 3: Use indexing (instead of creating new variables) to loop through the thermal subsidence calculations
function subsidence_strat_ll(E₀, τ, St, St_sigma, ages, heights, converted_strat, top, bottom, beta)
    bottom_t = ages[bottom]
    bottom_St = St[findclosest(heights[bottom], converted_strat)]
    ll = zero(float(eltype(St)))
    #@assert eachindex(model_St) == eachindex(model_St_sigma) == eachindex(model_t)
    for i in bottom+1:top
        # Calculate subsidence_model_heights given this unique smooth thermal subsidence model at each age in model_t
        x = (E₀*beta/pi)*sin(pi/beta)*(1 -exp(-(bottom_t -ages[i])/τ))
        s_idx = findclosest(heights[i], converted_strat)
        mu = St[s_idx]-bottom_St
        sigma = St_sigma[s_idx]
        # Turn that into a log likelihood using some age uncertainty of the curve
        ll -= (x-mu)*(x-mu) / (2*sigma*sigma)
    end
    return ll
end

# Part 2a: Modified StratMetropolis for extensional basins - without hiatus
# Method 1: Strat position for rift-drift transition is known
"""
```julia
SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_depths, Sμ, Sσ, beta_ip, t0_ip;
    subsidencebotom=minimum(smpl.Height),
    subsidencetop=maximum(smpl.Height),
    lithosphere = Normal(125000, 100),
)
```
Runs the main SubsidenceChron.jl age-depth model routine given
a set of Gaussian age constraints specified in the `smpl` struct,
an age-depth model configuration specified in the `config` struct,
thermal subsidence parameters defined in the `therm` struct,
decompation and backstripping outputs in the form of `subsidence_strat_depths`, `Sμ`, and `Sσ`,
and prior estimates for stretching factor Beta (`beta_ip`) and time of thermal subsedence onset (`t0_ip`).

### Examples:
```julia
(subsmdl, agedist, lldist, beta_t0dist, lldist_burnin) = SubsidenceStratMetropolis(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ_corr, Beta_sigma/10, T0_sigma/10)
```
"""
function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_depths, Sμ, Sσ, beta_ip, t0_ip;
        subsidencebottom=minimum(smpl.Height),
        subsidencetop=maximum(smpl.Height),
        lithosphere = Normal(125000, 100), # Meters!
    )

    # Run stratigraphic MCMC model
    print("Generating stratigraphic age-depth model...\n")

    # Define thermal subsidence model parameters
        ρ_mantle = 3330
        ρ_water = 1000
        αᵥ = 3.28*10^(-5)
        T_mantle = 1333
        κ = 0.00804/100^2*60*60*24*365.25 # m^2/yr
        z = zₚ = rand(lithosphere)
        τ = z^2/(pi^2*κ*1e6) # Myr
        E₀ = (4*z*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water)) # Also meters!

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
        if count(t) < 2
            (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        else
            (a,b) = hcat(fill!(similar(Height[t]), 1), Height[t]) \ Age[t]
        end
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
        ll += normpdf_ll(mean(lithosphere), std(lithosphere), z)
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
        lldist_burnin = Array{Float64}(undef,burnin÷sieve)

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
            # Prevent beta from going below 1
            while subs_parametersₚ[1] < 1
                subs_parametersₚ[1] = subs_parameters[1] + randn() * ideal_subs_parameters_sigma[1]
            end
            zₚ = z .+ randn() * std(lithosphere)

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

            copyat!(ts_model_agesₚ, model_agesₚ, subsidence_height_t)
            τ = zₚ^2/(pi^2*κ*1e6) # Myr
            E₀ = (4*zₚ*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water)) # Also meters!
            llₚ += normpdf_ll(mean(lithosphere), std(lithosphere), zₚ)
            llₚ += subsidence_ll(E₀, τ, ts_Sμ, ts_Sσ, ts_model_agesₚ, subs_parametersₚ)/length(ts_Sμ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                z = zₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
                # acceptancedist[i] = true
            end
            if mod(n,sieve) == 0
                lldist_burnin[n÷sieve] = ll
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

    # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        lldist = Array{Float64}(undef,nsteps)
        zdist = Array{Float64}(undef,nsteps)
        agedist = Array{Float64}(undef,npoints,nsteps)
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
            # Prevent beta from going below 1
            while subs_parametersₚ[1] < 1
                subs_parametersₚ[1] = subs_parameters[1] + randn() * ideal_subs_parameters_sigma[1]
            end
            zₚ = z .+ randn() * std(lithosphere)

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

            copyat!(ts_model_agesₚ, model_agesₚ, subsidence_height_t)
            τ = zₚ^2/(pi^2*κ*1e6) # Myr
            E₀ = (4*zₚ*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water)) # Also meters!
            llₚ += normpdf_ll(mean(lithosphere), std(lithosphere), zₚ)
            llₚ += subsidence_ll(E₀, τ, ts_Sμ, ts_Sσ, ts_model_agesₚ, subs_parametersₚ)/length(ts_Sμ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                z = zₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                zdist[n÷sieve] = z
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
            nanpctile(t0dist,97.5,dim=1), # 97.5th percentile
            nanmean(zdist,dim=1), # Mean lithospheric thickness
            nanstd(zdist,dim=1), # Standard deviation
            nanmedian(zdist,dim=1), # Median lithospheric thickness
            nanpctile(zdist,2.5,dim=1), # 2.5th percentile
            nanpctile(zdist,97.5,dim=1), # 97.5th percentile
        )
    return subsmdl, agedist, lldist, beta_t0dist, lldist_burnin, zdist
end

# Part 2a - Method 2: Strat position for rift-drift transition is unknown
"""
```julia
SubsidenceStratMetropolis_Height(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_depths, Sμ, Sσ, beta_ip, t0_ip;
    subsidencetop=maximum(smpl.Height),
    lithosphere = Normal(125000, 100),
)
```
Runs the main SubsidenceChron.jl age-depth model routine given
a set of Gaussian age constraints specified in the `smpl` struct,
an age-depth model configuration specified in the `config` struct,
thermal subsidence parameters defined in the `therm` struct,
decompation and backstripping outputs in the form of `subsidence_strat_depths`, `Sμ`, and `Sσ`,
and prior estimates for stretching factor Beta (`beta_ip`) and time of thermal subsedence onset (`t0_ip`).

### Examples:
```julia
(subsmdl, agedist, lldist, beta_tsdist, lldist_burnin, zdist) = SubsidenceStratMetropolis_Height(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ_corr, Beta_sigma/10, T0_sigma/10)
```
"""
function SubsidenceStratMetropolis_Height(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, subsidence_strat_depths, Sμ, Sσ, beta_ip, ts_ip;
    subsidencetop=maximum(smpl.Height),
    lithosphere = Normal(125000, 100), # Meters!
    )

    # Run stratigraphic MCMC model
    print("Generating stratigraphic age-depth model...\n")

    # Define thermal subsidence model parameters
        ρ_mantle = 3330
        ρ_water = 1000
        αᵥ = 3.28*10^(-5)
        T_mantle = 1333
        κ = 0.00804/100^2*60*60*24*365.25 # m^2/yr
        z = zₚ = rand(lithosphere)
        τ = z^2/(pi^2*κ*1e6) # Myr
        E₀ = (4*z*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water)) # Also meters!

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
        if iszero(t) == true
            (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        else
            (a,b) = hcat(fill!(similar(Height[t]), 1), Height[t]) \ Age[t]
        end
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
        ll += normpdf_ll(mean(lithosphere), std(lithosphere), z)

    # STEP 2: calculate log likelihood of the subsidence model parameters in the initial proposal
        # Define thermal subsidence model parameters
        ideal_subs_parameters = therm.Param
        ideal_subs_parameters_sigma = therm.Sigma

        # Initial proposal for subsidence parameters - randomly pick a set of values from the distribution
        subs_parameters = [ideal_subs_parameters[1]+beta_ip, ideal_subs_parameters[2]+ts_ip]
        # Add subsidence parameter likelihood to the initial proposal ll
        ll += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)

    # STEP 3: calculate log likelihood for the fit of the thermal subsidence curve in the initial proposal
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
        
        # Attempt 1: tried to only modified the subsidence_ll function slightly; majority of additional calculations are placed outside of this function
        # This didn't work - would take too long 
        #=
        subsidence_height_t = subs_parameters[2] .<= model_heights .<= subsidencetop
        ts_model_ages = model_ages[subsidence_height_t] # N.B., subsidence stuff goes top down
        # Find the subsidence model horizons that match each age model horizon within the range where subsidence modelling is active
        closest_subsidence = findclosest(model_heights[subsidence_height_t], equivalent_strat_height)
        #closest_subsidencebottom = findclosest(subs_parameters[2], equivalent_strat_height)
        # N.B. this will reverse ts_Sμ and ts_Sσ to match the (bottom-up) order in model_heights 
        ts_Sμ = Sμ[closest_subsidence] .- Sμ[closest_subsidence][1] # We only consider subsidence that happens after `subsidencebottom`
        ts_Sσ = Sσ[closest_subsidence]
        # Add subsidence likelihood to initial proposal ll
        ll += subsidence_strat_ll(E₀, τ, ts_Sμ, ts_Sσ, ts_model_ages, subs_parameters[1])/length(ts_Sμ)
        =#

        # Attempt 2:
        # This didn't work - would take too long 
        # ll += subsidence_ll_strat(E₀, τ, Sμ, Sσ, model_ages, model_heights, subsidencetop, equivalent_strat_height, subs_parameters)

        # Attempt 3: try to speed things up by using indexing to loop through the thermal subsidence calculations
        top_idx = findclosest(subsidencetop, model_heights)
        bottom_idx = findclosest(subs_parameters[2], model_heights)
        ll += subsidence_strat_ll(E₀, τ, Sμ, Sσ, model_ages, model_heights, equivalent_strat_height, top_idx, bottom_idx, subs_parameters[1])/(top_idx-bottom_idx+1)

    # Preallocate variables for MCMC proposals
        llₚ = ll
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)
        subs_parametersₚ = copy(subs_parameters)
        #ts_model_agesₚ = copy(ts_model_ages)

    # Run burnin
        # acceptancedist = fill(false,burnin)
        lldist_burnin = Array{Float64}(undef,burnin÷sieve)

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
            # Prevent beta from going below 1
            while subs_parametersₚ[1] < 1
                subs_parametersₚ[1] = subs_parameters[1] + randn() * ideal_subs_parameters_sigma[1]
            end
            zₚ = z .+ randn() * std(lithosphere)

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

            τ = zₚ^2/(pi^2*κ*1e6) # Myr
            E₀ = (4*zₚ*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water)) # Also meters!
            llₚ += normpdf_ll(mean(lithosphere), std(lithosphere), zₚ)
            
            #=
            # Both Attempt 1 and 2 take too long to run
            # Attempt 1
            subsidence_height_t = subs_parametersₚ[2] .<= model_heights .<= subsidencetop
            ts_model_agesₚ = model_agesₚ[subsidence_height_t]
            closest_subsidenceₚ = findclosest(model_heights[subsidence_height_t], equivalent_strat_height)
            #closest_subsidencebottomₚ = findclosest(subs_parametersₚ[2], equivalent_strat_height)
            # N.B. this will reverse ts_Sμ and ts_Sσ to match the (bottom-up) order in model_heights 
            ts_Sμₚ = Sμ[closest_subsidenceₚ] .- Sμ[closest_subsidenceₚ][1] # We only consider subsidence that happens after `subsidencebottom`
            ts_Sσₚ = Sσ[closest_subsidenceₚ]
            llₚ += subsidence_strat_ll(E₀, τ, ts_Sμₚ, ts_Sσₚ, ts_model_agesₚ, subs_parametersₚ[1])/length(ts_Sμₚ)

            # Attempt 2:
            # This didn't work - would take too long 
            # llₚ += subsidence_ll_strat(E₀, τ, Sμ, Sσ, model_agesₚ, model_heights, subsidencetop, equivalent_strat_height, subs_parametersₚ)
            =#

            # Attempt 3: 
            bottom_idx = findclosest(subs_parametersₚ[2], model_heights)
            llₚ += subsidence_strat_ll(E₀, τ, Sμ, Sσ, model_agesₚ, model_heights, equivalent_strat_height, top_idx, bottom_idx, subs_parametersₚ[1])/(top_idx-bottom_idx+1)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                z = zₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
                # acceptancedist[i] = true
            end
            if mod(n,sieve) == 0
                lldist_burnin[n÷sieve] = ll
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

    # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        lldist = Array{Float64}(undef,nsteps)
        zdist = Array{Float64}(undef,nsteps)
        agedist = Array{Float64}(undef,npoints,nsteps)
        beta_tsdist = Array{Float64}(undef,2,nsteps)

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
            # Prevent beta from going below 1
            while subs_parametersₚ[1] < 1
                subs_parametersₚ[1] = subs_parameters[1] + randn() * ideal_subs_parameters_sigma[1]
            end
            zₚ = z .+ randn() * std(lithosphere)

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

            τ = zₚ^2/(pi^2*κ*1e6) # Myr
            E₀ = (4*zₚ*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water)) # Also meters!
            llₚ += normpdf_ll(mean(lithosphere), std(lithosphere), zₚ)

            #=
            # Attempt 1
            subsidence_height_t = subs_parametersₚ[2] .<= model_heights .<= subsidencetop
            ts_model_agesₚ = model_agesₚ[subsidence_height_t]
            closest_subsidenceₚ = findclosest(model_heights[subsidence_height_t], equivalent_strat_height)
            #closest_subsidencebottomₚ = findclosest(subs_parametersₚ[2], equivalent_strat_height)
            # N.B. this will reverse ts_Sμ and ts_Sσ to match the (bottom-up) order in model_heights 
            ts_Sμₚ = Sμ[closest_subsidenceₚ] .- Sμ[closest_subsidenceₚ][1] # We only consider subsidence that happens after `subsidencebottom`
            ts_Sσₚ = Sσ[closest_subsidenceₚ]
            llₚ += subsidence_strat_ll(E₀, τ, ts_Sμₚ, ts_Sσₚ, ts_model_agesₚ, subs_parametersₚ[1])/length(ts_Sμₚ)

            # Attempt 2:
            # This didn't work - would take too long 
            # llₚ += subsidence_ll_strat(E₀, τ, Sμ, Sσ, model_agesₚ, model_heights, subsidencetop, equivalent_strat_height, subs_parametersₚ)
            
            #copyat!(ts_model_agesₚ, subsidence_height_t, model_agesₚ)
            #llₚ += subsidence_ll(E₀, τ, ts_Sμ, ts_Sσ, ts_model_agesₚ, subs_parametersₚ)/length(ts_Sμ)
            =#
            
            # Attempt 3: 
            bottom_idx = findclosest(subs_parametersₚ[2], model_heights)
            llₚ += subsidence_strat_ll(E₀, τ, Sμ, Sσ, model_agesₚ, model_heights, equivalent_strat_height, top_idx, bottom_idx, subs_parametersₚ[1])/(top_idx-bottom_idx+1)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                z = zₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                zdist[n÷sieve] = z
                agedist[:,n÷sieve] .= model_ages
                beta_tsdist[:,n÷sieve] .= subs_parameters
                #predicted_ages[n÷sieve] = beta_t0dist[2]+τ*log(1-(Sμ[?]*pi)/(E₀*beta_t0dist[1]*sin(pi/beta_t0dist[2])))
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

    # Crop the result
        agedist = agedist[active_height_t,:]
        betadist = beta_tsdist[1,:]
        tsdist = beta_tsdist[2,:]
        subsmdl = SubsidenceStratAgeModel_Height(
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
            nanmean(tsdist,dim=1), # Mean height for active rifting to thermal subsidence transition
            nanstd(tsdist,dim=1), # Standard deviation
            nanmedian(tsdist,dim=1), # Median height for active rifting to thermal subsidence transition
            nanpctile(tsdist,2.5,dim=1), # 2.5th percentile
            nanpctile(tsdist,97.5,dim=1), # 97.5th percentile
            nanmean(zdist,dim=1), # Mean lithospheric thickness
            nanstd(zdist,dim=1), # Standard deviation
            nanmedian(zdist,dim=1), # Median lithospheric thickness
            nanpctile(zdist,2.5,dim=1), # 2.5th percentile
            nanpctile(zdist,97.5,dim=1), # 97.5th percentile
        )
    return subsmdl, agedist, lldist, beta_tsdist, lldist_burnin, zdist
end



## --- Stratigraphic MCMC model with hiatus # # # # # # # # # # # # # # # # #

