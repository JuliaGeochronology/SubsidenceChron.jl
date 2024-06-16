## For subsidence modeling in extensional basins
## --- Load required pacages, install SubsidenceChron if required

    using SubsidenceChron
    using Plots

## --- Part 1: Decompaction and Backstripping

    # # # # # # # # # # # Enter stratigraphic information here! # # # # # # # # # # # #
    # Import the data file (.csv)
    data_csv = importdataset("examples/Svalbard_highres.csv",',', importas=:tuple)
    # Obtain stratigraphic info from the data file
    nLayers = length(data_csv.Thickness)
    strat = StratData(nLayers)
    strat.Lithology          = data_csv.Lithology
    strat.Thickness         .= data_csv.Thickness
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # OPTIONAL - Enter paleo water depth information here! # # # # # # # #
    # Import the data file (.csv)
    wd_csv = importdataset("examples/Svalbard_SequenceStrat.csv", ',', importas=:tuple)
    # Obtain paleo water depth info from the data file
    wd_nLayers = length(wd_csv.Thickness)
    wd = WaterDepth(wd_nLayers)
    wd.DepthID    = wd_csv.Type
    wd.Thickness .= wd_csv.Thickness
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # Configure MC model here! # # # # # # # # # #
    # Number of MCMC simulations
    nsims = 5000
    # Resolution for model horizons (in meters)
    res = 10.0 # meters
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Run the decompaction and backstripping MC model
    # (wd input is optional)
    @time (Sₜ, Sμ, Sσ, subsidence_strat_depths) = DecompactBackstrip(strat, wd, nsims, res)

    #= Code for storing and reading decompaction + backstripping results - will be useful when testing the age-depth modeling part of the model
    # Store results
    using DelimitedFiles
    writedlm("St_highres_5000sims.csv", Sₜ, ",")
    writedlm("St_mu_highres_5000sims.csv", Sμ)
    writedlm("St_sigma_highres_5000sims.csv", Sσ)
    writedlm("St_highres.txt", Sₜ)
    # Read results
    Sₜ = readdlm("St_highres_5000sims.csv", ',', Float64)
    Sμ = readdlm("St_mu_highres_5000sims.csv", ',', Float64)
    Sσ = readdlm("St_sigma_highres_5000sims.csv", ',', Float64)
    =#

    # Plot results - tectonic subsidence in comparison with present day stratigraphic heights
    p1 = plot(framestyle=:box, 
        xlabel="Stratigraphy removed [m]", 
        ylabel="Subsidence [m]", 
        fg_color_legend=:white, 
        xflip = true, 
        yflip = true,
    )
    plot!(p1, subsidence_strat_depths, Sₜ[:,2:end], alpha = 0.01, label = "", color = "blue")
    plot!(p1, subsidence_strat_depths, Sμ, alpha = 1, label = "Tectonic subsidence", color = "blue")
    savefig(p1, "DecompactBackstrip.pdf")
    display(p1)


## --- Part 2a: Age-depth modeling

    # # # # # # # # # # # Enter age constraint (sample) information here! # # # # # # # # # # # #
    # Input the number of samples we wish to model (must match below)
    nSamples = 4
    # Make an instance of a Chron Section object for nSamples
    smpl = ChronAgeData(nSamples)
    smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
    smpl.Age          .= [879.91,   791.1,    737.5,   717] # Measured ages
    smpl.Age_sigma    .= [  0.63,    2.45,      4.8,   0.4] # Measured 1-σ uncertainties
    smpl.Height       .= [ -2138,   -1223,     -100,     0] # Depths below surface should be negative
    smpl.Height_sigma .= [   10,      10,       10,    10] # Usually assume little or no sample height uncertainty
    smpl.Age_Sidedness .= [  -1,       0,        0,    +1] # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Age_Unit = "Ma" # Unit of measurement for ages
    smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

    # IMPORTANT: smpl.Height must increase with increasing stratigraphic height
    # -- i.e., stratigraphically younger samples must be more positive. For this
    # reason, it is convenient to represent depths below surface as negative
    # numbers.
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # Enter thermal subsidence parameter priors here! # # # # # # # # # #
    # Enter initial guesses for the beta factor and thermal subsidence onset age and their uncertainties
    Beta = 1.3
    Beta_sigma = 0.2
    T0 = 816
    T0_sigma = 50

    therm = ThermalSubsidenceParameters()
    therm.Param = [Beta, T0]
    therm.Sigma = [Beta_sigma, T0_sigma]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # Configure MCMC model here! # # # # # # # # # #
    # Configure the stratigraphic MCMC model
    config = StratAgeModelConfiguration()
    # If you in doubt, you can probably leave these parameters as-is
    config.resolution = res # Same units as sample height. Smaller is slower!
    config.bounding = 1.0 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
    (bottom, top) = extrema(smpl.Height)
    npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
    config.nsteps = 100000 # Number of steps to run in distribution MCMC
    config.burnin = 10000*npoints_approx # Number to discard
    config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Let's see if adding just a little bit of systematic error on tectonic subsidence solves the ll = -Inf problem
    # Currently the systematic error is set to be = 1% of the sigma for the tectonic subsidence of the first (oldest) layer
    Sσ_corr = Sσ .+ ((Sμ[end-1]-Sμ[end])/100)

## --- Stratigraphic MCMC model without hiatus

    #Run the model
    (subsmdl, agedist, lldist, beta_t0dist, lldist_burnin) = SubsidenceStratMetropolis(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ_corr, Beta_sigma/10, T0_sigma/10)

    #= Code for storing and reading age-depth model results (optional)
    # Store and read results
    writedlm("agedist.csv", agedist, ",")
    writedlm("lldist.csv", lldist)
    writedlm("beta_t0dist.csv", beta_t0dist, ",")
    writedlm("lldist_burnin.csv", lldist_burnin)

    agedist = readdlm("agedist.csv", ',', Float64)
    lldist = readdlm("lldist3.csv", ',', Float64)
    beta_t0dist = readdlm("beta_t0dist.csv", ',', Float64)
    lldist_burnin = readdlm("lldist_burnin.csv", ',', Float64)
    =#

    # Plot 1: Age-depth model (mean and 95% confidence interval for both model and data)
    hdl =  plot(framestyle=:box, xlabel="Age [$(smpl.Age_Unit)]", ylabel="Height ($(smpl.Height_Unit))",)
    plot!(hdl, [subsmdl.Age_025CI; reverse(subsmdl.Age_975CI)],[subsmdl.Height; reverse(subsmdl.Height)], fill=(round(Int,minimum(subsmdl.Height)),0.5,:blue), label="model")
    plot!(hdl, subsmdl.Age, subsmdl.Height, linecolor=:blue, label="", fg_color_legend=:white) # Center line
    t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(2*smpl.Age_sigma[t]),yerror=(2*smpl.Height_sigma[t]),label="data",seriestype=:scatter,color=:black)
    t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(2*smpl.Age_sigma[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
    any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),2*smpl.Age_sigma[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
    any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    savefig(hdl,"Svalbard_AgeDepth.pdf")
    display(hdl)

    # Add stratigraphy
    formation, top, bottom = find_formation_depths(data_csv.Formation, data_csv.Thickness)
    xl = collect(xlims())
    xmin = minimum(xl)
    xminshifted = xmin - nanrange(xl)/40
    ymax = maximum(subsmdl.Height)
    for i in eachindex(formation)
        x = [xmin, xminshifted]
        plot!(x, fill(ymax-top[i],2), fillto=(ymax-bottom[i]), label="")
        annotate!((xmin+xminshifted)/2, ymax-(top[i]+bottom[i])/2, text(formation[i], 6, :left, :bottom, rotation=-45))
    end
    xlims!(xminshifted, maximum(xl))
    savefig(hdl,"Svalbard_AgeDepth_strat.pdf")
    display(hdl)

    # Plot 2: Posterior distributions of beta
    post_beta = plot(framestyle=:box, xlabel="Beta [unitless]", ylabel="Probability density",)
    histogram!(post_beta, beta_t0dist[1,:], color="black", linecolor=nothing, alpha = 0.5, nbins=50, normalized=true, label="")
    vline!(post_beta, [subsmdl.Beta_Median], linecolor = "black", linestyle=:dot, linewidth = 3, label = "median")
    savefig(post_beta, "Svalbard_PosteriorBeta.pdf")
    display(post_beta)

    # Plot 3: Posterior distributions of t₀
    post_t0 = plot(framestyle=:box, xlabel="t0 Age [$(smpl.Age_Unit)]", ylabel="Probability density",)
    histogram!(post_t0, beta_t0dist[2,:], color="black", linecolor=nothing, alpha = 0.5, nbins=50, normalized=true, label="")
    vline!(post_t0, [subsmdl.T0_Median], linecolor = "black", linestyle=:dot, linewidth = 3, label = "median")
    savefig(post_t0, "Svalbard_PosteriorT0.pdf")
    display(post_t0)

    # Plot 4: Interpolate results for target horizons
        #= For Svalbard:
            -2138 = base of section (base of Akademikerbreen Gp)
            -1694 = onset of Bitter Springs CIE
            -1275 = terminatio of Bitter Springs CIE
            -1013 = base of Draken Fm
            -176 = base of Russoya Mb
            -56 = onset of Russoya CIE
            -1 = top of section (right below diamictite)
        =#
    target_height = [-2138, -1694, -1275, -1013, -176, -56, -1]
    interpolated_distribution = Array{Float64}(undef, length(target_height), size(agedist,2))
    # Interpolate between model horizons
    for i=1:size(agedist,2)
        interpolated_distribution[:,i] = linterp1s(subsmdl.Height,agedist[:,i],target_height)
    end
    # Calculate summary statistics
    predicted_means = nanmean(interpolated_distribution, dims=2)
    predicted_medians = nanmedian(interpolated_distribution, dims=2)
    predicted_025CI = nanpctile(interpolated_distribution, 2.5, dims=2)
    predicted_975CI = nanpctile(interpolated_distribution, 97.5, dims=2)
    # Plotting (full distribution with mean as a dashed line)
    age_interp = plot(framestyle=:box,
        xlabel="Age [$(smpl.Age_Unit)] for the onset of Bitter Springs CIE", 
        ylabel="Probability density",
    )
    histogram!(age_interp, interpolated_distribution[2,:], nbins=50, normalized=true, label="")
    vline!(age_interp, predicted_medians[2:2], linecolor = "black", linestyle=:dot, linewidth = 3, label="median: $(round(predicted_medians[2], digits=1))")
    savefig(age_interp, "Svalbard_InterpolatedAge_BSA_Onset.pdf")
    display(age_interp)

## --- Plot subsidence versus age

    # Tectonic subsidence as a functon of age
    p1 = plot(framestyle=:box,
        xlabel="Age [$(smpl.Age_Unit)]",
        ylabel="Subsidence [m]", 
        xflip = true,
        yflip = true, 
        legend = :bottomleft,
    )
    mean_age = linterp1s(subsmdl.Height, subsmdl.Age, maximum(subsmdl.Height).-subsidence_strat_depths)
    plot!(p1, mean_age, Sμ, alpha = 1, label = "Tectonic subsidence", color = "blue")
    plot!(p1, mean_age, Sₜ, alpha = 0.01, label = "", color = "blue", fg_color_legend=:white)

    # Add stratigraphy
    formation, top, bottom = find_formation_depths(data_csv.Formation, data_csv.Thickness)
    yl = ylims()
    ymin = min(yl...,)
    for i in eachindex(formation)
        x = linterp1s(subsmdl.Height, subsmdl.Age, maximum(subsmdl.Height).-[top[i], bottom[i]])
        plot!(x, fill(ymin, 2), fillto=0.0, label="")
        annotate!(sum(x)/2, 0, text(formation[i], 7, :left, :top, rotation=-45))
    end
    ylims!(yl)

    savefig(p1, "DecompactBackstripAge.pdf")
    display(p1)


## --- Plot subsidence versus age, including age uncertainty

    # Tectonic subsidence as a functon of age
    p1 = plot(framestyle=:box,
        xlabel="Age [$(smpl.Age_Unit)]",
        ylabel="Subsidence [m]", 
        xflip = true,
        yflip = true, 
    )
    xq = maximum(subsmdl.Height).-subsidence_strat_depths
    plot!(p1, linterp1s(subsmdl.Height, subsmdl.Age, xq), Sμ, alpha = 1, label = "Tectonic subsidence", color = "blue")
    ages = zeros(length(subsidence_strat_depths), nsims)
    for i in axes(ages,2)
        linterp1s!(view(ages, :, i), subsmdl.Height, view(agedist, :, rand(1:size(agedist,2))), xq)
    end
    plot!(p1, ages, Sₜ, alpha = 0.01, label = "", color = "blue", fg_color_legend=:white)

    # # Add stratigraphy
    # formation, top, bottom = find_formation_depths(data_csv.Formation, data_csv.Thickness)
    # yl = ylims()
    # ymin = min(yl...,)
    # for i in eachindex(formation)
    #     x = linterp1s(subsmdl.Height, subsmdl.Age, maximum(subsmdl.Height).-[top[i], bottom[i]])
    #     plot!(x, fill(ymin, 2), fillto=0.0, label="")
    #     annotate!(sum(x)/2, 0, text(formation[i], 7, :left, :top, rotation=-45))
    # end
    # ylims!(yl)

    savefig(p1, "DecompactBackstripAgeUncert.pdf")
    display(p1)

## --- Subsidence rate

    # Set bin width and spacing
    binoverlap = 10
    binwidth = round(Int, round(nanrange(subsmdl.Age)/10,sigdigits=1)) # Can also set manually, commented out below
    # binwidth = 1 # Same units as smpl.Age

    agebinedges = collect(minimum(subsmdl.Age):binwidth/binoverlap:maximum(subsmdl.Age))
    agebincenters = (agebinedges[1:end-binoverlap] + agebinedges[1+binoverlap:end])/2

    # Calculate rates for the stratigraphy of each markov chain step
    dsdt_dist = zeros(length(agebincenters), nsims)
    @time for i=1:nsims
        subsidence = linterp1s(view(ages, :, i), view(Sₜ, :, i), agebinedges, extrapolate=NaN)
        dsdt_dist[:,i] .= (subsidence[1:end-binoverlap] - subsidence[1+binoverlap:end]) ./ binwidth
    end

    # Find mean and 1-sigma (68%) CI
    dsdt = nanmean(dsdt_dist,dim=2)
    dsdt_50p = nanmedian(dsdt_dist,dim=2)
    dsdt_16p = nanpctile(dsdt_dist,15.865,dim=2) # Minus 1-sigma (15.865th percentile)
    dsdt_84p = nanpctile(dsdt_dist,84.135,dim=2) # Plus 1-sigma (84.135th percentile)

    # Plot results
    hdl = plot(framestyle=:box,
        xlabel="Age [$(smpl.Age_Unit)]", 
        ylabel="Subsidence [$(smpl.Height_Unit) / $(smpl.Age_Unit) over $binwidth $(smpl.Age_Unit)]", 
        fg_color_legend=:white,
        legend=:bottomleft,
        xflip=true,
    )
    plot!(hdl, agebincenters, dsdt, label="Mean", color=:black, linewidth=2)
    plot!(hdl, [agebincenters; reverse(agebincenters)],[dsdt_16p; reverse(dsdt_84p)], fill=(0,0.2,:darkblue), linealpha=0, label="68% CI")
    for lci in 20:5:45
        dsdt_lp = nanpctile(dsdt_dist,lci,dim=2)
        dsdt_up = nanpctile(dsdt_dist,100-lci,dim=2)
        plot!(hdl,[agebincenters; reverse(agebincenters)],[dsdt_lp; reverse(dsdt_up)], fill=(0,0.2,:darkblue), linealpha=0, label="")
    end
    plot!(hdl,agebincenters,dsdt_50p, label="Median", color=:grey, linewidth=1)

    # # Add stratigraphy
    # formation, top, bottom = find_formation_depths(data_csv.Formation, data_csv.Thickness)
    # yl = ylims()
    # ymax = max(yl...,)
    # for i in eachindex(formation)
    #     x = linterp1s(subsmdl.Height, subsmdl.Age, maximum(subsmdl.Height).-[top[i], bottom[i]])
    #     plot!(x, fill(ymax, 2), fillto=ymax+nanrange(yl)/25, label="")
    #     annotate!(sum(x)/2, ymax, text(formation[i], 7, :left, :top, rotation=-45))
    # end
    # ylims!(min(yl...),ymax+nanrange(yl)/25)
    
    savefig(hdl,"SubsidenceRateModelCI.pdf")
    display(hdl)

## --- End of File
