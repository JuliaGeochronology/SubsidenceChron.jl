## For subsidence modeling in extensional basins
## --- Load required pacages, install SubsidenceChron if required

    using SubsidenceChron
    using StatGeochem, Distributions, Plots

## There are currently two ways to deal with the rift-drift transition being in the middle of a succession
## RUN EITHER METHOD 1 OR METHOD 2

# # # # # # # # # # # # METHOD 1 # # # # # # # # # # # #
## --- Part 1: Decompaction and Backstripping

    # # # # # # # # # # # Enter stratigraphic information here! # # # # # # # # # # # #
    # Import the data file (.csv)
    data_csv = importdataset("examples/Test_DB_AllenAllenEx58.csv",',')
    # Obtain stratigraphic info from the data file
    nLayers = length(data_csv["Thickness"])
    strat = StratData(nLayers)
    strat.Lithology          = data_csv["Lithology"]
    strat.Thickness         .= data_csv["Thickness"]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # OPTIONAL - Enter paleo water depth information here! # # # # # # # #
    # Import the data file (.csv)
    #wd_csv = importdataset("examples/Svalbard_SequenceStrat.csv", ',')
    # Obtain paleo water depth info from the data file
    #wd_nLayers = length(wd_csv["Thickness"])
    #wd = WaterDepth(wd_nLayers)
    #wd.DepthID    = wd_csv["Type"]
    #wd.Thickness .= wd_csv["Thickness"]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # Configure MC model here! # # # # # # # # # #
    # Number of MCMC simulations
    nsims = 500
    # Resolution for model horizons (in meters)
    res = 1.0 # meters
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Run the decompaction and backstripping MC model
    # (wd input is optional)
    @time (Sₜ, Sμ, Sσ, subsidence_strat_depths) = DecompactBackstrip(strat, #=wd, =#nsims, res)

    #=
    # Code for storing and reading decompaction + backstripping results - will be useful when testing the age-depth modeling part of the model
    # Store results
    using DelimitedFiles
    writedlm("St_highres_StratTest.csv", Sₜ, ",")
    writedlm("St_mu_highres_StratTest.csv", Sμ)
    writedlm("St_sigma_highres_StratTest.csv", Sσ)
    writedlm("St_highres.txt", Sₜ)
    # Read results
    Sₜ = readdlm("St_highres_StratTest.csv", ',', Float64)
    Sμ = readdlm("St_mu_highres_StratTest.csv", ',', Float64)
    Sσ = readdlm("St_sigma_highres_StratTest.csv", ',', Float64)
    =#

    # Plot results - tectonic subsidence in comparison with present day stratigraphic heights
    p1 = plot(Sμ, alpha = 1, yflip = true, xflip = true, label = "Tectonic subsidence", color = "blue")
    #plot!(p1, reverse(subsidence_strat_depths), yflip = true, label = "Present-day thickness", color = "red")
    plot!(p1, Sₜ[:,2:end], alpha = 0.01, label = "", yflip = true, color = "blue", fg_color_legend=:white)
    #savefig(p1, "DecompactBackstrip_higherres.pdf")


## --- Part 2: Age-depth modeling

    # # # # # # # # # # # Enter age constraint (sample) information here! # # # # # # # # # # # #
    # Input the number of samples we wish to model (must match below)
    nSamples = 11
    # Make an instance of a Chron Section object for nSamples
    smpl = ChronAgeData(nSamples)
    smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", "Sample 6", "Sample 7", "Sample 8", "Sample 9", "Sample 10", "Sample 11") # Et cetera
    smpl.Age          .= [   260,     245,     210,     160,     145,     125,     100,      80,      55,      45,       0] # Measured ages
    smpl.Age_sigma    .= [     2,       2,       2,       2,       2,       2,       2,       2,       2,       2,       2] # Measured 1-σ uncertainties
    smpl.Height       .= [ -5400,   -5000,   -4250,   -4000,   -3600,   -3400,   -2500,   -1200,    -450,    -200,       0] # Depths below surface should be negative
    smpl.Height_sigma .= [     1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1] # Usually assume little or no sample height uncertainty
    smpl.Age_Sidedness .= [    0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0] # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Age_Unit = "Ma" # Unit of measurement for ages
    smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

    # IMPORTANT: smpl.Height must increase with increasing stratigraphic height
    # -- i.e., stratigraphically younger samples must be more positive. For this
    # reason, it is convenient to represent depths below surface as negative
    # numbers.
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # Enter thermal subsidence parameter priors here! # # # # # # # # # #
    # Enter initial guesses for the beta factor and stratigraphic height for the onset of thermal subsidence, as well as their uncertainties
    Beta = 2
    Beta_sigma = 0.5
    TS_Height = -3400
    TS_Height_sigma = 200

    therm = ThermalSubsidenceParameters()
    therm.Param = [Beta, TS_Height]
    therm.Sigma = [Beta_sigma, TS_Height_sigma]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # Configure MCMC model here! # # # # # # # # # #
    # Configure the stratigraphic MCMC model
    config = StratAgeModelConfiguration()
    # If you in doubt, you can probably leave these parameters as-is
    config.resolution = res*10 # Same units as sample height. Smaller is slower!
    config.bounding = 1.0 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
    (bottom, top) = extrema(smpl.Height)
    npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
    config.nsteps = 500 # Number of steps to run in distribution MCMC
    config.burnin = 100*npoints_approx # Number to discard
    config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Adding a little bit of systematic error on tectonic subsidence solves the ll = -Inf problem
    # Currently the systematic error is set to be = 1% of the sigma for the tectonic subsidence of the first (oldest) layer
    Sσ_corr = Sσ .+ ((Sμ[end-1]-Sμ[end])/100)

## --- Option a: Stratigraphic MCMC model without hiatus
    #Run the model
    (subsmdl, agedist, lldist, beta_tsdist, lldist_burnin, zdist) = SubsidenceStratMetropolis_Height(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ_corr, Beta_sigma/10, TS_Height_sigma/10)

    #= Code for storing and reading age-depth model results
    # Store and read results
    writedlm("agedist2.csv", agedist2, ",")
    writedlm("lldist2.csv", lldist2)
    writedlm("beta_t0dist2.csv", beta_t0dist2, ",")
    writedlm("lldist_burnin2.csv", lldist_burnin2)

    agedist3 = readdlm("agedist3.csv", ',', Float64)
    lldist3 = readdlm("lldist3.csv", ',', Float64)
    beta_t0dist3 = readdlm("beta_t0dist3.csv", ',', Float64)
    lldist_burnin3 = readdlm("lldist_burnin3.csv", ',', Float64)
    =#

    # Plot 1: Age-depth model (mean and 95% confidence interval for both model and data)
    hdl = plot([subsmdl.Age_025CI; reverse(subsmdl.Age_975CI)],[subsmdl.Height; reverse(subsmdl.Height)], fill=(round(Int,minimum(subsmdl.Height)),0.5,:blue), label="model", yflip = true, xflip = true)
    plot!(hdl, subsmdl.Age, subsmdl.Height, linecolor=:blue, label="", fg_color_legend=:white) # Center line
    t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(2*smpl.Age_sigma[t]),yerror=(2*smpl.Height_sigma[t]),label="data",seriestype=:scatter,color=:black)
    t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(2*smpl.Age_sigma[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
    any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),2*smpl.Age_sigma[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
    any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
    #plot!(hdl, [smpl.Age[1], smpl.Height[1]],[smpl.Age[3], smpl.Height[3]])
    savefig(hdl,"AAEx58_AgeDepth_UnknownHeight.pdf")
    #display(hdl)

    # Plot 2: Posterior distributions of beta
    post_beta = histogram(beta_tsdist[1,:], linecolor=nothing, alpha = 0.5, nbins=50)
    vline!([subsmdl.Beta_Median], linecolor=:blue, linestyle=:dot, linewidth = 3)
    vline!([1.71], linecolor=:black, linewidth = 3)
    savefig(post_beta, "AAEx58_PosteriorBeta_UnknownHeight.pdf")

    # Plot 3: Posterior distributions for the stratigraphic height of the start of thermal subsidence
    post_ts = histogram(beta_tsdist[2,:], linecolor=nothing, alpha = 0.5, nbins=50)
    vline!([subsmdl.TSHeight_Median], linecolor=:blue, linestyle=:dot, linewidth = 3)
    vline!([-3400], linecolor=:black, linewidth = 3)
    savefig(post_ts, "AAEx58_PosteriorTSHeight_UnknownHeight.pdf")

    # Plot 4: Age prediction for the active rifting-thermal subsidence transition
    interpolated_distribution = Array{Float64}(undef, size(beta_tsdist,2), size(agedist,2))
    tsdist_plot = repeat(beta_tsdist[2,:]', 500)
    ts_age = plot()
    # Interpolate between model horizons and plot
    # Right now I'm plotting them as a scatter plot, but maybe a contour/heatmap is more appropriate?
    for i=1:size(agedist,2)
        interpolated_distribution[:,i] = linterp1s(subsmdl.Height,agedist[:,i],beta_tsdist[2,:])
        scatter!(ts_age, interpolated_distribution[:,i], tsdist_plot[:,i], legend = false, alpha = 0.01, color = "blue")
    end
    display(ts_age)

    # Plot 5: Interpolate results for target horizons
    target_height = [-1000, -2000, -3000]
    interpolated_distribution = Array{Float64}(undef, length(target_height), size(agedist,2))
    # Interpolate between model horizons
    for i=1:size(agedist,2)
        interpolated_distribution[:,i] = linterp1s(subsmdl.Height,agedist[:,i],target_height)
    end
    # Calculate summary statistics
    predicted_medians = nanmedian(interpolated_distribution, dims=2)
    predicted_025CI = nanpctile(interpolated_distribution, 2.5, dims=2)
    predicted_975CI = nanpctile(interpolated_distribution, 97.5, dims=2)
    # Plotting (full distribution with mean as a dashed line)
    age_interp = histogram(interpolated_distribution[3,:], nbins=50, label="")
    vline!([predicted_medians[3]], linecolor = "black", linestyle=:dot, linewidth = 3)
    plot!(age_interp, xlabel="Age ($(smpl.Age_Unit)) for XXX", ylabel="Likelihood (unnormalized)")
    #savefig(age_interp, "InterpolatedAge_XXX.pdf")
    display(age_interp)


# # # # # # # # # # # # METHOD 2 # # # # # # # # # # # #
# TBC - Method 2 would follow the Allen and Allen treatment for a rift-drift transition in the middle of the succession
# In their workflow, I think we'll need to propage the entire Y matrix (from the DecompactBackstrip.jl script) to the MCMC functions, which I imagine would be too computationally expensive...
   
## --- End of File
