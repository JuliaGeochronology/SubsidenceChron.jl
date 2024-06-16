## For subsidence modeling in extensional basins - synthetic test case

## --- Load required pacages, install SubsidenceChron if required

    using SubsidenceChron
    using Distributions, Plots, Statistics, StatsBase

## --- Part 1: Decompaction and Backstripping

    # # # # # # # # # # # Enter stratigraphic information here! # # # # # # # # # # # #
    # Import the data file (.csv)
    data_csv = importdataset("examples/Test_DB_PerfectSubsidence.csv",',')
    # Obtain stratigraphic info from the data file
    nLayers = length(data_csv["Thickness"])
    strat = StratData(nLayers)
    strat.Lithology          = data_csv["Lithology"]
    strat.Thickness         .= data_csv["Thickness"]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # OPTIONAL - Enter paleo water depth information here! # # # # # # # #
    # Import the data file (.csv)
    wd_csv = importdataset("examples/PerfectSubsidence_SequenceStrat.csv", ',')
    # Obtain paleo water depth info from the data file
    wd_nLayers = length(wd_csv["Thickness"])
    wd = WaterDepth(wd_nLayers)
    wd.DepthID    = wd_csv["Type"]
    wd.Thickness .= wd_csv["Thickness"]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # Configure MC model here! # # # # # # # # # #
    # Number of MC simulations
    nsims = 1000
    # Resolution for model horizons (in meters)
    res = 20.0 # meters
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Run the decompaction and backstripping MC model
    # (wd input is optional)
    @time (Sₜ, Sμ, Sσ, subsidence_strat_depths) = DecompactBackstrip(strat, wd, nsims, res)

    #= Code for storing and reading decompaction + backstripping results - will be useful when testing the age-depth modeling part of the model
    # Store results
    using DelimitedFiles
    writedlm("St_highres.csv", Sₜ, ",")
    writedlm("St_mu_highres.csv", Sμ)
    writedlm("St_sigma_highres.csv", Sσ)
    writedlm("St_highres.txt", Sₜ)
    # Read results
    Sₜ = readdlm("St_highres.csv", ',', Float64)
    Sμ = readdlm("St_mu_highres.csv", ',', Float64)
    Sσ = readdlm("St_sigma_highres.csv", ',', Float64)
    =#

    # Plot results - tectonic subsidence in comparison with present day stratigraphic heights
    p1 = plot(Sμ, alpha = 1, yflip = true, xflip = true, label = "Tectonic subsidence", color = "blue")
    #plot!(p1, reverse(subsidence_strat_depths), yflip = true, label = "Present-day thickness", color = "red")
    plot!(p1, Sₜ[:,2:end], alpha = 0.01, label = "", yflip = true, color = "blue", fg_color_legend=:white)
    savefig(p1, "Test_DecompactBackstrip_higherres.pdf")


## --- Part 2a: Age-depth modeling

    # # # # # # # # # # # Enter age constraint (sample) information here! # # # # # # # # # # # #
    # Input the number of samples we wish to model (must match below)
    nSamples = 4
    # Make an instance of a Chron Section object for nSamples
    smpl = ChronAgeData(nSamples)
    smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
    smpl.Age          .= [   388.67,   362.63,    328.16,   220.08] # Measured ages
    smpl.Age_sigma    .= [        2,        2,         2,        2] # Measured 1-σ uncertainties
    smpl.Height       .= [    -1700,    -1100,      -700,     -100] # Depths below surface should be negative
    smpl.Height_sigma .= [     0.01,     0.01,      0.01,     0.01] # Usually assume little or no sample height uncertainty
    smpl.Age_Sidedness .= [       0,        0,         0,        0] # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Age_Unit = "Ma" # Unit of measurement for ages
    smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

    # IMPORTANT: smpl.Height must increase with increasing stratigraphic height
    # -- i.e., stratigraphically younger samples must be more positive. For this
    # reason, it is convenient to represent depths below surface as negative
    # numbers.
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # Enter thermal subsidence parameter priors here! # # # # # # # # # #
    # Enter initial guesses for the beta factor and thermal subsidence onset age and their uncertainties
    Beta = 1.8
    Beta_sigma = 0.2
    T0 = 410
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
    config.nsteps = 5000 # Number of steps to run in distribution MCMC
    config.burnin = 10000*npoints_approx # Number to discard
    config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## --- Option a: Stratigraphic MCMC model without hiatus
    #Run the model
    (subsmdl, agedist, lldist, beta_t0dist, lldist_burnin) = SubsidenceStratMetropolis(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ, 0.05, 10)

    #= Code for storing and reading age-depth model results
    # Store and read results
    writedlm("agedist.csv", agedist, ",")
    writedlm("lldist.csv", lldist)
    writedlm("beta_t0dist.csv", beta_t0dist, ",")
    writedlm("lldist_burnin.csv", lldist_burnin)

    agedist = readdlm("agedist.csv", ',', Float64)
    lldist = readdlm("lldist.csv", ',', Float64)
    beta_t0dist = readdlm("beta_t0dist.csv", ',', Float64)
    lldist_burnin = readdlm("lldist_burnin.csv", ',', Float64)
    =#

    # Plot 1: Age-depth model (mean and 95% confidence interval for both model and data)
    hdl = plot([subsmdl.Age_025CI; reverse(subsmdl.Age_975CI)],[subsmdl.Height; reverse(subsmdl.Height)], fill=(round(Int,minimum(subsmdl.Height)),0.5,:blue), label="model")
    plot!(hdl, subsmdl.Age, subsmdl.Height, linecolor=:blue, label="", fg_color_legend=:white) # Center line
    t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]),yerror=(2*smpl.Height_sigma[t]),label="data",seriestype=:scatter,color=:black)
    t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
    any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
    any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
    #plot!(hdl, [smpl.Age[1], smpl.Height[1]],[smpl.Age[3], smpl.Height[3]])
    savefig(hdl,"Test_AgeModel.pdf")
    #display(hdl)

    # Plot 2: Posterior distributions of beta
    post_beta = histogram(beta_t0dist[1,:], color="black", linecolor=nothing, alpha = 0.5, nbins=50)
    vline!([subsmdl.Beta_Median], linecolor = "black", linestyle=:dot, linewidth = 3)
    savefig(post_beta, "Test_PosteriorBeta.pdf")

    # Plot 3: Posterior distributions of t₀
    post_t0 = histogram(beta_t0dist[2,:], color="black", linecolor=nothing, alpha = 0.5, nbins=50)
    vline!([subsmdl.T0_Median], linecolor = "black", linestyle=:dot, linewidth = 3)
    savefig(post_t0, "Test_PosteriorT0.pdf")

    # Plot 4: Interpolate results for target horizons
    target_height = [-1840, -1320, -980, -360]
    interpolated_distribution = Array{Float64}(undef, length(target_height), size(agedist,2))
    # Interpolate between model horizons
    for i=1:size(agedist,2)
        interpolated_distribution[:,i] = linterp1s(subsmdl.Height,agedist[:,i],target_height)
    end
    # Calculate summary statistics
    predicted_medians = nanmedian(interpolated_distribution, dims=1)
    predicted_025CI = nanpctile(interpolated_distribution, 2.5, dims=1)
    predicted_975CI = nanpctile(interpolated_distribution, 97.5, dims=1)
    # Plotting (full distribution with mean as a dashed line)
    age_interp = histogram(interpolated_distribution[:,2], nbins=50, label="")
    vline!([predicted_medians[2]], linecolor = "black", linestyle=:dot, linewidth = 3)
    plot!(age_interp, xlabel="Age ($(smpl.Age_Unit))")
    savefig(age_interp, "Test_InterpolatedAge.pdf")
    display(age_interp)

## --- Code for test plots
    #=
    #Test Plot 2: see how well the model predicted beta and t0 matches the actual values:
    beta_range = 1:0.005:therm.Param[1]+therm.Sigma[1]*3 #three sigmas
    t0_range = therm.Param[2]-therm.Sigma[2]*3:1:therm.Param[2]+therm.Sigma[2]*3 #three sigmas
    beta_pdf = pdf.(Normal(therm.Param[1], therm.Sigma[1]), beta_range)
    t0_pdf = pdf.(Normal(therm.Param[2], therm.Sigma[2]), t0_range)

    h1 = fit(Histogram, beta_t0dist[1,:], nbins = 50)

    g1 = fit(Histogram, beta_t0dist[2,:], 372:2:536)

    testplot2_1 = plot(beta_range, beta_pdf, linecolor = "black", label = "prior", alpha = 0.5, xlims = (1,2), legend=:topright)
    vline!([therm.Param[1]], label = "actual beta", linecolor = "black", alpha = 0.5, linestyle=:dot, linewidth = 2, xlims = (1,2)) #1.4+/-0.2
    histogram!(twinx(), beta_t0dist[1,:], label = "posterior", color ="black", linecolor=nothing, xlims = (1,2), alpha = 0.5)

    testplot2_1_v2 = plot(beta_range, beta_pdf, linecolor = "black", label = "prior", legend=:topright, xlims = (1,2))
    vline!([therm.Param[1]], label = "actual beta", linecolor = "black", linestyle=:dot, linewidth = 2, xlims = (1,2)) #1.4+/-0.2

    testplot2_1_v3 = plot([therm1.Param[1]], [2.5], markerstrokecolor = "black", markerstrokewidth = 3, xlims = (1,2), xerror = therm1.Sigma[1], alpha = 0.5, ylims = (0,3), legend = false)
    vline!([therm1.Param[1]], label = "actual beta", linecolor = "black", linestyle=:dot, linewidth = 2, xlims = (1,2), alpha = 0.5) #1.4+/-0.2
    plot!([therm2.Param[1]], [2.6], markerstrokecolor = palette(:RdBu)[2], markerstrokewidth = 3, xlims = (1,2), xerror = therm2.Sigma[1], alpha = 0.5)
    plot!([therm3.Param[1]], [2.7], markerstrokecolor = palette(:RdBu)[4], markerstrokewidth = 3, xlims = (1,2), xerror = therm3.Sigma[1], alpha = 0.5)
    histogram!(twinx(), [beta_t0dist1[1,:],beta_t0dist2[1,:],beta_t0dist3[1,:],beta_t0dist4[1,:],beta_t0dist5[1,:]], seriescolor =["black",palette(:RdBu)[2],palette(:RdBu)[4],palette(:RdBu)[10],palette(:RdBu)[8]], linecolor=nothing, xlims = (1,2), alpha = 0.5)

    testplot2_1_v4 = scatter([therm1.Param[1]], [1000], ms = 3, markershape=:circle, markercolor = "black", markerstrokecolor = "black", markerstrokewidth = 2, xlims = (0.9,2), xerror = therm1.Sigma[1], ylims = (0,1100), legend = false)
    vline!([therm1.Param[1]], label = "actual beta", linecolor = "black", linestyle=:dot, linewidth = 2, xlims = (0.9,2)) #1.4+/-0.2
    plot!([therm2.Param[1]], [1040], ms = 3, markershape=:circle, markercolor = palette(:RdBu)[2], markerstrokecolor = palette(:RdBu)[2], markerstrokewidth = 2, xlims = (0.9,2), xerror = therm2.Sigma[1])
    plot!([therm3.Param[1]], [1080], ms = 3, markershape=:circle, markercolor = palette(:RdBu)[4], markerstrokecolor = palette(:RdBu)[4], markerstrokewidth = 2, xlims = (0.9,2), xerror = therm3.Sigma[1])
    plot!(h1, seriestype=:steps, label = "posterior", color ="black", xlims = (0.9,2))
    plot!(h2, seriestype=:steps, label = "posterior", color =palette(:RdBu)[2], xlims = (0.9,2))
    plot!(h3, seriestype=:steps, label = "posterior", color =palette(:RdBu)[4], xlims = (0.9,2))
    plot!(h4, seriestype=:steps, label = "posterior", color =palette(:RdBu)[10], xlims = (0.9,2))
    plot!(h5, seriestype=:steps, label = "posterior", color =palette(:RdBu)[8], xlims = (0.9,2))
    savefig(testplot2_1_v4, "Fig2a_beta_t0perturbs.pdf")

    vline!([nanmean(beta_t0dist1,dims=2)[1]], label = "posterior mean", linecolor = "blue", linewidth = 2)
    vline!([1.4], label = "actual beta", linecolor = "black", linewidth =2)
    savefig(testplot2_1,"beta_test6-3_0721.pdf")

    testplot2_2 = histogram(beta_t0dist1[2,:], label = "posterior", color = "black", linecolor = nothing, alpha = 0.5, legend=:topright)
    vline!([nanmean(beta_t0dist,dims=2)[2]], label = "posterior mean", linecolor = "blue", linewidth = 2) #120+/-20
    plot!(twinx(), t0_range, t0_pdf, linecolor = "red", label = "prior", legend=:topleft)
    vline!([T0], label = "prior mean", linecolor = "red", linewidth = 2)
    vline!([400], label = "actual beta", linecolor = "black", linewidth =2)
    savefig(testplot2_2,"t0_test6-3_0721.pdf")

    testplot2_2_v4 = scatter([therm1.Param[2]], [1200], ms = 3, markershape=:circle, markercolor = "black", markerstrokecolor = "black", markerstrokewidth = 2, xerror = therm1.Sigma[2], ylims = (0,1300), legend = false)
    vline!([therm1.Param[2]], label = "actual beta", linecolor = "black", linestyle=:dot, linewidth = 2) #1.4+/-0.2
    plot!([therm4.Param[2]], [1240], ms = 3, markershape=:circle, markercolor = palette(:RdBu)[10], markerstrokecolor = palette(:RdBu)[10], markerstrokewidth = 2, xerror = therm4.Sigma[2])
    plot!([therm5.Param[2]], [1280], ms = 3, markershape=:circle, markercolor = palette(:RdBu)[8], markerstrokecolor = palette(:RdBu)[8], markerstrokewidth = 2, xerror = therm5.Sigma[2])
    plot!(g1, seriestype=:steps, label = "posterior", color ="black", xlims = (280,550))
    plot!(g2, seriestype=:steps, label = "posterior", color =palette(:RdBu)[2], xlims = (280,550))
    plot!(g3, seriestype=:steps, label = "posterior", color =palette(:RdBu)[4], xlims = (280,550))
    plot!(g4, seriestype=:steps, label = "posterior", color =palette(:RdBu)[10], xlims = (280,550))
    plot!(g5, seriestype=:steps, label = "posterior", color =palette(:RdBu)[8], xlims = (280,550))
    savefig(testplot2_2_v4, "Fig2b_t0_betaperturbs.pdf")

    #Test Plot 4: see how well is the model doing matching the ideal subsidence curve:
    Sₜ_025CI = nanpctile(Sₜ, 2.5, dims = 2)[6:86]
    Sₜ_975CI = nanpctile(Sₜ, 97.5, dims = 2)[6:86]
    Sμ_crop = Sμ[6:86]
    Sμ_sample = [Sμ[6],Sμ[36],Sμ[56],Sμ[86]]
    testplot4 = plot([subsmdl5.Age; reverse(subsmdl5.Age)],[reverse(Sₜ_025CI); Sₜ_975CI], fill=(round(Int,minimum(Sₜ_025CI)),0.4,:blue), label="model")
    plot!(testplot4, subsmdl5.Age, reverse(Sμ_crop), linecolor=:blue, label="", fg_color_legend=:white) # Center line
    t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(testplot4, smpl.Age[t], reverse(Sμ_sample)[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]), label="data",seriestype=:scatter,color=:black)
    t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
    any(t) && plot!(testplot4, smpl.Age[t], reverse(Sμ_sample)[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
    any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, reverse(Sμ_sample[t])) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
    any(t) && plot!(testplot4, smpl.Age[t], reverse(Sμ_sample)[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
    any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, reverse(Sμ_sample)[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    plot!(testplot4, xlabel="Age ($(smpl.Age_Unit))", ylabel="Tectonic Subsidence ($(smpl.Height_Unit))")
    E₀ = 3165.6475782289444
    beta = nanmean(beta_t0dist5, dims=2)[1]
    t0 = nanmean(beta_t0dist5, dims=2)[2]
    Sμ_calculated = (E₀*beta/pi)*sin(pi/beta).*(1 .-exp.(-(t0 .-subsmdl5.Age)./50))
    plot!(testplot4, subsmdl5.Age, Sμ_calculated, linecolor=:purple, label="TS curve based on posterior")
    Sμ_ideal = (E₀*1.4/pi)*sin(pi/1.4).*(1 .-exp.(-(400 .-subsmdl5.Age)./50))
    plot!(testplot4, subsmdl5.Age, Sμ_ideal, linecolor=:red, label="actual TS curve")
    savefig(testplot4, "SubsidenceCurveComparison_SensAna_Test5.pdf")

    #Test Plot 5: see how the predicted ages for horizons-of-interest match the true values:
    testplot5 = histogram(predicted_ages[4,:], label = "predicted age for horizon 1", color ="blue", alpha = 0.5, legend=:topleft)
    vline!([393.74], label = "actual age for horizon 1", linecolor = "black", linewidth = 2)
    vline!([nanmedian(predicted_ages,dims=2)[4]], label = "posterior median", linecolor = "blue", linewidth = 2)
    savefig(testplot5,"PredictedAge_SensAna_Test1_h4.pdf")

    =#
    
## --- End of File
