## For subsidence modeling in extensional basins
## --- Load required pacages, install SubsidenceChron if required

    using SubsidenceChron
    using StatGeochem, Distributions, Plots

## --- Part 1: Decompaction and Backstripping

    # # # # # # # # # # # Enter stratigraphic information here! # # # # # # # # # # # #
    # Import the data file (.csv)
    data_csv = importdataset("examples/Svalbard_highres.csv",',')
    # Obtain stratigraphic info from the data file
    nLayers = length(data_csv["Thickness"])
    strat = StratData(nLayers)
    strat.Lithology          = data_csv["Lithology"]
    strat.Thickness         .= data_csv["Thickness"]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # OPTIONAL - Enter paleo water depth information here! # # # # # # # #
    # Import the data file (.csv)
    wd_csv = importdataset("examples/Svalbard_SequenceStrat.csv", ',')
    # Obtain paleo water depth info from the data file
    wd_nLayers = length(wd_csv["Thickness"])
    wd = WaterDepth(wd_nLayers)
    wd.DepthID    = wd_csv["Type"]
    wd.Thickness .= wd_csv["Thickness"]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # Configure MC model here! # # # # # # # # # #
    # Number of MCMC simulations
    nsims = 5000
    # Resolution for model horizons (in meters)
    res = 1.0 # meters
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
    p1 = plot(Sμ, alpha = 1, yflip = true, xflip = true, label = "Tectonic subsidence", color = "blue")
    #plot!(p1, reverse(subsidence_strat_depths), yflip = true, label = "Present-day thickness", color = "red")
    plot!(p1, Sₜ[:,2:end], alpha = 0.01, label = "", yflip = true, color = "blue", fg_color_legend=:white)
    savefig(p1, "DecompactBackstrip_higherres.pdf")


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
    config.nsteps = 5000 # Number of steps to run in distribution MCMC
    config.burnin = 1000*npoints_approx # Number to discard
    config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Let's see if adding just a little bit of systematic error on tectonic subsidence solves the ll = -Inf problem
    # Currently the systematic error is set to be = 1% of the sigma for the tectonic subsidence of the first (oldest) layer
    Sσ_corr = Sσ .+ ((Sμ[end-1]-Sμ[end])/100)

## --- Option a: Stratigraphic MCMC model without hiatus
    #Run the model
    (subsmdl, agedist, lldist, beta_t0dist, lldist_burnin) = SubsidenceStratMetropolis(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ_corr, Beta_sigma/10, T0_sigma/10)

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
    hdl = plot([subsmdl.Age_025CI; reverse(subsmdl.Age_975CI)],[subsmdl.Height; reverse(subsmdl.Height)], fill=(round(Int,minimum(subsmdl.Height)),0.5,:blue), label="model")
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
    savefig(hdl,"Svalbard_Result.pdf")
    #display(hdl)

    # Plot 2: Posterior distributions of beta
    post_beta = histogram(beta_t0dist[1,:], color="black", linecolor=nothing, alpha = 0.5, nbins=50)
    vline!([subsmdl.Beta_Median], linecolor = "black", linestyle=:dot, linewidth = 3)
    savefig(post_beta, "Svalbard_PosteriorBeta.pdf")

    # Plot 3: Posterior distributions of t₀
    post_t0 = histogram(beta_t0dist[2,:], color="black", linecolor=nothing, alpha = 0.5, nbins=50)
    vline!([subsmdl.T0_Median], linecolor = "black", linestyle=:dot, linewidth = 3)
    savefig(post_t0, "Svalbard_PosteriorT0.pdf")

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
    predicted_medians = nanmedian(interpolated_distribution, dims=1)
    predicted_025CI = nanpctile(interpolated_distribution, 2.5, dims=1)
    predicted_975CI = nanpctile(interpolated_distribution, 97.5, dims=1)
    # Plotting (full distribution with mean as a dashed line)
    age_interp = histogram(interpolated_distribution[:,2], nbins=50, label="")
    vline!([predicted_medians[2]], linecolor = "black", linestyle=:dot, linewidth = 3)
    plot!(age_interp, xlabel="Age ($(smpl.Age_Unit)) for the onset of Bitter Springs CIE", ylabel="Likelihood (unnormalized)")
    savefig(age_interp, "Svalbard_InterpolatedAge_BSA_Onset.pdf")
    display(age_interp)


## --- Plot 4 alternative ----
    #=
    # If target horizon is not in the range of the age-depth model...
    #Calculate ages for target horizons based on the ideal subsidence curve
    target_height = [0, 0.444, 0.863, 1.125, 1.962, 2.082, 2.137]
    target_index = findclosest(target_height, subsidence_strat_depths)
    target_subs = copy(Sₜ[target_index,:])
    beta_t0dist_filter = Array{Float64,2}(undef, 2, config.nsteps)
    idx = 1
    for i = 1:config.nsteps
        if beta_t0dist[1,i]<=subsmdl.Beta_975CI[1] && beta_t0dist[1,i]>=subsmdl.Beta_025CI[1] && beta_t0dist[2,i]<=subsmdl.T0_975CI[1] && beta_t0dist[2,i]>=subsmdl.T0_025CI[1]
            beta_t0dist_filter[1,idx] = beta_t0dist[1,i]
            beta_t0dist_filter[2,idx] = beta_t0dist[2,i]
            idx += 1
        end
    end
    beta_t0dist_95CI = beta_t0dist_filter[:,1:(idx-1)]

    beta_t0_sample_size = div((idx-1),30)
    Sₜ_sample_size = 30
    beta_t0_sampled = Array{Float64,2}(undef, 2, beta_t0_sample_size)
    predicted_ages = Array{Float64,2}(undef, length(target_index), beta_t0_sample_size*Sₜ_sample_size)
    τ = 50
    E₀ = 3165.6475782289444
    using ProgressMeter
    pgrs = Progress(beta_t0_sample_size, desc="Progress...")
    pgrs_interval = ceil(Int,sqrt(beta_t0_sample_size))
    for i = 1:beta_t0_sample_size
        beta_t0_sampled[:,i] = beta_t0dist_95CI[:,(rand(1:(idx-1)))]
        print(beta_t0_sampled[:,i])
        for j = 1:length(target_index)
            k = 1
            while k <= Sₜ_sample_size
                random_draw = rand(target_subs[j,:])
                temp_calc = (random_draw*pi)/(E₀*beta_t0_sampled[1,i]*sin(pi/beta_t0_sampled[1,i]))
                if temp_calc<1
                    predicted_ages[j,k+(i-1)*Sₜ_sample_size] = τ*log(1-temp_calc)+beta_t0_sampled[2,i]
                    k += 1
                end
            end
        end
        mod(i,pgrs_interval)==0 && update!(pgrs, i)
    end
    update!(pgrs, beta_t0_sample_size)

    predicted = histogram(predicted_ages[1,:], color="black", linecolor=nothing, alpha = 0.5)
    vline!([nanmedian(predicted_ages, dims=2)[1]], linecolor = "black", linestyle=:dot, linewidth = 3)
    savefig(predicted1, "Svalbard_BetaT0PredictedAge1.pdf")
    =#


    #= Code for test plots
    # Test Plot 1: see if younger strata can be accommodated by this one rifting event
    p6 = plot([631.2], [1100], xerror = 1.9, yerror = 15, yflip = true, xflip = true)
    plot!([541],[1323], yerror = 8)
    τ = 50
    E₀ = 3165.6475782289444
    for i = 1:5000
        xrange = 500:((beta_t0dist[2,i]-500)/1000):(beta_t0dist[2,i])
        tec_subs = (E₀*beta_t0dist[1,i]/pi)*sin(pi/beta_t0dist[1,i]).*(1 .-exp.(-(beta_t0dist[2,i] .-xrange)./τ))
        plot!(p6, xrange, tec_subs, alpha = 0.01, color = "blue", legend = false)
    end

    #Test Plot 2: see how well the model predicted beta and t0 matches the actual values:
    beta_range1 = 1:0.005:therm1.Param[1]+therm1.Sigma[1]*3 #three sigmas
    t0_range1 = therm1.Param[2]-therm1.Sigma[2]*3:1:therm1.Param[2]+therm1.Sigma[2]*3 #three sigmas
    beta_pdf1 = pdf.(Normal(therm1.Param[1], therm1.Sigma[1]), beta_range1)
    t0_pdf1 = pdf.(Normal(therm1.Param[2], therm1.Sigma[2]), t0_range1)
    beta_range2 = 1:0.005:therm2.Param[1]+therm2.Sigma[1]*3 #three sigmas
    beta_pdf2 = pdf.(Normal(therm2.Param[1], therm2.Sigma[1]), beta_range2)
    beta_range3 = 1:0.005:therm3.Param[1]+therm3.Sigma[1]*3 #three sigmas
    beta_pdf3 = pdf.(Normal(therm3.Param[1], therm3.Sigma[1]), beta_range3)

    h1 = fit(Histogram, beta_t0dist1[1,:], nbins = 50)
    h2 = fit(Histogram, beta_t0dist2[1,:], nbins = 50)
    h3 = fit(Histogram, beta_t0dist3[1,:], nbins = 50)
    h4 = fit(Histogram, beta_t0dist4[1,:], nbins = 50)
    h5 = fit(Histogram, beta_t0dist5[1,:], nbins = 50)

    g1 = fit(Histogram, beta_t0dist1[2,:], 372:2:536)
    g2 = fit(Histogram, beta_t0dist2[2,:], 372:2:556)
    g3 = fit(Histogram, beta_t0dist3[2,:], 370:2:546)
    g4 = fit(Histogram, beta_t0dist4[2,:], 370:2:630)
    g5 = fit(Histogram, beta_t0dist5[2,:], 370:2:726)

    testplot2_1 = plot(beta_range1, beta_pdf1, linecolor = "black", label = "prior", alpha = 0.5, xlims = (1,2), legend=:topright)
    vline!([therm1.Param[1]], label = "actual beta", linecolor = "black", alpha = 0.5, linestyle=:dot, linewidth = 2, xlims = (1,2)) #1.4+/-0.2
    histogram!(twinx(), beta_t0dist1[1,:], label = "posterior", color ="black", linecolor=nothing, xlims = (1,2), alpha = 0.5)
    plot!(beta_range2, beta_pdf2, linecolor = palette(:RdBu)[2], label = "prior", alpha = 0.5, xlims = (1,2), legend=:topright)
    histogram!(twinx(), beta_t0dist2[1,:], label = "posterior", color=palette(:RdBu)[2], linecolor=nothing, xlims = (1,2), alpha = 0.5)
    plot!(beta_range3, beta_pdf3, linecolor = palette(:RdBu)[4], label = "prior", alpha = 0.5, xlims = (1,2), legend=:topright)
    histogram!(twinx(), beta_t0dist3[1,:], label = "posterior", color=palette(:RdBu)[4], linecolor=nothing, xlims = (1,2), alpha = 0.5)
    histogram!(twinx(), beta_t0dist4[1,:], label = "posterior", color=palette(:RdBu)[10], linecolor=nothing, xlims = (1,2), alpha = 0.5)
    histogram!(twinx(), beta_t0dist5[1,:], label = "posterior", color=palette(:RdBu)[8], linecolor=nothing, xlims = (1,2), alpha = 0.5)

    testplot2_1_v2 = plot(beta_range1, beta_pdf1, linecolor = "black", label = "prior", legend=:topright, xlims = (1,2))
    vline!([therm1.Param[1]], label = "actual beta", linecolor = "black", linestyle=:dot, linewidth = 2, xlims = (1,2)) #1.4+/-0.2
    plot!(beta_range2, beta_pdf2, linecolor = palette(:RdBu)[2], label = "prior", legend=:topright, xlims = (1,2))
    plot!(beta_range3, beta_pdf3, linecolor = palette(:RdBu)[4], label = "prior", legend=:topright, xlims = (1,2))
    plot!(twinx(), [h1,h2,h3,h4,h5], seriestype=:steps, label = "posterior", color =["black",palette(:RdBu)[2],palette(:RdBu)[4],palette(:RdBu)[10],palette(:RdBu)[8]], xlims = (1,2))

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


    #Test Plot 3: see how ll and the different ll components change throughout the whole run:
    testplot3_1 = plot(lldist_subs, label = "Subsidence" , legend=:bottomright)
    plot!(testplot3_1, lldist_age, label = "Age")
    plot!(testplot3_1, lldist_tsparam, label = "TS Parameters")
    plot!(testplot3_1, lldist_height, label = "Height")
    savefig(testplot3_1, "LL_components_SensAna.pdf")

    testplot3_2 = plot(lldist1, label = "ll of accepted proposals" , legend=:bottomright)
    savefig(testplot3_2, "LL_SensAna_Test1.pdf")

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
## --- Option 2: Stratigraphic MCMC model including a hiatus with unknown duration - NOT WORKING YET!
    #=
    #Input strat height of hiatus
    hiatus_height = -500

    #Run the model with the additional term of "hiatus_height"
    (subsmdl, agedist, lldist, beta_t0dist, lldist_burnin) = SubsidenceStratMetropolis(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ, hiatus_height, 0.05, -5)

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl_hiatus = plot([subsmdl.Age_025CI; reverse(subsmdl.Age_975CI)],[subsmdl.Height; reverse(subsmdl.Height)], fill=(round(Int,minimum(subsmdl.Height)),0.5,:blue), label="model")
    plot!(hdl_hiatus, subsmdl.Age, subsmdl.Height, linecolor=:blue, label="", fg_color_legend=:white) # Center line
    t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl_hiatus, smpl.Age[t], smpl.Height[t], xerror=smpl.Age_sigma[t],label="data",seriestype=:scatter,color=:black)
    t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
    any(t) && plot!(hdl_hiatus, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age_sigma[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
    any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
    any(t) && plot!(hdl_hiatus, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),smpl.Age_sigma[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
    any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    plot!(hdl_hiatus, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
    savefig(hdl_hiatus,"AgeDepthModel_hiatus_failed_moreages.pdf")
    display(hdl_hiatus)
    =#

## --- Option 3: Stratigraphic MCMC model including hiata with known durations
    #=
    # Data about hiata
    nHiatuses = 1 # The number of hiata you have data for
    hiatus = HiatusData(nHiatuses) # Struct to hold data
    hiatus.Height         = [  -500]
    hiatus.Height_sigma   = [   0.0]
    hiatus.Duration       = [ 100.0]
    hiatus.Duration_sigma = [  20.0]

    # Run the model. Note the additional `hiatus` arguments
    @time (subsmdl, agedist, hiatusdist, lldist, beta_t0dist) = StratMetropolis(smpl, config, therm, subsidence_strat_depths, Sμ, Sσ, hiatus); sleep(0.5)

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([subsmdl.Age_025CI; reverse(subsmdl.Age_975CI)],[subsmdl.Height; reverse(subsmdl.Height)], fill=(minimum(subsmdl.Height),0.5,:blue), label="model")
    plot!(hdl, subsmdl.Age, subsmdl.Height, linecolor=:blue, label="", fg_color_legend=:white)
    plot!(hdl, smpl.Age, smpl.Height, xerror=smpl.Age_sigma*2,label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
    =#

    #=
    # More plots - NOT READY TO RUN - still need to readjust the format

    curve_ages = similar(Sμ)
    curve_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta*sin(pi/beta)))).+t0
    subsidence_strat_depths_m = copy(-subsidence_strat_depths[2:end]).*1000

    hiatus_age_depth_curve = plot([curve_1_ages; reverse(curve_4_ages)], [subsidence_strat_depths_m; reverse(subsidence_strat_depths_m)], fill = (round(Int,minimum(mdl_height)),0.5,:white), label = "Interpolated model results - 2σ error envelope")
    plot!(Age[2:4], Height[2:4], xerror=Age_sigma[2:4],label="Age data",seriestype=:scatter,color=:black, legend =:bottomleft, grid = false, dpi = 300)
    plot!(hiatus_age_depth_curve, curve_2_ages, subsidence_strat_depths_m, linecolor =:red, label = "lowest beta and t0 at 975CI")
    plot!(hiatus_age_depth_curve, curve_3_ages, subsidence_strat_depths_m, linecolor =:orange, label = "beta at 975CI and t0 at 025CI")
    plot!(hiatus_age_depth_curve, curve_4_ages, subsidence_strat_depths_m, linecolor =:yellow, label = "beta at 975CI and t0 at 975CI")
    #plot!(age_depth_curve, mdl_age, mdl_height, linecolor=:blue, label="", fg_color_legend=:white)

    Sᵣ = reverse(Sμ)
    target_ages = [815.29, 811.51, 810.7, 788.7, 752.7, 739.9]
    target_age_errors = [0.32, 0.25, 6.3, 0.24, 5.5, 6.1]
    target_heights = [-(2.232-0.25), -(2.232-0.425), -(2.232-0.435), -(2.232-1.14), -(2.232-1.93), -(2.232-2.125)]
    target_heights_m = copy(target_heights)*1000
    target_height_errors = [0.212, 0.1, 0.05, 0.2, 0.035, 0.035]
    target_height_errors_m = copy(target_height_errors)*1000
    target_subs_heights = [-Sᵣ[250], -Sᵣ[425], -Sᵣ[435], -Sᵣ[1140], -Sᵣ[1930], -Sᵣ[2125]]
    target_subs_height_errors = [(Sᵣ[250]-Sᵣ[38]), (Sᵣ[425]-Sᵣ[325]), (Sᵣ[435]-Sᵣ[385]), (Sᵣ[1140]-Sᵣ[940]), (Sᵣ[1930]-Sᵣ[1895]), (Sᵣ[2125]-Sᵣ[2090])]

    scatter!(hiatus_age_depth_curve, target_ages, target_heights_m, xerr = target_age_errors, yerr = target_height_errors_m, label = "correlated ages")
    scatter!(hiatus_thermal_subs_curve, target_ages, target_subs_heights, xerr = target_age_errors, yerr = target_subs_height_errors, label = "correlated ages")
    =#

## --- End of File
