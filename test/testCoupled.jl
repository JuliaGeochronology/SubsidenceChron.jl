## --- Standard coupled eruption-deposition modelling

nSamples = 5 # The number of samples you have data for
smpl = ChronAgeData(nSamples)
smpl.Name      =  ("KJ08-157", "KJ04-75", "KJ09-66",  "KJ04-72", "KJ04-70",)
smpl.Height   .=  [     -52.0,      44.0,      54.0,      82.0,      93.0,]
smpl.Height_sigma .= [    3.0,       1.0,       3.0,       3.0,       3.0,]
smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Path = abspath("../examples/DenverUPbExampleData/") # Where are the data files?
smpl.inputSigmaLevel = 2 # i.e., are the data files 1-sigma or 2-sigma. Integer.
smpl.Age_Unit = "Ma" # Unit of measurement for ages and errors in the data files
smpl.Height_Unit = "cm" # Unit of measurement for Height and Height_sigma

# Remove outliers (if any)
smpl = screen_outliers(smpl, maxgap=50; make_plots)

# Distribution boostrapping from chron strat object
BootstrappedDistribution = BootstrapCrystDistributionKDE(smpl)
@test BootstrappedDistribution isa Vector{Float64}
@test BootstrappedDistribution[1] ≈ 1.0511396290122654

# Estimate the eruption age distributions for each sample  - - - - - - - -

# Configure distribution model here
distSteps = 10^5 # Number of steps to run in distribution MCMC
distBurnin = floor(Int,distSteps/2) # Number to discard

# Run MCMC to estimate saturation and eruption/deposition age distributions
@time tMinDistMetropolis(smpl,distSteps,distBurnin,BootstrappedDistribution; make_plots)

# Configure the stratigraphic Monte Carlo model
config = StratAgeModelConfiguration()
config.resolution = 10.0 # Same units as sample height. Smaller is slower!
config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 100000 # Number of steps to run in distribution MCMC
config.burnin = 100000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps

# Run the stratigraphic MCMC model
println("StratMetropolisDist:")
@time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config)

# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [66.07, 66.06, 66.05, 66.03, 66.02, 66.01, 65.99, 65.98, 65.97, 65.96, 65.94, 65.94, 65.93, 65.93, 65.89] atol=0.1
@test mdl.Age_025CI ≈ [66.04, 66.0, 65.98, 65.96, 65.95, 65.94, 65.94, 65.93, 65.93, 65.92, 65.92, 65.91, 65.9, 65.89, 65.82] atol=0.15
@test mdl.Age_975CI ≈ [66.1, 66.09, 66.09, 66.08, 66.08, 66.07, 66.06, 66.05, 66.03, 66.01, 65.97, 65.97, 65.96, 65.96, 65.95] atol=0.15
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test size(lldist) == (config.nsteps,)
@test !isnan(mean(lldist))

# Try adding systematic uncertainties too
smpl.Chronometer = (:UPb, :UPb, :ArAr, :UPb, :UPb)
systematic=SystematicUncertainty()
systematic.ArAr = 0.005/2 # One-sigma
systematic.UPb = 0.005/2 # One-sigma

# Run the stratigraphic MCMC model
println("StratMetropolisDist with systematic uncertainties:")
@time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config, systematic)

# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [66.07, 66.06, 66.05, 66.03, 66.02, 66.01, 65.99, 65.98, 65.97, 65.95, 65.94, 65.94, 65.93, 65.93, 65.89] atol=0.1
@test mdl.Age_025CI ≈ [66.04, 66.0, 65.98, 65.96, 65.95, 65.94, 65.93, 65.93, 65.92, 65.92, 65.92, 65.91, 65.9, 65.88, 65.82] atol=0.15
@test mdl.Age_975CI ≈ [66.1, 66.09, 66.09, 66.08, 66.08, 66.07, 66.06, 66.05, 66.03, 66.01, 65.97, 65.97, 65.96, 65.96, 65.95] atol=0.15
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test size(lldist) == (config.nsteps,)
@test !isnan(mean(lldist))

## --- As above, but with hiata

nHiatuses = 2 # The number of hiatuses you have data for
hiatus = HiatusData(nHiatuses) # Struct to hold data
hiatus.Height         = [-7.0, 35.0 ]
hiatus.Height_sigma   = [ 0.0,  0.0 ]
hiatus.Duration       = [ 0.3,  0.3 ]
hiatus.Duration_sigma = [ 0.05, 0.05]

# Run the model. Note the additional `hiatus` arguments
println("StratMetropolisDist with hiata:")
@time (mdl, agedist, hiatusdist, lldist) = StratMetropolisDist(smpl, hiatus, config)

# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [66.08, 66.08, 66.07, 66.07, 66.07, 66.02, 66.02, 66.01, 66.01, 65.94, 65.94, 65.93, 65.93, 65.92, 65.89] atol=0.1
@test mdl.Age_025CI ≈ [66.06, 66.05, 66.04, 66.03, 66.02, 65.94, 65.94, 65.94, 65.93, 65.91, 65.91, 65.9, 65.89, 65.88, 65.82] atol=0.15
@test mdl.Age_975CI ≈ [66.11, 66.1, 66.1, 66.1, 66.1, 66.08, 66.08, 66.08, 66.08, 65.98, 65.96, 65.96, 65.96, 65.96, 65.95] atol=0.15
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test size(hiatusdist) == (nHiatuses, config.nsteps)
@test mean(hiatusdist, dims=2) ≈ [0.064; 0.061;;] atol=0.3
@test size(lldist) == (config.nsteps,)
@test !isnan(mean(lldist))

## --- As above, but treat everything as a gaussian/weighted mean

# Tel tMinDistMetropolis to treat these as gaussians, using the first row of data file
smpl.Age_DistType .= 1
smpl.inputSigmaLevel = 1

# Run MCMC to estimate saturation and eruption/deposition age distributions
@time tMinDistMetropolis(smpl,distSteps,distBurnin,BootstrappedDistribution; make_plots)

# Run the stratigraphic MCMC model
println("StratMetropolisDist with fitted Gaussians:")
@time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config)
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [65.99, 65.98, 65.97, 65.96, 65.96, 65.95, 65.94, 65.93, 65.92, 65.91, 65.91, 65.9, 65.88, 65.86, 65.84] atol=0.1
@test mdl.Age_025CI ≈ [65.87, 65.86, 65.86, 65.85, 65.85, 65.84, 65.84, 65.83, 65.83, 65.83, 65.82, 65.8, 65.75, 65.73, 65.7] atol=0.15
@test mdl.Age_975CI ≈ [66.09, 66.09, 66.08, 66.08, 66.08, 66.07, 66.06, 66.05, 66.04, 66.02, 66.0, 65.99, 65.98, 65.97, 65.96] atol=0.15
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test size(lldist) == (config.nsteps,)
@test !isnan(mean(lldist))

## -- Test coupled with Isoplot.jl Pb-loss-aware eruption ages

nSamples = 3 # The number of samples you have data for
smpl = ChronAgeData(nSamples)
smpl.Name      =  ("KR18-04", "KR18-01", "KR18-05")
smpl.Height   .=  [      0.0,     100.0,     200.0] # Arbitrary heights
smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Path = abspath("../examples/ConcordiaExampleData/") # Where are the data files?
smpl.inputSigmaLevel = 1 # i.e., are the data files 1-sigma or 2-sigma. Integer.
smpl.Age_Unit = "Ma" # Unit of measurement for ages and errors in the data files
smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

smpl = screen_outliers(smpl, maxgap=100; make_plots, discordancemin=-5)

BootstrappedDistribution = BootstrapCrystDistributionKDE(smpl, tpbloss=0)
@test BootstrappedDistribution isa Vector{Float64}
@test BootstrappedDistribution[1] ≈ 0.3697041339884247

# Configure distribution model here
distSteps = 2*10^5 # Number of steps to run in distribution MCMC
distBurnin = floor(Int,distSteps/2) # Number to discard

# Run MCMC to estimate saturation and eruption/deposition age distributions
@time tMinDistMetropolis(smpl,distSteps,distBurnin,HalfNormalDistribution; make_plots)


# Configure the stratigraphic Monte Carlo model
config = StratAgeModelConfiguration()
config.resolution = 25.0 # Same units as sample height. Smaller is slower!
config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 100000 # Number of steps to run in distribution MCMC
config.burnin = 100000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps

println("StratMetropolisDist, Pb-loss-aware:")
@time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config)

# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [752.2, 752.14, 752.09, 752.03, 751.97, 751.68, 751.39, 751.09, 750.78] atol=0.3
@test mdl.Age_025CI ≈ [751.86, 751.8, 751.75, 751.7, 751.65, 750.93, 750.72, 750.61, 750.52] atol=0.7
@test mdl.Age_975CI ≈ [752.51, 752.47, 752.43, 752.37, 752.29, 752.2, 752.08, 751.85, 750.99] atol=0.7
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test size(lldist) == (config.nsteps,)
@test !isnan(mean(lldist))
