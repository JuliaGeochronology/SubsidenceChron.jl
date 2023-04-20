# SubsidenceChron.jl
[![Docs][docs-dev-img]][docs-dev-url]
[![CI][ci-img]][ci-url]
[![CI (julia nightly)][ci-nightly-img]][ci-nightly-url]
[![codecov.io][codecov-img]][codecov-url]

A two-part package that incorporates and propagates uncertainties in the calculations of(1) decompaction and backstripping, and (2) thermal subsidence age-depth modeling in post-rift sedimentary basins. The first step uses Monte Carlo (MC) methods, and the second step uses Markov chain Monte Carlo (MCMC) methods.

The decompaction and backstripping MC model (part 1) is based on the workflow demonstrated in [Allen and Allen (2005)](https://www.wiley.com/en-ie/Basin+Analysis:+Principles+and+Application+to+Petroleum+Play+Assessment,+3rd+Edition-p-9780470673768). Users first input stratigraphic information (lithologies and thicknesses of measured section) and, optionally, water depth estimations (as distributions). Then the users need to define the resolution (that decompaction and backstripping calculations perform at) and the number of MC simulations. For each MC simulation, lithology-dependent parameters and water depth estimations are drawn from their respective distributions. The output of this first MC model describes how tectonic subsidence (with uncertainties) changes throughout deposition of the target sedimentary succession.

The second part of this package is an MCMC model that probes the age-depth relationship for strata deposited in a post-rift thermally subsiding basin. It builds upon the stratigraphic MCMC model of [Chron.jl](https://github.com/brenhinkeller/Chron.jl), which only uses stratigraphic superposition as the underlying principle (i.e. without making any assumptions about how sedimentation rate changes). Because this package is specifically developed for the thermal subsidence phase of the [Mckenzie](https://doi.org/10.1016/0012-821X(78)90071-7)-type basins, our MCMC model adds two extra log-likelihood terms to examine 1) the shape of the proposed tectonic subsidence-age curve, and 2) the rifting conditions (stretching factor &beta; and onset time for thermal subsidence t<sub>0</sub>). Inputs for this MCMC model include: 1) results propagated from the first MC model (i.e. tectonic subsidence through time), 2) existing age constraints, 3) prior distributions of &beta; and t<sub>0</sub>, and 4) user-defined MCMC model configurations. The outputs of this second MCMC model include both age estimations at all model horizons throughout the section, and updated distributions for &beta; and t<sub>0</sub>.

## Installation
SubsidenceChron.jl is written in the Julia programming language and registered on the General registry. To use this package, first download and install [Julia](https://julialang.org/downloads/). Then, SubsidenceChron.jl can be installed by typing ] in the REPL to enter the Julia package manager and typing:
```
pkg> add SubsidenceChron
```

## Standard Usage
To see the standard useage of this program, run [example/SubsidenceChron1.0Synthetic.jl](), which should follow the following steps:

### Load necessary packages
```julia
using SubsidenceChron

using StatGeochem, Distributions, Plots, Statistics, StatsBase
```
### Enter stratigraphic information
Here we use a synthetic stratigraphic succession as an example. This synthetic succession consists of, from the youngest to the oldest:

| Thickness (km) |  Lithology  | Water depth                         |
|:--------------:|:-----------:|:-----------------------------------:|
| 1              | Shale       | Between storm wave base and fair weather wave base|
| 0.6            | Sandstone   | Above fair weather wave base        |
| 0.4            | Limestone   | Above fair weather wave base        |

Stratigraphic information and (optional) water depth information should be stored in `.csv` files following the same format as those provided in the [examples/](examples/) folder. For example, the stratigraphic information `.csv` file should look like [examples/Test_DB_PerfectSubsidence.csv](examples/Test_DB_PerfectSubsidence.csv):
```
layer,Thickness,Lithology
3,1,Shale
2,0.6,Sandstone
1,0.4,Limestone
```
while the water depth information `.csv` file should look like [examples/PerfectSubsidence_SequenceStrat.csv](examples/PerfectSubsidence_SequenceStrat.csv):
```
layer,Thickness,Type,Notes
3,1,SWB,
2,0.6,FWWB,
1,0.4,FWWB,
```
Below is a list of currently acceptable lithology and water depth identifiers:

<table>
<tr><th>Lithology </th><th>Water depth </th></tr>
<tr><td>

| | Identifiers | |
|:---------:|:---------:|:---------:|
| Shale | Siltstone | Sandstone |
| Chalk | Limestone | Dolostone |
| Anhydrite | Quartzite | Diabase |
| Rhyolite | Diamictite | Conglomerate|
| Breccia  | Basalt | Andesite |

</td><td>

|Identifier | Meaning                                          | 
|:---------:|:------------------------------------------------:|
|Exposure   |Subaerial exposure surface                        |
|FWWB       |Above fair weather wave base                      |
|SWB        |Between storm wave base and fair weather wave base|
|BWB        |Below storm wave base                             |

</td></tr> </table>

Once the `.csv` files is in the correct format, the user can import them into the model:
```julia
# # # # # # # # # # # Enter stratigraphic information here! # # # # # # # # # # # #
# Import the data file (.csv)
data_csv = importdataset("examples/Test_DB_PerfectSubsidence.csv",',')
# Obtain stratigraphic info from the data file
nLayers = length(data_csv["Thickness"])
strat = NewStratData(nLayers)
strat.Lithology          = data_csv["Lithology"]
strat.Thickness         .= data_csv["Thickness"]
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # OPTIONAL - Enter paleo water depth information here! # # # # # # # # 
# Import the data file (.csv)
wd_csv = importdataset("examples/PerfectSubsidence_SequenceStrat.csv", ',')
# Obtain paleo water depth info from the data file
wd_nLayers = length(wd_csv["Thickness"])
wd = NewWaterDepth(wd_nLayers)
wd.DepthID    = wd_csv["Type"]
wd.Thickness .= wd_csv["Thickness"] 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
```
If water depth information is not applicable/does not exist for the user's case study, the code block for entering water depth information should be commented out.

#### Configure and run the decompaction and backstripping MC model
After importing data, the user need to configure the MC model (i.e. specify the number of MC simulations and the resolution):
```julia
# # # # # # # # # # Configure MC model here! # # # # # # # # # #
# Number of MC simulations
nsims = 1000 
# Resolution for model horizons (in km)
res = 0.02
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Run the decompaction and backstripping MC model
# (wd input is optional)
@time (Sₜ, Sμ, Sσ, model_strat_heights) = DecompactBackstrip(strat, wd, nsims, res)
```
If water depth information is not applicable/does not exist for the user's case study, the user should exclude the input `wd` when calling the `DecompactBackstrip` function:
```julia
    @time (Sₜ, Sμ, Sσ, model_strat_heights) = DecompactBackstrip(strat, nsims, res)
```
The outputs of this first MC model are

1. tectonic subsidence (throughout the deposition of the target succession) from all MC simulations (`Sₜ`);
2. tectonic subsidence mean, averaging across all MC simulations (`Sμ`);
3. tectonic subsidence standard deviation (`Sσ`);
4. model stratigraphic heights (spacing = `res`) (`model_strat_heights`).

### Enter age constraints and priors for thermal subsidence parameters
For the synthetic test case, we used the "real" ages for four randomly-picked horizons as age constraints. The "real" ages are calculated using &beta; = 1.85 and t<sub>0</sub> = 400 Ma, and assuming that this succession was deposited in a perfect McKenzie-type thermally subsiding basin. 
```julia
# # # # # # # # # # # Enter age constraint (sample) information here! # # # # # # # # # # # #
# Input the number of samples we wish to model (must match below)
nSamples = 4
# Make an instance of a Chron Section object for nSamples
smpl = NewChronAgeData(nSamples)
smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
smpl.Age          .= [   388.67,   362.63,    328.16,   220.08] # Measured ages
smpl.Age_sigma    .= [        2,        2,         2,        2] # Measured 1-σ uncertainties
smpl.Height       .= [    -1700,    -1100,      -700,     -100] # Depths below surface should be negative
smpl.Height_sigma .= [     0.01,     0.01,      0.01,     0.01] # Usually assume little or no sample height uncertainty
smpl.Age_Sidedness .= [       0,        0,         0,        0] # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Age_Unit = "Ma" # Unit of measurement for ages
smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

# # # # # # # # # # Enter thermal subsidence parameter priors here! # # # # # # # # # #
# Enter initial guesses for the beta factor and thermal subsidence onset age and their uncertainties
Beta = 1.8
Beta_sigma = 0.2
T0 = 410
T0_sigma = 50

therm = NewThermalSubsidenceParameters()
therm.Param = [Beta, T0]
therm.Sigma = [Beta_sigma, T0_sigma]
```
Note that when entering the age constraints, `smpl.Height` must increase with increasing stratigraphic height.

#### Configure and run the thermal subsidence MCMC model
To run the MCMC model, the user should first configure the model by setting appropiate values for `config.bounding`, `config.nsteps`, and `config.burnin`.
```julia
 # # # # # # # # # # Configure MCMC model here! # # # # # # # # # #
# Configure the stratigraphic MCMC model
config = NewStratAgeModelConfiguration()
# If you in doubt, you can probably leave these parameters as-is
config.resolution = res*1000 # Same units as sample height. Smaller is slower!
config.bounding = 1.0 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 5000 # Number of steps to run in distribution MCMC 
config.burnin = 10000*npoints_approx # Number to discard 
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
```
Then we call the `SubsidenceStratMetropolis` function. The user is welcomed to change the last two inputs of this function, which represents the perturbation to &beta; and t<sub>0</sub> for the initial proposal - changing these two values shouldn't affect the posterior distribution.
```julia
#Run the model
(subsmdl, agedist, lldist, beta_t0dist, lldist_burnin) = SubsidenceStratMetropolis(smpl, config, therm, model_strat_heights, Sμ, Sσ, 0.05, 10)
```
The outputs of this second MCMC model are

1. summary statistics for the age predictions of all model horizons and for the thermal subsidence parameters (`subsmdl`);
2. age predictions of all model horizons from all recorded MCMC runs (`agedist`);
3. log likelihood from all recorded MCMC runs (`lldist`);
4. subsidence parameters from all recorded MCMC runs (`beta_t0dist`);
5. log likelihood from all runs during the burn-in stage (`lldist_burnin`).

### Plot results
1. Plot the age-depth model;
```julia
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
```
2. Plot the posterior distribution of &beta;;
```julia
# Plot 2: Posterior distributions of beta
post_beta = histogram(beta_t0dist[1,:], color="black", linecolor=nothing, alpha = 0.5, nbins=50)
vline!([subsmdl.Beta_Median], linecolor = "black", linestyle=:dot, linewidth = 3)
savefig(post_beta, "Test_PosteriorBeta.pdf")
```
3. Plot the posterior distribution of t<sub>0</sub>;
```julia
# Plot 3: Posterior distributions of t₀
post_t0 = histogram(beta_t0dist[2,:], color="black", linecolor=nothing, alpha = 0.5, nbins=50)
vline!([subsmdl.T0_Median], linecolor = "black", linestyle=:dot, linewidth = 3)
savefig(post_t0, "Test_PosteriorT0.pdf")
```
4. Calculate and plot the age predictions (as distributions) for undated horizons of interest.
```julia
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
#display(age_interp)
```



[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JuliaGeochronology.github.io/SubsidenceChron.jl/dev/
[ci-img]: https://github.com/JuliaGeochronology/SubsidenceChron.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/JuliaGeochronology/SubsidenceChron.jl/actions/workflows/CI.yml
[ci-nightly-img]:https://github.com/JuliaGeochronology/SubsidenceChron.jl/workflows/CI%20(Julia%20nightly)/badge.svg
[ci-nightly-url]:https://github.com/JuliaGeochronology/SubsidenceChron.jl/actions/workflows/CI-julia-nightly.yml
[codecov-img]: http://codecov.io/github/JuliaGeochronology/SubsidenceChron.jl/coverage.svg?branch=main
[codecov-url]: http://codecov.io/github/JuliaGeochronology/SubsidenceChron.jl?branch=main
