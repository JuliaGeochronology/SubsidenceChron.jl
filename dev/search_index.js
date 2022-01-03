var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Chron","category":"page"},{"location":"#Chron","page":"Home","title":"Chron","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Chron.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Chron]","category":"page"},{"location":"#Chron.BootstrapCrystDistributionKDE-Tuple{AbstractArray{<:Number}}","page":"Home","title":"Chron.BootstrapCrystDistributionKDE","text":"BootstrapCrystDistributionKDE(data::AbstractArray, [sigma::AbstractArray]; cutoff=-0.05)\n\nBootstrap an estimate of the pre-eruptive (or pre-depositional) mineral crystallization distribution shape from a 1- or 2-d array of sample ages (one row per sample, one column per datum, padded with NaNs as needed) and an equivalent-size array of one-sigma uncertainties, using a kernel density estimate of stacked sample data.\n\nExamples\n\n# Bootstrap crystallization distribution for a synthetic dataset with ages\n# [1,2,3,...10] Ma and uncertainties of 1 Ma each\nBootstrappedDistribution = BootstrapCrystDistributionKDE(1:10, ones(10))\n\n\n\n\n\n","category":"method"},{"location":"#Chron.BootstrapCrystDistributionKDE-Tuple{ChronAgeData}","page":"Home","title":"Chron.BootstrapCrystDistributionKDE","text":"BootstrapCrystDistributionKDE(smpl::ChronAgeData; cutoff=-0.05)\n\nBootstrap an estimate of the pre-eruptive (or pre-depositional) mineral crystallization distribution shape from a Chron.ChronAgeData object containing data for several samples, using a kernel density estimate of stacked sample data.\n\nExamples\n\nBootstrappedDistribution = BootstrapCrystDistributionKDE(smpl)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.StratMetropolis-Tuple{ChronAgeData, StratAgeModelConfiguration}","page":"Home","title":"Chron.StratMetropolis","text":"StratMetropolis(smpl::ChronAgeData, [hiatus::HiatusData,] config::StratAgeModelConfiguration)\n\nRuns the main Chron.jl age-depth model routine for a stratigraphic set of samples defined by sample heights and simple Gaussian age constraints in the smpl struct, and an age-depth model configuration defined by the config struct.\n\nOptionally, if a hiatus struct is provided, the model will additionally incorporate information about the durations of known hiatuses at the specified model heights.\n\nExamples:\n\n(mdl, agedist, lldist) = StratMetropolis(smpl, config)\n\n(mdl, agedist, hiatusdist, lldist) = StratMetropolis(smpl, hiatus, config)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.StratMetropolis14C-Tuple{ChronAgeData, StratAgeModelConfiguration}","page":"Home","title":"Chron.StratMetropolis14C","text":"StratMetropolis14C(smpl::ChronAgeData, [hiatus::HiatusData,] config::StratAgeModelConfiguration)\n\nRuns the main Chron.jl age-depth model routine for a stratigraphic set of samples defined by sample heights and interpolated radiocarbon age constraints in the smpl struct, and an age-depth model configuration defined by the config struct.\n\nOptionally, if a hiatus struct is provided, the model will additionally incorporate information about the durations of known hiatuses at the specified model heights.\n\nExamples:\n\n(mdl, agedist, lldist) = StratMetropolis14C(smpl, config)\n\n(mdl, agedist, hiatusdist, lldist) = StratMetropolis14C(smpl, hiatus, config)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.StratMetropolisDist-Tuple{ChronAgeData, StratAgeModelConfiguration}","page":"Home","title":"Chron.StratMetropolisDist","text":"StratMetropolisDist(smpl::ChronAgeData, [hiatus::HiatusData,] config::StratAgeModelConfiguration)\n\nRuns the main Chron.jl age-depth model routine for a stratigraphic set of samples defined by sample heights and fitted asymmetric age distributions (bilinear_exponential) in the smpl struct, and an age-depth model configuration defined by the config struct.\n\nOptionally, if a hiatus struct is provided, the model will additionally incorporate information about the durations of known hiatuses at the specified model heights.\n\nExamples:\n\n(mdl, agedist, lldist) = StratMetropolisDist(smpl, config)\n\n(mdl, agedist, hiatusdist, lldist) = StratMetropolisDist(smpl, hiatus, config)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.bilinear_exponential-Tuple{Number, AbstractVector}","page":"Home","title":"Chron.bilinear_exponential","text":"bilinear_exponential(x::Number, p::AbstractVector)\n\nEvaluate the value of a \"bilinear exponential\" function defined by the parameters p at point x. This function, which can be used to approximate various asymmetric probability distributions found in geochronology, is defined by an exponential function with two log-linear segments joined by an atan sigmoid:\n\nℯ^p_1 + v x_s p_4^2 p_5^2 - (1-v) x_s p_4^2p_5^2\n\nwhere\n\nv = 12 - atan(x_s)π\n\nis a sigmoid, positive on the left-hand side, and\n\nx_s = (x - p_2)p_3^2\n\nis x scaled by mean and variance.\n\nThe elements of the parameter array p may be considered to approximately represent\n\np[1] # pre-exponential (log normaliation constant)\np[2] # mean (central moment)\np[3] # standard deviation\np[4] # sharpness\np[5] # skew\n\n\n\n\n\n","category":"method"},{"location":"#Chron.bilinear_exponential_ll-Tuple{Number, AbstractVector}","page":"Home","title":"Chron.bilinear_exponential_ll","text":"bilinear_exponential_ll(x, p)\n\nReturn the log likelihood corresponding to a bilinear_exponential distribution defined by the parameters p evaluate at point x.\n\nIf x is provided as a number and p as a vector, a single evaluation will be returned; if x is provided as a vector and p as a matrix, the sum over all i for each distribution p[:,i] evaluated at x[i] will be returned instead.\n\nSee also bilinear_exponential\n\n\n\n\n\n","category":"method"},{"location":"#Chron.check_dist_ll-Tuple{AbstractArray, AbstractArray, AbstractArray, Number, Number}","page":"Home","title":"Chron.check_dist_ll","text":"check_dist_ll(dist::AbstractArray, mu::AbstractArray, sigma::AbstractArray, tmin::Number, tmax::Number)\n\nReturn the log-likelihood of a set of mineral ages with means mu and uncertianty sigma being drawn from a given source (i.e., crystallization / closure) distribution dist, with terms to prevent runaway at low N.\n\nExamples\n\nmu, sigma = collect(100:0.1:101), 0.01*ones(11)\nll = check_dist_ll(MeltsVolcanicZirconDistribution, mu, sigma, 100, 101)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.metropolis_min!-Tuple{AbstractArray, Int64, AbstractArray, AbstractArray, AbstractArray}","page":"Home","title":"Chron.metropolis_min!","text":"metropolis_min!(tminDist::Array, nsteps::Int, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nIn-place (non-allocating) version of metropolis_min, fills existing array tminDist.\n\nRun a Metropolis sampler to estimate the minimum of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon eruption ages from a distribution of zircon crystallization ages.\n\nExamples\n\nmetropolis_min!(tminDist, 2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.metropolis_min-Tuple{Int64, AbstractArray, AbstractArray, AbstractArray}","page":"Home","title":"Chron.metropolis_min","text":"tminDist = metropolis_min(nsteps::Int, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nRun a Metropolis sampler to estimate the minimum of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon eruption ages from a distribution of zircon crystallization ages.\n\nExamples\n\ntmindist = metropolis_min(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.metropolis_minmax!-Tuple{AbstractArray, AbstractArray, AbstractArray, AbstractArray, Int64, AbstractArray, AbstractArray, AbstractArray}","page":"Home","title":"Chron.metropolis_minmax!","text":"metropolis_minmax!(tminDist, tmaxDist, llDist, acceptanceDist, nsteps::Int, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nIn-place (non-allocating) version of metropolis_minmax, filling existing arrays\n\nRun a Metropolis sampler to estimate the extrema of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon saturation and eruption ages from a distribution of zircon crystallization ages.\n\nExamples\n\nmetropolis_minmax!(tmindist, tmaxdist, lldist, acceptancedist, 2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.metropolis_minmax-Tuple{Int64, AbstractArray, AbstractArray, AbstractArray}","page":"Home","title":"Chron.metropolis_minmax","text":"metropolis_minmax(nsteps::Int, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nRun a Metropolis sampler to estimate the extrema of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon saturation and eruption ages from a distribution of zircon crystallization ages.\n\nExamples\n\ntmindist, tmaxdist, lldist, acceptancedist = metropolis_minmax(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.tMinDistMetropolis-Tuple{ChronAgeData, Int64, Int64, Array{Float64}}","page":"Home","title":"Chron.tMinDistMetropolis","text":"tMinDistMetropolis(smpl::ChronAgeData, nsteps::Int, burnin::Int, dist::Array{Float64})\n\nCalculate the minimum limiting (eruption/deposition) age of each sample defined in the smpl struct, using the metropolis_min function, assuming mineral ages for each sample are drawn from the source distribution dist. Fits a bilinear_exponential function to the resulting stationary distribution for each sample and stores the results in smpl.Params for use by the StratMetropolisDist function.\n\nExamples\n\nsmpl = tMinDistMetropolis(smpl, 5*10^5, 2*10^5, TriangularDistribution)\n\n\n\n\n\n","category":"method"}]
}
