var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Chron","category":"page"},{"location":"#Chron","page":"Home","title":"Chron","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Chron.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Chron]","category":"page"},{"location":"#Chron.BootstrapCrystDistributionKDE-Tuple{AbstractArray{var\"#s76\", N} where {var\"#s76\"<:Number, N}}","page":"Home","title":"Chron.BootstrapCrystDistributionKDE","text":"BootstrapCrystDistributionKDE(data::AbstractArray, [sigma::AbstractArray])\n\nBootstrap an estimate of thq;e pre-eruptive (or pre-deposition) mineral crystallization distribution shape from a 2-d array of sample ages (one row per sample, one column per datum, padded with NaNs as needed) and an equivalent-size array of one-sigma uncertainties, using a kernel density estimate of stacked sample data.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.BootstrapCrystDistributionKDE-Tuple{ChronAgeData}","page":"Home","title":"Chron.BootstrapCrystDistributionKDE","text":"BootstrapCrystDistributionKDE(smpl::ChronAgeData)\n\nBootstrap an estimate of the pre-eruptive (or pre-deposition) mineral crystallization distribution shape from a Chron.ChronAgeData object containing data for several samples, using a kernel density estimate of stacked sample data.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.MSWD-Tuple{Any, Any}","page":"Home","title":"Chron.MSWD","text":"MSWD(x, σ)\n\nReturn the Mean Square of Weighted Deviates (AKA the reduced chi-squared statistic) of a dataset with values x and one-sigma uncertainties σ\n\n\n\n\n\n","category":"method"},{"location":"#Chron.awmean-Tuple{Any, Any}","page":"Home","title":"Chron.awmean","text":"(wx, wσ, mswd) = awmean(x, σ)\n\nWeighted mean, absent the geochonologist's MSWD correction to uncertainty.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.bilinear_exponential-Tuple{Number, AbstractVector{T} where T}","page":"Home","title":"Chron.bilinear_exponential","text":"bilinear_exponential(x::Number, p::AbstractVector)\n\nEvaluate the value of a \"bilinear exponential\" function defined by the parameters p at point x. This function, which can be used to approximate various asymmetric probability distributions found in geochronology, is defined by an exponential function with two log-linear segments joined by an atan sigmoid:\n\nℯ^p_1 + v x_s p_4^2 p_5^2 - (1-v) x_s p_4^2p_5^2\n\nwhere\n\nv = 12 - atan(x_s)π\n\nis a sigmoid, positive on the left-hand side, and\n\nx_s = (x - p_2)p_3^2\n\nis x scaled by mean and variance.\n\nThe elements of the parameter array p may be considered to approximately represent\n\np[1] # pre-exponential (log normaliation constant)\np[2] # mean (central moment)\np[3] # standard deviation\np[4] # sharpness\np[5] # skew\n\n\n\n\n\n","category":"method"},{"location":"#Chron.bilinear_exponential_ll-Tuple{Number, AbstractVector{T} where T}","page":"Home","title":"Chron.bilinear_exponential_ll","text":"bilinear_exponential_ll(x, p)\n\nReturn the log likelihood corresponding to a bilinear_exponential distribution defined by the parameters p evaluate at point x.\n\nIf x is provided as a number and p as a vector, a single evaluation will be returned; if x is provided as a vector and p as a matrix, the sum over all i for each distribution p[:,i] evaluated at x[i] will be returned instead.\n\nSee also bilinear_exponential\n\n\n\n\n\n","category":"method"},{"location":"#Chron.check_dist_ll-Tuple{AbstractArray, AbstractArray, AbstractArray, Number, Number}","page":"Home","title":"Chron.check_dist_ll","text":"check_dist_ll(dist::AbstractArray, mu::AbstractArray, sigma::AbstractArray, tmin::Number, tmax::Number)\n\nReturn the log-likelihood of a set of mineral ages with means mu and uncertianty sigma being drawn from a given source (i.e., crystallization / closure) distribution dist, with terms to prevent runaway at low N.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.cntr-Tuple{AbstractArray}","page":"Home","title":"Chron.cntr","text":"cntr(edges::AbstractArray)\n\nGiven an array of bin edges, return a corresponding vector of bin centers\n\n\n\n\n\n","category":"method"},{"location":"#Chron.draw_from_distribution!-Tuple{Array{var\"#s19\", N} where {var\"#s19\"<:AbstractFloat, N}, AbstractArray{var\"#s20\", N} where {var\"#s20\"<:AbstractFloat, N}}","page":"Home","title":"Chron.draw_from_distribution!","text":"draw_from_distribution!(dist::AbstractArray{<:AbstractFloat}, x::Array{<:AbstractFloat})\n\nFill an existing variable x with random floating point numbers drawn from a continuous probability distribution specified by a vector dist defining the PDF curve thereof.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.draw_from_distribution-Tuple{AbstractArray{var\"#s3\", N} where {var\"#s3\"<:AbstractFloat, N}, Integer}","page":"Home","title":"Chron.draw_from_distribution","text":"x = draw_from_distribution(dist::AbstractArray{<:AbstractFloat}, n::Integer)\n\nDraw n random floating point numbers from a continuous probability distribution specified by a vector dist defining the PDF curve thereof.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.findclosest-Tuple{Any, Any}","page":"Home","title":"Chron.findclosest","text":"findclosest(source, target)\n\nReturn the index of the numerically closest value in the indexable collection target for each value in source. If muliple values are equally close, the first one is used\n\n\n\n\n\n","category":"method"},{"location":"#Chron.findclosestabove-Tuple{Any, Any}","page":"Home","title":"Chron.findclosestabove","text":"findclosestabove(source, target)\n\nReturn the index of the nearest value of the indexable collection target that is greater than (i.e., \"above\") each value in source. If no such values exist in target, returns an index of 0.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.findclosestbelow-Tuple{Any, Any}","page":"Home","title":"Chron.findclosestbelow","text":"findclosestbelow(source, target)\n\nReturn the index of the nearest value of the indexable collection target that is less than (i.e., \"below\") each value in source. If no such target values exist in target, returns an index of 0.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.findmatches-Tuple{Any, Any}","page":"Home","title":"Chron.findmatches","text":"findmatches(source, target)\n\nReturn the index of the first value in target (if any) that is equal to a given value in source for each value in source; else 0.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.gwmean-Tuple{Any, Any}","page":"Home","title":"Chron.gwmean","text":"(wx, wσ, mswd) = gwmean(x, σ)\n\nGeochronologist's weighted mean, with \"MSWD correction\" to uncertainty, i.e., wσ is increased by a factor of sqrt(mswd)\n\n\n\n\n\n","category":"method"},{"location":"#Chron.metropolis_min!-Tuple{AbstractArray, Int64, AbstractArray, AbstractArray, AbstractArray}","page":"Home","title":"Chron.metropolis_min!","text":"metropolis_min!(tminDist::Array, nsteps::Int, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nIn-place (non-allocating) version of metropolis_min, fills existing array tminDist.\n\nRun a Metropolis sampler to estimate the minimum of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon eruption ages from a distribution of zircon crystallization ages.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.metropolis_min-Tuple{Int64, AbstractArray, AbstractArray, AbstractArray}","page":"Home","title":"Chron.metropolis_min","text":"tminDist = metropolis_min(nsteps::Int, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nRun a Metropolis sampler to estimate the minimum of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon eruption ages from a distribution of zircon crystallization ages.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.metropolis_minmax!-Tuple{AbstractArray, AbstractArray, AbstractArray, AbstractArray, Int64, AbstractArray, AbstractArray, AbstractArray}","page":"Home","title":"Chron.metropolis_minmax!","text":"metropolis_minmax!(tminDist, tmaxDist, llDist, acceptanceDist, nsteps::Int, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nIn-place (non-allocating) version of metropolis_minmax, filling existing arrays\n\nRun a Metropolis sampler to estimate the extrema of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon saturation and eruption ages from a distribution of zircon crystallization ages.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.metropolis_minmax-Tuple{Int64, AbstractArray, AbstractArray, AbstractArray}","page":"Home","title":"Chron.metropolis_minmax","text":"(tminDist, tmaxDist, llDist, acceptanceDist) = metropolis_minmax(nsteps::Int, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nRun a Metropolis sampler to estimate the extrema of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon saturation and eruption ages from a distribution of zircon crystallization ages.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.midpointintegrate-Tuple{AbstractRange, AbstractArray}","page":"Home","title":"Chron.midpointintegrate","text":"midpointintegrate(bincenters, values)\n\nAdd up the area under a curve with y positions specified by a vector of values and x positions specfied by a vector of bincenters using midpoint integration.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanextrema-Tuple{Any}","page":"Home","title":"Chron.nanextrema","text":"nanextrema(A; dims)\n\nFind the extrema (maximum & minimum) of an indexable collection A, ignoring NaNs, optionally along a dimension specified by dims.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanmask!-Tuple{Any, Any}","page":"Home","title":"Chron.nanmask!","text":"nanmask!(mask, A)\n\nFill a Boolean mask of dimensions size(A) that is false wherever A is NaN\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanmask-Tuple{Any}","page":"Home","title":"Chron.nanmask","text":"nanmask(A)\n\nCreate a Boolean mask of dimensions size(A) that is false wherever A is NaN\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanmax-Tuple{Any, Any}","page":"Home","title":"Chron.nanmax","text":"nanmax(a,b)\n\nAs max(a,b), but if either argument is NaN, return the other one\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanmaximum-Tuple{Any}","page":"Home","title":"Chron.nanmaximum","text":"nanmaximum(A; dims)\n\nFind the largest non-NaN value of an indexable collection A, optionally along a dimension specified by dims.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanmean-Tuple{Any}","page":"Home","title":"Chron.nanmean","text":"nanmean(A, [W]; dims)\n\nIgnoring NaNs, calculate the mean (optionally weighted) of an indexable collection A, optionally along dimensions specified by dims.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanmedian-Tuple{Any}","page":"Home","title":"Chron.nanmedian","text":"nanmedian(A; dims)\n\nCalculate the median, ignoring NaNs, of an indexable collection A, optionally along a dimension specified by dims.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanmin-Tuple{Any, Any}","page":"Home","title":"Chron.nanmin","text":"nanmin(a,b)\n\nAs min(a,b), but if either argument is NaN, return the other one\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanminimum-Tuple{Any}","page":"Home","title":"Chron.nanminimum","text":"nanminimum(A; dims)\n\nAs minimum but ignoring NaNs: Find the smallest non-NaN value of an indexable collection A, optionally along a dimension specified by dims.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanrange-Tuple{Any}","page":"Home","title":"Chron.nanrange","text":"nanrange(A; dims)\n\nCalculate the range (maximum - minimum) of an indexable collection A, ignoring NaNs, optionally along a dimension specified by dims.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.nanstd-Tuple{Any}","page":"Home","title":"Chron.nanstd","text":"nanstd(A, [W]; dims)\n\nCalculate the standard deviation (optionaly weighted), ignoring NaNs, of an indexable collection A, optionally along a dimension specified by dims.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.norm_quantile-Tuple{Any}","page":"Home","title":"Chron.norm_quantile","text":"norm_quantile(F::Number)\n\nHow far away from the mean (in units of sigma) should we expect proportion F of the samples to fall in a standard Gaussian (Normal[0,1]) distribution\n\n\n\n\n\n","category":"method"},{"location":"#Chron.norm_width-Tuple{Any}","page":"Home","title":"Chron.norm_width","text":"norm_width(N::Number)\n\nHow dispersed (in units of sigma) should we expect a sample of N numbers drawn from a standard Gaussian (Normal[0,1]) distribution to be?\n\n\n\n\n\n","category":"method"},{"location":"#Chron.normcdf-Tuple{Any, Any, Any}","page":"Home","title":"Chron.normcdf","text":"normcdf(mu,sigma,x)\n\nCumulative density function of the Normal (Gaussian) distribution\n\n12 + erf(\fracx-μσ2)2\n\nwith mean mu and standard deviation sigma, evaluated at x.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.normpdf-Tuple{Any, Any, Any}","page":"Home","title":"Chron.normpdf","text":"normpdf(mu,sigma,x)\n\nProbability density function of the Normal (Gaussian) distribution\n\nℯ^-(x-μ)^2  (2σ^2)  σ2π\n\nwith mean mu and standard deviation sigma, evaluated at x\n\n\n\n\n\n","category":"method"},{"location":"#Chron.normpdf_ll-Tuple{Any, Any, Any}","page":"Home","title":"Chron.normpdf_ll","text":"normpdf_ll(mu, sigma, x)\n\nFast log likelihood corresponding to a Normal (Gaussian) distribution with mean mu and standard deviation sigma, evaluated at x.\n\nIf x, [mu, and sigma] are given as arrays, the sum of the log likelihood over all x will be returned.\n\nSee also normpdf\n\n\n\n\n\n","category":"method"},{"location":"#Chron.normproduct-NTuple{4, Any}","page":"Home","title":"Chron.normproduct","text":"normproduct(μ1, σ1, μ2, σ2)\n\nThe integral of the product of two normal distributions N[μ1,σ1] * N[μ2,σ2]. This is itself just another Normal distribution! Specifically, one with variance σ1^2 + σ2^2, evaluated at distance |μ1-μ2| from the mean\n\n\n\n\n\n","category":"method"},{"location":"#Chron.normproduct_ll-NTuple{4, Any}","page":"Home","title":"Chron.normproduct_ll","text":"normproduct_ll(μ1, σ1, μ2, σ2)\n\nLog likelihood corresponding to the integral of N[μ1,σ1] * N[μ2,σ2] As normproduct, but using the fast log likelihood of a Normal distribution\n\n\n\n\n\n","category":"method"},{"location":"#Chron.pctile-Tuple{Any, Any}","page":"Home","title":"Chron.pctile","text":"pctile(A, p; dims)\n\nFind the pth percentile of an indexable collection A, ignoring NaNs, optionally along a dimension specified by dims.\n\nA valid percentile value must satisfy 0 <= p <= 100.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.tMinDistMetropolis-Tuple{ChronAgeData, Int64, Int64, Array{Float64, N} where N}","page":"Home","title":"Chron.tMinDistMetropolis","text":"smpl = tMinDistMetropolis(smpl::ChronAgeData,nsteps::Int,burnin::Int,dist::Array{Float64})\n\nCalculate the minimum limiting (eruption/deposition) age of each sample in smpl using the metropolis_min function, assuming mineral ages for each sample are drawn from the source distribution dist. Fits a bilinearexponential function to the resulting stationary distribution for each sample.\n\n\n\n\n\n","category":"method"},{"location":"#Chron.trapz-Tuple{AbstractRange, AbstractArray}","page":"Home","title":"Chron.trapz","text":"trapz(edges, values)\n\nAdd up the area under a curve with y positions specified by a vector of values and x positions specfied by a vector of edges using trapezoidal integration. Bins need not be evenly spaced, though it helps.\n\n\n\n\n\n","category":"method"}]
}