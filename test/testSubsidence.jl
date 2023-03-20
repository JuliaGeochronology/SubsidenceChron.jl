## -- Test decompact function
y = [0, 0.2, 0.45, 1.2, 2.5, 3.4, 3.6, 4, 4.25, 5, 5.4]
ϕ₀ = [0.63, 0.49, 0.63, 0.70, 0.49, 0.40, 0.20, 0.49, 0.05, 0.20]
c = [0.51, 0.27, 0.51, 0.71, 0.27, 0.60, 0.60, 0.27, 0.20, 0.30]
Y = fill(1E-18, (11, 11))
Y[:,1] .= y
Y[1,:] .= 0
for i = 2:10
    for j = i:10
        SubsidenceChron.decompact!(view(Y,:,i), y, ϕ₀[j], c[j], j, i)
    end
end
Y_target = vcat(
    [0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0],
    [0.2  0.261474  0.887641  1.74732  1.20451  0.300026  0.474364  0.392066  0.770834  0.471032  1.0e-18],
    [0.45  1.06158  2.26118  2.70626  1.43681  0.758822  0.830523  1.1601  1.22082  1.0e-18  1.0e-18],
    [1.2  2.39165  3.17795  2.9117  1.86018  1.09798  1.59568  1.60168  1.0e-18  1.0e-18  1.0e-18],
    [2.5  3.29906  3.37943  3.31601  2.15361  1.86152  2.02924  1.0e-18  1.0e-18  1.0e-18  1.0e-18],
    [3.4  3.4997  3.78061  3.57605  2.91158  2.29078  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18],
    [3.6  3.90022  4.03355  4.32824  3.32741  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18],
    [4.0  4.15153  4.78422  4.73219  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18],
    [4.25  4.90183  5.1854  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18],
    [5.0  5.30235  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18],
    [5.4  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18]
)
@test isapprox(Y, Y_target, atol=0.001)

## --- Test DecompactBackstrip function
strat_test = NewStratData(3)
strat_test.Lithology = ["Shale", "Sandstone", "Limestone"]
strat_test.Thickness .= [1, 0.6, 0.4]
nsims_test = 5000
res_test = 0.02
(Sₜ_test, Sμ_test, Sσ_test, model_strat_heights_test) = SubsidenceChron.DecompactBackstrip(strat_test, nsims_test, res_test)

model_strat_heights_target = 0:0.02:2
@test isapprox(model_strat_heights_test, model_strat_heights_target, atol=0.0001)

Sμ_target = [1098.927, 1093.926, 1088.810, 1083.576, 1078.223, 1072.751, 1067.158, 1061.443, 1055.603, 1049.638, 1043.546, 1037.324, 1030.971, 1024.484, 1017.862, 1011.102, 1004.200, 997.155, 989.963, 982.621, 975.125, 967.473, 959.658, 951.678, 943.528, 935.202, 926.695, 918.001, 909.114, 900.027, 890.731, 881.218, 871.480, 861.506, 851.284, 840.802, 830.046, 819.000, 807.646, 795.964, 783.930, 771.517, 758.693, 745.420, 731.652, 717.332, 702.390, 686.730, 670.225, 652.686, 633.802, 623.305, 612.731, 602.079, 591.350, 580.543, 569.657, 558.691, 547.646, 536.520, 525.313, 514.025, 502.655, 491.203, 479.668, 468.051, 456.349, 444.564, 432.694, 420.740, 408.701, 396.576, 384.367, 372.071, 359.689, 347.222, 334.667, 322.026, 309.298, 296.481, 283.577, 271.977, 260.189, 248.205, 236.016, 223.613, 210.984, 198.119, 185.005, 171.627, 157.970, 144.015, 129.744, 115.131, 100.150, 84.768, 68.948, 52.638, 35.778, 18.277, 0.000]
@test all(isapprox.(Sμ_target, Sμ_test, rtol=0.02))

#Still need a test for Sσ - will run the model for longer/a few times and figure out the mean and variance of Sσ

## --- Test subsidence_ll function
y_litho= 125000
ρ_mantle = 3330
ρ_water = 1000
αᵥ = 3.28*10^(-5)
T_mantle = 1333
τ = 50 #Myr
E₀ = (4*y_litho*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water))

calc_ages = [119.0972886, 159.6099062, 182.1044563, 197.8668511, 210.0732074, 220.075277, 228.5755609, 235.9863097, 242.570271, 248.5052774, 253.9173321, 258.8990017, 263.520312, 267.8355388, 271.8876202, 275.7111255, 279.3343123, 282.7805892, 286.0695784, 289.2179037, 292.2397845, 295.1474908, 297.9516962, 300.6617565, 303.2859309, 305.8315604, 308.3052128, 310.7128024, 313.0596895, 315.3507635, 317.5905147, 319.7830954, 321.9323735, 324.0419813, 326.115359, 328.1557961, 330.1664717, 332.1504946, 334.1109458, 336.050925, 337.9736047, 339.8822946, 341.7805242, 343.6721515, 345.5615154, 347.4536569, 349.3546633, 351.2722381, 353.2167383, 355.2033331, 357.2577057, 358.3641581, 359.454486, 360.5292613, 361.5890253, 362.6342915, 363.665547, 364.683254, 365.6878516, 366.6797572, 367.6593674, 368.6270595, 369.5831926, 370.5281083, 371.462132, 372.3855734, 373.2987277, 374.2018758, 375.0952859, 375.9792132, 376.8539016, 377.7195836, 378.5764815, 379.4248081, 380.2647678, 381.0965574, 381.9203676, 382.7363844, 383.5447916, 384.3457736, 385.1395192, 385.8424306, 386.5467493, 387.252755, 387.9607487, 388.6710551, 389.3840265, 390.100048, 390.8195425, 391.5429789, 392.2708803, 393.0038368, 393.742521, 394.48771, 395.2403155, 396.0014282, 396.7723861, 397.5548851, 398.351178, 399.1644792, 400]
beta_t0_test = [1.4, 400]
subs_ll_test = SubsidenceChron.subsidence_ll(E₀, τ, Sμ_test[1:end-1], Sσ_test[1:end-1], calc_ages[1:end-1], beta_t0_test)
@test isapprox(subs_ll_test, 0, atol=0.1)

## --- Test SubsidenceStratMetropolis function
# Make an instance of a ChronAgeData object for nSamples
nSamples = 6
smpl = NewChronAgeData(nSamples)
smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", "Sample 6") # Et cetera
smpl.Age          .= [ 220.075277,  328.1557961,  362.6342915,  388.6710551, 396.0014282, 398.351178]# Measured ages
smpl.Age_sigma    .= [   1.0,    1.0,    1.0,    1.0,  1.0,  1.0] # Measured 1-σ uncertainties
smpl.Height       .= [ -100,  -700, -1100, -1700,  -1900,  -1960] # Depths below surface should be negative
smpl.Height_sigma .= fill(0.01, nSamples) # Usually assume little or no sample height uncertainty
smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Age_Unit = "Ma" # Unit of measurement for ages
smpl.Height_Unit = "m"

# Enter initial guesses for the beta factor and thermal subsidence onset age and their uncertainties
Beta = 1.42
Beta_sigma = 0.2
T0 = 420
T0_sigma = 50
therm = NewThermalSubsidenceParameters()
therm.Param = [Beta, T0]
therm.Sigma = [Beta_sigma, T0_sigma]

# Configure the stratigraphic Monte Carlo model
config = NewStratAgeModelConfiguration()
config.resolution = 20 # Same units as sample height. Smaller is slower!
config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 20000 # Number of steps to run in distribution MCMC
config.burnin = 20000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps

@time (subsmdl_test, agedist_test, lldist_test, beta_t0dist_test, lldist_burnin_test) = SubsidenceStratMetropolis(smpl, config, therm, model_strat_heights_test[1:end-1], Sμ_test[1:end-1], Sσ_test[1:end-1], 0.05, 10)

# Test that results match expectation, within some tolerance
@test subsmdl_test.Age isa Vector{Float64}
@test subsmdl_test.Beta isa Vector{Float64}
@test subsmdl_test.T0 isa Vector{Float64}
@test agedist_test isa Matrix{Float64}
@test beta_t0dist_test isa Matrix{Float64}
@test lldist_test isa Vector{Float64}
@test lldist_burnin_test isa Vector{Float64}

# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist_test)])

@test isapprox(only(subsmdl_test.Beta), 1.385317084366247, atol=0.2)
@test isapprox(only(subsmdl_test.Beta_025CI), 1.256171601851893, atol=0.2)
@test isapprox(only(subsmdl_test.Beta_975CI), 1.5277726246698682, atol=0.2)
@test isapprox(only(subsmdl_test.T0), 402.94742883910374, atol=50)
@test isapprox(only(subsmdl_test.T0_025CI), 391.5831873363192, atol=50)
@test isapprox(only(subsmdl_test.T0_975CI), 419.939403901756, atol=50)
