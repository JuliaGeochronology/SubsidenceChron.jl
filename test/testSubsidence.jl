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
Y_target = [0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.2  0.261474  0.887641  1.74732  1.20451  0.300026  0.474364  0.392066  0.770834  0.471032  1.0e-18
            0.45  1.06158  2.26118  2.70626  1.43681  0.758822  0.830523  1.1601  1.22082  1.0e-18  1.0e-18
            1.2  2.39165  3.17795  2.9117  1.86018  1.09798  1.59568  1.60168  1.0e-18  1.0e-18  1.0e-18
            2.5  3.29906  3.37943  3.31601  2.15361  1.86152  2.02924  1.0e-18  1.0e-18  1.0e-18  1.0e-18
            3.4  3.4997  3.78061  3.57605  2.91158  2.29078  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18
            3.6  3.90022  4.03355  4.32824  3.32741  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18
            4.0  4.15153  4.78422  4.73219  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18
            4.25  4.90183  5.1854  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18
            5.0  5.30235  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18
            5.4  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18  1.0e-18]

@test isapprox(Y, Y_target, atol=0.001)

## --- Test DecompactBackstrip function
strat_test = StratData(3)
strat_test.Lithology = ["Shale", "Sandstone", "Limestone"]
strat_test.Thickness .= [1000, 600, 400]
nsims_test = 5000
res_test = 20.0 # meters
(Sₜ_test, Sμ_test, Sσ_test, subsidence_strat_depths_test) = SubsidenceChron.DecompactBackstrip(strat_test, nsims_test, res_test)

subsidence_strat_depths_target = 0:res_test:2000
@test isapprox(subsidence_strat_depths_test, subsidence_strat_depths_target, atol=0.0001)

Sμ_target = [1098.927, 1093.926, 1088.810, 1083.576, 1078.223, 1072.751, 1067.158, 1061.443, 1055.603, 1049.638, 1043.546, 1037.324, 1030.971, 1024.484, 1017.862, 1011.102, 1004.200, 997.155, 989.963, 982.621, 975.125, 967.473, 959.658, 951.678, 943.528, 935.202, 926.695, 918.001, 909.114, 900.027, 890.731, 881.218, 871.480, 861.506, 851.284, 840.802, 830.046, 819.000, 807.646, 795.964, 783.930, 771.517, 758.693, 745.420, 731.652, 717.332, 702.390, 686.730, 670.225, 652.686, 633.802, 623.305, 612.731, 602.079, 591.350, 580.543, 569.657, 558.691, 547.646, 536.520, 525.313, 514.025, 502.655, 491.203, 479.668, 468.051, 456.349, 444.564, 432.694, 420.740, 408.701, 396.576, 384.367, 372.071, 359.689, 347.222, 334.667, 322.026, 309.298, 296.481, 283.577, 271.977, 260.189, 248.205, 236.016, 223.613, 210.984, 198.119, 185.005, 171.627, 157.970, 144.015, 129.744, 115.131, 100.150, 84.768, 68.948, 52.638, 35.778, 18.277, 0.000]
@test all(isapprox.(Sμ_target, Sμ_test, rtol=0.03))

Sσ_target = [98.7, 98.9, 99.1, 99.3, 99.5, 99.6, 99.7, 99.8, 99.9, 100.0, 100.0, 100.0, 100.0, 99.9, 99.8, 99.7, 99.6, 99.4, 99.2, 98.9, 98.6, 98.3, 97.9, 97.5, 97.1, 96.6, 96.1, 95.5, 94.9, 94.2, 93.5, 92.7, 91.9, 91.1, 90.2, 89.2, 88.2, 87.2, 86.1, 85.0, 83.9, 82.7, 81.5, 80.4, 79.3, 78.3, 77.4, 76.8, 76.5, 76.9, 78.4, 77.8, 77.3, 76.7, 76.2, 75.7, 75.3, 74.8, 74.4, 74.0, 73.7, 73.5, 73.3, 73.1, 73.1, 73.2, 73.3, 73.6, 74.0, 74.5, 75.2, 76.1, 77.1, 78.4, 79.8, 81.5, 83.4, 85.6, 88.0, 90.8, 93.8, 91.2, 88.5, 85.7, 82.7, 79.6, 76.4, 73.0, 69.4, 65.6, 61.6, 57.3, 52.8, 48.0, 42.9, 37.4, 31.4, 24.8, 17.6, 9.5, 0.0]
@test all(isapprox.(Sσ_target, Sσ_target, atol=10))

#Still need a test for Sσ - will run the model for longer/a few times and figure out the mean and variance of Sσ

## --- Test subsidence_ll function
a= 125000 # meters!
ρ_mantle = 3330
ρ_water = 1000
αᵥ = 3.28*10^(-5)
T_mantle = 1333
τ = 50 #Myr
E₀ = (4*a*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water)) # Also meters!

calc_ages = [119.0972886, 159.6099062, 182.1044563, 197.8668511, 210.0732074, 220.075277, 228.5755609, 235.9863097, 242.570271, 248.5052774, 253.9173321, 258.8990017, 263.520312, 267.8355388, 271.8876202, 275.7111255, 279.3343123, 282.7805892, 286.0695784, 289.2179037, 292.2397845, 295.1474908, 297.9516962, 300.6617565, 303.2859309, 305.8315604, 308.3052128, 310.7128024, 313.0596895, 315.3507635, 317.5905147, 319.7830954, 321.9323735, 324.0419813, 326.115359, 328.1557961, 330.1664717, 332.1504946, 334.1109458, 336.050925, 337.9736047, 339.8822946, 341.7805242, 343.6721515, 345.5615154, 347.4536569, 349.3546633, 351.2722381, 353.2167383, 355.2033331, 357.2577057, 358.3641581, 359.454486, 360.5292613, 361.5890253, 362.6342915, 363.665547, 364.683254, 365.6878516, 366.6797572, 367.6593674, 368.6270595, 369.5831926, 370.5281083, 371.462132, 372.3855734, 373.2987277, 374.2018758, 375.0952859, 375.9792132, 376.8539016, 377.7195836, 378.5764815, 379.4248081, 380.2647678, 381.0965574, 381.9203676, 382.7363844, 383.5447916, 384.3457736, 385.1395192, 385.8424306, 386.5467493, 387.252755, 387.9607487, 388.6710551, 389.3840265, 390.100048, 390.8195425, 391.5429789, 392.2708803, 393.0038368, 393.742521, 394.48771, 395.2403155, 396.0014282, 396.7723861, 397.5548851, 398.351178, 399.1644792, 400]
beta_t0_test = [1.4, 400]
subs_ll_test = SubsidenceChron.subsidence_ll(E₀, τ, Sμ_test[1:end-1], Sσ_test[1:end-1], calc_ages[1:end-1], beta_t0_test)
@test isapprox(subs_ll_test, 0, atol=0.1)

## --- Test SubsidenceStratMetropolis function
# Make an instance of a ChronAgeData object for nSamples
nSamples = 6
smpl = ChronAgeData(nSamples)
smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", "Sample 6") # Et cetera
smpl.Age          .= [ 220.075277,  328.1557961,  362.6342915,  388.6710551, 396.0014282, 398.351178] # Measured ages
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
therm = ThermalSubsidenceParameters()
therm.Param = [Beta, T0]
therm.Sigma = [Beta_sigma, T0_sigma]

# Configure the stratigraphic Monte Carlo model
config = StratAgeModelConfiguration()
config.resolution = res_test # Same units as sample height. Smaller is slower!
config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 40000  # Number of steps to run in distribution MCMC
config.burnin = 20000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps

@time (subsmdl_test, agedist_test, lldist_test, beta_t0dist_test, lldist_burnin_test) = SubsidenceStratMetropolis(smpl, config, therm, subsidence_strat_depths_test[1:end-1], Sμ_test[1:end-1], Sσ_test[1:end-1], 0.05, 10)

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

# Test mean age-depth model
mean_ages = nanmean(agedist_test, dim=2)
expected_mean_ages = [398.51, 397.72, 396.96, 396.19, 395.43, 394.67, 393.92, 393.17, 392.41, 391.65, 390.9, 390.15, 389.4, 388.65, 387.8, 386.92, 386.06, 385.19, 384.3, 383.42, 382.53, 381.65, 380.77, 379.89, 379.02, 378.14, 377.26, 376.37, 375.49, 374.61, 373.73, 372.85, 371.97, 371.08, 370.21, 369.33, 368.45, 367.56, 366.68, 365.8, 364.93, 364.05, 363.18, 362.31, 360.6, 358.87, 357.17, 355.45, 353.74, 352.0, 350.31, 348.61, 346.93, 345.24, 343.53, 341.79, 340.07, 338.35, 336.61, 334.87, 333.13, 331.4, 329.66, 327.94, 324.34, 320.83, 317.26, 313.63, 310.12, 306.46, 302.86, 299.29, 295.63, 291.97, 288.41, 284.76, 281.0, 277.38, 273.77, 270.23, 266.64, 262.96, 259.28, 255.54, 251.75, 248.12, 244.46, 240.89, 237.42, 233.89, 230.45, 227.02, 223.57, 220.26]
@test all(isapprox.(mean_ages, expected_mean_ages, atol=25))
# println(round.(mean_ages, digits=2))

# Test subsidence parameters
@test isapprox(only(subsmdl_test.Beta), 1.385317084366247, atol=0.1)
@test isapprox(only(subsmdl_test.Beta_025CI), 1.256171601851893, atol=0.1)
@test isapprox(only(subsmdl_test.Beta_975CI), 1.5277726246698682, atol=0.1)
@test isapprox(only(subsmdl_test.T0), 402.94742883910374, atol=5)
@test isapprox(only(subsmdl_test.T0_025CI), 390.44768044716994, atol=10)
@test isapprox(only(subsmdl_test.T0_975CI), 414.28186494452217, atol=10)

## Specify subsidence bottom and top
@time (subsmdl_test, agedist_test, lldist_test, beta_t0dist_test, lldist_burnin_test) = SubsidenceStratMetropolis(smpl, config, therm, subsidence_strat_depths_test[1:end-1], Sμ_test[1:end-1], Sσ_test[1:end-1], 0.05, 10, subsidencebottom=-1000, subsidencetop=-500)

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

# Test mean age-depth model
mean_ages = nanmean(agedist_test, dim=2)
expected_mean_ages = [398.5, 397.72, 396.95, 396.2, 395.46, 394.71, 393.95, 393.19, 392.43, 391.66, 390.89, 390.11, 389.36, 388.65, 387.8, 386.92, 386.08, 385.21, 384.34, 383.44, 382.56, 381.69, 380.82, 379.97, 379.11, 378.25, 377.35, 376.47, 375.58, 374.69, 373.78, 372.87, 371.97, 371.06, 370.17, 369.29, 368.4, 367.5, 366.6, 365.72, 364.85, 363.99, 363.15, 362.31, 360.63, 358.91, 357.16, 355.46, 353.7, 352.04, 350.24, 348.51, 346.78, 344.97, 343.24, 341.49, 339.76, 338.02, 336.23, 334.61, 332.92, 331.29, 329.63, 327.98, 325.37, 322.71, 319.84, 316.61, 313.29, 310.11, 306.75, 303.37, 299.87, 296.55, 292.74, 289.19, 285.59, 281.84, 277.8, 273.63, 269.47, 265.38, 261.12, 257.3, 253.1, 249.34, 245.41, 241.86, 238.09, 234.59, 230.73, 226.96, 223.6, 220.28]
@test all(isapprox.(mean_ages, expected_mean_ages, atol=25))
# println(round.(mean_ages, digits=2))

# Test subsidence parameters
@test isapprox(only(subsmdl_test.Beta), 1.1247054117222945, atol=0.1)
@test isapprox(only(subsmdl_test.Beta_025CI), 1.0148599420923872, atol=0.1)
@test isapprox(only(subsmdl_test.Beta_975CI), 1.4426669092339366, atol=0.1)
@test isapprox(only(subsmdl_test.T0), 384.9922091880086, atol=5)
@test isapprox(only(subsmdl_test.T0_025CI), 337.6653192413728, atol=10)
@test isapprox(only(subsmdl_test.T0_975CI), 488.08880450361073, atol=10)

## --- As above, but specify heights as positive numbers above bottom of section

# Make an instance of a ChronAgeData object for nSamples
nSamples = 6
smpl = ChronAgeData(nSamples)
smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", "Sample 6") # Et cetera
smpl.Age          .= [ 220.075277,  328.1557961,  362.6342915,  388.6710551, 396.0014282, 398.351178] # Measured ages
smpl.Age_sigma    .= [   1.0,    1.0,    1.0,    1.0,  1.0,  1.0] # Measured 1-σ uncertainties
smpl.Height       .= [  1900,   1300,    900,    300,  100,   40] # depths above bottom of section as positive
smpl.Height_sigma .= fill(0.01, nSamples) # Usually assume little or no sample height uncertainty
smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Age_Unit = "Ma" # Unit of measurement for ages
smpl.Height_Unit = "m"

@time (subsmdl_test, agedist_test, lldist_test, beta_t0dist_test, lldist_burnin_test) = SubsidenceStratMetropolis(smpl, config, therm, subsidence_strat_depths_test[1:end-1], Sμ_test[1:end-1], Sσ_test[1:end-1], 0.05, 10, subsidencebottom=1000, subsidencetop=1500)

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

# Test mean age-depth model
mean_ages = nanmean(agedist_test, dim=2)
expected_mean_ages = [398.61, 397.78, 396.99, 396.26, 395.55, 394.87, 394.18, 393.38, 392.74, 392.04, 391.27, 390.53, 389.65, 388.74, 387.93, 387.15, 386.38, 385.6, 384.81, 384.06, 383.33, 382.53, 381.73, 380.95, 380.14, 379.3, 378.42, 377.62, 376.79, 375.94, 375.09, 374.1, 373.26, 372.27, 371.32, 370.25, 369.26, 368.4, 367.37, 366.26, 365.2, 364.1, 363.21, 362.3, 360.66, 358.53, 356.82, 354.9, 352.73, 350.94, 348.94, 347.02, 345.05, 343.42, 341.77, 340.15, 338.74, 337.24, 335.72, 334.24, 332.81, 331.04, 329.31, 327.89, 324.6, 320.52, 316.51, 313.65, 309.71, 306.69, 303.72, 300.81, 297.91, 294.51, 291.23, 288.37, 285.97, 282.59, 278.41, 275.78, 273.22, 270.78, 267.98, 264.7, 261.29, 258.09, 254.56, 250.22, 246.32, 242.24, 237.7, 232.23, 226.76, 220.52]
@test all(isapprox.(mean_ages, expected_mean_ages, atol=25))
# println(round.(mean_ages, digits=2))

# Test subsidence parameters
@test isapprox(only(subsmdl_test.Beta), 1.1247054117222945, atol=0.1)
@test isapprox(only(subsmdl_test.Beta_025CI), 1.0148599420923872, atol=0.15)
@test isapprox(only(subsmdl_test.Beta_975CI), 1.3509359022349554, atol=0.15)
@test isapprox(only(subsmdl_test.T0), 384.7424612127314, atol=5)
@test isapprox(only(subsmdl_test.T0_025CI), 337.7980819566939, atol=10)
@test isapprox(only(subsmdl_test.T0_975CI), 489.37262155648244, atol=10)
