# Dataset to use -- imported from examples folder
ds = importdataset(joinpath("..", "examples", "Svalbard_highres.csv"), importas=:Tuple)

# Depths (from top of section) to the top and bottom of each unique formation
formations, tops, bottoms = find_formation_depths(ds.Formation, ds.Thickness)
@test formations == ["Russøya Fm", "Backlundtoppen Fm", "Draken Fm", "Svanbergfjellet Fm", "Grusdievbreen Fm"]
@test tops == [0.0, 212.0, 690.0, 989.0, 1482.0]
@test bottoms == [212.0, 690.0, 989.0, 1482.0, 2138.0]

# Heights (from bottom of section) to the top and bottom of each unique formation
formations, tops, bottoms = find_formation_heights(ds.Formation, ds.Thickness)
@test formations == ["Russøya Fm", "Backlundtoppen Fm", "Draken Fm", "Svanbergfjellet Fm", "Grusdievbreen Fm"]
@test tops == [2138.0, 1926.0, 1448.0, 1149.0, 656.0]
@test bottoms == [1926.0, 1448.0, 1149.0, 656.0, 0.0]