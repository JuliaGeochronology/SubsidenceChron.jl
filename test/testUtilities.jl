
ds = importdataset(joinpath("..", "examples", "Svalbard_highres.csv"), importas=:Tuple)
formations, tops, bottoms = find_formation_depths(ds.Formation, ds.Thickness)
@test formations == ["Russøya Fm", "Backlundtoppen Fm", "Draken Fm", "Svanbergfjellet Fm", "Grusdievbreen Fm"]
@test tops == [0.0, 212.0, 690.0, 989.0, 1482.0]
@test bottoms == [212.0, 690.0, 989.0, 1482.0, 2138.0]