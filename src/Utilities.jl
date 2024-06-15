
# TODO: docstring
"""
```julia
find_formation_depths(formation, thickness)
```
Find the top and bottom depth of each unique formation given a list of stratigraphic units and thicknesses

### Examples:
```julia
ds = importdataset(joinpath("..", "examples", "Svalbard_highres.csv"), importas=:Tuple)
formations, tops bottoms = find_formation_depths(ds.Formation, ds.Thickness)
```
"""
function find_formation_depths(formation, thickness)
    depth = cumsum([0; thickness])
    unique_formations = unique(formation)
    unique_formation_tops = [depth[findfirst(x->x==n, formation)] for n in unique_formations]
    unique_formation_bottoms = [depth[findlast(x->x==n, formation)+1] for n in unique_formations]
    return unique_formations, unique_formation_tops, unique_formation_bottoms
end

"""
```julia
find_formation_heights(formation, thickness)
```
Find the top and bottom height of each unique formation given a list of stratigraphic units and thicknesses

### Examples:
```julia
ds = importdataset(joinpath("..", "examples", "Svalbard_highres.csv"), importas=:Tuple)
formations, tops bottoms = find_formation_heights(ds.Formation, ds.Thickness)
```
"""
function find_formation_heights(formation, thickness)
    height = reverse!(cumsum(reverse!([thickness; 0])))
    unique_formations = unique(formation)
    unique_formation_tops = [height[findfirst(x->x==n, formation)] for n in unique_formations]
    unique_formation_bottoms = [height[findlast(x->x==n, formation)+1] for n in unique_formations]
    return unique_formations, unique_formation_tops, unique_formation_bottoms
end