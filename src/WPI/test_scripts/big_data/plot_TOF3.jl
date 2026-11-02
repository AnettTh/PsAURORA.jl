using AURORA
using CairoMakie

#===

Make figures with one ax per L and E(ω), for given α

===#

## Choosing my favorite values
function nearest(grid, target)
    i = argmin(abs.(grid .- target))
    return grid[i], i
end

favorite_L, i_favorite_L = nearest(data["L_grid"], 6.0)
favorite_μ, i_favorite_μ = nearest(data["μ_grid"], -cos(deg2rad(2)))
favorite_E, i_favorite_E = nearest(data["E_grid"], 25*1e3)
favorite_ω, i_favorite_ω = nearest(data["ω_frac"], 0.3)
