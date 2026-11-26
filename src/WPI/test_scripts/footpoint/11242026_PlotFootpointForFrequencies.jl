using AURORA
using Dates
using ProgressMeter
using CairoMakie
using JLD2
##
data = load("src/WPI/data/11-26-2026/FootpointForFrequencies2.jld2")
FP_enu = data["FP_enu"]                      # [xyz, E, ω, Ψ₀]
E_grid = data["E_grid"]                      # [eV]
ω_frac = data["ω_frac"]
params = data["params"]

## Choose the energy: the grid point closest to this value
E_show = 40e3
iE = argmin(abs.(E_grid .- E_show))

## Slice the energy axis: both arrays are now [ω, Ψ₀]
east_km  = FP_enu[1, iE, :, :] ./ 1e3
north_km = FP_enu[2, iE, :, :] ./ 1e3

nω, nΨ = size(east_km)                       # two axes now, not three
ω_all  = [ω_frac[iω] for iω in 1:nω, iΨ in 1:nΨ]

valid = @. isfinite(east_km) & isfinite(north_km)
@show count(valid), length(valid)

## One axis, one energy
fig_E = Figure(size=(500, 300))
ax = Axis(fig_E[1, 1];
    xlabel = "East [km]",
    ylabel = "North [km]",
    title  = "Footpoints relative to $(params.station) at $(params.time)\n" *
             "E = $(round(E_grid[iE] / 1e3, digits=1)) keV, α₀ = $(params.α0)°",
    aspect = DataAspect(),
)

sc = scatter!(ax, east_km[valid], north_km[valid];
    color      = ω_all[valid],
    colormap   = :turbo,
    colorrange = extrema(ω_frac),
    markersize = 8,
)
Colorbar(fig_E[1, 2], sc; label="ω / Ωₑ")
#xlims!(ax, 1.07, 1.4)
#ylims!(ax, 4.75, 5.0)
#save("src/WPI/test_scripts/footpoint/figures/FP_and_frequency/11242026_FrequencyAndFP_E$(E_show/1e3)keV.png", fig_E)
