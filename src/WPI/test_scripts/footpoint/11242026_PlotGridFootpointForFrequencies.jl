using AURORA
using Dates
using ProgressMeter
using CairoMakie
using JLD2

@load "src/WPI/data/11-24-2026/FootpointForFrequencies2.jld2"

east_km  = FP_enu[1, :, :, :] ./ 1e3
north_km = FP_enu[2, :, :, :] ./ 1e3

nE, nω, nΨ = size(east_km)
ω_all = [ω_frac[iω] for iE in 1:nE, iω in 1:nω, iΨ in 1:nΨ]

valid = @. isfinite(east_km) & isfinite(north_km)

##
n_rows, n_cols = 3, 2
fig = Figure(size=(900, 500))
panel_axes = Matrix{Axis}(undef, n_rows, n_cols)
sc = nothing

for iE in 1:nE
    row, col = fldmod1(iE, n_cols)
    ax = Axis(fig[row, col]; xlabel="East [km]", ylabel="North [km]", aspect=DataAspect())
    panel_axes[row, col] = ax

    m = valid[iE, :, :]
    sc = scatter!(ax, east_km[iE, :, :][m], north_km[iE, :, :][m];
        color      = ω_all[iE, :, :][m],
        colormap   = :turbo,
        colorrange = extrema(ω_frac),
        markersize = 7,
    )

    text!(ax, 0.03, 0.95; text="E = $(round(E_grid[iE] / 1e3, digits=1)) keV",
        space=:relative, align=(:left, :top), fontsize=12)
end

Colorbar(fig[1:n_rows, n_cols + 1], sc; label="ω / Ωₑ")

linkaxes!(vec(panel_axes)...)
for row in 1:n_rows, col in 1:n_cols
    row < n_rows && hidexdecorations!(panel_axes[row, col]; grid=false)
    col > 1      && hideydecorations!(panel_axes[row, col]; grid=false)
end

Label(fig[0, 1:n_cols + 1], "Footpoints relative to $(params.station) at $(params.time), α₀ = $(params.α0)°";
    fontsize=16, font=:bold, tellwidth=false)


save("src/WPI/test_scripts/footpoint/figures/11242026_FrequencyAndFP.png", fig)
