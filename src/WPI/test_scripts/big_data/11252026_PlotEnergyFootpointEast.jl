using AURORA
using JLD2
using Dates
using CairoMakie

@load "src/WPI/data/11-13-2026/energy-time-position.jld2"

##
ref = LTPReference(params.station)
fp_enu = fill(NaN, size(FP))

for iΨ in Base.axes(FP, 5), iω in Base.axes(FP, 4), iμ in Base.axes(FP, 3), iE in Base.axes(FP, 2)
    fp_enu[:, iE, iμ, iω, iΨ] .= FP_to_ENU(FP[:, iE, iμ, iω, iΨ], ref, params.time)
end

## Everything below uses east and north as [E, μ, ω, Ψ₀], in km
east  = fp_enu[1, :, :, :, :] ./ 1e3
north = fp_enu[2, :, :, :, :] ./ 1e3
nE, nμ, nω, nΨ = size(east)

E_keV = collect(E_grid) ./ 1e3
ωs    = collect(ω_frac)
α_deg = collect(rad2deg.(α_grid))                      # use rad2deg.(α_grid) if it was saved in radians

## Coordinates of every entry, same shape as east
E_all = [E_keV[iE] for iE in 1:nE, iμ in 1:nμ, iω in 1:nω, iΨ in 1:nΨ]
ω_all = [ωs[iω] for iE in 1:nE, iμ in 1:nμ, iω in 1:nω, iΨ in 1:nΨ]

valid = @. isfinite(east) & isfinite(north)  # precipitated runs only

## One colour range per direction, shared by both pitch angles
range_east = extrema(east[valid])
range_north = (4.5, 4.804579595669551) #extrema(north[valid])

## Figure: one row per direction, one column per pitch angle
fig = Figure(size=(300 * nμ + 200, 600))
panel_axes = Matrix{Axis}(undef, 2, nμ)          # [direction, pitch angle]
sc_east = sc_north = nothing

for iμ in 1:nμ
    m = valid[:, iμ, :, :]
    x = E_all[:, iμ, :, :][m]
    y = ω_all[:, iμ, :, :][m]

    ax_e = Axis(fig[1, iμ]; title="α₀ = $(round(α_deg[iμ], digits=2))°")
    ax_n = Axis(fig[2, iμ]; xlabel="E [keV]")
    panel_axes[1, iμ], panel_axes[2, iμ] = ax_e, ax_n

    sc_east = scatter!(ax_e, x, y; color=east[:, iμ, :, :][m],
        colormap=:viridis, colorrange=range_east, markersize=12, alpha=0.8)
    sc_north = scatter!(ax_n, x, y; color=north[:, iμ, :, :][m],
        colormap=:plasma, colorrange=range_north, markersize=12, alpha=0.8)
end

# One colorbar per row, to the right of the last column
Colorbar(fig[1, nμ + 1], sc_east;  label="East [km]")
Colorbar(fig[2, nμ + 1], sc_north; label="North [km]")

# The y-label only on the first column, one label for each row
panel_axes[1, 1].ylabel = "ω / Ωₑ"
panel_axes[2, 1].ylabel = "ω / Ωₑ"

linkaxes!(vec(panel_axes)...)
for iμ in 1:nμ
    hidexdecorations!(panel_axes[1, iμ]; grid=false)          # x labels only on the bottom row
    iμ > 1 && hideydecorations!(panel_axes[1, iμ]; grid=false) # y labels only on the first column
    iμ > 1 && hideydecorations!(panel_axes[2, iμ]; grid=false)
end

Label(fig[0, 1:nμ + 1], "Footpoint offset from $(params.station) at $(params.time)";
    fontsize=16, font=:bold, tellwidth=false)

save("src/WPI/test_scripts/big_data/figures/11252026_EnergyFootpointEastNorth.png", fig)
