using AURORA
using JLD2
using ProgressMeter
using Dates
using CairoMakie

@load "src/WPI/data/11-16-2026/energy-time-position.jld2"

##
ref = LTPReference(station)
fp_enu = fill(NaN, size(FP))

for iΨ in Base.axes(FP, 5), iω in Base.axes(FP, 4), iμ in Base.axes(FP, 3), iE in Base.axes(FP, 2)
    fp_enu[:, iE, iμ, iω, iΨ] .= FP_to_ENU(FP[:, iE, iμ, iω, iΨ], ref, t)
end

E_keV = collect(E_grid) ./ 1e3

east  = fp_enu[1, :, :, :, :] ./ 1e3      # n_E × n_μ × n_ω × n_Ψ, in km
north = fp_enu[2, :, :, :, :] ./ 1e3

# Energy of every entry, same shape as east and north
E_all = repeat(E_keV, 1, size(east, 2), size(east, 3), size(east, 4))

# Flatten and keep only precipitated cases (NaN = no resonance or bounced)
percipitated = findall(isfinite, vec(east))
Δeast = vec(east)[percipitated]
Δnorth = vec(north)[percipitated]
E = vec(E_all)[percipitated]

@show length(percipitated), length(east)      # number of footpoints out of the total

## Plot footpoints only
fig_FP = Figure(size = (700, 300))
ax_FP = Axis(fig_FP[1, 1];
    xlabel = "East [km]",
    ylabel = "North [km]",
    title  = "Footpoints relative to $(station) at $t",
    aspect = DataAspect(),
)

sc_FP = scatter!(
    ax_FP,
    Δeast,
    Δnorth;
    color      = E,
    colormap   = :turbo,
    colorrange = extrema(E_keV),
    markersize = 6,
)
Colorbar(fig_FP[1, 2], sc_FP; label = "E [keV]")
#ylims!(4.4, 4.9)
display(fig_FP)
#save("src/WPI/test_scripts/big_data/figures/11132026-FP/FP_only.png", fig_FP)
