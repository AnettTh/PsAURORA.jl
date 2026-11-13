using AURORA
using CairoMakie
using JLD2

@load "src/WPI/data/11-03-2026/03112026-ducted_TOF_energy_heatmap.jld2"

## Make Figure
fig = Figure(size=(700,500))
ax = Axis(
    fig[1,1],
    xlabel="ω/Ωₑ",
    ylabel="α [°]",
    title="Time-of-flight using "
    * if use_boris "boris mover \n" else "integration \n " end
    * "E=$(E / 1e3) keV"
    * ", MLT $(round(longitude_to_MLT(ϕ_eq; degrees=false), digits = 1))"
    * ", L=$(L)",
)

hm = heatmap!(ax, ωs ./ Ω_e_eq, rad2deg.(acos.(-μ0_grid)), TOF'; colormap=:turbo, colorrange(2.0, 8.0))
Colorbar(fig[1,2], hm; label="TOF [s]")
##
#xlims!(ax, 0.8, 1.0)
#ylims!(ax, 100, 1000)

save("src/WPI/test_scripts/TOF/figures/03112026_ducted_TOF_alpha_heatmap.png", fig)
