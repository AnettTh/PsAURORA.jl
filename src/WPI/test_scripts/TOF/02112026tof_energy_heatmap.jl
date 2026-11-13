using AURORA
using CairoMakie
using ProgressMeter
using JLD2


@load "src/WPI/data/02112026TOF_heatmap.jld2"

jldsave(outfile;
    TOF,                              # [E, ω], size (60, 60)
    E_grid  = collect(E_grid),        # [eV]
    ωs      = collect(ωs),            # [rad/s]
    ω_frac  = collect(ωs) ./ Ω_e_eq,  # ω/Ωₑ,eq, the axis you plot against
    Ω_e_eq,
    μ0, L, ϕ_eq,
    Bw      = amplitude.Bw,           # [T]
    Ψ0      = wave.Ψ0,
    t_chirp = modulation.t_duration,
    use_boris,
)

## Make Figure
data = load(outfile)

fig = Figure(size=(700, 500))
ax  = Axis(fig[1, 1];
    xlabel = "ω/Ωₑ",
    ylabel = "E [keV]",
    title  = "Time-of-flight using " * (data["use_boris"] ? "boris mover" : "integration") * "\n" *
             "α = $(round(rad2deg(acos(abs(data["μ0"]))), digits=1))°" *
             ", MLT $(round(longitude_to_MLT(data["ϕ_eq"]; degrees=false), digits=1))" *
             ", L = $(data["L"])",
)

hm = heatmap!(ax, data["ω_frac"], data["E_grid"][10:end] ./ 1e3, data["TOF"][10:end, :]';
    colormap  = :turbo,
    nan_color = :white,
)
Colorbar(fig[1, 2], hm; label="TOF [s]")

#save("src/WPI/test_scripts/TOF/figures/02112026TOF_energy_heatmap.png", fig)
fig
